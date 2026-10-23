"""Versioned design inputs, codon metrics and explicit vector size accounting."""

import hashlib
import json
import math
from pathlib import Path

from Bio.Data import CodonTable
from Bio.SeqFeature import Location
from pydantic import ValidationError

from payload_compiler.exceptions import InputError
from payload_compiler.models import (
    CodonMetrics,
    CodonProfile,
    OptimizeSettings,
    Payload,
    ProtectedRegion,
    VectorArchitecture,
)
from payload_compiler.profiles import (
    HOST_PROFILE_VERSIONS,
    HOST_PROFILES,
    RETIRED_PROFILES,
)
from payload_compiler.validation import coding_region_positions, extract_cds


def read_json_model(path: Path, model):
    try:
        return model.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError, ValidationError) as exc:
        raise InputError(f"Invalid {model.__name__} file {path.name}: {exc}") from exc


def profile_digest(profile: CodonProfile) -> str:
    canonical = json.dumps(profile.model_dump(mode="json"), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def get_profile(host: str, custom: CodonProfile | None = None) -> CodonProfile:
    if custom is not None:
        return custom
    if host in RETIRED_PROFILES:
        raise InputError(RETIRED_PROFILES[host])
    if host not in HOST_PROFILES:
        raise InputError(
            f"Unknown host profile '{host}'; available: {', '.join(sorted(HOST_PROFILES))}"
        )
    return CodonProfile(
        name=host,
        version=HOST_PROFILE_VERSIONS[host],
        source="Kazusa CUTG Homo sapiens, GenBank release 160.0, 2007-06-15",
        citation="https://zearth.kazusa.or.jp/codon/readme_codon.html",
        description="Species-level historical proxy; not a cell-line or expression model",
        frequencies=HOST_PROFILES[host],
    )


def codon_metrics(payload: Payload, profile: CodonProfile, rare_threshold: float):
    metrics = []
    for region in payload.coding_regions:
        if region.translation_table != profile.translation_table:
            raise InputError(f"Profile translation table does not match CDS '{region.name}'")
        coding = extract_cds(payload, region)
        table = CodonTable.unambiguous_dna_by_id[region.translation_table]
        codons = [coding[i:i + 3] for i in range(0, len(coding), 3)]
        sense = [c for c in codons if c in table.forward_table]
        weights = []
        for codon in sense:
            aa = table.forward_table[codon]
            family = [c for c, a in table.forward_table.items() if a == aa]
            # Met/Trp (single-codon families) contribute no codon adaptation information.
            if len(family) > 1:
                weights.append(profile.frequencies[codon] / max(
                    profile.frequencies[c] for c in family
                ))
        cai = None if not weights else (
            0.0 if 0.0 in weights else math.exp(sum(map(math.log, weights)) / len(weights))
        )
        positions = coding_region_positions(region)
        pairs = [a + b for a, b in zip(codons, codons[1:], strict=False)
                 if a in table.forward_table and b in table.forward_table]
        pair_values = [profile.pair_scores[p] for p in pairs if p in profile.pair_scores]
        metrics.append(CodonMetrics(
            region=region.name,
            cai=None if cai is None else round(cai, 6),
            gc3_percent=round(100 * sum(c[2] in "GC" for c in sense) / max(1, len(sense)), 6),
            distinct_codons=len(set(sense)),
            rare_codon_positions=tuple(
                positions[i * 3] + 1 for i, c in enumerate(codons)
                if c in table.forward_table and profile.frequencies[c] < rare_threshold
            ),
            mean_pair_score=round(sum(pair_values) / len(pair_values), 6) if pair_values else None,
            scored_pairs=len(pair_values), total_pairs=len(pairs),
        ))
    return tuple(metrics)


def codon_constraints_pass(metrics, settings: OptimizeSettings) -> bool:
    return all(
        settings.gc3_min <= m.gc3_percent <= settings.gc3_max
        and m.distinct_codons >= settings.min_distinct_codons
        and (settings.max_rare_codons is None
             or len(m.rare_codon_positions) <= settings.max_rare_codons)
        for m in metrics
    )


def resolve_protection(payload: Payload, settings: OptimizeSettings, vector=None):
    regions = list(settings.protected_regions)
    if vector:
        regions.extend(ProtectedRegion(start=c.start, end=c.end, name=c.name)
                       for c in vector.components if c.protected or c.kind == "itr")
    for feature in payload.features:
        labels = " ".join(v for values in feature.qualifiers.values() for v in values).lower()
        is_itr = "itr" in labels.split() or "inverted terminal repeat" in labels
        if feature.kind not in settings.protected_feature_types and not is_itr:
            continue
        if feature.location is None:
            raise InputError("Protected feature has no coordinates")
        try:
            loc = Location.fromstring(feature.location, len(payload.sequence))
            for part in loc.parts:
                regions.append(ProtectedRegion(
                    start=int(part.start) + 1, end=int(part.end), name=feature.kind
                ))
        except (ValueError, TypeError) as exc:
            raise InputError(f"Cannot resolve protected feature: {exc}") from exc
    if any(r.end > len(payload.sequence) for r in regions):
        raise InputError("Protected interval extends beyond the sequence")
    regions = tuple(sorted(set((r.start, r.end, r.name) for r in regions)))
    resolved = tuple(ProtectedRegion(start=s, end=e, name=n) for s, e, n in regions)
    return settings.model_copy(update={"protected_regions": resolved})


def protected_positions(settings: OptimizeSettings) -> frozenset[int]:
    return frozenset(p for r in settings.protected_regions for p in range(r.start - 1, r.end))


def vector_budget(payload: Payload, architecture: VectorArchitecture) -> dict:
    length = len(payload.sequence)
    if any(c.end > length for c in architecture.components):
        raise InputError("Vector component extends beyond the input sequence")
    itrs = sorted((c for c in architecture.components if c.kind == "itr"), key=lambda c: c.start)
    if architecture.representation == "itr_to_itr_unit":
        if len(itrs) != 2 or itrs[0].start != 1 or itrs[1].end != length:
            raise InputError("ITR-to-ITR input requires two terminal ITR components")
    elif itrs:
        raise InputError("Cassette excludes ITRs; do not also annotate ITR components")
    overhead = sum(architecture.omitted_itr_lengths or ())
    limit = architecture.limit or (4700 if architecture.genome_type == "ssAAV" else 2400)
    total = length + overhead
    components = []
    # Dynamic programming retains one deterministic minimum-swap plan for each size.
    # The planning space is capped to prevent untrusted catalogs exhausting memory.
    plans: dict[int, tuple] = {total: ()}
    for component in architecture.components:
        size = component.end - component.start + 1
        alternatives = sorted(component.alternatives, key=lambda a: (a.length, a.name))
        components.append({
            **component.model_dump(mode="json"), "length": size,
            "shorter_alternatives": [
                {**a.model_dump(), "saving": size - a.length} for a in alternatives
                if a.length < size and not component.protected and component.kind != "itr"
            ],
        })
        if component.protected or component.kind == "itr":
            continue
        next_plans = dict(plans)
        for current, choices in plans.items():
            for alternative in alternatives:
                if alternative.length >= size:
                    continue
                new_size = current - size + alternative.length
                choice = choices + ((component.name, alternative.name, alternative.length),)
                old = next_plans.get(new_size)
                if old is None or (len(choice), choice) < (len(old), old):
                    next_plans[new_size] = choice
        if len(next_plans) > 100000:
            raise InputError("Vector alternative catalog exceeds the 100000-size planning limit")
        plans = next_plans
    feasible = [(size, choices) for size, choices in plans.items() if size <= limit]
    chosen_size, choices = min(
        feasible, key=lambda p: (len(p[1]), -p[0], p[1])
    ) if feasible else min(plans.items())
    return {
        "genome_type": architecture.genome_type,
        "representation": architecture.representation,
        "input_length": length, "omitted_itr_overhead": overhead,
        "accounted_unit_length": total, "limit": limit,
        "within_limit": total <= limit, "headroom": limit - total,
        "unassigned_input_bases": length - sum(
            c.end - c.start + 1 for c in architecture.components
        ),
        "components": components,
        "plan": {"estimated_length": chosen_size, "fits": chosen_size <= limit,
                 "changes": [{"component": c, "alternative": a, "length": n}
                             for c, a, n in choices]},
        "limitations": [
            "Component lengths and identities are user-declared; ITR sequence integrity "
            "is not verified.",
            "Alternative plans are size estimates, not validated biological replacements "
            "or edited DNA.",
            "scAAV uses a conservative 2400-nt ITR-to-ITR unit budget, not a physical "
            "dimer model.",
            "Packaging limits are configurable screening thresholds, not packaging-yield "
            "predictions.",
        ],
        "source": "https://www.addgene.org/guides/aav/",
    }
