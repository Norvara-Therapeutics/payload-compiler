"""Deterministic bounded exhaustive synonymous search with overlap constraints."""

from Bio.Data import CodonTable
from Bio.Seq import Seq

from payload_compiler.models import SearchResult
from payload_compiler.validation import coding_region_positions


def global_search(payload, profile, settings, feasible, seed):
    """Enumerate complete synonymous assignments, pruning fixed/new-CpG conflicts.

    A state budget counts attempted codon assignments (not elapsed time). Exhaustion
    proves the declared lexicographic optimum over all feasible assignments; a cutoff
    returns the best feasible incumbent and never claims global optimality.
    """
    from payload_compiler.optimize import _codon_in_sequence, _cpg_positions, _translate_codon

    original = payload.sequence
    original_cpg = set(_cpg_positions(original))
    fixed = {p: original[p] for r in settings.protected_regions
             for p in range(r.start - 1, r.end)}
    variables = []
    covered = set()
    for region in sorted(payload.coding_regions, key=lambda r: (r.start, r.end, r.name)):
        positions = coding_region_positions(region)
        table = CodonTable.unambiguous_dna_by_id[region.translation_table]
        for offset in range(0, len(positions), 3):
            triplet = positions[offset:offset + 3]
            covered.update(triplet)
            codon = _codon_in_sequence(original, triplet, region.strand)
            aa = _translate_codon(codon, table, first=offset == 0)
            # Start and stop identity are protected, including coding fragments.
            synonyms = [codon] if offset == 0 or aa == "*" else [
                c for c, amino in table.forward_table.items() if amino == aa
            ]
            domains = []
            ordered = sorted(synonyms, key=lambda c: (c != codon, -profile.frequencies[c], c))
            for replacement in ordered:
                bases = replacement if region.strand == 1 else str(Seq(replacement).complement())
                domains.append(tuple(zip(triplet, bases, strict=True)))
            if len(domains) == 1:
                fixed.update(domains[0])
            else:
                variables.append(domains)
    fixed.update({p: b for p, b in enumerate(original) if p not in covered})
    owner = {}
    for depth, domains in enumerate(variables):
        for position, _ in domains[0]:
            owner.setdefault(position, depth)
    # Fixed bases are known from the root, including single-choice overlapping CDSs.
    def known(p, depth):
        return p in fixed or owner.get(p, len(variables)) < depth

    def objective(sequence):
        preference = 0.0
        for region in payload.coding_regions:
            positions = coding_region_positions(region)
            for offset in range(0, len(positions), 3):
                preference += profile.frequencies[_codon_in_sequence(
                    sequence, positions[offset:offset + 3], region.strand
                )]
        changed = sum(a != b for a, b in zip(original, sequence, strict=True))
        return (len(_cpg_positions(sequence)), changed, -round(preference, 9), sequence)

    best = seed
    best_key = objective(seed) if seed is not None else None
    explored = 0
    exhausted = True
    # Frames hold the next domain index; siblings are generated lazily. This avoids
    # recursion limits and bounds memory to O(CDS length * sequence length).
    stack = [(0, original, 0)]
    while stack:
        depth, sequence, choice_index = stack.pop()
        if depth == len(variables):
            key = objective(sequence)
            if (best_key is None or key < best_key) and feasible(sequence):
                best, best_key = sequence, key
            continue
        if choice_index >= len(variables[depth]):
            continue
        if explored >= settings.search_budget:
            exhausted = False
            break
        stack.append((depth, sequence, choice_index + 1))
        explored += 1
        choice = variables[depth][choice_index]
        if any((p in fixed and fixed[p] != b)
               or (known(p, depth) and sequence[p] != b) for p, b in choice):
            continue
        chars = list(sequence)
        for p, base in choice:
            chars[p] = base
        updated = "".join(chars)
        determined_cpgs = {
            p for p in _cpg_positions(updated) if known(p - 1, depth + 1) and known(p, depth + 1)
        }
        # Under "strict_subset" a settled CpG outside the original set can never be
        # undone by deeper assignments, so the branch is pruned. "count_only" allows
        # relocation, so only the incumbent count bound below applies.
        if settings.cpg_policy == "strict_subset" and not determined_cpgs.issubset(original_cpg):
            continue
        if best_key is not None and len(determined_cpgs) > best_key[0]:
            continue
        stack.append((depth + 1, updated, 0))
    status = ("proven_optimal" if best is not None else "infeasible") if exhausted else (
        "best_found" if best is not None else "no_feasible_found"
    )
    return best, SearchResult(
        solver="global", status=status, explored_states=explored, budget=settings.search_budget,
        proof_scope=(
            "All synonymous assignments under declared translation, fixed start/stop, "
            "protected-base, "
            "no-new-CpG, motif, splice, codon and enabled screening/folding constraints; "
            "linear stored sequence; no claim of biological optimality"
        ),
    )
