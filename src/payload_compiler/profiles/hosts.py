"""Relative codon frequencies transcribed from the Homo sapiens CUTG table.

Source: Kazusa Codon Usage Database, Homo sapiens [gbpri], GenBank Flat File
Release 160.0 (2007-06-15), distributed as RNA codons in relative-frequency
format. Values are used only to break ties between equally CpG-reducing,
synonymous candidates; this is a general human proxy, not an HEK293-specific
expression model.
"""

HUMAN_CODON_FREQUENCIES = {
    "TTT": 0.46,
    "TTC": 0.54,
    "TTA": 0.08,
    "TTG": 0.13,
    "CTT": 0.13,
    "CTC": 0.20,
    "CTA": 0.07,
    "CTG": 0.40,
    "ATT": 0.36,
    "ATC": 0.47,
    "ATA": 0.17,
    "ATG": 1.00,
    "GTT": 0.18,
    "GTC": 0.24,
    "GTA": 0.12,
    "GTG": 0.46,
    "TCT": 0.19,
    "TCC": 0.22,
    "TCA": 0.15,
    "TCG": 0.05,
    "AGT": 0.15,
    "AGC": 0.24,
    "CCT": 0.29,
    "CCC": 0.32,
    "CCA": 0.28,
    "CCG": 0.11,
    "ACT": 0.25,
    "ACC": 0.36,
    "ACA": 0.28,
    "ACG": 0.11,
    "GCT": 0.27,
    "GCC": 0.40,
    "GCA": 0.23,
    "GCG": 0.11,
    "TAT": 0.44,
    "TAC": 0.56,
    "TAA": 0.30,
    "TAG": 0.24,
    "TGA": 0.47,
    "CAT": 0.42,
    "CAC": 0.58,
    "CAA": 0.27,
    "CAG": 0.73,
    "AAT": 0.47,
    "AAC": 0.53,
    "AAA": 0.43,
    "AAG": 0.57,
    "GAT": 0.46,
    "GAC": 0.54,
    "GAA": 0.42,
    "GAG": 0.58,
    "TGT": 0.46,
    "TGC": 0.54,
    "TGG": 1.00,
    "CGT": 0.08,
    "CGC": 0.18,
    "CGA": 0.11,
    "CGG": 0.20,
    "AGA": 0.21,
    "AGG": 0.21,
    "GGT": 0.16,
    "GGC": 0.34,
    "GGA": 0.25,
    "GGG": 0.25,
}

HOST_PROFILES = {
    "human": HUMAN_CODON_FREQUENCIES,
}
HOST_PROFILE_VERSIONS = {
    "human": "Homo sapiens CUTG, GenBank Release 160.0, 2007-06-15",
}

#: Profile names withdrawn because the bundled data did not support the claim the
#: name made. Requesting one is an error that explains the substitute, rather than
#: silently returning different data than the name implies.
RETIRED_PROFILES = {
    "human-hek293": (
        "'human-hek293' was an alias for the general Homo sapiens CUTG table and was "
        "never HEK293-derived. It has been withdrawn rather than left to imply "
        "cell-line measurements. Use --host human for the species-level proxy, or "
        "supply a measured cell-line table with --profile FILE."
    ),
}
