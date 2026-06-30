"""Amino acid alphabets, scales, and residue groupings shared across feature modules.

Charge convention follows localCIDER: K and R are positive, D and E are negative, and
histidine is treated as neutral. The disorder-promoting set is the positive end of the
TOP-IDP scale (Campen et al. 2008), the residues enriched in disordered sequence.
"""

from __future__ import annotations

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = frozenset(AA)

POSITIVE = frozenset("KR")
NEGATIVE = frozenset("DE")
CHARGED = POSITIVE | NEGATIVE
AROMATIC = frozenset("FWY")
POLAR = frozenset("STNQCYH")  # polar uncharged plus histidine
PROLINE = frozenset("P")
GLYCINE = frozenset("G")
DISORDER_PROMOTING = frozenset("ARGQSEKP")

# Kyte-Doolittle hydropathy. Normalized to [0, 1] as (kd + 4.5) / 9.0, the localCIDER
# convention, so 1.0 is the most hydrophobic residue and 0.0 the most hydrophilic.
KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}
KD_NORM = {aa: (kd + 4.5) / 9.0 for aa, kd in KYTE_DOOLITTLE.items()}
