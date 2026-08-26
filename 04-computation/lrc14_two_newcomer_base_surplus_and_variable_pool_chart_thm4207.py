#!/usr/bin/env python3
"""Exact THM-4207 base-surplus, depth-four, and chart-number controls."""

from fractions import Fraction as F
from math import comb, lcm
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "04-computation"))

from lrc14_certificates import L_exact
from lrc14_two_scale_tail_continuation_thm2184 import safe_measure

P = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
A = {120, 126, 143}
O = set(P) - A
ALPHA = F(4, 63)
POOL_LATTICE = 18_241_159_416_480
L50 = lcm(POOL_LATTICE, 14 * 50)
L51 = lcm(POOL_LATTICE, 14 * 51)
L5051 = lcm(POOL_LATTICE, 14 * 50, 14 * 51)

# Complete E_4^A(50) list proved/audited in THM-4174, equation (19).
EDGES = (
    (88, 95, 176, 193),
    (88, 145, 176, 193),
    (88, 145, 193, 290),
    (145, 168, 193, 290),
)


def exact_mass(row):
    first = L_exact(row)
    second = safe_measure(row)
    assert first == second
    return first


def text(x):
    return f"{x.numerator}/{x.denominator}"


for edge in EDGES:
    assert set(edge) <= O
    base = tuple(v for v in P if v not in edge)
    m50 = exact_mass(base + (50,))
    m51 = exact_mass(base + (51,))
    m5051 = exact_mass(base + (50, 51))
    d50 = m50 - ALPHA
    d51 = m51 - ALPHA
    d5051 = m5051 - ALPHA
    # These are the exact threshold-comparison integers 63*N-4*L on each
    # row's natural pool/joint-wall denominator L.
    delta63 = (
        int(63 * L50 * d50),
        int(63 * L51 * d51),
        int(63 * L5051 * d5051),
    )
    assert d50 > 0 and d51 > 0 and d5051 < 0
    print(
        f"R={edge}; m50={text(m50)}; m51={text(m51)}; "
        f"m5051={text(m5051)}; margins=({text(d50)},{text(d51)},"
        f"{text(d5051)}); wall_denominators=({L50},{L51},{L5051}); "
        f"delta63={delta63}"
    )

# The first edge also freezes the exact information discarded by ordinary
# deck membership.  If U=G_(P\R), inclusion-exclusion inside U gives
# m_(50,51) >= m_50+m_51-u.  Relative to alpha, its retained margin is
# sigma_50+sigma_51-delta_U.
edge = EDGES[0]
base = tuple(v for v in P if v not in edge)
u = exact_mass(base)
m50 = exact_mass(base + (50,))
m51 = exact_mass(base + (51,))
m5051 = exact_mass(base + (50, 51))
bonferroni = m50 + m51 - u
bonferroni_margin = bonferroni - ALPHA
assert m5051 >= bonferroni
assert bonferroni_margin == F(-320_371, 24_674_650)
print(
    "BASE_SURPLUS R=(88,95,176,193); "
    f"u={text(u)}; pair_lower={text(bonferroni)}; "
    f"lower_margin={text(bonferroni_margin)}; "
    f"actual_pair={text(m5051)}"
)


def uncovered_k_sets(n, k, chart_count):
    """Number of k-sets containing every omitted chart label."""
    if chart_count > k:
        return 0
    return comb(n - chart_count, k - chart_count)


# For an n-label base and k-subset targets, charts indexed by omitted labels
# p cover B precisely when p is not in B.  A chart set O leaves exactly the
# k-sets containing O uncovered, hence the minimum chart number is k+1.
n = len(P)
k = 9
uncovered = tuple(uncovered_k_sets(n, k, d) for d in range(k + 2))
assert uncovered[k] == 1 and uncovered[k + 1] == 0
print(
    f"CHART_NUMBER n={n} k={k} minimum={k + 1} "
    f"uncovered_by_chart_count={uncovered}"
)

print("PASS: interval engines, base-surplus control, and chart count agree")
