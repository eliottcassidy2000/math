"""Exact referee for THM-2054's seven-sector Fejer atom budget."""

from fractions import Fraction
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


ATOM_WEIGHT = sum(comb(6, a) * a for a in range(7))
require(ATOM_WEIGHT == 192, "six-sector atom weight")

MARGINS = {
    8: Fraction(1087, 5880),
    9: Fraction(129643, 980980),
    10: Fraction(5583, 35672),
    11: Fraction(311453, 1605240),
    12: Fraction(6295, 24696),
}

# Since 2^20+1<2^21 and log(2)<1, log(2^20+1)<21.
N20 = 2**20 + 1
require(N20 < 2**21, "H=2^20 logarithm bound")
UNIFORM_H20 = Fraction(192 * 11 * 43, N20)
require(UNIFORM_H20 == Fraction(90816, 1048577), "uniform H=2^20 error")
UNIFORM_RESIDUAL = min(MARGINS.values()) - UNIFORM_H20
require(
    UNIFORM_RESIDUAL == Fraction(46851988331, 1028633065460),
    "uniform H=2^20 residual",
)
require(UNIFORM_RESIDUAL > 0, "positive uniform residual")

# Since 2^19+1<2^20 and log(2)<1, log(2^19+1)<20.
N19 = 2**19 + 1
require(N19 < 2**20, "H=2^19 logarithm bound")
ROW_ERRORS = {
    k: Fraction(192 * (k - 1) * 41, N19) for k in range(8, 13)
}
ROW_RESIDUALS = {k: MARGINS[k] - ROW_ERRORS[k] for k in MARGINS}
require(all(value > 0 for value in ROW_RESIDUALS.values()), "rowwise residuals")
require(min(ROW_RESIDUALS, key=ROW_RESIDUALS.get) == 9, "tightest row is k=9")
require(
    ROW_RESIDUALS[9] == Fraction(2064067449, 171439007740),
    "exact k=9 residual",
)

# Boundary alias showing that |M|>H sum |b_i| must be strict.
chi_1 = (1, 0)
chi_2 = (0, 1)
lam = (1, 1)
n = (1, -1)
scalar = n[0] * sum(x * y for x, y in zip(chi_1, lam)) + n[1] * sum(
    x * y for x, y in zip(chi_2, lam)
)
vector = tuple(n[0] * chi_1[j] + n[1] * chi_2[j] for j in range(2))
require(scalar == 0 and vector != (0, 0), "strict cutoff alias")

print("THM-2054 RELATIVE FEJER ATOM-BUDGET REFEREE")
print(f"sum_A |A|={ATOM_WEIGHT}")
print(f"H=2^20 uniform_error<{UNIFORM_H20}")
print(f"H=2^20 residual>{UNIFORM_RESIDUAL}")
for k in range(8, 13):
    print(
        f"H=2^19 k={k} error<{ROW_ERRORS[k]} "
        f"residual>{ROW_RESIDUALS[k]}"
    )
print(f"tightest_rowwise_residual=k9:{ROW_RESIDUALS[9]}")
print(f"strict_cutoff_alias=scalar:{scalar},vector:{vector}")
print("TOURNAMENT ANALYSIS=not applicable: signed atom hyperedges are essential")
print("RESULT=PASS")
