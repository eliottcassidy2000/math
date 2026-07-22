"""Exact referee for THM-2053's transverse residue deck and AP conductor."""

from fractions import Fraction
from math import floor


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def abs_residue(x: int, modulus: int) -> int:
    residue = x % modulus
    return min(residue, modulus - residue)


def deck(template: tuple[int, ...], modulus: int) -> Fraction:
    numerator = max(
        min(abs_residue(ell * speed, modulus) for speed in template)
        for ell in range(modulus)
    )
    return Fraction(numerator, modulus)


ap_checks = 0
for k in range(1, 13):
    template = tuple(range(1, k + 1))
    for modulus in range(1, 301):
        expected = Fraction(floor(modulus / (k + 1)), modulus)
        require(deck(template, modulus) == expected, f"AP deck k={k}, N={modulus}")
        ap_checks += 1

ap12 = (1, -1, *range(2, 13))
bad_ap12 = [modulus for modulus in range(1, 301) if deck(ap12, modulus) < Fraction(1, 14)]
require(max(bad_ap12) == 155, "AP12 conductor")
require(deck(ap12, 13) == Fraction(1, 13), "D_13")
require(deck(ap12, 15) == Fraction(1, 15), "D_15 nonmonotonicity")

templates = [
    ap12,
    (1, -1, 3, 5, 8),
    (2, -2, 3, 7, 11),
    (1, -1, 4, 6, 9, 13),
]
divisor_checks = 0
for template in templates:
    for modulus in range(1, 121):
        for divisor in range(1, modulus + 1):
            if modulus % divisor == 0:
                require(
                    deck(template, modulus) >= deck(template, divisor),
                    f"divisor monotonicity template={template}, d={divisor}, N={modulus}",
                )
                divisor_checks += 1

small_divisor_checks = 0
for template in templates:
    for divisor in range(1, 15):
        if all(speed % divisor != 0 for speed in template):
            require(deck(template, divisor) >= Fraction(1, divisor), "small-divisor exit")
            small_divisor_checks += 1

# Fixed chirotope, unbounded transverse magnitude.
for K in (12, 20, 100):
    columns = [(1, r) for r in range(12)] + [(1, K)]
    positive_row = [x + y for x, y in columns]
    require(len(set(positive_row)) == 13 and min(positive_row) > 0, f"positive row K={K}")
    normalized_coordinate = [Fraction(y, x + y) for x, y in columns]
    require(
        normalized_coordinate == sorted(normalized_coordinate)
        and normalized_coordinate[0] < normalized_coordinate[1],
        f"adjacent normalized pair K={K}",
    )
    signs = []
    for i in range(len(columns)):
        for j in range(i + 1, len(columns)):
            determinant = columns[i][0] * columns[j][1] - columns[i][1] * columns[j][0]
            signs.append((determinant > 0) - (determinant < 0))
    require(all(sign == 1 for sign in signs), f"fixed chirotope K={K}")
    transverse = [2 * r - 1 for r in range(12)] + [2 * K - 1]
    require(max(abs(value) for value in transverse) == 2 * K - 1, f"unbounded R K={K}")

print("THM-2053 TRANSVERSE RESIDUE-DECK REFEREE")
print(f"AP_formula_checks={ap_checks}")
print(f"divisor_monotonicity_checks={divisor_checks}")
print(f"small_divisor_exit_checks={small_divisor_checks}")
print(f"AP12_D13={deck(ap12, 13)} AP12_D15={deck(ap12, 15)}")
print(f"AP12_largest_bad_modulus={max(bad_ap12)} conductor={max(bad_ap12)+1}")
print("fixed_chirotope_R_values=23,39,199")
print("TOURNAMENT ANALYSIS=not applicable: determinant magnitudes are load-bearing")
print("RESULT=PASS")
