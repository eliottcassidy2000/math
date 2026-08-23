#!/usr/bin/env python3
"""Independent arithmetic referee for the proposed THM-3743 reduction.

Literature inputs (BHS closed-zonotope equivalence and the flatness bound) are
audited in the accompanying report.  This file independently recomputes every
integer constant without importing the incoming companion.
"""

from fractions import Fraction
from math import comb, gcd, isqrt


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


d = 12
flatness_square = Fraction((d + 1) * (2 * d + 1), 6) * d**3
require(flatness_square == 93600, flatness_square)
relation_square = Fraction(7, 6) ** 2 * flatness_square
require(relation_square == 127400, relation_square)
relation_cap = isqrt(relation_square.numerator)
require(relation_cap == 356, relation_cap)
require(relation_cap**2 < relation_square < (relation_cap + 1) ** 2, relation_square)

# Independent totient-by-definition count of unordered reduced p<q,
# p+q<=356.  For each sum s>=3 there are phi(s)/2 such pairs.
totient_sum = 0
for total in range(3, relation_cap + 1):
    phi = sum(gcd(residue, total) == 1 for residue in range(1, total))
    require(phi % 2 == 0, (total, phi))
    totient_sum += phi // 2
require(totient_sum == 19314, totient_sum)

# Independent dynamic-program count of all nonzero vectors in Z^13 with
# l1 norm <=356.  One coordinate contributes 1+2x+2x^2+... .
coefficient_counts = [0] * (relation_cap + 1)
coefficient_counts[0] = 1
for _ in range(13):
    updated = [0] * (relation_cap + 1)
    for old_mass, old_count in enumerate(coefficient_counts):
        if not old_count:
            continue
        updated[old_mass] += old_count
        for magnitude in range(1, relation_cap - old_mass + 1):
            updated[old_mass + magnitude] += 2 * old_count
    coefficient_counts = updated
ambient_count = sum(coefficient_counts) - 1
require(ambient_count == 1978967793896659449022201064, ambient_count)

# Check the closed AP boundary hostile directly.
ap_speeds = tuple(range(1, 14))
ap_time = Fraction(1, 14)
ap_minimum = min(
    min((ap_time * speed) % 1, 1 - ((ap_time * speed) % 1))
    for speed in ap_speeds
)
require(ap_minimum == Fraction(1, 14), ap_minimum)
ap_relation = (1, -2, 1) + (0,) * 10
require(sum(a * n for a, n in zip(ap_relation, ap_speeds)) == 0, ap_relation)
require(sum(abs(a) for a in ap_relation) == 4, ap_relation)
require(Fraction(6, 7) * 4 == Fraction(24, 7), "AP width")

# Exhaust the short (3,4,5) box independently.
triple_minimum = None
triple_minimizers = []
for a in range(-8, 9):
    for b in range(-8, 9):
        for c in range(-8, 9):
            if (a, b, c) == (0, 0, 0) or 3 * a + 4 * b + 5 * c:
                continue
            norm = abs(a) + abs(b) + abs(c)
            if triple_minimum is None or norm < triple_minimum:
                triple_minimum = norm
                triple_minimizers = [(a, b, c)]
            elif norm == triple_minimum:
                triple_minimizers.append((a, b, c))
require(triple_minimum == 4, triple_minimum)
require(set(triple_minimizers) == {(1, -2, 1), (-1, 2, -1)}, triple_minimizers)

Q = 91**6
hadamard_square = relation_cap**2 * 3**11 * Q**22
hadamard_cap = isqrt(hadamard_square)
expected_hadamard_cap = int(
    "296721347184071259951513500572385227832063299530012809732642802018052529730741866571256016315918699080945038248181476620474893833141605"
)
require(hadamard_cap == expected_hadamard_cap, hadamard_cap)
require(hadamard_cap**2 <= hadamard_square < (hadamard_cap + 1) ** 2, "Hadamard floor")

# Cross-check the DP count against the closed support formula.
support_formula = sum(
    2**support * comb(13, support) * comb(relation_cap, support)
    for support in range(1, 14)
)
require(support_formula == ambient_count, (support_formula, ambient_count))

print("THM-3743 independent arithmetic referee")
print(f"flatness_square={flatness_square}")
print(f"relation_square={relation_square}")
print(f"relation_cap={relation_cap}")
print(f"reduced_pair_ratio_count={totient_sum}")
print(f"coordinate_support_times_unordered_ratio_count={comb(13,2)*totient_sum}")
print(f"ambient_l1_vector_count={ambient_count}")
print(f"triple_minimizers={tuple(triple_minimizers)}")
print(f"ap_boundary_minimum={ap_minimum}")
print(f"hadamard_cap={hadamard_cap}")
print("ALL CHECKS PASSED;ARTIFACT=THM-3743_INDEPENDENT_AUDIT")
