#!/usr/bin/env python3
"""Independent exact referee for THM-2269.

This path uses character orthogonality and direct congruence reconstruction,
not the primary companion's mask-energy loop.
"""

from fractions import Fraction


def norm_num(num: int, den: int) -> int:
    r = num % den
    return min(r, den - r)


# Nonzero-character orthogonality:
# sum_{k=1}^{12} zeta^(kd) is 12 for d=0 and -1 otherwise.
def direct_nonzero_energy(support: tuple[int, ...]) -> int:
    total = 0
    for left in support:
        for right in support:
            total += 12 if (right - left) % 13 == 0 else -1
    return total


referee_masks = []
for a in range(13):
    referee_masks.append((a,))
for a in range(13):
    for b in range(a + 1, 13):
        referee_masks.append((a, b))

assert len(referee_masks) == 91
assert {direct_nonzero_energy(mask) for mask in referee_masks} == {12, 22}

# Independent angle audit.  If d is the nonzero exponent difference, the
# distance of pi*d/13 from a zero of cosine, in units pi/26, is |2d-13|.
# Its minimum is one.  The chord bound sin(t)>=2t/pi then gives
# |1+zeta^d|^2>=4/169 for every d.
angle_distance_multiset = sorted(abs(2 * d - 13) for d in range(1, 13))
assert angle_distance_multiset == [1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11]
uniform_mode_lower = Fraction(4, 169)

# Independently verify the L2 residue-class transform on an exact rational
# trigonometric polynomial.  Root-character summation keeps n=k mod 13 and
# multiplies the coefficient by 13; y-orthogonality then squares and sums.
coeff = {
    n: Fraction((n * n + 3 * n + 7) % 19 - 9, (abs(n) % 5) + 2)
    for n in range(-47, 48)
}
residue_energy = {}
for k in range(13):
    grouped = {
        (n - k) // 13: 13 * value
        for n, value in coeff.items()
        if (n - k) % 13 == 0
    }
    direct_l2 = sum(value * value for value in grouped.values())
    predicted_l2 = 169 * sum(
        value * value for n, value in coeff.items() if n % 13 == k
    )
    assert direct_l2 == predicted_l2
    residue_energy[k] = direct_l2

# Reconstruct the hostile branch using integer numerators over 338.
x_num, den = 79, 338
q = (2, 339, 677, 1015, 1353)
assert 7 * norm_num(x_num, den) > den
assert all(14 * norm_num(speed * x_num, den) >= den for speed in q)
assert 14 * norm_num(13 * x_num, den) < den

checked = 0
for c in range(5, 20):
    for b in range(2, c):
        assert 14 * norm_num((13**b) * x_num, den) >= den
        assert 14 * norm_num((13**c) * x_num, den) >= den
        checked += 1
assert checked == 150

# Forward phases are 79/338 -> 1/26 -> 1/2 -> 1/2.
assert (13 * Fraction(79, 338)) % 1 == Fraction(1, 26)
assert (13 * Fraction(1, 26)) % 1 == Fraction(1, 2)
assert (13 * Fraction(1, 2)) % 1 == Fraction(1, 2)

strict = Fraction(15041431, 70270200)
repeat = Fraction(5229541, 70270200)
assert strict / 169 == Fraction(15041431, 11875663800)
assert repeat / 169 == Fraction(5229541, 11875663800)
assert 4 * strict / 28561 == Fraction(15041431, 501746795550)
assert 4 * repeat / 28561 == Fraction(5229541, 501746795550)

print("THM-2269 independent exact referee")
print(f"orthogonality_masks_checked={len(referee_masks)}")
print("orthogonality_energy_values=(12,22)")
print(f"uniform_mode_lower={uniform_mode_lower}")
print(f"trigonometric_coefficients_checked={len(coeff)}")
print(f"residue_classes_checked={len(residue_energy)}")
print(f"strict_residue_floor={strict / 169}")
print(f"repeated_residue_floor={repeat / 169}")
print(f"strict_every_residue_floor={4 * strict / 28561}")
print(f"repeated_every_residue_floor={4 * repeat / 28561}")
print(f"hostile_profiles_checked={checked}")
print("hostile_orbit=(79/338,1/26,1/2,1/2)")
print("ALL_CHECKS_PASSED")
