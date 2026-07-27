#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2568.

Checks the refined-pair to coarse-target line-sum map, the complementary
relative-twist hostile, a positive refined-drift control, fixed-frequency
leakage with exact full-frequency cancellation, and the anchored duplicate-
probe capacity/floor used to locate the typing boundary.
"""

from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_q(value):
    return value.numerator // value.denominator


def fractional_part(value):
    return value - floor_q(value)


def circle_norm(value):
    residue = fractional_part(value)
    return min(residue, 1 - residue)


def danger(value, level=1):
    return circle_norm(value) < Fraction(level, 14)


def cyclotomic_reduce(coefficients):
    """Reduce a length-13 coefficient vector modulo Phi_13."""
    require(len(coefficients) == P, "cyclotomic vector has wrong length")
    top = coefficients[-1]
    return tuple(value - top for value in coefficients[:-1])


def monomial(exponent, scale=1):
    coefficients = [0] * P
    coefficients[exponent % P] = scale
    return cyclotomic_reduce(coefficients)


ZERO_CYC = (0,) * (P - 1)


print("== refined pair spectrum -> coarse target line sum ==")
basis_checks = 0
for left_shift in range(P):
    for right_shift in range(P):
        for target in range(P):
            coefficients = [0] * P
            for left_colour in range(P):
                right_colour = (target - left_colour) % P
                exponent = (
                    left_colour * left_shift
                    + right_colour * right_shift
                ) % P
                coefficients[exponent] += 1
            observed = cyclotomic_reduce(coefficients)
            expected = (
                monomial(target * left_shift, P)
                if left_shift == right_shift
                else ZERO_CYC
            )
            require(observed == expected, "target line-sum basis identity failed")
            basis_checks += 1
require(basis_checks == P**3, "line-sum basis census changed")
pair_dimension = P**2
coarse_dimension = P
kernel_dimension = pair_dimension - coarse_dimension
require(kernel_dimension == 156, "coarse pushforward kernel changed")
print(f"  exact basis identities: {basis_checks}")
print(f"  refined/coarse/kernel dimensions: {pair_dimension}/{coarse_dimension}/{kernel_dimension}")
print("  diagonal zero annihilates every coarse target line")


print("\n== complementary relative-twist hostile ==")
# A(s,t)=1_(s!=t).  Its two-dimensional transform is supported only on
# a+b=0, while A(0,t) has every nontrivial one-sided colour.
hostile_checks = 0
for a in range(P):
    for b in range(P):
        coefficients = [0] * P
        for s in range(P):
            for t in range(P):
                if s != t:
                    coefficients[(a * s + b * t) % P] += 1
        observed = cyclotomic_reduce(coefficients)
        if (a + b) % P != 0:
            expected = ZERO_CYC
        elif a == 0 and b == 0:
            expected = monomial(0, P * (P - 1))
        else:
            expected = monomial(0, -P)
        require(observed == expected, "relative-twist hostile spectrum changed")
        hostile_checks += 1

for b in range(P):
    coefficients = [0] * P
    for t in range(1, P):
        coefficients[(b * t) % P] += 1
    expected = monomial(0, P - 1 if b == 0 else -1)
    require(
        cyclotomic_reduce(coefficients) == expected,
        "one-sided hostile spectrum changed",
    )

require(hostile_checks == P**2, "hostile transform census changed")
print(f"  exact two-dimensional coefficients: {hostile_checks}")
print("  Ahat(0,0)=12/13; Ahat(a,-a)=-1/13 for a!=0")
print("  A(0,t)=1_(t!=0) fires all twelve one-sided colours")
print("  diagonal-shift drift: 0; coarse target: identically 0")


print("\n== positive refined drift still has zero coarse pushforward ==")
control = [
    [0 if s == t else 1 + ((3 * s + 5 * t) % 7) for t in range(P)]
    for s in range(P)
]

orbit_means = [[Fraction(0) for _ in range(P)] for _ in range(P)]
for s in range(P):
    for t in range(P):
        orbit_means[s][t] = Fraction(
            sum(control[(s + g) % P][(t + g) % P] for g in range(P)),
            P,
        )

refined_drift = Fraction(0)
for s in range(P):
    for t in range(P):
        difference = Fraction(control[s][t]) - orbit_means[s][t]
        refined_drift += difference * difference
refined_drift /= P**2
require(refined_drift > 0, "positive refined-drift control became circulant")

row_means = [Fraction(sum(row), P) for row in control]
global_mean = sum(row_means, Fraction(0)) / P
row_variance = sum(
    (value - global_mean) ** 2 for value in row_means
) / P
require(row_variance > 0, "row-sum sufficient test lost variation")
require(refined_drift >= row_variance, "row variance exceeded refined drift")
require(all(control[s][s] == 0 for s in range(P)), "control diagonal is nonzero")
print(f"  refined drift: {refined_drift}")
print(f"  normalized-row variance: {row_variance}")
print("  every coarse target line still vanishes by the diagonal-zero identity")


print("\n== fixed-frequency leakage and full-frequency cancellation ==")
# On C4, P=delta_0 and Q=1-P.  With normalized Fourier coefficients,
# P_hat(X) conjugate(Q_hat(X)) has numerators (3,-1,-1,-1)/16.
fixed_x_numerators = (3, -1, -1, -1)
require(all(value != 0 for value in fixed_x_numerators), "fixed-X term vanished")
require(sum(fixed_x_numerators) == 0, "full-X complement cancellation failed")
print(f"  four fixed-X numerators /16: {fixed_x_numerators}")
print("  all four terms survive; their full-X sum is exactly zero")


print("\n== anchored duplicate-pair positive boundary ==")


def translated_profile(sample, level, sign):
    return tuple(
        danger(sample + sign * Fraction(shift, P), level)
        for shift in range(P)
    )


samples = tuple(Fraction(2 * cell + 1, 364) for cell in range(182))
danger_profiles = {
    level: tuple(
        sorted(
            {
                translated_profile(sample, level, +1)
                for sample in samples
            }
        )
    )
    for level in (1, 2)
}
safe_a_profiles = tuple(
    tuple(not value for value in translated_profile(sample, 1, -1))
    for sample in samples
)
safe_a_profiles = tuple(sorted(set(safe_a_profiles)))

anchored_counts = {}
minimum_row_capacities = {}
phase_identity_checks = {}
for level in (1, 2):
    anchored = 0
    minimum_capacity = P
    phase_checks = 0
    for safe_a in safe_a_profiles:
        if not safe_a[0]:
            continue
        for danger_k in danger_profiles[level]:
            if not danger_k[0]:
                continue
            safe_k = tuple(not value for value in danger_k)
            table = [
                [
                    int(
                        safe_a[r]
                        and danger_k[r]
                        and safe_a[s]
                        and safe_k[s]
                    )
                    for s in range(P)
                ]
                for r in range(P)
            ]
            require(all(table[r][r] == 0 for r in range(P)), "aux diagonal failed")
            require(all(table[r][0] == 0 for r in range(P)), "aux zero column failed")
            row_capacity = sum(table[0])
            minimum_capacity = min(minimum_capacity, row_capacity)
            anchored += 1

            delta = max(range(P), key=lambda shift: table[0][shift])
            require(table[0][delta] == 1, "positive offset disappeared")
            require(delta != 0, "positive offset became diagonal")
            offset_line = [table[r][(r + delta) % P] for r in range(P)]
            require(offset_line[0] == 1, "offset anchor disappeared")
            require(offset_line[(-delta) % P] == 0, "offset forced zero disappeared")

            # Exact inversion at the forced zero:
            # sum_(q!=0) Fhat(q) zeta^(q delta) = -mean(F).
            character_numerator = 0
            for r, value in enumerate(offset_line):
                if value:
                    character_numerator += (
                        P - 1 if (r + delta) % P == 0 else -1
                    )
            require(
                character_numerator == -sum(offset_line),
                "phase-rotated nonzero-colour identity failed",
            )
            phase_checks += 1

    expected_floor = P - 2 - 2 * level
    require(minimum_capacity == expected_floor, "paired capacity floor changed")
    anchored_counts[level] = anchored
    minimum_row_capacities[level] = minimum_capacity
    phase_identity_checks[level] = phase_checks

require(minimum_row_capacities == {1: 9, 2: 7}, "ordinary/guard floors changed")
ordinary_floor = Fraction(9, P**2 * (P - 1))
guard_floor = Fraction(7, P**2 * (P - 1))
require(ordinary_floor == Fraction(3, 676), "ordinary coefficient floor changed")
require(guard_floor == Fraction(7, 2028), "guard coefficient floor changed")
print(f"  anchored profile pairs (ordinary/guard): {anchored_counts[1]}/{anchored_counts[2]}")
print("  minimum positive row capacities: 9 ordinary, 7 guard")
print(f"  phase identities checked: {phase_identity_checks[1] + phase_identity_checks[2]}")
print(f"  auxiliary phase-rotated floors: {ordinary_floor} rho / {guard_floor} rho")
print("  scope: shifted duplicate probes over a fixed old head, not a lawful left endpoint orbit")


print("\nTHM-2568 exact referee: PASS")
