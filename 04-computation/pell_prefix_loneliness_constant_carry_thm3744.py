#!/usr/bin/env python3
"""Exact arithmetic companion for THM-3744.

For Pell numbers P_0=0, P_1=1, P_{n+2}=2P_{n+1}+P_n, this checks the
all-length certificate

    M({P_1,...,P_N}) = A_N/(P_N+1),
    A_N=(P_N-P_{N-1}+1)/2.

The global upper bound is checked through the exact interval-collapse
inequalities used in the proof.  A logically independent lower-envelope
candidate scan recomputes the maximum for N<=9.  Further checks cover the
constant-carry residue formula, odd square-triangular factorizations, profile
formulas, and hostile perturbations.  No floating point is used.
"""

from __future__ import annotations

from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def distance(numerator: int, denominator: int) -> Fraction:
    residue = numerator % denominator
    return Fraction(min(residue, denominator - residue), denominator)


def pell_numbers(limit: int) -> list[int]:
    pell = [0, 1]
    while len(pell) <= limit:
        pell.append(2 * pell[-1] + pell[-2])
    return pell


def direct_lower_envelope(speeds: tuple[int, ...]) -> tuple[Fraction, tuple[Fraction, ...]]:
    """Exact independent maximum of min_p ||p t|| on [0,1/2].

    Every distance function is affine between its half-integral cusps.  A
    maximum of their lower envelope is therefore an endpoint, a cusp, or an
    equality point of two affine pieces.  Such points have denominator 2p,
    p+q, or |p-q|.  No unenumerated flat maximizing interval exists because
    every individual affine branch has nonzero slope.  Enumerating all
    corresponding rationals is exhaustive.
    """

    denominators = {2 * speed for speed in speeds}
    for first_index, first in enumerate(speeds):
        for second in speeds[first_index + 1 :]:
            denominators.add(first + second)
            denominators.add(abs(first - second))
    denominators.discard(0)

    candidates = {Fraction(0), Fraction(1, 2)}
    for denominator in denominators:
        candidates.update(
            Fraction(numerator, denominator)
            for numerator in range(denominator // 2 + 1)
        )

    best = Fraction(-1)
    owners: list[Fraction] = []
    for time in candidates:
        margin = min(
            distance(speed * time.numerator, time.denominator) for speed in speeds
        )
        if margin > best:
            best = margin
            owners = [time]
        elif margin == best:
            owners.append(time)
    return best, tuple(sorted(set(owners)))


MAX_N = 200
pell = pell_numbers(MAX_N + 2)
auxiliary = [0] * (MAX_N + 3)
checks = 0

print("=== Pell and affine companion recurrences ===")
for index in range(1, MAX_N + 2):
    numerator = pell[index] - pell[index - 1] + 1
    require(numerator % 2 == 0, "A_n is integral")
    auxiliary[index] = numerator // 2
    checks += 1

for index in range(1, MAX_N):
    require(
        auxiliary[index + 2]
        == 2 * auxiliary[index + 1] + auxiliary[index] - 1,
        "affine A recurrence",
    )
    require(
        pell[index + 1] * pell[index - 1] - pell[index] ** 2 == (-1) ** index,
        "Pell Cassini identity",
    )
    checks += 2
print("P_n and A_n checked through n=201; A_(n+2)=2A_(n+1)+A_n-1.")

print("\n=== Monotone candidate phases ===")
for index in range(1, MAX_N + 1):
    left = auxiliary[index] * (pell[index + 1] + 1)
    right = auxiliary[index + 1] * (pell[index] + 1)
    difference = left - right
    predicted = (pell[index] - pell[index - 1] + (-1) ** (index + 1)) // 2
    require(difference == predicted, "successive U_n determinant")
    require(difference >= 0, "U_n is nonincreasing")
    require((difference == 0) == (index == 2), "only U_2=U_3 ties")
    checks += 3
print("U_n=A_n/(P_n+1) decreases, with the sole adjacent tie U_2=U_3=1/3.")

print("\n=== All-length safe witness and global interval collapse ===")
for length in range(1, MAX_N + 1):
    denominator = pell[length] + 1
    delta = Fraction(auxiliary[length], denominator)

    for index in range(1, length + 1):
        value = Fraction(pell[index] * auxiliary[length], denominator)
        lower = auxiliary[index] - 1 + delta
        upper = auxiliary[index] - delta
        require(lower <= value <= upper, "candidate lies in every closed safe band")
        require(value.numerator // value.denominator == auxiliary[index] - 1, "exact floor")
        require(distance(pell[index] * auxiliary[length], denominator) >= delta, "safe margin")

        residue = pell[index] * auxiliary[length] - (auxiliary[index] - 1) * denominator
        residue_twice = (
            denominator
            + pell[index - 1]
            + (-1) ** index * pell[length - index]
        )
        require(2 * residue == residue_twice, "closed residue formula")
        profile_distance = Fraction(
            denominator
            - abs(pell[index - 1] + (-1) ** index * pell[length - index]),
            2 * denominator,
        )
        require(
            profile_distance
            == distance(pell[index] * auxiliary[length], denominator),
            "closed distance formula",
        )
        checks += 5

    for index in range(1, length - 1):
        first = Fraction(
            pell[index] * auxiliary[length], denominator
        ) - (auxiliary[index] - 1)
        second = Fraction(
            pell[index + 1] * auxiliary[length], denominator
        ) - (auxiliary[index + 1] - 1)
        third = Fraction(
            pell[index + 2] * auxiliary[length], denominator
        ) - (auxiliary[index + 2] - 1)
        require(third == 2 * second + first - 1, "constant carry equals one")
        checks += 1

    if length >= 2:
        require(delta > Fraction(auxiliary[length] - 1, pell[length] - 1), "terminal lower flank")
        require(Fraction(1, 2) * (1 - delta) == Fraction(auxiliary[2] - delta, pell[2]), "P_2 base bound")
        for index in range(2, length):
            require(
                pell[index + 1] * delta > auxiliary[index + 1] - 1 + delta,
                "inductive interval enters above the previous safe band",
            )
            require(
                Fraction(pell[index + 1], pell[index])
                * (auxiliary[index] - delta)
                < auxiliary[index + 1] + delta,
                "inductive interval stays below the next forbidden-band exit",
            )
            checks += 2
        require(
            Fraction(auxiliary[length] - delta, pell[length]) == delta,
            "terminal upper bound collapses to delta",
        )
        checks += 3

print("For every N<=200 the phase delta_N is safe and the proof interval collapses to t=delta_N.")
print("The same displayed inequalities are symbolic in N; the computation is a hostile replay, not induction by sampling.")

print("\n=== Independent exact lower-envelope scan ===")
for length in range(1, 10):
    expected = Fraction(auxiliary[length], pell[length] + 1)
    measured, owners = direct_lower_envelope(tuple(pell[1 : length + 1]))
    require(measured == expected, "direct lower envelope matches formula")
    require(owners == (expected,), "unique maximizer on [0,1/2]")
    checks += 2
    print(f"N={length:2d}  M={measured!s:>7}  unique_half_circle_time={owners[0]}")

print("\n=== Odd square-triangular factorization ===")
print(" k   N=4k-1       M       N=4k+1       M       (a,s,x,b)")
for depth in range(1, 9):
    a = pell[2 * depth - 1]
    s = pell[2 * depth]
    x = s + a
    b = pell[2 * depth + 1]
    left_length = 4 * depth - 1
    right_length = 4 * depth + 1
    left_margin = Fraction(auxiliary[left_length], pell[left_length] + 1)
    right_margin = Fraction(auxiliary[right_length], pell[right_length] + 1)

    require(x * x - 2 * s * s == 1, "square-triangular Pell equation")
    require(a * b == s * s + 1, "cannonball neighbor product")
    require(auxiliary[left_length] == 2 * a * a, "left A is twice a square")
    require(pell[left_length] + 1 == 2 * a * x, "left denominator factorization")
    require(auxiliary[right_length] == x * x, "right A is a square")
    require(pell[right_length] + 1 == 2 * x * b, "right denominator factorization")
    require(left_margin == Fraction(a, x), "left reduced margin")
    require(right_margin == Fraction(x, 2 * b), "right reduced margin")
    checks += 8
    print(
        f"{depth:2d}  {left_length:3d} {str(left_margin):>10}"
        f"  {right_length:3d} {str(right_margin):>10}"
        f"  ({a},{s},{x},{b})"
    )

print("\n=== Closed symmetric distance profiles ===")
for depth in range(1, 9):
    length = 4 * depth + 1
    center = 2 * depth + 1
    denominator = pell[length] + 1
    delta_numerator = auxiliary[length]
    for index in range(1, length + 1):
        measured = distance(pell[index] * delta_numerator, denominator)
        predicted = Fraction(
            pell[center] - pell[abs(center - index)], 2 * pell[center]
        )
        require(measured == predicted, "4k+1 symmetric profile")
        checks += 1

    length = 4 * depth - 1
    center = 2 * depth
    denominator = pell[length] + 1
    delta_numerator = auxiliary[length]
    x = pell[2 * depth] + pell[2 * depth - 1]
    for index in range(1, length + 1):
        profile_index = max(abs(index - center) - 1, 0)
        companion = pell[profile_index + 1] + pell[profile_index]
        measured = distance(pell[index] * delta_numerator, denominator)
        predicted = Fraction(x - companion, 2 * x)
        require(measured == predicted, "4k-1 three-center profile")
        checks += 1
print("Both odd-index profile formulas checked for k=1,...,8 at every runner.")

print("\n=== The 13-speed square-triangular row and hostiles ===")
length = 13
delta = Fraction(auxiliary[length], pell[length] + 1)
profile_numerators = [
    distance(pell[index] * delta.numerator, delta.denominator) * delta.denominator
    for index in range(1, length + 1)
]
require(delta == Fraction(99, 338), "Pell prefix N=13 margin")
require(
    profile_numerators == [99, 140, 157, 164, 167, 168, 169, 168, 167, 164, 157, 140, 99],
    "N=13 symmetric numerator profile",
)
modified_distance = distance(170 * delta.numerator, delta.denominator)
require(modified_distance == Fraction(35, 169) < delta, "single-speed perturbation hostile")
require(distance(0, 1) == 0, "including the zero speed kills loneliness")
checks += 4
print("S_13=(P_1,...,P_13) has M=99/338 and profile numerators over 338:")
print(" ", ",".join(str(value) for value in profile_numerators))
print("Replacing the central speed P_7=169 by 170 makes the proposed phase margin 35/169.")
print("Including P_0=0 forces the maximum to zero; neither hostile inherits the theorem.")

print("\n=== Scope ===")
print("The Pell word is the deterministic silver continued fraction with digit 2, not Khinchin-typical data.")
print("At N=13 the margin 99/338 is far above 1/14, so this is a structured safe family, not an LRC(14) frontier closure.")
print("No planar-Jacobian carrier or variable polynomial Cohn factor is produced.")
print(f"CHECKS={checks}")
