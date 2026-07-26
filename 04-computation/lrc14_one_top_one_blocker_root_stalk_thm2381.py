#!/usr/bin/env python3
"""Exact companion for THM-2381.

The proof is analytic.  This dependency-free Fraction/integer script checks
the two finite root-count mechanisms and an independent wall-cell instance of
the sharp mixed-depth anti-shield constant.
"""

from fractions import Fraction as F
from math import lcm


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def danger_rat(num: int, den: int, width: int = 1) -> bool:
    """Return ||num/den|| < width/14, exactly."""
    rem = num % den
    dist = min(rem, den - rem)
    return 14 * dist < width * den


def endpoint_rat(num: int, den: int, width: int = 1) -> bool:
    rem = num % den
    dist = min(rem, den - rem)
    return 14 * dist == width * den


def root_word_13(width: int, z_num: int, z_den: int) -> frozenset[int]:
    """Grid word {h: ||(z+h)/13|| < width/14}."""
    den = 13 * z_den
    return frozenset(
        h for h in range(13)
        if danger_rat(z_num + h * z_den, den, width)
    )


def fibre_word(v: int, z_num: int, z_den: int, N: int) -> frozenset[int]:
    """Word {j: ||v(z+j)/N||<1/14}, for z=z_num/z_den."""
    den = z_den * N
    return frozenset(
        j for j in range(N)
        if danger_rat(v * (z_num + z_den * j), den)
    )


def wall_cell_measure(speeds_in: tuple[int, ...], speeds_out: tuple[int, ...]) -> F:
    """Measure of points dangerous for every in-speed and safe for out-speeds.

    The common wall denominator 14*lcm(speeds) resolves every open cell.
    Endpoints have measure zero, so midpoint membership is exact.
    """
    speeds = speeds_in + speeds_out
    D = 14 * lcm(*speeds)
    total = 0
    for cell in range(D):
        # midpoint=(2*cell+1)/(2D)
        num = 2 * cell + 1
        den = 2 * D
        if (
            all(danger_rat(v * num, den) for v in speeds_in)
            and all(not danger_rat(v * num, den) for v in speeds_out)
        ):
            total += 1
    return F(total, D)


# The thirteen-root support squeeze.
root_phases = 0
ordinary_sizes: set[int] = set()
guard_sizes: set[int] = set()
for z_num in range(337):
    z_den = 337
    ordinary = root_word_13(1, z_num, z_den)
    guard = root_word_13(2, z_num, z_den)
    # 337 is coprime to 2*7*13, so no sampled phase is an endpoint.
    for h in range(13):
        check(
            not endpoint_rat(z_num + h * z_den, 13 * z_den, 1),
            "ordinary thirteen-root sample hit an endpoint",
        )
        check(
            not endpoint_rat(z_num + h * z_den, 13 * z_den, 2),
            "guard thirteen-root sample hit an endpoint",
        )
    ordinary_sizes.add(len(ordinary))
    guard_sizes.add(len(guard))
    root_phases += 1
check(ordinary_sizes == {1, 2}, "ordinary root-word size classification failed")
check(guard_sizes == {3, 4}, "guard root-word size classification failed")
check(
    min(guard_sizes) > max(ordinary_sizes),
    "guard/top-unit support squeeze failed",
)


# Seven-bin counts and mixed-depth intersections.
unit_reps = tuple(u for u in range(1, 14) if u % 7)
phase_reps = (1, 19, 73, 127, 191)
phase_den = 211
bin_cases = 0
mixed_cases = 0

for M in range(1, 4):
    N = 7 ** (M + 1)

    # A depth-M width-one factor occupies one whole residue bin.
    for u_top in unit_reps:
        top = 7**M * u_top
        for z_num in phase_reps:
            word_top = fibre_word(top, z_num, phase_den, N)
            check(len(word_top) == N // 7, "top word has wrong size")
            occupied = {j % 7 for j in word_top}
            check(len(occupied) == 1, "top word is not one complete bin")

    for r in range(M):
        for u_low in unit_reps:
            low = 7**r * u_low
            for z_num in phase_reps:
                word_low = fibre_word(low, z_num, phase_den, N)
                check(len(word_low) == N // 7, "low word has wrong size")
                for residue in range(7):
                    count = sum(j % 7 == residue for j in word_low)
                    check(count == N // 49, "low per-bin count failed")
                    bin_cases += 1

                # Width two has exactly twice the same per-bin count.
                den = phase_den * N
                word_guard = frozenset(
                    j for j in range(N)
                    if danger_rat(
                        low * (z_num + phase_den * j), den, width=2
                    )
                )
                check(
                    len(word_guard) == 2 * N // 7,
                    "width-two low word has wrong size",
                )
                for residue in range(7):
                    count = sum(j % 7 == residue for j in word_guard)
                    check(
                        count == 2 * N // 49,
                        "width-two low per-bin count failed",
                    )
                    bin_cases += 1

        for u_low in unit_reps:
            low = 7**r * u_low
            for u_top in unit_reps:
                top = 7**M * u_top
                for z_num in phase_reps:
                    word_low = fibre_word(low, z_num, phase_den, N)
                    word_top = fibre_word(top, z_num, phase_den, N)
                    check(len(word_low) == N // 7, "mixed low word size failed")
                    check(len(word_top) == N // 7, "mixed top word size failed")
                    check(
                        len(word_low & word_top) == N // 49,
                        "mixed low/top intersection failed",
                    )
                    check(
                        len(word_low - word_top) == 6 * N // 49,
                        "mixed low-minus-top count failed",
                    )
                    mixed_cases += 1


anti_shield = F(6, 7) * F(6, 49)
check(anti_shield == F(36, 343), "anti-shield arithmetic failed")
cell_control = wall_cell_measure((1,), (7, 49))
check(cell_control == anti_shield, "wall-cell anti-shield control failed")


print("THM-2381 one-top-one-blocker septimal root-stalk closure")
print(f"thirteen-root generic phases: {root_phases}")
print(f"ordinary root-word sizes: {sorted(ordinary_sizes)}")
print(f"guard root-word sizes: {sorted(guard_sizes)}")
print(f"seven-bin lower-factor counts checked: {bin_cases}")
print(f"mixed-depth fibres checked: {mixed_cases}")
print("mixed counts: low=N/7, top=N/7, intersection=N/49")
print(f"anti-shield floor: {anti_shield}")
print(f"independent D_1 \\\\ (D_7 union D_49) wall-cell measure: {cell_control}")
print("closed septimal alternative: k=2, (t,b)=(1,1)")
print("remaining k=2 alternatives: (1,0), (2,0), (5,2)")
print("ledger=165; LRC(14) remains open")
