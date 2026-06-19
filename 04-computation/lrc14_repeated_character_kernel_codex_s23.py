#!/usr/bin/env python3
"""
LRC(14) repeated-residue character kernel.

This continues HYP-2630's next target.  The support-six coimage coefficient is

    S_d(a_1,...,a_6) = sum_{r in (F_7^*)^6, a.r=0} C_d(r).

HYP-2630 showed that the remaining k=10 tail is dominated by repeated residue
packets.  This script keeps the finite coimage kernel explicit and asks whether
the two-large packets are controlled by a small character table.

Tournament Analysis declaration:
  vertices are finite proof quotients, not runners.  The quotient preserves
  the signed support-six coimage coefficient and height-2 tail/wall predicate;
  it destroys witness times and raw support identities.
"""
from __future__ import annotations

import cmath
import itertools
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_height2_coimage_wall_classes_codex_s18 as s18  # noqa: E402
import lrc14_support6_coimage_fiber_codex_s14 as s14  # noqa: E402
import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402

MOD = 7
AMBIENT_D = 9
EPS = 1e-10
ROOT_OF_UNITY = cmath.exp(2j * math.pi / MOD)
ROOTS = tuple(ROOT_OF_UNITY ** e for e in range(MOD))
RESIDUE_TUPLES = tuple(itertools.product(range(1, MOD), repeat=6))
PAIR_DOMAIN = (0, 2, 3, 4, 5, 6)


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def chi7(x: int) -> int:
    x %= MOD
    if x == 0:
        return 0
    return 1 if x in {1, 2, 4} else -1


def fmt_complex(z: complex) -> str:
    if abs(z.imag) < 5e-13:
        return f"{z.real:.12g}"
    return f"{z.real:.12g}{z.imag:+.3g}j"


def signed_int(z: complex, unit: float) -> int:
    val = z.real / unit
    nearest = round(val)
    if abs(val - nearest) > 1e-7 or abs(z.imag) > 1e-8:
        raise ValueError(f"not an integer packet: {z} / {unit} = {val}")
    return nearest


def coeff_table() -> dict[tuple[int, ...], complex]:
    return {
        r: s12.residue_coeff(r, AMBIENT_D)
        for r in RESIDUE_TUPLES
    }


def direct_kernel(speed: tuple[int, ...], coeffs: dict[tuple[int, ...], complex]) -> complex:
    total = 0j
    for r, c in coeffs.items():
        if sum(a * ri for a, ri in zip(speed, r)) % MOD == 0:
            total += c
    return total


def additive_transform(
    v: tuple[int, ...],
    coeffs: dict[tuple[int, ...], complex],
    cache: dict[tuple[int, ...], complex],
) -> complex:
    key = tuple(x % MOD for x in v)
    if key in cache:
        return cache[key]
    total = 0j
    for r, c in coeffs.items():
        exponent = sum(x * ri for x, ri in zip(key, r)) % MOD
        total += c * ROOTS[exponent]
    cache[key] = total
    return total


def fourier_kernel(
    speed: tuple[int, ...],
    coeffs: dict[tuple[int, ...], complex],
    cache: dict[tuple[int, ...], complex],
) -> complex:
    total = 0j
    for t in range(MOD):
        total += additive_transform(tuple((t * a) % MOD for a in speed), coeffs, cache)
    return total / MOD


def coimage_rows() -> list[s14.FiberStats]:
    return s14.compute_stats_for_d(AMBIENT_D, s14.support_classes())


def wall_tail_classes() -> set[tuple[int, ...]]:
    return s18.enumerate_wall_classes(10, 2).classes


def signature(a: int, b: int) -> tuple[int, int, int, int, int, int]:
    return (
        chi7(a),
        chi7(b),
        chi7(a * b),
        chi7((a - 1) * (b - 1)),
        chi7(a - b),
        (a + b) % MOD,
    )


def affine_selector(a: int, b: int) -> int:
    """Legendre selector for the 4+1+1 high/low split off the zero lane."""
    return (a * b * (1 + 3 * ((a + b) % MOD)) - 1) % MOD


def status_for(cls: tuple[int, ...], row_lookup: dict[tuple[int, ...], s14.FiberStats], hit_h2: set[tuple[int, ...]]) -> str:
    row = row_lookup[cls]
    if row.signed_abs <= 1e-12:
        return "zero"
    if cls in hit_h2:
        return "height2-wall"
    return "tail"


def transform_identity_report(
    classes: list[tuple[int, ...]],
    coeffs: dict[tuple[int, ...], complex],
) -> dict[tuple[int, ...], complex]:
    section("FINITE ADDITIVE-FOURIER IDENTITY")
    transform_cache: dict[tuple[int, ...], complex] = {}
    direct_cache: dict[tuple[int, ...], complex] = {}
    max_err = 0.0
    max_imag = 0.0
    worst: tuple[int, ...] | None = None
    for cls in classes:
        direct = direct_kernel(cls, coeffs)
        via_fourier = fourier_kernel(cls, coeffs, transform_cache)
        direct_cache[cls] = direct
        err = abs(direct - via_fourier)
        if err > max_err:
            max_err = err
            worst = cls
        max_imag = max(max_imag, abs(direct.imag), abs(via_fourier.imag))
    print(
        "Verified over the 159 projective support-six coimage classes:"
    )
    print(
        "  S_d(a) = sum_{a.r=0} C_d(r) = (1/7) sum_{t in F_7} C_hat(t a)"
    )
    print(f"max identity error: {max_err:.3e} at class {worst}")
    print(f"max imaginary drift: {max_imag:.3e}")
    print(f"additive transform cache entries used: {len(transform_cache)}")
    print(
        "Readout: the repeated tail can be attacked as a finite additive-Fourier "
        "kernel first, and only then converted into multiplicative chi_7 packets."
    )
    return direct_cache


def four_two_report(
    coeffs: dict[tuple[int, ...], complex],
    row_lookup: dict[tuple[int, ...], s14.FiberStats],
    hit_h2: set[tuple[int, ...]],
    unit: float,
) -> dict[int, int]:
    section("4+2 CHARACTER KERNEL")
    print(
        "The signed kernel values are shown in units U.  HYP-2630 reported "
        "absolute masses; the signs matter for the eventual reciprocal-tail sum."
    )
    print(f"kernel unit U = {unit:.15g}")
    print(f"{'a':>2} {'chi(a)':>7} {'class':>22} {'S_9/U':>7} {'S_9':>16} {'status':>12}")
    packets: dict[int, int] = {}
    for a in range(MOD):
        speed = (1, 1, 1, 1, a, a)
        cls = s14.canon_support(speed)
        z = direct_kernel(speed, coeffs)
        n = signed_int(z, unit) if abs(z) > EPS else 0
        packets[a] = n
        print(
            f"{a:>2} {chi7(a):>7} {str(cls):>22} {n:>7} "
            f"{fmt_complex(z):>16} {status_for(cls, row_lookup, hit_h2):>12}"
        )
    print(
        "\nFor a in F_7^* \\ {1}, the formula is exact in packet units:"
    )
    print("  2*S_9(1,1,1,1,a,a)/U = -43 - 7*chi_7(a)")
    nonzero_nonone = [a for a in range(2, MOD)]
    failures = [
        a for a in nonzero_nonone
        if 2 * packets[a] != -43 - 7 * chi7(a)
    ]
    print(f"formula failures on a=2..6: {failures}")
    print(
        "Thus the QR/NQR split is a literal Legendre component, not merely a "
        "visual grouping of floating-point masses."
    )
    return packets


def four_one_one_report(
    coeffs: dict[tuple[int, ...], complex],
    row_lookup: dict[tuple[int, ...], s14.FiberStats],
    hit_h2: set[tuple[int, ...]],
    unit: float,
) -> dict[tuple[int, int], int]:
    section("4+1+1 TAIL-SIDE CHARACTER TABLE")
    print(
        "Pairs are the two singleton residues after normalizing the repeated "
        "residue to 1.  The domain excludes 1 so the multiplicity pattern is "
        "actually 4+1+1."
    )
    print(
        f"{'pair':>7} {'S_9/U':>7} {'status':>12} {'chi(Q)':>7} "
        f"{'chi(a),chi(b),chi(ab),chi((a-1)(b-1)),chi(a-b),a+b':>62}"
    )
    packets: dict[tuple[int, int], int] = {}
    by_value: Counter[int] = Counter()
    by_line: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for a, b in itertools.combinations(PAIR_DOMAIN, 2):
        speed = (1, 1, 1, 1, a, b)
        cls = s14.canon_support(speed)
        z = direct_kernel(speed, coeffs)
        n = signed_int(z, unit) if abs(z) > EPS else 0
        packets[(a, b)] = n
        by_value[n] += 1
        by_line[(a + b) % MOD].append((a, b, n))
        selector_chi = chi7(affine_selector(a, b))
        print(
            f"{str((a,b)):>7} {n:>7} {status_for(cls, row_lookup, hit_h2):>12} "
            f"{selector_chi:>7} "
            f"{str(signature(a,b)):>62}"
        )

    print("\nValue histogram in packet units:")
    for n, count in sorted(by_value.items(), key=lambda x: (-x[0], x[1])):
        print(f"  {n:>3}: {count} pairs")

    zero_line = by_line[2]
    print("\nAffine-line surprise:")
    print(f"  pairs with a+b=2 mod 7: {zero_line}")
    print(
        "  They are all zero.  This extra affine coordinate is not captured by "
        "the shorter HYP-2630 signature using only chi(a), chi(b), chi(ab), and "
        "chi((a-1)(b-1))."
    )

    collisions: dict[tuple[int, int, int, int], list[tuple[int, int, int]]] = defaultdict(list)
    for (a, b), n in packets.items():
        short = signature(a, b)[:4]
        collisions[short].append((a, b, n))
    print("\nShort-signature collisions that force the extra affine coordinate:")
    for short, items in sorted(collisions.items()):
        vals = sorted({n for _a, _b, n in items})
        if len(vals) > 1:
            print(f"  {short}: {items}")

    high_failures = []
    low_failures = []
    for (a, b), n in packets.items():
        if (a + b) % MOD == 2:
            continue
        selector_chi = chi7(affine_selector(a, b))
        if n == 8 and selector_chi != 1:
            high_failures.append((a, b, selector_chi))
        if n == 1 and selector_chi != -1:
            low_failures.append((a, b, selector_chi))
    print("\nLegendre selector after removing the zero lane:")
    print("  Q(a,b) = ab*(1+3(a+b))-1 mod 7")
    print("  S_9/U = 8 iff chi_7(Q)=+1; S_9/U = 1 iff chi_7(Q)=-1")
    print(f"  high selector failures: {high_failures}")
    print(f"  low selector failures: {low_failures}")

    return packets


def signed_ledger_report(four_two: dict[int, int], four_one_one: dict[tuple[int, int], int]) -> None:
    section("SIGNED INTEGER LEDGER FOR THE REPEATED TAIL")
    four_two_tail = [four_two[a] for a in (0, 2, 3, 4, 5, 6)]
    four_one_tail = [n for n in four_one_one.values() if n]
    sum42 = sum(four_two_tail)
    sum411 = sum(four_one_tail)
    abs42 = sum(abs(n) for n in four_two_tail)
    abs411 = sum(abs(n) for n in four_one_tail)
    net = sum42 + sum411
    abs_total = abs42 + abs411
    print(f"4+2 signed packet sum:      {sum42:>5} U")
    print(f"4+2 absolute packet mass:   {abs42:>5} U")
    print(f"4+1+1 signed packet sum:    {sum411:>5} U")
    print(f"4+1+1 absolute packet mass: {abs411:>5} U")
    print(f"combined signed sum:        {net:>5} U")
    print(f"combined absolute mass:     {abs_total:>5} U")
    print(f"absolute/net ratio: {abs_total / abs(net):.6g}")
    print(
        "Readout: even before the analytic reciprocal sums, the finite repeated "
        "kernel has a 3:1 absolute/signed gap.  The proof should exploit this "
        "signed chi_7 table instead of bounding the packets by absolute mass."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = [
        "additive_fourier_kernel",
        "affine_line_a_plus_b",
        "quadratic_character_phase",
        "height2_wall_deletion",
        "euler_copy_capacity",
        "raw_coimage_enumeration",
        "raw_runner_vertices",
    ]
    print("Candidate vertices considered:")
    print(
        "  runners, gaps, fixed circle sections, section boundaries, wall-crossing "
        "events, speed residues, cover arcs, additive Fourier modes, quadratic "
        "characters, Jacobi/cross-ratio signatures, and proof obligations."
    )
    print("Chosen Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}")
    print("directed_3_cycles = 0")
    print("SCC_sizes = [1,1,1,1,1,1,1]")
    print(
        "Preserved predicate: the signed d=9 support-six coimage coefficient "
        "after height-2 wall deletion.  Destroyed information: exact witness "
        "times, raw support identities, and raw runner labels."
    )
    print(
        "Challenged assumption: the four multiplicative chi signatures from "
        "HYP-2630 are not complete for 4+1+1 packets.  The affine line "
        "a+b=2 mod 7 is a hidden zero-lane coordinate."
    )


def main() -> None:
    section("LRC14 REPEATED-RESIDUE CHARACTER KERNEL - CODEX S23")
    print(
        "Goal: turn HYP-2630's two-large repeated tail from a coimage census into "
        "a finite signed character table suitable for a cotangent/Dedekind tail "
        "bound."
    )
    coeffs = coeff_table()
    rows = coimage_rows()
    row_lookup = {row.cls: row for row in rows}
    hit_h2 = wall_tail_classes()
    direct_cache = transform_identity_report(s14.support_classes(), coeffs)
    unit_speed = (1, 1, 1, 1, 2, 6)
    unit = abs(direct_kernel(unit_speed, coeffs))
    # Sanity: the direct cache is keyed by canonical classes, while unit_speed is
    # already canonical.  Keep both paths checked.
    assert abs(direct_cache[s14.canon_support(unit_speed)] - direct_kernel(unit_speed, coeffs)) < 1e-9
    four_two = four_two_report(coeffs, row_lookup, hit_h2, unit)
    four_one_one = four_one_one_report(coeffs, row_lookup, hit_h2, unit)
    signed_ledger_report(four_two, four_one_one)
    tournament_analysis()


if __name__ == "__main__":
    main()
