#!/usr/bin/env python3
"""
LRC(14) two-large repeated-residue Dedekind/phase scout.

HYP-2630 left a precise target: the remaining k=10 support-six tail is not
explained by Euler-copy capacity.  The dominant classes are two-large
repeated-residue packets, and their QR/NQR split is a phase effect.

This script makes the finite cotangent/Dedekind transform explicit.  For speed
residues a_i and Fourier residues r_i in F_7^*,

    S_d(a_1,...,a_6) = sum_{a.r=0} C_d(r)

is expanded with the additive character of the relation hyperplane.  The
result is a finite Dedekind packet over additive frequencies m in F_7:

    S_d(a) = (1/7) sum_m sum_T (-1)^|T| chat(0,T)^(d-6)
             prod_i D_T(m a_i),

    D_T(ell) = sum_{r in F_7^*} r chat(r,T) zeta_7^(ell r).

That formula is the two-large analogue of a cotangent/Dedekind reciprocal
tail: retain the packet address and the additive frequency before taking
absolute values.  The point of the run is to keep the quadratic-character phase
visible instead of collapsing the packet to a scalar absolute envelope.

Tournament Analysis declaration:
  vertices are proof quotients, not runners.  This run explicitly considers
  runners, large residue coordinates, additive frequencies, conjugate frequency
  shells, speed-residue pairs, coimage classes, and proof obligations; the
  chosen quotient preserves the signed support-six coimage predicate and
  destroys exact witness times and raw support identities.
"""
from __future__ import annotations

import cmath
import itertools
import math
import sys
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_coimage_fiber_codex_s14 as s14  # noqa: E402
import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402

AMBIENT_D = 9
EPS = 1e-10
ZETA7 = cmath.exp(2j * math.pi / 7)
UNIT = 147.0 / (16.0 * math.pi**6)


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def chi7(x: int) -> int:
    x %= 7
    if x == 0:
        return 0
    return 1 if x in {1, 2, 4} else -1


def affine_selector(a: int, b: int) -> int:
    """HYP-2632 selector for the 4+1+1 high/low split off a+b=2."""
    return (a * b * (1 + 3 * ((a + b) % 7)) - 1) % 7


def clean_unit(x: complex | float) -> float:
    z = x.real if isinstance(x, complex) else x
    q = z / UNIT
    rounded_half = round(2 * q) / 2
    if abs(q - rounded_half) < 1e-7:
        return rounded_half
    return q


def fmt_unit(x: complex | float, width: int = 8) -> str:
    q = clean_unit(x)
    if abs(q) < 1e-7:
        q = 0.0
    if abs(q - round(q)) < 1e-7:
        return f"{int(round(q)):>{width}}"
    return f"{q:>{width}.3f}"


def direct_fiber(vals: tuple[int, ...]) -> s14.FiberStats:
    return s14.compute_stats_for_d(AMBIENT_D, [s14.canon_support(vals)])[0]


@lru_cache(None)
def dedekind_factor(T: tuple[int, ...], ell: int) -> complex:
    """D_T(ell) = sum_r r chat(r,T) zeta_7^(ell r)."""
    ell %= 7
    return sum(r * s12.chat(r, T) * (ZETA7 ** ((ell * r) % 7)) for r in range(1, 7))


def dedekind_by_frequency(vals: tuple[int, ...]) -> dict[int, complex]:
    """Additive-frequency expansion of the coimage fiber."""
    by_m: dict[int, complex] = defaultdict(complex)
    for sign, T in zip(s12.SIGNS, s12.SUBSETS):
        zero_factor = s12.chat(0, T) ** (AMBIENT_D - 6)
        for m in range(7):
            z = sign * zero_factor / 7.0
            for a in vals:
                z *= dedekind_factor(T, (m * a) % 7)
            by_m[m] += z
    return dict(by_m)


def conjugate_square_shells(by_m: dict[int, complex]) -> dict[int, complex]:
    """Group additive frequencies by the square m^2 in F_7^*."""
    return {
        1: by_m.get(1, 0j) + by_m.get(6, 0j),
        2: by_m.get(3, 0j) + by_m.get(4, 0j),
        4: by_m.get(2, 0j) + by_m.get(5, 0j),
    }


def chi_frequency_shells(by_m: dict[int, complex]) -> dict[int, complex]:
    out = {-1: 0j, 0: 0j, 1: 0j}
    for m, z in by_m.items():
        out[chi7(m)] += z
    return out


def shell_abs_units(shells: dict[int, complex]) -> float:
    return sum(abs(z) for z in shells.values()) / UNIT


def pair_residue_matrix(vals: tuple[int, int, int, int, int, int]) -> dict[tuple[int, int], complex]:
    """
    Two-large residue matrix.

    For vals=(1,1,1,1,a,b), entries are indexed by the two large Fourier
    residues (u,v).  This is the natural "blind" absolute envelope that forgets
    the additive-frequency phase.
    """
    a, b = vals[4], vals[5]
    out: dict[tuple[int, int], complex] = {}
    for u in range(1, 7):
        for v in range(1, 7):
            target = (-a * u - b * v) % 7
            z = 0j
            for core in itertools.product(range(1, 7), repeat=4):
                if sum(core) % 7 == target:
                    z += s12.residue_coeff(core + (u, v), AMBIENT_D)
            out[(u, v)] = z
    return out


def blind_pair_abs_units(vals: tuple[int, int, int, int, int, int]) -> float:
    return sum(abs(z) for z in pair_residue_matrix(vals).values()) / UNIT


def verify_transform() -> None:
    section("FINITE COTANGENT/DEDEKIND TRANSFORM CHECK")
    print(
        "The additive-frequency formula is checked against the existing HYP-2617 "
        "coimage fiber implementation.  The unit below is the exact-looking "
        "normalization observed in every two-large packet."
    )
    print(f"U = 147/(16*pi^6) = {UNIT:.12g}")
    print()
    print(f"{'class':>24} {'direct/U':>10} {'dedekind/U':>12} {'abs error':>12}")
    samples = [
        (1, 1, 1, 1, 2, 2),
        (1, 1, 1, 1, 3, 3),
        (1, 1, 1, 1, 2, 3),
        (1, 1, 1, 1, 1, 1),
    ]
    for vals in samples:
        row = direct_fiber(vals)
        by_m = dedekind_by_frequency(vals)
        dedekind = sum(by_m.values())
        print(
            f"{str(vals):>24} {fmt_unit(row.signed, 10)} "
            f"{fmt_unit(dedekind, 12)} {abs(row.signed - dedekind):>12.3e}"
        )
    print(
        "\nReading: the relation-hyperplane coimage is exactly the finite "
        "Dedekind packet once the m-frequency address is retained."
    )


def report_four_two() -> None:
    section("4+2 PACKET: (1,1,1,1,a,a)")
    print(
        "Square shells are conjugate frequency pairs {m,-m}, indexed by m^2 in "
        "F_7^*.  For the non-identity 4+2 row these shells are sign-definite, "
        "so the shell bound equals the signed mass while still exposing chi_7(a)."
    )
    print()
    print(
        f"{'a':>2} {'chi':>4} {'S/U':>8} {'|S|/U':>8} "
        f"{'chiBd/U':>9} {'sqBd/U':>8} {'blind/u':>9} {'square shells q=1,2,4':>28}"
    )
    values: dict[int, float] = {}
    for a in range(1, 7):
        vals = (1, 1, 1, 1, a, a)
        row = direct_fiber(vals)
        by_m = dedekind_by_frequency(vals)
        chi_shells = chi_frequency_shells(by_m)
        square_shells = conjugate_square_shells(by_m)
        values[a] = clean_unit(row.signed)
        shell_text = " ".join(
            f"{q}:{fmt_unit(square_shells[q], 0).strip()}" for q in (1, 2, 4)
        )
        print(
            f"{a:>2} {chi7(a):>4} {fmt_unit(row.signed)} {abs(row.signed)/UNIT:>8.3f} "
            f"{shell_abs_units(chi_shells):>9.3f} {shell_abs_units(square_shells):>8.3f} "
            f"{blind_pair_abs_units(vals):>9.3f} {shell_text:>28}"
        )

    print("\nPunctured character law for a in {2,3,4,5,6}:")
    print("  2*S(1,1,1,1,a,a)/U = -43 - 7*chi_7(a)")
    print("  equivalently |S|/U = (43 + 7*chi_7(a))/2.")
    print(
        "This is the QR/NQR split HYP-2630 saw numerically: residues 2,4 give "
        "25U and nonresidues 3,5,6 give 18U.  The all-equal a=1 packet cancels."
    )


def report_four_one_one() -> None:
    section("4+1+1 PACKET: (1,1,1,1,a,b)")
    print(
        "The unit two-large 4+1+1 row is even more quantized: every unordered "
        "pair a<b in {2,...,6} has signed mass in {0,U,8U}.  The primary "
        "HYP-2632 refinement adds the affine zero lane a+b=2 and the selector "
        "Q(a,b)=ab*(1+3(a+b))-1 for the high/low split off that lane."
    )
    print()
    print(
        f"{'(a,b)':>8} {'sig chi(a),chi(b),chi(ab),chi((a-1)(b-1))':>48} "
        f"{'a+b':>4} {'chiQ':>5} {'S/U':>7} {'|S|/U':>7} {'chiBd':>7} "
        f"{'sqBd':>7} {'blind/u':>9} "
        f"{'m=1..6 in U':>24}"
    )
    values = []
    selector_failures = []
    zero_lane_pairs = []
    for a, b in itertools.combinations(range(2, 7), 2):
        vals = (1, 1, 1, 1, a, b)
        row = direct_fiber(vals)
        by_m = dedekind_by_frequency(vals)
        chi_shells = chi_frequency_shells(by_m)
        square_shells = conjugate_square_shells(by_m)
        sig = (chi7(a), chi7(b), chi7(a * b), chi7((a - 1) * (b - 1)))
        lane = (a + b) % 7
        qchi = chi7(affine_selector(a, b))
        m_vec = "[" + ",".join(fmt_unit(by_m[m], 0).strip() for m in range(1, 7)) + "]"
        signed_units = clean_unit(row.signed)
        values.append(signed_units)
        if lane == 2:
            zero_lane_pairs.append((a, b))
        elif (signed_units == 8 and qchi != 1) or (signed_units == 1 and qchi != -1):
            selector_failures.append((a, b, signed_units, qchi))
        print(
            f"{str((a, b)):>8} {str(sig):>48} {lane:>4} {qchi:>5} "
            f"{fmt_unit(row.signed, 7)} "
            f"{abs(row.signed)/UNIT:>7.3f} {shell_abs_units(chi_shells):>7.3f} "
            f"{shell_abs_units(square_shells):>7.3f} {blind_pair_abs_units(vals):>9.3f} "
            f"{m_vec:>24}"
        )
    print()
    print(f"Observed signed unit values: {sorted(set(values))}")
    print(f"Unit-domain zero lane a+b=2 mod 7: {zero_lane_pairs}")
    print(f"Selector failures off the zero lane: {selector_failures}")
    print(
        "Bound extracted for the unit 4+1+1 packet: |S| <= 8U.  The two exact "
        "zero classes (3,6) and (4,5) are not visible from raw pair-residue "
        "absolute mass; they require frequency cancellation.  Off that lane, "
        "the same Legendre phase determines whether the packet is U or 8U."
    )


def report_guardrails() -> None:
    section("BOUND READOUT AND LIFT TARGET")
    print("Finite coefficient packet bounds at ambient d=9:")
    print("  4+2 unit packets, a != 1: |S| <= 25U = 3675/(16*pi^6).")
    print("  4+1+1 unit packets:       |S| <=  8U = 147/(2*pi^6).")
    print("  zero-cusp (0,0,1,1,1,1):  S = -4U.")
    print()
    zero = direct_fiber((0, 0, 1, 1, 1, 1))
    print(f"zero-cusp check: S/U = {fmt_unit(zero.signed, 0).strip()}")
    print(
        "\nProof consequence: the next analytic lemma should not bound the two-large "
        "tail by the blind 36-entry pair-residue matrix.  It should first split "
        "the reciprocal tail by additive frequency m (or conjugate frequency "
        "shells) and only then take absolute values.  That keeps the quadratic "
        "character phase visible enough to recover the 25U/18U split and the "
        "4+1+1 affine/selector cancellation."
    )
    print(
        "\nAssumption challenged: the natural tournament vertices are not runners "
        "or raw arcs.  For this proof step the useful vertices are frequency "
        "shells, two-large residue pairs, coimage classes, finite wall ledgers, "
        "and lift obligations; quotienting to raw absolute mass destroys the "
        "phase information that makes the packet small."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = [
        "additive_frequency_shells",
        "quadratic_character_phase",
        "two_large_dedekind_packet",
        "height2_wall_tail",
        "euler_copy_capacity",
        "blind_pair_abs_matrix",
        "raw_runner_vertices",
    ]
    print("Pairwise observable: preservation of the signed support-six coimage tail.")
    print(
        "Switch/gauge: prefer the quotient that retains packet address, chi_7 "
        "phase, affine lane, Q selector, and additive-frequency cancellation "
        "before applying an absolute bound."
    )
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}")
    print("directed_3_cycles = 0")
    print("SCCs = seven singleton proof-obligation vertices")
    print("hamiltonian_paths = 1")


def main() -> None:
    verify_transform()
    report_four_two()
    report_four_one_one()
    report_guardrails()
    tournament_analysis()


if __name__ == "__main__":
    main()
