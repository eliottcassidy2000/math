#!/usr/bin/env python3
"""
codex-S110: coprime density, Euler-totient packets, and the three tiling
recursion modes.

The point of this audit is not another raw LRC search.  It makes precise a
common algebraic skeleton:

  * Euler/totient:   phi = mu * id on the divisor lattice.
  * Full tiling:     +++---+ is Boolean Mobius on three generators.
  * Even half:       ++- is Boolean Mobius on two generators.
  * Odd half:        ++-+--+ is again Boolean Mobius on three generators after
                     the corner letters are ordered correctly.

The LRC proof lesson is that exact-period denominator packets and tiling/far
packets should be retained as a product-lattice address before scalarizing.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, gcd, pi


F = Fraction


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def divisors(n: int) -> list[int]:
    ds = [1]
    for p, a in factor(n).items():
        ds = [d * p**e for d in ds for e in range(a + 1)]
    return sorted(ds)


def mobius(n: int) -> int:
    fac = factor(n)
    if any(a > 1 for a in fac.values()):
        return 0
    return -1 if len(fac) % 2 else 1


def phi(n: int) -> int:
    out = n
    for p in factor(n):
        out = out // p * (p - 1)
    return out


def phi_by_mobius(n: int) -> int:
    return sum(mobius(d) * (n // d) for d in divisors(n))


def fmt_factor(n: int) -> str:
    fac = factor(n)
    if not fac:
        return "1"
    return "*".join(f"{p}^{a}" if a > 1 else str(p) for p, a in fac.items())


def full_cells(n: int) -> int:
    return comb(n - 1, 2)


def half_cells(n: int) -> int:
    return ((n - 1) * (n - 1)) // 4


def sign_string(vals: list[int]) -> str:
    return "".join("+" if v > 0 else "-" if v < 0 else "0" for v in vals)


def legendre7(a: int) -> int:
    r = a % 7
    if r == 0:
        return 0
    return 1 if r in {1, 2, 4} else -1


def prompt_chi7_slot(a: int) -> int:
    return 1 if a % 7 == 0 else legendre7(a)


def print_totient_layer() -> None:
    print("=" * 96)
    print("A. Euler-copy / coprime-density layer")
    print("=" * 96)
    rows = [14, 30, 84, 89, 173, 210, 425, 1260, 2520]
    print("phi is the divisor-Mobius inverse of the copy rule sum_{d|n} c(d)=n:")
    for n in rows:
        sigma_phi = sum(phi(d) for d in divisors(n))
        density = F(phi(n), n)
        print(
            f"  n={n:4d}={fmt_factor(n):<14s} phi={phi(n):4d}, "
            f"sum_phi_div={sigma_phi:4d}, mu*id={phi_by_mobius(n):4d}, "
            f"phi/n={density}"
        )
    print()
    print("Mertens/Farey coprime-density approximations:")
    target1 = 3 / (pi * pi)
    target2 = 6 / (pi * pi)
    for q in [14, 41, 83, 210, 1000]:
        s_phi = sum(phi(b) for b in range(1, q + 1))
        s_phi_over_b = sum(F(phi(b), b) for b in range(1, q + 1))
        print(
            f"  Q={q:4d}: sum_phi/Q^2={s_phi/(q*q):.6f} "
            f"(target 3/pi^2={target1:.6f}); "
            f"sum(phi/b)/Q={float(s_phi_over_b/q):.6f} "
            f"(target 6/pi^2={target2:.6f})"
        )
    print()


def print_recursion_layer() -> None:
    print("=" * 96)
    print("B. The three tiling recurrences as Boolean-Mobius kernels")
    print("=" * 96)
    full_letters = [1, 1, 1, -1, -1, -1, 1]
    even_letters = [1, 1, -1]
    odd_prompt_letters = [1, 1, -1, 1, -1, -1, 1]
    odd_corner_reorder = [1, 1, 1, -1, -1, -1, 1]
    print("Letter-level sign packets:")
    print(f"  full B3 A+B+C-D-E-F+G:          {sign_string(full_letters)}")
    print(f"  even half B2 A+B-C:             {sign_string(even_letters)}")
    print(f"  odd half prompt A+B-C+D-E-F+G:  {sign_string(odd_prompt_letters)}")
    print(f"  odd half reordered corners A,B,D: {sign_string(odd_corner_reorder)}")
    print("Readout: odd half has two useful addresses.  In prompt order over")
    print("A..G = 1..7 it is chi_7 with the zero/triple slot positive; in corner order A,B,D")
    print("it is B3 inclusion-exclusion with pair-overlaps C,E,F.")
    print("  chi_7 plus zero/triple slot:       " + sign_string([prompt_chi7_slot(a) for a in range(1, 8)]))
    print()

    print("Pushforward by size/depth offset:")
    print("  full staircase:  offsets n-1,n-2,n-3 ->  3, -3, +1")
    print("  even half:       offsets n-1,n-2     ->  2, -1")
    print("  odd half:        offsets n-1..n-4    ->  2,  0, -2, +1")
    print("The zero at n-2 is geometric cancellation (-C + D), not absence of data.")
    print()

    full_ok = []
    even_ok = []
    odd_ok = []
    global_half_ok = []
    for n in range(5, 21):
        full_ok.append(full_cells(n) == 3 * full_cells(n - 1) - 3 * full_cells(n - 2) + full_cells(n - 3))
        if n % 2 == 0:
            even_ok.append(half_cells(n) == 2 * half_cells(n - 1) - half_cells(n - 2))
        else:
            odd_ok.append(
                half_cells(n)
                == 2 * half_cells(n - 1) - half_cells(n - 2)
                + half_cells(n - 2) - 2 * half_cells(n - 3) + half_cells(n - 4)
            )
        global_half_ok.append(
            half_cells(n) - 2 * half_cells(n - 1) + 2 * half_cells(n - 3) - half_cells(n - 4) == 0
        )
    print(f"Verification n=5..20: full Delta^3={all(full_ok)}, even half={all(even_ok)}, odd half={all(odd_ok)}")
    print(f"Global half recurrence (x-1)^3(x+1): {all(global_half_ok)}")
    print("Half sequence h_n=floor((n-1)^2/4), n=2..14:")
    print("  " + ", ".join(str(half_cells(n)) for n in range(2, 15)))
    print()


def print_product_lattice_layer() -> None:
    print("=" * 96)
    print("C. Product-lattice packet bridge")
    print("=" * 96)
    x = F(1, 7)
    b2 = 1 - (1 - x) ** 2
    b3 = 1 - (1 - x) ** 3
    print("Independent danger probability x=1/7 through Boolean IE:")
    print(f"  B2 nonempty packet ++-       = 1-(6/7)^2 = {b2}")
    print(f"  B3 nonempty packet +++---+   = 1-(6/7)^3 = {b3}")
    print("This is the same Mobius mechanism as phi=mu*id, but on a Boolean")
    print("far/tiling lattice instead of a divisor lattice.")
    print()
    print("Main-term packet sizes for a 13-speed exact-period witness model:")
    for d in [41, 83, 89, 173, 210, 1260]:
        main = F(phi(d), 1) * F(6, 7) ** 13
        three_packet = F(phi(d), 1) * b3
        print(
            f"  D={d:4d}: phi(D)={phi(d):4d}, "
            f"(6/7)^13*phi={float(main):8.4f}, B3_packet_capacity={float(three_packet):8.3f}"
        )
    print("Proof readout: the denominator divisor label and the far/tiling subset")
    print("label should live on the product lattice Div(D) x B_r until the residual")
    print("error is bounded.  Scalarizing to phi(D) or to a seven-letter sign too")
    print("early loses the multiplicativity defect and the geometry of cancellations.")
    print()


def print_one_tail_denominator_layer() -> None:
    print("=" * 96)
    print("D. The S109 covering one-tail denominators as coprime packets")
    print("=" * 96)
    print("For C={1..11,13}, the q-covering one-tail branch is w=84m and")
    print("M=7m/(84m+5).  The witness denominator D=84m+5 is always coprime")
    print("to the killed primes 2,3,7, so the proof escapes through unit packets.")
    print()
    for m in range(1, 13):
        d = 84 * m + 5
        val = F(7 * m, d)
        margin = val - F(1, 14)
        print(
            f"  m={m:2d}, D={d:4d}={fmt_factor(d):<14s}, "
            f"gcd(D,84)={gcd(d,84)}, phi/D={F(phi(d), d)}, "
            f"M-1/14={margin}"
        )
    print()


def print_tournament_analysis() -> None:
    print("=" * 96)
    print("Tournament Analysis")
    print("=" * 96)
    vertices = [
        "divisor_mobius_phi",
        "exact_period_unit_packets",
        "full_B3_tiling_recursion",
        "even_B2_half_recursion",
        "odd_B3_half_recursion",
        "product_lattice_packet_ledger",
        "scalar_projection_guardrail",
    ]
    print("Vertices are proof carriers, not runners or arcs.")
    print("Pairwise observable: how much primitive packet/address data survives before")
    print("scalarization.  Switch: carriers preserving both divisor and Boolean")
    print("Mobius labels rank above carriers that collapse them.")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("Challenged assumption: the three tiling recurrences are three unrelated")
    print("formulas.  They are all Mobius kernels; the differences come from which")
    print("quotient pushes the kernel to sizes/depths.")


def main() -> None:
    print("codex-S110 totient/coprime density and tiling-recursion packet bridge")
    print_totient_layer()
    print_recursion_layer()
    print_product_lattice_layer()
    print_one_tail_denominator_layer()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
