#!/usr/bin/env python3
"""HYP-3135 scout: the De Moivre resolvent numbers 120/320 as packet payloads.

This is not a new LRC proof.  It records a narrow arithmetic check and uses it
as a controlled-forgetting guardrail for the current LRC14 endgame:

* the example quintic's resolvent quartic has roots 2,-4,8,-16;
* the coefficients 120 and 320 are elementary-symmetric middle layers, not
  roots themselves;
* the live LRC14 proof has the same shape: the useful invariant is a packet of
  low/middle interaction data (Q-floor, signed SPEC low/tail, Lee-Yang radius,
  edge tail/tip sectors, global-consistency quotient, assembled multi-far
  factorization), not any single scalar shadow.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import comb, factorial


def elem_sym(values: list[int], k: int) -> int:
    total = 0
    for combo in combinations(values, k):
        prod = 1
        for value in combo:
            prod *= value
        total += prod
    return total


def real_fifth_root(value: int) -> float:
    sign = -1.0 if value < 0 else 1.0
    return sign * (abs(value) ** (1.0 / 5.0))


def quintic(x: float) -> float:
    return x**5 + 20 * x**3 + 20 * x**2 + 30 * x + 10


def print_resolvent_check() -> None:
    roots = [2, -4, 8, -16]
    e1 = elem_sym(roots, 1)
    e2 = elem_sym(roots, 2)
    e3 = elem_sym(roots, 3)
    e4 = elem_sym(roots, 4)
    coeffs = [1, -e1, e2, -e3, e4]

    print("DE MOIVRE RESOLVENT CHECK")
    print(f"resolvent roots: {roots}")
    print(f"elementary symmetric layers: e1={e1}, e2={e2}, e3={e3}, e4={e4}")
    print(f"quartic coefficients from roots: {coeffs}")
    print("quartic: x^4 + 10*x^3 - 120*x^2 - 320*x + 1024")
    print("reading: 120 = -e2 is the pair layer; 320 = e3 is the triple layer.")

    x = sum(real_fifth_root(root) for root in roots)
    print()
    print("REAL QUINTIC ROOT CHECK")
    print("branch expression: fifthroot(2) - fifthroot(4) + fifthroot(8) - fifthroot(16)")
    print(f"x ~= {x:.15f}")
    print(f"x^5 + 20*x^3 + 20*x^2 + 30*x + 10 ~= {quintic(x):.3e}")


def print_a000568_sandwich() -> None:
    # A000568 values in the apex-7 window used by incoming HYP-3131.
    a000568 = {4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
    print()
    print("A000568 SANDWICH WINDOW FROM HYP-3131 / HYP-3134")
    print(" n   C(n,3)   A000568(n)   2*(n-1)!/3   status")
    for n, value in a000568.items():
        lower = comb(n, 3)
        upper = Fraction(2 * factorial(n - 1), 3)
        status = "inside" if lower <= value <= upper else "breaks"
        print(f"{n:2d} {lower:8d} {value:12d} {str(upper):>13s}   {status}")


def print_repo_motifs() -> None:
    motifs = [
        (
            "pair/triple cycle correction",
            "05-knowledge/variables/cycle-counts.md",
            "c_2 contains -120*bc33: a middle interaction correction, not a raw count.",
        ),
        (
            "bounded resonant-center audit",
            "00-navigation/SESSION-LOG.md",
            "HYP-2852 recorded a provable 0/320 bounded-V resonant-center floor check.",
        ),
        (
            "support-six signed cancellation",
            "00-navigation/SESSION-LOG.md",
            "HYP-2608 used 320 random six-support vectors to expose residue-addressed reciprocal sums.",
        ),
        (
            "cap/slack signal",
            "00-navigation/SESSION-LOG.md",
            "k=13 cap rows reported +0.320 slack: margin, not theorem by itself.",
        ),
        (
            "tournament/fixed-path scale",
            "07-reflections/the-blue-line-morse-theory.md",
            "120 appears as 5! labelled tournament mass and fixed-path/profile scale.",
        ),
        (
            "torsion moat",
            "07-reflections/lonely-runner-recursive-metrics-s376.md",
            "n=15 order-3 puncture moat has missed-cell count 120.",
        ),
        (
            "modular/support horizon",
            "07-reflections/eta-powers-type-II-and-matroid-support-gates.md",
            "q^120 is a finite support horizon where sparsity becomes modular control.",
        ),
    ]
    print()
    print("SELECTED 120/320 MOTIFS ALREADY IN THE REPO")
    for label, path, reading in motifs:
        print(f"- {label}: {path}")
        print(f"  {reading}")


def print_lrc_packet_map() -> None:
    print()
    print("CURRENT LRC14 REMAINING-PROOF MAP")
    rows = [
        (
            "Q/apex block",
            "closed as a positive factor",
            "HYP-3130 exact c_r floors; HYP-3128 Lee-Yang Q-block; HYP-3129 union/SPEC agreement.",
        ),
        (
            "far additions",
            "help rather than obstruct",
            "HYP-3131: far elements push miss-PGF zeros outward; multi-far reduces to bounded core.",
        ),
        (
            "global quotient discipline",
            "A000568 wedge is a proof-legality signal",
            "HYP-3134: A000568 sits between raw edge sectors and paired tail/tip child packets as a global-consistency quotient.",
        ),
        (
            "assembled multi-far floor",
            "factorized into three named proof fields",
            "HYP-3136: L(S)=Rprime*meas(R-safe)*meas(Q-lonely), with Q/apex minorant, R-safe wide-V/bounded-core floor, and signed SPEC coupling.",
        ),
        (
            "Rprime coupling",
            "certified but not fully theorem-packaged",
            "HYP-3129/HYP-3136: signed SPEC gives Rprime >= 0.64178 on the tested family; needs closed-form constant chase.",
        ),
        (
            "absolute minorant envelope",
            "ruled out",
            "HYP-3130: B1/h0 > 1, so the proof must use sign cancellation, not absolute Schur mass.",
        ),
        (
            "bounded core",
            "main proof bottleneck",
            "HYP-3122/HYP-3132(S70)/HYP-3131: the remaining core is the k=8 phi4 dip, now a solvable biquadratic resolvent plus the rho>1 => Rprime>=c bridge.",
        ),
        (
            "packet legality",
            "still required",
            "HYP-3116/HYP-3124/HYP-3125: retain edge tail/tip sectors, finite address, observer-gluing, Phi/P, or name debt.",
        ),
    ]
    for carrier, status, note in rows:
        print(f"- {carrier}: {status}")
        print(f"  {note}")

    print()
    print("SYNTHESIS")
    print(
        "The example quartic teaches the useful lesson: the branch roots are simple, "
        "but the load-bearing data live in the middle elementary-symmetric payload "
        "(-120, 320).  The LRC14 analogue is to stop promoting raw p0, raw root "
        "radius, or raw tournament counts alone.  The current proof packet should "
        "carry Q-floor constants, signed SPEC low/tail, Lee-Yang radius/far-push "
        "status, edge tail/tip deletion sectors, a named global-consistency "
        "quotient, and the HYP-3136 assembled multi-far factorization.  The "
        "remaining work is a closed-form middle-layer packet "
        "theorem, sharpened by S70 to a k=8 biquadratic/phi4 coefficient bound."
    )


def main() -> None:
    print_resolvent_check()
    print_a000568_sandwich()
    print_repo_motifs()
    print_lrc_packet_map()


if __name__ == "__main__":
    main()
