#!/usr/bin/env python3
"""
lrc14_ap_boundary_core_certificate_codex_s90.py

Node-1 follow-up: attack the coordinated-growth AP boundary core

    S(t,V) = {t,2t,...,12t,V}.

The prompt suggested that the finite-V Part-A cash-out might require a
Diophantine condition such as V/t -> infinity on this core.  This script checks
the opposite proof route: S(t,V) has an explicit one-point certificate, so this
AP family is not a genuine slow-fast obstruction.

After dividing by gcd(t,V), either:

  * A denominator 13*t residue certificate gives tau=a/(13*t).  Then
    j*t*tau = j*a/13 is 1/13 away from Z for all j=1..12 because a is not
    divisible by 13; the residue a*V mod 13t is chosen away from 0 by at
    least 1/14.

  * In the exceptional reduced case t=1 and V=13m, the tail-law certificate
    tau=m/(13m+1) gives margin m/(13m+1), with equality 1/14 only at m=1.

The residue existence is just counting.  For gcd(t,V)=1 and 13 not dividing V,
multiplication by V permutes Z/(13t); the middle residues have cardinality
13t - 2*ceil(13t/14) + 1 > t, so not all are images of a divisible by 13.
If 13 divides V and t>1, the same search lives on the image 13*Z/tZ, and the
middle residues modulo t are nonempty; a representative can avoid divisibility
by 13.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import ceil, gcd


HALF14 = Fraction(1, 14)


def fmt(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def dist_mod(r: int, n: int) -> int:
    r %= n
    return min(r, n - r)


def norm_frac(x: Fraction) -> Fraction:
    r = x % 1
    return min(r, 1 - r)


@dataclass(frozen=True)
class Cert:
    t: int
    V: int
    g: int
    t0: int
    V0: int
    tau_reduced: Fraction
    tau_original: Fraction
    case: str
    margin: Fraction
    ap_margin: Fraction
    V_margin: Fraction


def margin_for_reduced(t: int, V: int, tau: Fraction) -> tuple[Fraction, Fraction, Fraction]:
    ap = min(norm_frac(j * t * tau) for j in range(1, 13))
    vm = norm_frac(V * tau)
    return min(ap, vm), ap, vm


def denominator_13t_certificate(t: int, V: int) -> int | None:
    """Find a modulo 13t with AP margin 1/13 and V-margin at least 1/14."""
    n = 13 * t
    for a in range(1, n):
        if a % 13 == 0:
            continue
        if 14 * dist_mod(a * V, n) >= n:
            return a
    return None


def certificate(t: int, V: int) -> Cert:
    if t <= 0 or V <= 0:
        raise ValueError("t,V must be positive")
    g = gcd(t, V)
    t0, V0 = t // g, V // g

    a = denominator_13t_certificate(t0, V0)
    if a is not None:
        tau0 = Fraction(a, 13 * t0)
        margin, ap, vm = margin_for_reduced(t0, V0, tau0)
        return Cert(t, V, g, t0, V0, tau0, tau0 / g, "denominator_13t", margin, ap, vm)

    # The only expected no-13t-certificate case after gcd reduction is t0=1,
    # V0=13m.  Then {1,...,12,13m} has the explicit tail-law witness below.
    if t0 == 1 and V0 % 13 == 0:
        m = V0 // 13
        tau0 = Fraction(m, 13 * m + 1)
        margin, ap, vm = margin_for_reduced(t0, V0, tau0)
        return Cert(t, V, g, t0, V0, tau0, tau0 / g, "tail_13m", margin, ap, vm)

    raise AssertionError(f"unexpected uncovered case t={t} V={V} reduced=({t0},{V0})")


def verify_box(t_max: int, ratio_max: int) -> list[Cert]:
    rows: list[Cert] = []
    failures: list[tuple[int, int, Cert]] = []
    for t in range(1, t_max + 1):
        for V in range(1, ratio_max * t + 1):
            if V in {j * t for j in range(1, 13)}:
                continue
            cert = certificate(t, V)
            rows.append(cert)
            if cert.margin < HALF14:
                failures.append((t, V, cert))
    if failures:
        for t, V, cert in failures[:10]:
            print(f"FAIL t={t} V={V} case={cert.case} margin={fmt(cert.margin)} tau={fmt(cert.tau_original)}")
        raise SystemExit(1)
    return rows


def case_stats(rows: list[Cert]) -> None:
    counts: dict[str, int] = {}
    mins: dict[str, Cert] = {}
    for rec in rows:
        counts[rec.case] = counts.get(rec.case, 0) + 1
        if rec.case not in mins or rec.margin < mins[rec.case].margin:
            mins[rec.case] = rec
    print("Case counts and minimum margins")
    for case in sorted(counts):
        rec = mins[case]
        print(
            f"  {case:17s}: count={counts[case]:6d} min_margin={fmt(rec.margin):>12s} "
            f"at t={rec.t}, V={rec.V}, reduced=({rec.t0},{rec.V0}), "
            f"tau={fmt(rec.tau_original)}, ap={fmt(rec.ap_margin)}, V={fmt(rec.V_margin)}"
        )
    print()


def sample_rows() -> None:
    examples = [
        (1, 13),
        (1, 26),
        (1, 182),
        (7, 92),
        (20, 261),
        (20, 2007),
        (37, 13 * 37 + 1),
        (37, 78 * 37 + 5),
        (84, 13 * 97),
    ]
    print("Representative certificates")
    for t, V in examples:
        if V in {j * t for j in range(1, 13)}:
            continue
        rec = certificate(t, V)
        print(
            f"  S={{t..12t,V}} t={t:3d} V={V:5d}: case={rec.case:17s} "
            f"tau={fmt(rec.tau_original):>12s} margin={fmt(rec.margin):>10s} "
            f"(AP={fmt(rec.ap_margin)}, V={fmt(rec.V_margin)})"
        )
    print()


def proof_ledger() -> None:
    print("Proof ledger")
    print("  1. Divide by g=gcd(t,V).  Scale invariance sends a reduced witness tau0")
    print("     to tau=tau0/g for the original set.")
    print("  2. If a is not divisible by 13, tau0=a/(13t) makes the AP block")
    print("     {t,2t,...,12t} exactly 1/13-safe.")
    print("  3. It remains to choose a with ||aV/(13t)|| >= 1/14.")
    print("     When gcd(t,V)=1 and 13∤V, multiplication by V permutes Z/(13t),")
    print("     and the middle residues have size 13t-2ceil(13t/14)+1 > t,")
    print("     so one middle residue has preimage a not divisible by 13.")
    print("  4. When gcd(t,V)=1, 13|V, and t>1, the image is 13*Z/tZ; choose")
    print("     a middle residue modulo t and then a representative a mod 13t")
    print("     not divisible by 13.")
    print("  5. The exceptional reduced family {1,...,12,13m} has")
    print("     tau=m/(13m+1), giving every runner distance m/(13m+1) >= 1/14.")
    print("     Equality occurs only at m=1, the tight AP {1,...,13}.")
    print()


def tournament_analysis() -> None:
    print("Tournament Analysis")
    vertices = [
        ("mod_13t_residue_certificate", 7),
        ("tail_13m_exact_witness", 6),
        ("scale_invariance_gcd_reduction", 5),
        ("slow_fast_arc_error_budget", 3),
        ("raw_V_over_t_growth", 2),
        ("runner_vertices", 1),
    ]
    hist: dict[int, int] = {}
    for _name, score in vertices:
        hist[score] = hist.get(score, 0) + 1
    path = " > ".join(name for name, _score in sorted(vertices, key=lambda row: -row[1]))
    print("  pairwise observable: explicit certified margin for S(t,V)")
    print("  switch/gauge: a vertex wins if it proves the AP boundary core uniformly")
    print(f"  score_hist={hist}")
    print(f"  Hamiltonian path: {path}")
    print("  challenged assumption: the boundary AP core needs V/t->infinity.")
    print("  The residue certificate uses denominator 13t and works at finite ratio;")
    print("  only the non-AP wide clusters still need the slow-fast/decorrelation node.")


def main() -> None:
    print("LRC14 AP boundary-core certificate scout (Codex S90)")
    print("Family S(t,V)={t,2t,...,12t,V}; exact rational certificates")
    print()
    sample_rows()
    rows = verify_box(t_max=80, ratio_max=120)
    print(f"Verified certificate box: t<=80, V<=120*t, rows={len(rows)}, failures=0")
    print()
    case_stats(rows)
    proof_ledger()
    tournament_analysis()


if __name__ == "__main__":
    main()
