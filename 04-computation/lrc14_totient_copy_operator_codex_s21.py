#!/usr/bin/env python3
"""
LRC(14) totient-copy operator and squarefree divisor profiles.

User reframe:
  Create c(n) copies of each natural number n, with c(n) >= 1, so that

      sum_{d|n} c(d) = n.

The unique solution is Euler's totient phi(n).  Equivalently, on the q-grid
{0/q,1/q,...,(q-1)/q}, the number of points whose reduced denominator is exactly
d is phi(d), and these exact-period packets partition q.

This script asks what that does to the active LRC14 squarefree-profile route:
  * q=14: phi(14)=6 is the unit seam (Z/14Z)^* -> F_7^*.
  * q=210: the squarefree profile {2,3,5,7} has copy atoms prod(p-1).
  * q=1260: the raw K_14 Hill product keeps prime powers 2^2 and 3^2,
    so selected-prime mask weights become p^a-1 before projection to rad.
  * q=2520: the THM-523 half-clash denominator has phi(2520)=576, equal
    to the full-mask copy mass of q=1260.

Tournament Analysis declaration:
  Vertices are proof quotients, not raw runners:
    totient_copy_operator, exact_period_q_grid, raw_1260_profile,
    squarefree_210_profile, unit_seam_phi14, coimage_character_tail,
    raw_divisor_counts, raw_runner_vertices.
  Pairwise observable: preservation of the divisor-copy identity, the LRC14
    exact-period witness predicate, squarefree transfer value, and compatibility
    with HYP-2626's mod-7 coimage seam.
  Tie path: the listed order.

Assumption challenge:
  Candidate vertices considered: runners, speeds, q-grid residues, reduced
  denominators, divisor-lattice nodes, prime masks, exact-period packets, unit
  groups, coimage classes, and proof obligations.  The quotient preserved here
  is "which exact-period packets can carry witness/coimage mass"; it destroys
  the order of residues around the circle and individual witness times.
"""
from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from operator import mul

PRIMES = (2, 3, 5, 7)


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def prod(vals) -> int:
    return reduce(mul, vals, 1)


def divisors(n: int) -> list[int]:
    out = []
    for d in range(1, int(math.isqrt(n)) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


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


def fmt_factor(n: int) -> str:
    fs = factor(n)
    if not fs:
        return "1"
    return "*".join(str(p) if e == 1 else f"{p}^{e}" for p, e in fs.items())


def phi(n: int) -> int:
    if n == 1:
        return 1
    out = n
    for p in factor(n):
        out = out // p * (p - 1)
    return out


def mobius(n: int) -> int:
    fs = factor(n)
    if any(e > 1 for e in fs.values()):
        return 0
    return -1 if len(fs) % 2 else 1


def copy_by_recursion(nmax: int) -> dict[int, int]:
    c: dict[int, int] = {}
    for n in range(1, nmax + 1):
        c[n] = n - sum(c[d] for d in divisors(n) if d < n)
    return c


def mask_of(n: int, primes: tuple[int, ...] = PRIMES) -> int:
    m = 0
    for i, p in enumerate(primes):
        if n % p == 0:
            m |= 1 << i
    return m


def mask_name(mask: int, primes: tuple[int, ...] = PRIMES) -> str:
    vals = [str(p) for i, p in enumerate(primes) if mask & (1 << i)]
    return "{" + ",".join(vals) + "}" if vals else "{}"


def radical(n: int) -> int:
    return prod(factor(n).keys()) if n > 1 else 1


def exact_period_counts(q: int) -> Counter[int]:
    """Counts reduced denominators among residues a/q, 0 <= a < q."""
    counts: Counter[int] = Counter()
    for a in range(q):
        d = q // math.gcd(a, q)
        counts[d] += 1
    return counts


def mask_profile(q: int, primes: tuple[int, ...] = PRIMES) -> Counter[int]:
    """Compress exact-period counts by the prime mask of the reduced denominator."""
    profile: Counter[int] = Counter()
    for d in divisors(q):
        profile[mask_of(d, primes)] += phi(d)
    return profile


def theoretical_mask_profile(q: int, primes: tuple[int, ...] = PRIMES) -> Counter[int]:
    """Product formula: selecting p contributes p^a-1 across all positive exponents."""
    fs = factor(q)
    profile: Counter[int] = Counter({0: 1})
    for i, p in enumerate(primes):
        if p not in fs:
            continue
        selected_weight = p ** fs[p] - 1
        nxt: Counter[int] = Counter()
        bit = 1 << i
        for m, weight in profile.items():
            nxt[m] += weight
            nxt[m | bit] += weight * selected_weight
        profile = nxt
    return profile


def profile_rows(profile: Counter[int], primes: tuple[int, ...] = PRIMES) -> list[tuple[str, int]]:
    return [(mask_name(m, primes), profile[m]) for m in sorted(profile)]


def verify_totient_identity() -> None:
    section("COPY RULE IS EULER TOTIENT")
    c = copy_by_recursion(80)
    bad_phi = [n for n in range(1, 81) if c[n] != phi(n)]
    bad_sum = [n for n in range(1, 81) if sum(c[d] for d in divisors(n)) != n]
    print("recursive copy rule: c(n)=n-sum_{d|n,d<n} c(d)")
    print("Möbius inversion: c = mu * id, hence c(n)=phi(n)")
    print(f"checked n<=80 against phi: mismatches={bad_phi}")
    print(f"checked divisor sums n<=80: mismatches={bad_sum}")
    print()
    print(f"{'n':>3} {'factor':>12} {'c(n)':>5} {'sum_{d|n}c(d)':>14} {'mu*id':>6}")
    for n in range(1, 25):
        mob = sum(mobius(d) * (n // d) for d in divisors(n))
        print(f"{n:>3} {fmt_factor(n):>12} {c[n]:>5} {sum(c[d] for d in divisors(n)):>14} {mob:>6}")


def exact_period_report() -> None:
    section("EXACT-PERIOD Q-GRID READOUT")
    print(
        "Residues a/q are partitioned by reduced denominator d|q. "
        "The packet size is phi(d), so the copy rule is literally the q-grid exact-period census."
    )
    for q in (14, 30, 210):
        counts = exact_period_counts(q)
        ok = all(counts[d] == phi(d) for d in divisors(q))
        print()
        print(f"q={q}={fmt_factor(q)}  exact-period packets match phi(d): {ok}")
        print("  " + ", ".join(f"d={d}: {counts[d]}" for d in divisors(q)))


def profile_report() -> None:
    section("SQUAREFREE AND RAW DENOMINATOR COPY PROFILES")
    qs = (14, 210, 315, 1260, 2520, 36, 70, 504)
    print(
        f"{'q':>5} {'factor':>18} {'rad':>5} {'phi(q)':>7} {'phi/q':>9} "
        f"{'full-mask mass':>14} {'profile sum':>11}"
    )
    for q in qs:
        profile = mask_profile(q)
        theory = theoretical_mask_profile(q)
        if profile != theory:
            raise AssertionError((q, profile, theory))
        full_mask = mask_of(radical(q))
        print(
            f"{q:>5} {fmt_factor(q):>18} {radical(q):>5} {phi(q):>7} "
            f"{str(Fraction(phi(q), q)):>9} {profile[full_mask]:>14} {sum(profile.values()):>11}"
        )

    print()
    for q in (210, 1260):
        print(f"mask profile for q={q}={fmt_factor(q)}:")
        rows = profile_rows(mask_profile(q))
        line = []
        for name, weight in rows:
            line.append(f"{name}:{weight}")
        for i in range(0, len(line), 4):
            print("  " + "  ".join(line[i : i + 4]))
        print()

    p210 = mask_profile(210)
    p1260 = mask_profile(1260)
    print("raw-vs-squarefree amplification from 210 to 1260:")
    print("  repeated 2 contributes factor (2^2-1)/(2-1)=3")
    print("  repeated 3 contributes factor (3^2-1)/(3-1)=4")
    print(f"{'mask':>12} {'mass210':>8} {'mass1260':>9} {'ratio':>8}")
    for m in sorted(p210):
        if p210[m] == 0:
            continue
        print(f"{mask_name(m):>12} {p210[m]:>8} {p1260[m]:>9} {str(Fraction(p1260[m], p210[m])):>8}")


def lrc14_telescope_report() -> None:
    section("LRC14 1260/2520 TOTIENT READOUT")
    half = Fraction(15, 36) - Fraction(2, 5) - Fraction(1, 70) - Fraction(1, 504)
    print(f"THM-523 half-clash = {half}; denominator = {half.denominator}={fmt_factor(half.denominator)}")
    print(f"symmetry doubles to {2 * half}; denominator = {(2 * half).denominator}={fmt_factor((2 * half).denominator)}")
    print()
    print("Totient-copy resonance:")
    print(f"  phi(14)   = {phi(14)}   -> unit seam size (Z/14Z)^* and F_7^*")
    print(f"  phi(210)  = {phi(210)}  -> full squarefree mod-210 exact-period copies")
    print(f"  phi(1260) = {phi(1260)} -> exact-period copies at raw K14 product")
    print(f"  phi(2520) = {phi(2520)} -> exact-period copies at half-clash denominator")
    print(
        "  full-mask mass of q=1260 = "
        f"{mask_profile(1260)[mask_of(210)]}, equal to phi(2520)."
    )
    print()
    print("Denominators in the half-clash, with exact-period copy densities:")
    for q in (36, 5, 70, 504, 2520, 1260):
        print(
            f"  q={q:>4}={fmt_factor(q):>12}: phi={phi(q):>4}, "
            f"phi/q={Fraction(phi(q), q)}, rad={radical(q)}"
        )
    print()
    print(
        "Readout: the dyadic doubling from 1260 to 2520 is not cosmetic.  It turns "
        "the exact-period unit count into 576, the same mass obtained by selecting "
        "all four prime masks inside the raw 1260 profile.  That is the arithmetic "
        "shadow of the half-gap before tau<->1-tau symmetry halves the denominator."
    )


@dataclass(frozen=True)
class Tournament:
    vertices: tuple[str, ...]
    edges: dict[tuple[str, str], str]


def make_order_tournament(vertices: tuple[str, ...]) -> Tournament:
    edges: dict[tuple[str, str], str] = {}
    rank = {v: i for i, v in enumerate(vertices)}
    for a, b in itertools.combinations(vertices, 2):
        winner = a if rank[a] < rank[b] else b
        edges[(a, b)] = winner
    return Tournament(vertices, edges)


def tournament_fingerprints(t: Tournament) -> dict[str, object]:
    scores = Counter({v: 0 for v in t.vertices})
    for (a, b), winner in t.edges.items():
        scores[winner] += 1
    cycles = 0
    for tri in itertools.combinations(t.vertices, 3):
        wins = {v: 0 for v in tri}
        for a, b in itertools.combinations(tri, 2):
            winner = t.edges[(a, b)] if (a, b) in t.edges else t.edges[(b, a)]
            wins[winner] += 1
        if sorted(wins.values()) == [1, 1, 1]:
            cycles += 1
    hist = Counter(scores.values())
    return {
        "score_hist": dict(sorted(hist.items())),
        "directed_3_cycles": cycles,
        "hamiltonian_path_count": 1 if cycles == 0 and len(set(scores.values())) == len(scores) else None,
    }


def mask_mass_tournament(q: int) -> dict[str, object]:
    profile = mask_profile(q)
    masks = tuple(sorted(profile))
    # Gauge: larger exact-period mass wins; ties go to the fixed mask order.
    edges: dict[tuple[int, int], int] = {}
    for a, b in itertools.combinations(masks, 2):
        if profile[a] > profile[b]:
            winner = a
        elif profile[b] > profile[a]:
            winner = b
        else:
            winner = a if a < b else b
        edges[(a, b)] = winner
    scores = Counter({m: 0 for m in masks})
    for winner in edges.values():
        scores[winner] += 1
    cycles = 0
    for tri in itertools.combinations(masks, 3):
        wins = {m: 0 for m in tri}
        for a, b in itertools.combinations(tri, 2):
            winner = edges[(a, b)] if (a, b) in edges else edges[(b, a)]
            wins[winner] += 1
        if sorted(wins.values()) == [1, 1, 1]:
            cycles += 1
    top = sorted(masks, key=lambda m: (-scores[m], -profile[m], m))[:5]
    return {
        "q": q,
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "directed_3_cycles": cycles,
        "top_masks": [(mask_name(m), profile[m], scores[m]) for m in top],
    }


def tournament_report() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = (
        "totient_copy_operator",
        "exact_period_q_grid",
        "raw_1260_profile",
        "squarefree_210_profile",
        "unit_seam_phi14",
        "coimage_character_tail",
        "raw_divisor_counts",
        "raw_runner_vertices",
    )
    t = make_order_tournament(vertices)
    fp = tournament_fingerprints(t)
    print("proof-quotient Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("observable: preserve divisor-copy identity, exact-period LRC witness packets, squarefree transfer, coimage seam")
    print(f"fingerprints: {fp}")
    print()
    print("mask-mass tournaments use exact-period mass as the pairwise observable; ties follow mask order.")
    for q in (210, 1260):
        print(f"q={q}: {mask_mass_tournament(q)}")


def main() -> None:
    verify_totient_identity()
    exact_period_report()
    profile_report()
    lrc14_telescope_report()
    tournament_report()


if __name__ == "__main__":
    main()
