#!/usr/bin/env python3
"""
lrc_recursive_multisieve_lower_bound_s563.py

codex-2026-06-02 S563

Recursive multi-sieve lower-bound ledger for the Lonely Runner Conjecture.

Benchmark:
  Tao proves a global asymptotic improvement over the first-moment bound
  1/(2k), namely 1/(2k) plus a logarithmic surplus.  Bedert later improves the
  global surplus by Riesz products.  This script is not a new global asymptotic
  proof.  It tests a different route: recursive rational-denominator sieves
  that certify much larger lower bounds on structured cases, including parts
  of the sieve-covered core where the coarse THM-369 sieve is blind.

Tier model:
  T0 small_sieve       : q <= n = k+1.  This is THM-369.
  T1 first_fine_window : n < q <= 2n.  This catches many pinch descendants.
  T2 prime_power_lift  : q is a prime power with base prime <= n.
  T3 crt_smooth_lift   : every prime factor of q is <= n.
  T4 external_fine     : all remaining q <= Q.

For a denominator q, t=a/q is a closed theta-witness when
  min_i ||a v_i / q|| >= theta.

Tournament Analysis declaration:
  Vertices: speed-set proof obligations, not runners.
  Pairwise observable: (best certified margin up to Q, first tier reaching
    the conjectural threshold 1/(k+1), witness denominator).
  Switch/gauge: stronger certified margin wins; ties prefer earlier recursive
    tier and smaller denominator.
  Tie Hamiltonian path: displayed sample order.
  Fingerprints: score histogram, directed 3-cycles, SCCs, and Hamiltonian path
    count are reported.

Assumption challenge:
  Vertices could be runners, denominators, residues, CRT channels, endpoint
  owners, packets, Fourier modes, or proof obligations.  Here vertices are
  whole speed-set obligations because the quotient preserves the lower-bound
  certification problem, but destroys individual interval ownership.  The
  result is a proof-search ledger, not a global theorem beyond Tao/Bedert.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from math import gcd, isqrt, lcm, log
from itertools import combinations


@dataclass(frozen=True)
class Witness:
    margin: F
    a: int
    q: int
    tier: str
    tier_rank: int


@dataclass(frozen=True)
class Sample:
    label: str
    speeds: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class Audit:
    sample: Sample
    k: int
    n: int
    trivial: F
    conjecture: F
    coarse_missed_q: int | None
    first_conjecture: Witness | None
    best: Witness
    tier_bests: tuple[Witness, ...]


def normalize(speeds: tuple[int, ...]) -> tuple[int, ...]:
    cleaned = sorted({abs(s) for s in speeds if s})
    g = 0
    for s in cleaned:
        g = gcd(g, s)
    return tuple(s // g for s in cleaned)


def factor(value: int) -> tuple[tuple[int, int], ...]:
    x = value
    out: list[tuple[int, int]] = []
    p = 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            out.append((p, e))
        p += 1 if p == 2 else 2
    if x > 1:
        out.append((x, 1))
    return tuple(out)


def is_prime_power_with_small_base(q: int, n: int) -> bool:
    f = factor(q)
    return len(f) == 1 and f[0][0] <= n and q > n


def is_smooth_over_n(q: int, n: int) -> bool:
    return q > n and all(p <= n for p, _ in factor(q))


def denominator_tiers(n: int, qmax: int) -> list[tuple[str, list[int]]]:
    seen: set[int] = set()

    def take(name: str, qs: list[int]) -> tuple[str, list[int]]:
        fresh = [q for q in qs if 2 <= q <= qmax and q not in seen]
        seen.update(fresh)
        return name, fresh

    tiers = [
        take("T0_small_sieve", list(range(2, n + 1))),
        take("T1_first_fine_window", list(range(n + 1, 2 * n + 1))),
        take("T2_prime_power_lift", [q for q in range(2, qmax + 1) if is_prime_power_with_small_base(q, n)]),
        take("T3_crt_smooth_lift", [q for q in range(2, qmax + 1) if is_smooth_over_n(q, n)]),
    ]
    tiers.append(take("T4_external_fine", list(range(2, qmax + 1))))
    return tiers


def residue_margin(speeds: tuple[int, ...], a: int, q: int) -> F:
    return min(F(min((a * s) % q, q - ((a * s) % q)), q) for s in speeds)


def best_in_qs(speeds: tuple[int, ...], qs: list[int], tier: str, rank: int) -> Witness | None:
    best: Witness | None = None
    for q in qs:
        for a in range(1, q):
            margin = residue_margin(speeds, a, q)
            cand = Witness(margin=margin, a=a, q=q, tier=tier, tier_rank=rank)
            if best is None or witness_key(cand) > witness_key(best):
                best = cand
    return best


def witness_key(w: Witness) -> tuple[F, int, int]:
    return (w.margin, -w.tier_rank, -w.q)


def coarse_missed_q(speeds: tuple[int, ...], n: int) -> int | None:
    for q in range(2, n + 1):
        if all(s % q for s in speeds):
            return q
    return None


def audit_sample(sample: Sample, qmax: int) -> Audit:
    speeds = normalize(sample.speeds)
    k = len(speeds)
    n = k + 1
    conjecture = F(1, n)
    trivial = F(1, 2 * k)
    tier_bests: list[Witness] = []
    first_conjecture: Witness | None = None
    best = Witness(F(0), 0, 1, "none", 99)
    for rank, (tier, qs) in enumerate(denominator_tiers(n, qmax)):
        tier_best = best_in_qs(speeds, qs, tier, rank)
        if tier_best is None:
            continue
        tier_bests.append(tier_best)
        if witness_key(tier_best) > witness_key(best):
            best = tier_best
        if first_conjecture is None:
            for q in qs:
                for a in range(1, q):
                    margin = residue_margin(speeds, a, q)
                    if margin >= conjecture:
                        first_conjecture = Witness(margin, a, q, tier, rank)
                        break
                if first_conjecture is not None:
                    break
    return Audit(
        sample=Sample(sample.label, speeds, sample.note),
        k=k,
        n=n,
        trivial=trivial,
        conjecture=conjecture,
        coarse_missed_q=coarse_missed_q(speeds, n),
        first_conjecture=first_conjecture,
        best=best,
        tier_bests=tuple(tier_bests),
    )


def fmt(x: F | None) -> str:
    if x is None:
        return "-"
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def lcm_range(lo: int, hi: int) -> int:
    out = 1
    for q in range(lo, hi + 1):
        out = lcm(out, q)
    return out


def packet(n: int, scale: int, skip: int) -> tuple[int, ...]:
    return (1,) + tuple(scale * q for q in range(1, n) if q != skip)


def deterministic_samples() -> list[Sample]:
    n = 14
    blind14 = tuple(range(1, 13)) + (lcm_range(2, 14),)
    blind24 = tuple(range(1, 13)) + (lcm_range(2, 24),)
    return [
        Sample("AP_k13", tuple(range(1, n)), "wall; small sieve reaches q=n exactly"),
        Sample("noncovered_random_k13", (2, 5, 11, 17, 23, 31, 37, 41, 43, 47, 53, 59, 61), "misses q=3, so coarse sieve already beats Tao"),
        Sample("blind_lcm_2_14", blind14, "sieve-covered through q<=14; fine witness should be 2/27"),
        Sample("blind_lcm_2_24", blind24, "sieve-blind far beyond n; same fine witness survives"),
        Sample("S562_packet_n14", packet(14, 7, 6), "residual packet from HYP-2073"),
        Sample("S562_packet_n14_lift", packet(14, 14, 6), "dyadic lift of the n=14 residual packet"),
        Sample("S562_packet_n17", packet(17, 17, 8), "prime n=17 control with residual dyadic denominator"),
        Sample("S562_packet_n18", packet(18, 9, 8), "n=18 square-payload residual packet"),
    ]


def tao_shape_unit(k: int) -> float:
    # Constant-free shape only; Tao's theorem has an unspecified absolute c.
    ll = log(log(k)) if k > 3 else 1.0
    return log(k) / (k * k * ll * ll)


def bedert_shape_unit(k: int) -> float:
    # Constant-free shape from the later Riesz-product improvement.
    return k ** (-5 / 3)


def tournament_fingerprint(audits: list[Audit]) -> dict[str, object]:
    def key(a: Audit) -> tuple[F, int, int]:
        first_rank = a.first_conjecture.tier_rank if a.first_conjecture else 99
        first_q = a.first_conjecture.q if a.first_conjecture else 10**9
        return (a.best.margin, -first_rank, -first_q)

    n = len(audits)
    adj = [[False] * n for _ in range(n)]
    for i, left in enumerate(audits):
        for j, right in enumerate(audits):
            if i != j:
                adj[i][j] = key(left) > key(right) or (key(left) == key(right) and i < j)

    scores = [sum(row) for row in adj]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]):
            c3 += 1

    # Tiny SCC routine.
    def reach(start: int) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v, edge in enumerate(adj[u]):
                if edge and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        u = next(iter(remaining))
        ru = reach(u)
        comp = {v for v in remaining if v in ru and u in reach(v)}
        sccs.append(len(comp))
        remaining -= comp

    return {
        "vertices": [a.sample.label for a in audits],
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": 1 if c3 == 0 and len(set(scores)) == n else "not_counted",
    }


def main() -> None:
    qmax = 220
    print("S563 recursive multi-sieve lower-bound ledger")
    print("=" * 78)
    print("Benchmark frame")
    print("  Tao: global 1/(2k) + c log(k)/(k^2 (log log k)^2), unspecified c.")
    print("  Bedert: later Riesz-product global improvement of polynomial scale.")
    print("  This script: case/recursive certificates, often reaching 1/(k+1).")
    print(f"  denominator cutoff Q={qmax}\n")

    audits = [audit_sample(sample, qmax) for sample in deterministic_samples()]

    print("1. Samples")
    print("  label                  k  n  coarse_q  first>=1/n                 best<=Q                  surplus_vs_1/(2k)")
    for audit in audits:
        first = audit.first_conjecture
        first_s = "-" if first is None else f"{first.tier}:a={first.a},q={first.q},m={fmt(first.margin)}"
        best = audit.best
        surplus = best.margin - audit.trivial
        print(
            f"  {audit.sample.label:<22} {audit.k:2d} {audit.n:2d} "
            f"{str(audit.coarse_missed_q or '-'):>8}  "
            f"{first_s:<28} "
            f"{best.tier}:a={best.a},q={best.q},m={fmt(best.margin):<7} "
            f"{fmt(surplus)}"
        )
        print(f"    note: {audit.sample.note}")

    print("\n2. Tier bests")
    for audit in audits:
        print(f"  {audit.sample.label}:")
        for w in audit.tier_bests:
            print(f"    {w.tier:<22} best a/q={w.a}/{w.q:<3} margin={fmt(w.margin)}")

    print("\n3. Tao/Bedert-scale orientation")
    for k in sorted({a.k for a in audits}):
        trivial = F(1, 2 * k)
        conj = F(1, k + 1)
        print(
            f"  k={k}: 1/(2k)={fmt(trivial)}  1/(k+1)={fmt(conj)}  "
            f"gap_to_conjecture={fmt(conj - trivial)}  "
            f"Tao_shape_unit(c=1)~{tao_shape_unit(k):.6f}  "
            f"Bedert_shape_unit(c=1)~{bedert_shape_unit(k):.6f}"
        )
    print("  The constants in the global theorems are not asserted here; the point is")
    print("  that a recursive sieve certificate at 1/(k+1) is far above the Tao regime")
    print("  on any instance where it applies.")

    print("\n4. Tournament Analysis")
    print("  vertices: speed-set proof obligations")
    print("  observable: best certified margin, first conjecture-reaching tier, denominator")
    print("  switch: larger margin wins; ties prefer earlier tier and smaller q")
    print(f"  fingerprints: {tournament_fingerprint(audits)}")

    print("\n5. Synthesis")
    print("  Coarse THM-369 already proves the full 1/(k+1) bound off the sieve-covered core.")
    print("  The multi-sieve tiers are useful precisely inside that core:")
    print("    blind_lcm_* has no q<=14 witness but gets t=2/27 with margin 2/27.")
    print("    S562 n=14 residual packets get fine-window witnesses such as 6/23.")
    print("    S562 n=17 residual packets get a dyadic-residual witness 9/25.")
    print("  This does not beat Tao as a universal asymptotic theorem; it shrinks the")
    print("  Tao/Bedert hard regime to a recursively generated near-AP residual core.")


if __name__ == "__main__":
    main()
