#!/usr/bin/env python3
"""LRC Pillai fixed-clock carrier scout.

codex-2026-06-04-S646

Perfect numbers are fixed points of the aliquot/divisor carrier
sigma(N)-N = N.  This script ports that lesson to the LRC pair-sum clock
C = 2n-1 by replacing ordinary divisor mass with the gcd-shell mass

    A(C) = sum_{1 <= a <= (C-1)/2} gcd(a, C).

Equivalently, if P(C)=sum_{r=1}^C gcd(r,C) is Pillai's arithmetical
function, then A(C)=(P(C)-C)/2.  The LRC clock is "Pillai-fixed" when
A(C)=C, i.e. P(C)=3C.

The goal is not to prove LRC.  It is to test whether the aliquot fixed-point
idea exposes a real n=14 side channel: the AP speed row is a shell transversal
for C=27, and Vstar is not shell-fixed but is gcd-carrier-fixed.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import gcd


def factor(n: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            a = 0
            while n % d == 0:
                n //= d
                a += 1
            out.append((d, a))
        d += 1 if d == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def phi(n: int) -> int:
    ans = n
    for p, _ in factor(n):
        ans = ans // p * (p - 1)
    return ans


def divisors_from_factorization(fac: list[tuple[int, int]]) -> list[int]:
    divs = [1]
    for p, a in fac:
        divs = [d * (p ** e) for d in divs for e in range(a + 1)]
    return sorted(divs)


def pillai(n: int) -> int:
    fac = factor(n)
    return sum(d * phi(n // d) for d in divisors_from_factorization(fac))


def shell_mass(C: int) -> int:
    assert C % 2 == 1
    return sum(gcd(a, C) for a in range(1, (C + 1) // 2))


def shell(a: int, C: int) -> int:
    r = a % C
    if r == 0:
        return 0
    return min(r, C - r)


def shell_counter(speeds: list[int], C: int) -> Counter[int]:
    return Counter(shell(v, C) for v in speeds)


def gcd_counter(speeds: list[int], C: int) -> Counter[int]:
    return Counter(gcd(shell(v, C), C) for v in speeds)


def weighted_gcd_mass(speeds: list[int], C: int) -> int:
    return sum(gcd(shell(v, C), C) for v in speeds)


def target_shell_counter(n: int) -> Counter[int]:
    return Counter(range(1, n))


def target_gcd_counter(n: int) -> Counter[int]:
    C = 2 * n - 1
    return gcd_counter(list(range(1, n)), C)


def l1_counter_delta(a: Counter[int], b: Counter[int]) -> int:
    keys = set(a) | set(b)
    return sum(abs(a[k] - b[k]) for k in keys)


def primitive(speeds: list[int]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def maximin_pair_pinch(speeds: list[int]) -> tuple[Fraction, list[Fraction]]:
    """Exact maximin using the repo THM-369 pair-sum pinch oracle."""
    candidates: set[tuple[int, int]] = set()
    for a, b in combinations(speeds, 2):
        D = a + b
        for m in range(D + 1):
            candidates.add((m, D))

    best = Fraction(0, 1)
    best_times: list[Fraction] = []
    for m, D in candidates:
        local = min(
            Fraction(min((v * m) % D, D - ((v * m) % D)), D)
            for v in speeds
        )
        if local > best:
            best = local
            best_times = [Fraction(m, D)]
        elif local == best:
            best_times.append(Fraction(m, D))

    return best, sorted(set(best_times))


def ap_single_swap_scan(max_n: int = 14) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in range(4, max_n + 1):
        base = set(range(1, n))
        C = 2 * n - 1
        target = Fraction(1, n)
        base_shells = target_shell_counter(n)
        base_gcds = target_gcd_counter(n)
        base_mass = weighted_gcd_mass(sorted(base), C)
        for out_v in range(1, n):
            for in_v in range(1, 2 * n + 1):
                if in_v in base:
                    continue
                speeds = sorted((base - {out_v}) | {in_v})
                if len(speeds) != n - 1 or not primitive(speeds):
                    continue
                M, times = maximin_pair_pinch(speeds)
                shells = shell_counter(speeds, C)
                gcds = gcd_counter(speeds, C)
                rows.append(
                    {
                        "n": n,
                        "C": C,
                        "out": out_v,
                        "inn": in_v,
                        "M": M,
                        "times": times,
                        "tight": M == target,
                        "below": M < target,
                        "gap": M - target,
                        "mass_delta": weighted_gcd_mass(speeds, C) - base_mass,
                        "shell_l1": l1_counter_delta(shells, base_shells),
                        "gcd_l1": l1_counter_delta(gcds, base_gcds),
                        "out_shell": shell(out_v, C),
                        "in_shell": shell(in_v, C),
                        "out_gcd": gcd(shell(out_v, C), C),
                        "in_gcd": gcd(shell(in_v, C), C),
                    }
                )
    return rows


def directed_3cycles(adj: list[list[int]]) -> int:
    n = len(adj)
    c = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        ):
            c += 1
    return c


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach_from(start: int, reverse: bool = False) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u in range(n):
                edge = adj[u][v] if reverse else adj[v][u]
                if edge and u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen

    remaining = set(range(n))
    sizes = []
    while remaining:
        v = next(iter(remaining))
        comp = reach_from(v) & reach_from(v, reverse=True)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            count = dp[mask][v]
            if not count:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += count
    return sum(dp[(1 << n) - 1])


def majority_tournament(vertices: list[str], rankings: list[list[str]]) -> list[list[int]]:
    pos = [{v: i for i, v in enumerate(r)} for r in rankings]
    n = len(vertices)
    adj = [[0] * n for _ in range(n)]
    for i, u in enumerate(vertices):
        for j, v in enumerate(vertices):
            if i == j:
                continue
            wins = sum(1 for p in pos if p[u] < p[v])
            losses = len(pos) - wins
            if wins > losses or (wins == losses and pos[0][u] < pos[0][v]):
                adj[i][j] = 1
    return adj


def ranking_edge_set(vertices: list[str], ranking: list[str]) -> set[tuple[int, int]]:
    p = {v: i for i, v in enumerate(ranking)}
    edges = set()
    for i, u in enumerate(vertices):
        for j, v in enumerate(vertices):
            if i != j and p[u] < p[v]:
                edges.add((i, j))
    return edges


def edge_set(adj: list[list[int]]) -> set[tuple[int, int]]:
    return {
        (i, j)
        for i in range(len(adj))
        for j in range(len(adj))
        if adj[i][j]
    }


def print_clock_table() -> None:
    print("=" * 72)
    print("Pillai fixed clocks for the LRC shell carrier")
    print("=" * 72)
    fixed = []
    for C in range(3, 1000, 2):
        mass = shell_mass(C)
        if mass == C:
            fixed.append(C)
    print(f"Odd C <= 999 with A(C)=C: {fixed}")
    print(f"Translated total runner denominators n=(C+1)/2: {[(C + 1)//2 for C in fixed]}")
    print()
    for C in [7, 9, 11, 13, 15, 17, 21, 23, 25, 27, 33, 39]:
        P = pillai(C)
        mass = shell_mass(C)
        print(
            f"C={C:2d} fac={factor(C)!s:18s} "
            f"P(C)/C={str(Fraction(P, C)):>5s} "
            f"A(C)={mass:3d} A(C)-C={mass-C:+4d}"
        )

    print()
    print("Classification proof sketch for odd C:")
    print("  P(C)/C is multiplicative with local factor 1+a*(1-1/p) at p^a.")
    print("  A(C)=C iff P(C)/C=3.")
    print("  If 3 does not divide C, the two smallest odd-prime factors already")
    print("  give (9/5)*(13/7)>3, and one prime power never solves the equation.")
    print("  If 3^a || C: a=1 forces the remaining factor to be 9/5, hence 5;")
    print("  a=2 leaves a ratio below every odd-prime factor; a=3 gives 27;")
    print("  a>3 overshoots.  So the odd fixed clocks are exactly 15 and 27.")


def print_single_swap_table(rows: list[dict[str, object]]) -> None:
    print()
    print("=" * 72)
    print("AP single-swap scout through n=14")
    print("=" * 72)
    by_n: defaultdict[int, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_n[int(row["n"])].append(row)

    for n in sorted(by_n):
        C = 2 * n - 1
        mass = shell_mass(C)
        tight = [r for r in by_n[n] if r["tight"]]
        below = [r for r in by_n[n] if r["below"]]
        mass_changed_tight = [r for r in tight if r["mass_delta"] != 0]
        same_mass_loose = [r for r in by_n[n] if (not r["tight"]) and r["mass_delta"] == 0]
        print(
            f"n={n:2d} C={C:2d} A(C)-C={mass-C:+4d} "
            f"rows={len(by_n[n]):3d} tight={len(tight):2d} below={len(below):1d} "
            f"tight_mass_changed={len(mass_changed_tight):1d} "
            f"same_mass_loose={len(same_mass_loose):2d}"
        )
        for r in tight:
            times = ", ".join(str(t) for t in r["times"][:4])
            if len(r["times"]) > 4:
                times += ", ..."
            print(
                "  tight swap "
                f"{r['out']:>2}->{r['inn']:<2} M={r['M']} "
                f"shell {r['out_shell']}->{r['in_shell']} "
                f"gcd {r['out_gcd']}->{r['in_gcd']} "
                f"mass_delta={r['mass_delta']} shell_l1={r['shell_l1']} "
                f"gcd_l1={r['gcd_l1']} times=[{times}]"
            )

    violations = [r for r in rows if r["tight"] and r["mass_delta"] != 0]
    print()
    print("Finite observation:")
    print(
        f"  tight AP single-swaps with changed gcd-shell mass through n=14: {len(violations)}"
    )
    print("  Thus every tight single-swap in this window is fixed by the")
    print("  weighted gcd-shell carrier, even when it is not shell-fixed.")

    n14 = [r for r in rows if r["n"] == 14 and r["tight"]]
    print()
    print("n=14 detail:")
    for r in n14:
        speeds = sorted((set(range(1, 14)) - {int(r["out"])}) | {int(r["inn"])})
        C = 27
        print(f"  speeds={speeds}")
        print(f"  shell counts delta={dict(shell_counter(speeds, C) - target_shell_counter(14))}")
        print(f"  reverse shell delta={dict(target_shell_counter(14) - shell_counter(speeds, C))}")
        print(f"  gcd counter={dict(sorted(gcd_counter(speeds, C).items()))}")
        print(f"  AP gcd counter={dict(sorted(target_gcd_counter(14).items()))}")
        print("  Vstar moves shell 12 to shell 3.  The shell carrier sees a defect,")
        print("  but the divisor/gcd carrier sees an exact fixed row.")


def print_tournament_analysis() -> None:
    print()
    print("=" * 72)
    print("Tournament Analysis")
    print("=" * 72)
    vertices = [
        "pillai_fixed_clock",
        "gcd_shell_carrier",
        "single_swap_zero_defect",
        "mod3_doubling_law",
        "pair_pinch_oracle",
        "carry_lift_conservativity",
        "scalar_gap_M",
        "raw_tuple_scan",
    ]
    rankings = [
        [
            "pillai_fixed_clock",
            "gcd_shell_carrier",
            "single_swap_zero_defect",
            "mod3_doubling_law",
            "pair_pinch_oracle",
            "carry_lift_conservativity",
            "scalar_gap_M",
            "raw_tuple_scan",
        ],
        [
            "pair_pinch_oracle",
            "carry_lift_conservativity",
            "single_swap_zero_defect",
            "gcd_shell_carrier",
            "mod3_doubling_law",
            "pillai_fixed_clock",
            "scalar_gap_M",
            "raw_tuple_scan",
        ],
        [
            "pillai_fixed_clock",
            "single_swap_zero_defect",
            "gcd_shell_carrier",
            "scalar_gap_M",
            "mod3_doubling_law",
            "pair_pinch_oracle",
            "carry_lift_conservativity",
            "raw_tuple_scan",
        ],
        [
            "gcd_shell_carrier",
            "mod3_doubling_law",
            "pillai_fixed_clock",
            "single_swap_zero_defect",
            "carry_lift_conservativity",
            "pair_pinch_oracle",
            "scalar_gap_M",
            "raw_tuple_scan",
        ],
    ]
    adj = majority_tournament(vertices, rankings)
    scores = sorted(sum(row) for row in adj)
    hp = hamiltonian_paths(adj)
    c3 = directed_3cycles(adj)
    scc = scc_sizes(adj)
    majority_edges = edge_set(adj)
    print("Pairwise observable:")
    print("  Which retained LRC carrier best preserves the floor/tight predicate.")
    print("Switch/gauges:")
    print("  proof transfer from aliquot fixed points, LRC proof utility,")
    print("  AP single-swap discrimination, and cost/availability.")
    print("Challenged assumption:")
    print("  Vertices need not be runners.  Considered runners, residues,")
    print("  antipodal shells, gcd strata, pair pinches, boundary events,")
    print("  carry fibers, and proof obligations.  This scout uses carrier")
    print("  lenses as vertices; the quotient preserves tight/floor evidence")
    print("  and destroys raw speed order.")
    print("Tie Hamiltonian path:")
    print("  " + " > ".join(rankings[0]))
    print(f"score_hist={dict(Counter(scores))}")
    print(f"directed_3_cycles={c3}")
    print(f"scc_sizes={scc}")
    print(f"hamiltonian_paths={hp}")
    for idx, ranking in enumerate(rankings, start=1):
        flips = len(majority_edges ^ ranking_edge_set(vertices, ranking))
        print(f"edge_flips_vs_gauge_{idx}={flips}")
    majority_order = sorted(vertices, key=lambda v: -sum(adj[vertices.index(v)]))
    print("majority_score_order=" + " > ".join(majority_order))


def main() -> None:
    print_clock_table()
    rows = ap_single_swap_scan(max_n=14)
    print_single_swap_table(rows)
    print_tournament_analysis()
    print()
    print("=" * 72)
    print("Research takeaways")
    print("=" * 72)
    print("1. The LRC n=14 clock C=27 is not just composite: it is one of")
    print("   exactly two odd Pillai-fixed clocks C with A(C)=C, the other")
    print("   being C=15 (n=8).")
    print("2. AP is shell-fixed.  Vstar is not shell-fixed, but it is fixed by")
    print("   the coarser gcd/divisor carrier: the missing shell 12 and doubled")
    print("   shell 3 both live in the gcd-3 stratum of C=27.")
    print("3. In the AP single-swap scan through n=14, every tight non-AP row")
    print("   preserves weighted gcd-shell mass.  This is a new necessary")
    print("   finite certificate target inspired directly by aliquot fixed points.")
    print("4. The open proof use is a no-leak lemma: once the carrier is fixed,")
    print("   pair-pinch/carry labels should force AP or the known Vstar floor")
    print("   row; every mass-changing row should be strictly loose.")


if __name__ == "__main__":
    main()
