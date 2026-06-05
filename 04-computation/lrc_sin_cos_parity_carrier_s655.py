#!/usr/bin/env python3
"""S655: sine/cosine parity carrier for LRC n=14.

The user prompt points at the analytic split

    sin x = odd powers / odd factorials,
    cos x = even powers / even factorials,
    cot x = (log sin x)' = cos x / sin x.

This scout ports that split to the repo's LRC n=14 proof state.  The analogy is
typed: "sin" is the odd boundary/zero carrier, while "cos" is the derivative
or slack carrier attached to that boundary.  In LRC n=14 the AP row's even
speed layer is the doubled prime-7 problem with exact slack; the odd layer is
the binder.

The computation is intentionally finite and exact where it matters:
  * AP parity strata at t=1/n for several even n.
  * AP and Vstar n=14 wall signatures at t=1/14.
  * Exact THM-369 pair-sum-pinch scan of n=14 AP single swaps.
  * A small proof-route Tournament Analysis over possible uses of the carrier.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import cos, gcd, pi, sin


def v2(x: int) -> int:
    if x == 0:
        return 99
    k = 0
    while x % 2 == 0:
        x //= 2
        k += 1
    return k


def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def dist(x: Fraction) -> Fraction:
    y = frac(x)
    return min(y, 1 - y)


def slope_at(v: int, t: Fraction) -> int | str:
    """Derivative of ||v t|| away from cusp points."""
    y = frac(v * t)
    if y == 0:
        return "zero"
    if y == Fraction(1, 2):
        return "cusp"
    return v if y < Fraction(1, 2) else -v


def ap(n: int) -> list[int]:
    return list(range(1, n))


def vstar(n: int) -> list[int]:
    return sorted([x for x in range(1, n) if x != n - 2] + [2 * n - 4])


def primitive(speeds: list[int]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def maximin_pair_pinch(speeds: list[int]) -> tuple[Fraction, list[Fraction]]:
    """Exact maximin over pair-sum pinch candidates, as in S646/THM-369 use."""
    candidates: set[tuple[int, int]] = set()
    for a, b in combinations(speeds, 2):
        D = a + b
        for m in range(D + 1):
            candidates.add((m, D))

    best = Fraction(0, 1)
    best_times: list[Fraction] = []
    for m, D in candidates:
        local = min(dist(Fraction(v * m, D)) for v in speeds)
        t = Fraction(m, D)
        if local > best:
            best = local
            best_times = [t]
        elif local == best:
            best_times.append(t)
    return best, sorted(set(best_times))


def shell(v: int, C: int) -> int:
    r = v % C
    return min(r, C - r) if r else 0


def gcd_shell_counter(speeds: list[int], C: int) -> Counter[int]:
    return Counter(gcd(shell(v, C), C) for v in speeds)


def wall_signature(n: int, speeds: list[int], t: Fraction) -> dict[str, object]:
    dists = {v: dist(v * t) for v in speeds}
    m = min(dists.values())
    active = [v for v in speeds if dists[v] == m]
    slopes = [slope_at(v, t) for v in active]
    numeric_slopes = [s for s in slopes if isinstance(s, int)]
    by_depth: dict[int, Fraction] = {}
    for v in speeds:
        depth = v2(v)
        by_depth[depth] = min(by_depth.get(depth, Fraction(99, 1)), dists[v])
    return {
        "t": t,
        "M_at_t": m,
        "active": active,
        "active_parity": Counter(v % 2 for v in active),
        "active_slopes": slopes,
        "right_derivative": min(numeric_slopes) if numeric_slopes else None,
        "left_derivative": max(numeric_slopes) if numeric_slopes else None,
        "depth_min": by_depth,
    }


def active_wall_pairs(n: int, speeds: list[int]) -> list[tuple[int, int, Fraction]]:
    """Unique active odd complement pairs at the floor times of a tight row."""
    M, times = maximin_pair_pinch(speeds)
    floor = Fraction(1, n)
    if M != floor:
        return []
    pairs: set[tuple[int, int, Fraction]] = set()
    for t in times:
        sig = wall_signature(n, speeds, t)
        active = sorted(sig["active"])
        if len(active) == 2:
            pairs.add((active[0], active[1], t))
    return sorted(pairs, key=lambda x: (x[0] + x[1], x[0], x[2]))


def format_depths(depth_min: dict[int, Fraction], floor: Fraction) -> str:
    parts = []
    for k in sorted(depth_min):
        ratio = depth_min[k] / floor
        parts.append(f"v2={k}: {depth_min[k]} ({ratio}*delta)")
    return "; ".join(parts)


def single_swap_scan_n14() -> list[dict[str, object]]:
    n = 14
    base = set(ap(n))
    rows = []
    for out_v in range(1, n):
        for in_v in range(1, 2 * n + 1):
            if in_v in base:
                continue
            speeds = sorted((base - {out_v}) | {in_v})
            if len(speeds) != n - 1 or not primitive(speeds):
                continue
            M, times = maximin_pair_pinch(speeds)
            rows.append(
                {
                    "out": out_v,
                    "in": in_v,
                    "speeds": speeds,
                    "M": M,
                    "times": times,
                    "gcd_shells": gcd_shell_counter(speeds, 27),
                }
            )
    return rows


def cot_product_half(C: int) -> float:
    prod = 1.0
    for k in range(1, (C + 1) // 2):
        prod *= (cos(pi * k / C) / sin(pi * k / C)) ** 2
    return prod


def directed_3cycles(adj: list[list[int]]) -> int:
    c = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        ):
            c += 1
    return c


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, reverse: bool = False) -> set[int]:
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
        comp = reach(v) & reach(v, reverse=True)
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
            if not dp[mask][v]:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[-1])


def majority_tournament(vertices: list[str], rankings: list[list[str]]) -> list[list[int]]:
    pos = [{v: i for i, v in enumerate(ranking)} for ranking in rankings]
    n = len(vertices)
    adj = [[0] * n for _ in range(n)]
    for i, u in enumerate(vertices):
        for j, v in enumerate(vertices):
            if i == j:
                continue
            wins = sum(1 for p in pos if p[u] < p[v])
            if wins > len(pos) / 2 or (wins == len(pos) / 2 and pos[0][u] < pos[0][v]):
                adj[i][j] = 1
    return adj


def main() -> None:
    print("S655 LRC Sin/Cos Parity Carrier")
    print("=" * 72)
    print()
    print("Analytic dictionary:")
    print("  sin  = odd boundary carrier: zeros / active floor constraints.")
    print("  cos  = derivative/slack carrier: adjacent response at those zeros.")
    print("  cot  = log-derivative ledger: divide by the boundary but retain every pole.")
    print("  LRC  = keep parity, pair-sum, gcd-shell, and owner/carry side channels.")

    print("\nAP parity strata at t=1/n")
    print("-" * 72)
    for n in [8, 10, 12, 14, 16, 18, 28]:
        floor = Fraction(1, n)
        sig = wall_signature(n, ap(n), floor)
        print(
            f"n={n:2d}, active={sig['active']}, slopes={sig['active_slopes']}, "
            f"depths: {format_depths(sig['depth_min'], floor)}"
        )

    print("\nn=14 named rows")
    print("-" * 72)
    named = [("AP", ap(14)), ("Vstar", vstar(14))]
    for name, speeds in named:
        M, times = maximin_pair_pinch(speeds)
        sig = wall_signature(14, speeds, Fraction(1, 14))
        print(f"{name}: speeds={speeds}")
        print(f"  exact M={M}, best_times_count={len(times)}, first_times={times[:8]}")
        print(
            f"  at t=1/14: M={sig['M_at_t']}, active={sig['active']}, "
            f"parity={dict(sig['active_parity'])}, slopes={sig['active_slopes']}"
        )
        print(f"  depth minima: {format_depths(sig['depth_min'], Fraction(1,14))}")
        print(f"  C=27 gcd shells: {dict(sorted(gcd_shell_counter(speeds, 27).items()))}")

    print("\nEven/odd 14/21 bridge")
    print("-" * 72)
    print("Goldbach/Lemoine duplicate center: 7+7=14 and 7+2*7=21.")
    print("The LRC floor wall does not use the duplicate center directly; it uses")
    print("odd complement pairs around the central 7, all summing to the even row 14.")
    for name, speeds in named:
        pairs = active_wall_pairs(14, speeds)
        unique = sorted({(a, b) for a, b, _ in pairs})
        print(f"{name}: active odd complement pairs={unique}")
        for a, b in unique:
            left = 7 - a
            right = b - 7
            print(
                f"  {a}+{b}=14, centered offsets=(-{left},+{right}), "
                f"pair product={a*b}, diagonal defect={(7-a)*(b-7)}"
            )
    print("Reading: the even target 14 is held by odd boundary pairs; the odd")
    print("target 21 is the diagonal scalar shadow.  Progress for LRC 14 needs")
    print("the off-diagonal boundary pairs plus C=27 carry/gcd labels, not the")
    print("raw 14/21 scalar match.")

    print("\nn=14 AP single-swap pinch census")
    print("-" * 72)
    rows = single_swap_scan_n14()
    floor = Fraction(1, 14)
    tight = [r for r in rows if r["M"] == floor]
    below = [r for r in rows if r["M"] < floor]
    loose = [r for r in rows if r["M"] > floor]
    print(f"rows={len(rows)}, below={len(below)}, tight={len(tight)}, loose={len(loose)}")
    for r in tight:
        speeds = r["speeds"]
        times = r["times"]
        sigs = [wall_signature(14, speeds, t) for t in times[:6]]
        print(f"tight swap {r['out']} -> {r['in']}: M={r['M']}, gcd_shells={dict(sorted(r['gcd_shells'].items()))}")
        print(f"  first best times={times[:10]}")
        for sig in sigs:
            print(
                f"    t={sig['t']}: active={sig['active']}, "
                f"slopes={sig['active_slopes']}, parity={dict(sig['active_parity'])}"
            )
    if not below:
        print("No AP single-swap row falls below the floor; this rechecks the S646 guardrail.")

    print("\nCotangent polygon side check")
    print("-" * 72)
    for C in [15, 21, 27, 29]:
        prod = cot_product_half(C)
        print(f"C={C:2d}: prod_1^((C-1)/2) cot(pi*k/C)^2 = {prod:.12g}; C*prod={C*prod:.12g}")
    print("For odd C this half-product collapses to 1/C; composite side channels are invisible unless gcd shells are retained.")

    print("\nTournament Analysis")
    print("-" * 72)
    vertices = [
        "odd_boundary_sin_carrier",
        "even_derivative_cos_slack",
        "cot_log_derivative_pole_ledger",
        "pair_sum_pinch_oracle",
        "C27_gcd_shell_carrier",
        "owner_carry_conservativity",
        "raw_parity_numerology",
    ]
    rankings = [
        [
            "odd_boundary_sin_carrier",
            "even_derivative_cos_slack",
            "C27_gcd_shell_carrier",
            "pair_sum_pinch_oracle",
            "owner_carry_conservativity",
            "cot_log_derivative_pole_ledger",
            "raw_parity_numerology",
        ],
        [
            "C27_gcd_shell_carrier",
            "pair_sum_pinch_oracle",
            "owner_carry_conservativity",
            "odd_boundary_sin_carrier",
            "even_derivative_cos_slack",
            "cot_log_derivative_pole_ledger",
            "raw_parity_numerology",
        ],
        [
            "pair_sum_pinch_oracle",
            "odd_boundary_sin_carrier",
            "C27_gcd_shell_carrier",
            "even_derivative_cos_slack",
            "owner_carry_conservativity",
            "cot_log_derivative_pole_ledger",
            "raw_parity_numerology",
        ],
    ]
    adj = majority_tournament(vertices, rankings)
    scores = Counter(sum(row) for row in adj)
    print(f"vertices={vertices}")
    print(f"score_hist={dict(sorted(scores.items()))}")
    print(f"directed_3cycles={directed_3cycles(adj)}")
    print(f"scc_sizes={scc_sizes(adj)}")
    print(f"hamiltonian_paths={hamiltonian_paths(adj)}")
    print("ranking:")
    score_rows = sorted(((sum(adj[i]), v) for i, v in enumerate(vertices)), reverse=True)
    for score, vertex in score_rows:
        print(f"  {score}: {vertex}")
    print()
    print("Assumption challenge:")
    print("  Vertices are not runners.  They are proof carriers: odd boundary, even")
    print("  derivative/slack, pair-sum pinch clocks, C=27 gcd shells, and owner/carry")
    print("  obligations.  The quotient preserves floor/tight evidence and destroys")
    print("  raw runner order, which must be reattached only when needed.")

    print("\nSynthesis")
    print("-" * 72)
    print("1. The sine/cos split sharpens HYP-2116: for n=14 the odd layer is the")
    print("   active boundary carrier; the even layer is a derivative/slack certificate")
    print("   imported from solved prime 7.")
    print("2. AP and Vstar have the same odd active wall at t=1/14, but Vstar is")
    print("   distinguished by C=27 gcd-shell transport; hence the next proof must")
    print("   combine odd-wall slopes with the HYP-2222 gcd fixed pocket.")
    print("3. The single-swap census leaves no subfloor row and isolates the same")
    print("   AP/Vstar pocket, converting the user's trig analogy into a concrete")
    print("   no-leak lemma target: odd boundary + even slack + gcd shell + owner/carry.")


if __name__ == "__main__":
    main()
