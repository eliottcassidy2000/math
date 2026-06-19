#!/usr/bin/env python3
"""
lrc_invariant_atlas_codex_20260619.py

Creative invariant scout for Lonely Runner structure, focused on the LRC(14)
hard layer but written in ordinary speed-set language.

Question:
    which invariants actually determine the structure?

Method:
    Build a mixed bank of primitive 13-speed rows, compute exact gap data
    M(S)=max_t min_v ||v t|| and exact safe-set topology at level 1/14, then
    compare candidate invariant signatures by how much target variation remains
    inside their fibers.

Tournament Analysis declaration:
    Vertex set: safe-component boundary owner pairs, not runners.
    Pairwise observable: owner pair A -> B if A bounds more/larger safe
      components than B; ties follow lexicographic speed order.
    Fingerprints: owner score histogram, directed cycles, SCC sizes, and
      Hamiltonian path count when the owner-pair set is small.

Assumption challenge:
    I considered runners, residues, q-obligations, endpoint boundaries, safe
    components, additive anti-cosets, and proof obligations as vertices.
    Runner vertices preserve speed labels but collapse the endpoint geometry.
    Safe-component boundary owners preserve the exact level-1/14 obstruction
    and expose which speeds actually carve the surviving set.
"""

from __future__ import annotations

import itertools
import math
import random
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce
from math import gcd


H = F(1, 14)
RNG = random.Random(20260619)


def gcd_all(vals) -> int:
    return reduce(gcd, vals, 0)


def circ_norm(x: F) -> F:
    r = x - int(x)
    if r < 0:
        r += 1
    return r if r <= F(1, 2) else 1 - r


def candidate_taus(S: tuple[int, ...]) -> set[F]:
    out = {F(1, 2)}
    for v in S:
        k = 0
        while True:
            t = F(2 * k + 1, 2 * v)
            if t > F(1, 2):
                break
            out.add(t)
            k += 1
    for a, b in itertools.combinations(S, 2):
        for D in (a + b, abs(b - a)):
            if D <= 0:
                continue
            k = 1
            while True:
                t = F(k, D)
                if t > F(1, 2):
                    break
                out.add(t)
                k += 1
    return out


def exact_M(S: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    best = F(0)
    ats: list[F] = []
    for t in candidate_taus(S):
        val = min(circ_norm(v * t) for v in S)
        if val > best:
            best = val
            ats = [t]
        elif val == best:
            ats.append(t)
    return best, tuple(sorted(ats))


def danger_segments_with_owners(v: int, h: F = H):
    out = []
    for j in range(v):
        center = F(j, v)
        a = (center - h / v) % 1
        b = (center + h / v) % 1
        if a < b:
            out.append((a, b, v))
        else:
            out.append((a, F(1), v))
            out.append((F(0), b, v))
    return out


def safe_components(S: tuple[int, ...], h: F = H):
    """Exact safe intervals and the speed owners at each boundary."""
    bps = {F(0), F(1)}
    left_owners = defaultdict(set)   # intervals starting at boundary
    right_owners = defaultdict(set)  # intervals ending at boundary
    segments = []
    for v in S:
        for a, b, owner in danger_segments_with_owners(v, h):
            segments.append((a, b, owner))
            bps.add(a)
            bps.add(b)
            left_owners[a].add(owner)
            right_owners[b].add(owner)
    bps = sorted(bps)
    comps = []
    for a, b in zip(bps, bps[1:]):
        if a >= b:
            continue
        mid = (a + b) / 2
        if all(circ_norm(v * mid) > h for v in S):
            # If (a,b) is safe, then danger ended at a and begins at b.
            lo = tuple(sorted(right_owners[a]))
            hi = tuple(sorted(left_owners[b]))
            comps.append((a, b, lo, hi))
    return comps


def is_covering(S: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def residue_signature(S: tuple[int, ...]) -> tuple:
    out = []
    for m in (2, 3, 4, 5, 7, 13, 14):
        counts = Counter(v % m for v in S)
        out.append(tuple(sorted(counts.values(), reverse=True)))
    return tuple(out)


def q_coverage_signature(S: tuple[int, ...]) -> tuple[int, ...]:
    # Capped counts keep this a proof-facing coarse invariant.
    return tuple(min(3, sum(1 for v in S if v % q == 0)) for q in range(2, 15))


def first_unit_leak(S: tuple[int, ...], qmax: int = 60) -> tuple[int | None, int]:
    for q in range(2, qmax + 1):
        units = [a for a in range(1, q) if gcd(a, q) == 1]
        leaks = 0
        for a in units:
            if all(14 * min((v * a) % q, (-v * a) % q) >= q for v in S):
                leaks += 1
        if leaks:
            return q, leaks
    return None, 0


def additive_signature(S: tuple[int, ...]) -> tuple[int, int, int, int]:
    Sset = set(S)
    three = 0
    for a, b in itertools.combinations(S, 2):
        if a + b in Sset:
            three += 1
    sums = Counter(a + b for a, b in itertools.combinations(S, 2))
    four_collisions = sum(c * (c - 1) // 2 for c in sums.values() if c > 1)
    diffs = Counter(abs(a - b) for a, b in itertools.combinations(S, 2))
    diff_energy = sum(c * c for c in diffs.values())
    longest_run = 1
    cur = 1
    for a, b in zip(S, S[1:]):
        if b == a + 1:
            cur += 1
        else:
            longest_run = max(longest_run, cur)
            cur = 1
    longest_run = max(longest_run, cur)
    return three, four_collisions, diff_energy, longest_run


def subset_relation_signature(S: tuple[int, ...], max_support: int = 6) -> tuple:
    """Height-1 anti-coset proxy: disjoint equal subset sums."""
    by_sum = defaultdict(list)
    n = len(S)
    for mask in range(1, 1 << n):
        bits = mask.bit_count()
        if bits > max_support:
            continue
        sm = 0
        for i, v in enumerate(S):
            if mask & (1 << i):
                sm += v
        by_sum[sm].append((mask, bits))
    min_support = None
    support_counts = Counter()
    total = 0
    for masks in by_sum.values():
        if len(masks) < 2:
            continue
        for i in range(len(masks)):
            mi, bi = masks[i]
            for mj, bj in masks[i + 1:]:
                if mi & mj:
                    continue
                supp = bi + bj
                if supp < 2:
                    continue
                total += 1
                min_support = supp if min_support is None else min(min_support, supp)
                support_counts[min(supp, 8)] += 1
    return (
        min_support or 99,
        support_counts.get(4, 0),
        support_counts.get(5, 0),
        support_counts.get(6, 0),
        sum(c for s, c in support_counts.items() if s >= 6),
        total,
    )


def binding_signature(S: tuple[int, ...], M: F, tau: F) -> tuple:
    binders = tuple(v for v in S if circ_norm(v * tau) == M)
    rows = []
    for a, b in itertools.combinations(binders, 2):
        for kind, D in (("sum", a + b), ("diff", abs(a - b))):
            if D <= 0:
                continue
            if (D * tau).denominator != 1:
                continue
            j = M * D
            if j.denominator == 1:
                rows.append((14 * int(j) - D, D % 14, kind, a, b))
    rows.sort()
    if not rows:
        return (len(binders), 999, 999, "none")
    r, dmod, kind, a, b = rows[0]
    return (len(binders), r, dmod, kind)


def owner_pair_key(lo: tuple[int, ...], hi: tuple[int, ...]) -> tuple[int, int]:
    left = min(lo) if lo else -1
    right = min(hi) if hi else -1
    return tuple(sorted((left, right)))


def hamiltonian_paths(edge: list[list[int]]) -> int:
    n = len(edge)
    if n > 12:
        return -1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if edge[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def scc_sizes(edge: list[list[int]]) -> list[int]:
    n = len(edge)
    reach = [[bool(edge[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    seen = set()
    sizes = []
    for i in range(n):
        if i in seen:
            continue
        comp = {j for j in range(n) if reach[i][j] and reach[j][i]}
        seen |= comp
        sizes.append(len(comp))
    return sorted(sizes, reverse=True)


def owner_tournament_signature(comps) -> tuple:
    weights = Counter()
    for a, b, lo, hi in comps:
        weights[owner_pair_key(lo, hi)] += b - a
    owners = sorted(weights)
    n = len(owners)
    if n == 0:
        return (0, (), 0, (), 0)
    edge = [[0] * n for _ in range(n)]
    outdeg = [0] * n
    for i, A in enumerate(owners):
        for j, B in enumerate(owners):
            if i == j:
                continue
            if (weights[A], A) >= (weights[B], B):
                edge[i][j] = 1
                outdeg[i] += 1
    cycles = 0
    for i, j, k in itertools.combinations(range(n), 3):
        if edge[i][j] and edge[j][k] and edge[k][i]:
            cycles += 1
        if edge[i][k] and edge[k][j] and edge[j][i]:
            cycles += 1
    hist = tuple(sorted(Counter(outdeg).items()))
    return (n, hist, cycles, tuple(scc_sizes(edge)), hamiltonian_paths(edge))


def component_topology_signature(comps) -> tuple:
    owner_pairs = Counter(owner_pair_key(lo, hi) for _, _, lo, hi in comps)
    top_pair = owner_pairs.most_common(1)[0] if owner_pairs else (((-1, -1), 0))
    return (
        len(comps),
        len(owner_pairs),
        top_pair[0],
        top_pair[1],
    )


def component_geometry_signature(comps) -> tuple:
    lengths = [b - a for a, b, _, _ in comps]
    owner_pairs = Counter(owner_pair_key(lo, hi) for _, _, lo, hi in comps)
    top_pair = owner_pairs.most_common(1)[0] if owner_pairs else (((-1, -1), 0))
    return (
        len(comps),
        sum(lengths),
        max(lengths) if lengths else F(0),
        len(owner_pairs),
        top_pair[0],
        top_pair[1],
    )


def make_bank(limit: int = 520) -> dict[str, tuple[int, ...]]:
    bank: dict[str, tuple[int, ...]] = {}

    def add(label: str, vals):
        S = tuple(sorted(set(int(v) for v in vals)))
        if len(S) == 13 and min(S) > 0 and gcd_all(S) == 1:
            bank.setdefault(label, S)

    add("AP_1_13", range(1, 14))
    add("GW_24", list(range(1, 12)) + [13, 24])
    add("loose_36", list(range(1, 12)) + [13, 36])

    for w in range(14, 121):
        add(f"tail12_{w}", list(range(1, 12)) + [13, w])
        add(f"tail13_{w}", list(range(1, 13)) + [w])
        if len(bank) >= limit:
            break

    Ps = [(1, 2, 3, 12), (1, 2, 3, 5), (1, 2, 7, 12), (2, 3, 5, 7)]
    Es = [
        tuple(range(9)),
        (0, 1, 2, 3, 4, 5, 6, 21),
        (0, 1, 2, 3, 4, 5, 6, 7, 68),
        (0, 2, 3, 5, 8, 11, 13, 17, 19),
    ]
    for P in Ps:
        for E in Es:
            for V in (42, 56, 70, 84, 98, 112):
                S = set(P) | {V - e for e in E}
                add(f"S3_P{P}_V{V}_E{max(E)}", S)

    # Covering-biased random rows.  Force one multiple for each q, then fill.
    tries = 0
    while len(bank) < limit and tries < limit * 80:
        tries += 1
        vals = set()
        for q in range(2, 15):
            vals.add(q * RNG.randint(1, max(1, 90 // q)))
        while len(vals) < 13:
            vals.add(RNG.randint(1, 100))
        if len(vals) > 13:
            vals = set(RNG.sample(sorted(vals), 13))
        add(f"cover_rand_{tries}", vals)

    # Generic primitive rows, to keep the atlas from overfitting covering sets.
    tries = 0
    while len(bank) < limit + 120 and tries < 3000:
        tries += 1
        vals = tuple(sorted(RNG.sample(range(1, 75), 13)))
        add(f"gen_rand_{tries}", vals)

    return bank


def profile_row(label: str, S: tuple[int, ...]) -> dict:
    M, ats = exact_M(S)
    tau = ats[0]
    comps = safe_components(S)
    qleak, qleaks = first_unit_leak(S)
    return {
        "label": label,
        "S": S,
        "M": M,
        "slack": M - H,
        "safe_measure": sum(b - a for a, b, _, _ in comps),
        "component_count": len(comps),
        "widest_component": max((b - a for a, b, _, _ in comps), default=F(0)),
        "covering": is_covering(S),
        "qleak": qleak,
        "qleaks": qleaks,
        "residue": residue_signature(S),
        "qcov": q_coverage_signature(S),
        "add": additive_signature(S),
        "subset": subset_relation_signature(S),
        "component_topology": component_topology_signature(comps),
        "component_geometry": component_geometry_signature(comps),
        "owner_tournament": owner_tournament_signature(comps),
        "binding": binding_signature(S, M, tau),
    }


def variance(vals: list[float]) -> float:
    if not vals:
        return 0.0
    mu = sum(vals) / len(vals)
    return sum((x - mu) ** 2 for x in vals) / len(vals)


def fiber_report(rows: list[dict], key: str, target: str) -> tuple[float, tuple, tuple]:
    total = variance([float(r[target]) for r in rows])
    groups = defaultdict(list)
    for r in rows:
        groups[r[key]].append(r)
    weighted = 0.0
    worst = None
    for members in groups.values():
        vals = [float(r[target]) for r in members]
        weighted += len(vals) * variance(vals)
        spread = max(vals) - min(vals)
        if len(vals) >= 2 and (worst is None or spread > worst[0]):
            worst = (spread, members)
    ratio = weighted / (len(rows) * total) if total > 0 else 0.0
    sizes = [len(v) for v in groups.values()]
    meta = (len(groups), sum(1 for s in sizes if s == 1), max(sizes), sum(sizes) / len(sizes))
    return ratio, worst, meta


def corr(xs: list[float], ys: list[float]) -> float:
    if not xs or len(xs) != len(ys):
        return 0.0
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def numeric_invariants(row: dict) -> dict[str, float]:
    add = row["add"]
    subset = row["subset"]
    comp = row["component_geometry"]
    owner = row["owner_tournament"]
    return {
        "max_speed": max(row["S"]),
        "span": max(row["S"]) - min(row["S"]),
        "covering": float(row["covering"]),
        "first_unit_leak": float(row["qleak"] or 99),
        "qleaks": row["qleaks"],
        "three_term": add[0],
        "four_sum_collisions": add[1],
        "diff_energy": add[2],
        "longest_run": add[3],
        "subset_min_support": subset[0],
        "subset_support6plus": subset[4],
        "safe_components": comp[0],
        "safe_measure": float(comp[1]),
        "widest_component": float(comp[2]),
        "owner_pair_count": comp[3],
        "owner_tournament_vertices": owner[0],
        "owner_tournament_cycles": owner[2],
        "owner_tournament_hp": owner[4] if owner[4] >= 0 else 0,
        "binding_r": row["binding"][1],
    }


def main() -> None:
    print("LRC invariant atlas: what actually determines the structure?")
    print("=" * 78)
    bank = make_bank()
    print(f"bank size: {len(bank)} primitive 13-speed rows")
    rows = []
    for i, (label, S) in enumerate(bank.items(), 1):
        rows.append(profile_row(label, S))
        if i % 100 == 0:
            print(f"  profiled {i} rows...")
    rows.sort(key=lambda r: (r["M"], r["safe_measure"], r["label"]))

    print("\nHardest rows by exact M(S):")
    for r in rows[:12]:
        print(f"  {r['label']:<28} M={r['M']}={float(r['M']):.6f} "
              f"slack={float(r['slack']):.6f} L={float(r['safe_measure']):.6f} "
              f"comps={r['component_count']:>3} qleak={r['qleak']} "
              f"bind={r['binding']} S={r['S']}")

    print("\nNear-tight tail laws discovered in the bank:")
    tail13 = []
    tail12 = []
    by_label = {r["label"]: r for r in rows}
    for w in range(14, 121):
        r13 = by_label.get(f"tail13_{w}")
        if r13 and w % 13 == 0:
            a = w // 13
            tail13.append((w, r13["M"], F(a, 13 * a + 1), r13["binding"]))
        r12 = by_label.get(f"tail12_{w}")
        if r12 and w % 12 == 0 and w >= 36:
            a = w // 12
            tail12.append((w, r12["M"], F(a, 12 * a + 5), r12["binding"]))
    ok13 = all(m == pred for _, m, pred, _ in tail13)
    ok12 = all(m == pred for _, m, pred, _ in tail12)
    print(f"  family A: S={{1..12,13a}} has M=a/(13a+1) for tested a=2..9? {ok13}")
    for w, m, pred, bind in tail13[:5]:
        print(f"    w={w:<3} M={m} pred={pred} bind={bind}")
    print(f"  family B: S={{1..11,13,12a}} has M=a/(12a+5) for tested a>=3? {ok12}")
    for w, m, pred, bind in tail12[:5]:
        print(f"    w={w:<3} M={m} pred={pred} bind={bind}")

    print("\nBest single numeric correlations with hardness M(S):")
    target = [float(r["M"]) for r in rows]
    table = []
    for name in numeric_invariants(rows[0]):
        xs = [numeric_invariants(r)[name] for r in rows]
        table.append((abs(corr(xs, target)), corr(xs, target), name))
    for _, c, name in sorted(table, reverse=True)[:16]:
        print(f"  corr(M,{name}) = {c:+.3f}")

    print("\nFiber loss test: lower is more structure-determining.")
    print("  value = within-fiber variance / total variance for exact M.")
    keys = ("residue", "qcov", "add", "subset", "component_topology",
            "owner_tournament", "component_geometry", "binding")
    for key in keys:
        ratio, worst, meta = fiber_report(rows, key, "M")
        groups, singletons, maxfiber, avgfiber = meta
        print(f"  {key:<20} loss={ratio:.3f} groups={groups} "
              f"singletons={singletons} maxfiber={maxfiber} avgfiber={avgfiber:.2f}")
        if worst and worst[0] > 0:
            spread, members = worst
            a = min(members, key=lambda r: r["M"])
            b = max(members, key=lambda r: r["M"])
            print(f"    worst collision spread={spread:.6f}:")
            print(f"      low  {a['label']} M={float(a['M']):.6f} S={a['S']}")
            print(f"      high {b['label']} M={float(b['M']):.6f} S={b['S']}")

    print("\nFiber loss test for safe-measure L(S) at level 1/14.")
    for key in keys:
        ratio, worst, meta = fiber_report(rows, key, "safe_measure")
        groups, singletons, maxfiber, avgfiber = meta
        print(f"  {key:<20} loss={ratio:.3f} groups={groups} "
              f"singletons={singletons} maxfiber={maxfiber} avgfiber={avgfiber:.2f}")

    print("\nInterpretation:")
    print("  1. Residues and q-coverage are necessary filters, not structure.")
    print("     Their fibers still mix rows with visibly different M and L.")
    print("  2. Additive relations and height-1 anti-coset counts improve the picture,")
    print("     especially around AP/Goddyn-Wong families, but still do not determine")
    print("     the exact obstruction.")
    print("  3. Topology-only safe-component data helps, but exact safe geometry")
    print("     (lengths plus boundary owners) is where the structure becomes almost")
    print("     determined.  This is solution-facing, but less so than exact M: it")
    print("     records which speeds carve the witness reservoir.")
    print("  4. Binding denominator data nearly determines M, but it is even more")
    print("     solution-facing:")
    print("     it is best viewed as a readout, not a proof input.")
    print("  5. The recurring determinant is not one scalar invariant.  It is a stack:")
    print("       CRT obligation -> additive anti-coset -> safe-component owners ->")
    print("       binding denominator.")
    print("     The hard rows are exactly where a lower layer looks ambiguous and the")
    print("     next layer resolves the ambiguity.")


if __name__ == "__main__":
    main()
