#!/usr/bin/env python3
"""S591: sumset-support calculus for the LRC C=2n-1 ledger.

The point is not to enumerate a large box.  It is to separate the support
layers that THM-401/HYP-2132 put on the table:

  1. speed-shell support in Z/(2n-1), modulo antipodes;
  2. pair-sum shell support, i.e. where pinch denominators land mod C;
  3. actual low pair denominators and divisibility shields;
  4. unit-visible versus nonunit holes.

This is a small exact notebook over named rows from S571-S573 and n=14.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import gcd


def dist(x):
    x %= 1
    return min(x, 1 - x)


def shell(C, r):
    r %= C
    if r == 0:
        return 0
    return min(r, C - r)


def all_shells(n):
    return list(range(1, n))


def exact_M_and_arg(V):
    cands = {F(0), F(1, 2)}
    for a, b in combinations(V, 2):
        for D in (a + b, abs(a - b)):
            if D:
                for m in range(D + 1):
                    cands.add(F(m, D))
    for v in V:
        for m in range(2 * v + 1):
            cands.add(F(m, 2 * v))

    best = F(-1)
    best_t = F(0)
    best_active = ()
    for t in cands:
        vals = [dist(v * t) for v in V]
        mn = min(vals)
        if mn > best:
            best = mn
            best_t = t
            best_active = tuple(v for v, val in zip(V, vals) if val == mn)
    return best, best_t, best_active


def inv_mod(a, m):
    return pow(a % m, -1, m)


def format_counter(counter):
    return "{" + ",".join(f"{k}:{counter[k]}" for k in sorted(counter)) + "}"


def group_by_gcd(items, C):
    out = defaultdict(list)
    for a in items:
        out[gcd(a, C)].append(a)
    return {g: tuple(v) for g, v in sorted(out.items())}


def row_summary(name, n, V):
    C = 2 * n - 1
    V = tuple(V)
    shells = all_shells(n)
    unit_shells = {a for a in shells if gcd(a, C) == 1}

    speed_counts = Counter(shell(C, v) for v in V)
    missing = tuple(a for a in shells if speed_counts[a] == 0)
    doubled = tuple((a, speed_counts[a]) for a in shells if speed_counts[a] > 1)
    zero_hits = speed_counts[0]
    perfect = (not missing) and (not doubled) and zero_hits == 0 and len(V) == n - 1

    pair_shell_counts = Counter()
    pair_den_counts = Counter()
    for a, b in combinations(V, 2):
        pair_shell_counts[shell(C, a + b)] += 1
        pair_den_counts[a + b] += 1

    pair_missing = tuple(a for a in shells if pair_shell_counts[a] == 0)
    pair_private = tuple(a for a in shells if pair_shell_counts[a] == 1)
    pair_zero = pair_shell_counts[0]

    low_denominators = tuple(D for D in sorted(pair_den_counts) if D <= n)
    shielded_low = {}
    unshielded_low = []
    for D in low_denominators:
        shields = tuple(v for v in V if v % D == 0)
        if shields:
            shielded_low[D] = shields
        else:
            unshielded_low.append(D)

    M, t, active = exact_M_and_arg(V)
    unit_witnesses = []
    for a in missing:
        if gcd(a, C) == 1:
            k = inv_mod(a, C)
            vals = [dist(F(v * k, C)) for v in V]
            unit_witnesses.append((a, F(k, C), min(vals)))

    return {
        "name": name,
        "n": n,
        "C": C,
        "V": V,
        "M": M,
        "t": t,
        "active": active,
        "floor": F(1, n),
        "edge": F(2, C),
        "speed_counts": speed_counts,
        "missing": missing,
        "doubled": doubled,
        "zero_hits": zero_hits,
        "perfect": perfect,
        "missing_unit": tuple(a for a in missing if a in unit_shells),
        "missing_nonunit": tuple(a for a in missing if a not in unit_shells),
        "doubled_unit": tuple((a, c) for a, c in doubled if a in unit_shells),
        "doubled_nonunit": tuple((a, c) for a, c in doubled if a not in unit_shells),
        "missing_by_gcd": group_by_gcd(missing, C),
        "doubled_by_gcd": group_by_gcd([a for a, _ in doubled], C),
        "pair_shell_counts": pair_shell_counts,
        "pair_missing": pair_missing,
        "pair_private": pair_private,
        "pair_zero": pair_zero,
        "pair_missing_unit": tuple(a for a in pair_missing if a in unit_shells),
        "pair_missing_nonunit": tuple(a for a in pair_missing if a not in unit_shells),
        "pair_private_unit": tuple(a for a in pair_private if a in unit_shells),
        "pair_private_nonunit": tuple(a for a in pair_private if a not in unit_shells),
        "pair_den_counts": pair_den_counts,
        "low_denominators": low_denominators,
        "shielded_low": shielded_low,
        "unshielded_low": tuple(unshielded_low),
        "unit_witnesses": tuple(unit_witnesses),
    }


def ap_identity_report():
    rows = []
    for n in range(4, 21):
        C = 2 * n - 1
        V = tuple(range(1, n))
        summary = row_summary(f"AP_n{n}", n, V)
        speed_ok = summary["perfect"]
        pair_missing_ok = summary["pair_missing"] == (1,)
        low_unshielded_ok = summary["unshielded_low"] == (n,)
        rows.append((n, C, speed_ok, pair_missing_ok, low_unshielded_ok))
    return rows


def shell_strata_table():
    out = []
    for n in (11, 12, 13, 14):
        C = 2 * n - 1
        strata = defaultdict(list)
        for a in all_shells(n):
            strata[gcd(a, C)].append(a)
        out.append((n, C, {g: tuple(v) for g, v in sorted(strata.items())}))
    return out


def lens_value(summary, lens):
    if lens == "speed_shell_transversal":
        return (summary["perfect"], summary["missing"], summary["doubled"], summary["zero_hits"])
    if lens == "unit_visible_holes":
        return summary["missing_unit"]
    if lens == "nonunit_holes":
        return (summary["missing_nonunit"], summary["doubled_nonunit"])
    if lens == "pair_sum_shell_holes":
        return (summary["pair_missing"], summary["pair_private"], summary["pair_zero"])
    if lens == "low_denominator_shields":
        return (summary["low_denominators"], summary["unshielded_low"])
    if lens == "witness_denominator":
        return (summary["M"].denominator, summary["t"].denominator)
    if lens == "gcd_strata":
        return (tuple(summary["missing_by_gcd"].items()), tuple(summary["doubled_by_gcd"].items()))
    if lens == "raw_sumset_size":
        return (len(summary["pair_den_counts"]), len(summary["pair_shell_counts"]))
    raise KeyError(lens)


LENS_META = {
    "low_denominator_shields": {"payload": 5, "maturity": 5, "cost": 2},
    "unit_visible_holes": {"payload": 5, "maturity": 5, "cost": 1},
    "speed_shell_transversal": {"payload": 5, "maturity": 4, "cost": 1},
    "nonunit_holes": {"payload": 4, "maturity": 4, "cost": 2},
    "gcd_strata": {"payload": 4, "maturity": 3, "cost": 2},
    "witness_denominator": {"payload": 4, "maturity": 5, "cost": 5},
    "pair_sum_shell_holes": {"payload": 3, "maturity": 3, "cost": 2},
    "raw_sumset_size": {"payload": 1, "maturity": 1, "cost": 1},
}


def lens_metrics(summaries):
    metrics = {}
    for lens, meta in LENS_META.items():
        values = {lens_value(s, lens) for s in summaries}
        metrics[lens] = {
            **meta,
            "separates": len(values),
            "values": values,
        }
    return metrics


def proof_key(item):
    name, meta = item
    return (meta["maturity"], meta["payload"], meta["separates"], -meta["cost"], name)


def cost_key(item):
    name, meta = item
    return (-meta["cost"], meta["separates"], meta["payload"], meta["maturity"], name)


def tournament_from_key(names, metrics, key_fn):
    idx = {name: i for i, name in enumerate(names)}
    adj = {name: set() for name in names}
    for a, b in combinations(names, 2):
        ka = key_fn((a, metrics[a]))
        kb = key_fn((b, metrics[b]))
        if ka >= kb:
            adj[a].add(b)
        else:
            adj[b].add(a)
    return adj, idx


def directed_3_cycles(names, adj):
    total = 0
    for a, b, c in combinations(names, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            total += 1
        if c in adj[a] and b in adj[c] and a in adj[b]:
            total += 1
    return total


def scc_sizes(names, adj):
    def reach(start):
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    reaches = {u: reach(u) for u in names}
    unused = set(names)
    sizes = []
    while unused:
        u = next(iter(unused))
        comp = {v for v in unused if v in reaches[u] and u in reaches[v]}
        sizes.append(len(comp))
        unused -= comp
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count(names, adj):
    n = len(names)
    idx = {name: i for i, name in enumerate(names)}
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            u = names[last]
            for v in adj[u]:
                j = idx[v]
                if not (mask & (1 << j)):
                    dp[mask | (1 << j)][j] += count
    return sum(dp[-1])


def edge_flip_count(names, adj_a, adj_b):
    flips = 0
    for a, b in combinations(names, 2):
        if (b in adj_a[a]) != (b in adj_b[a]):
            flips += 1
    return flips


def print_row_summary(summary):
    M = summary["M"]
    print(f"{summary['name']}: n={summary['n']} C={summary['C']} V={summary['V']}")
    print(
        f"  M={M} at t={summary['t']} active={summary['active']} "
        f"floor={summary['floor']} edge=2/C={summary['edge']} "
        f"class={'floor' if M == summary['floor'] else ('open-gap' if M < summary['edge'] else 'loose/edge')}"
    )
    print(
        f"  speed shells: perfect={summary['perfect']} zero_hits={summary['zero_hits']} "
        f"missing={summary['missing']} doubled={summary['doubled']}"
    )
    print(
        f"    missing_unit={summary['missing_unit']} missing_nonunit={summary['missing_nonunit']} "
        f"doubled_nonunit={summary['doubled_nonunit']} gcd_missing={summary['missing_by_gcd']}"
    )
    print(
        f"  pair-sum shells mod C: missing={summary['pair_missing']} "
        f"private={summary['pair_private']} zero_pair_sums={summary['pair_zero']}"
    )
    print(
        f"  low pair denominators <=n: {summary['low_denominators']} "
        f"unshielded={summary['unshielded_low']}"
    )
    if summary["unit_witnesses"]:
        witnesses = ", ".join(f"miss {a}->t={t}, min={mn}" for a, t, mn in summary["unit_witnesses"])
        print(f"  unit-miss inverse witnesses: {witnesses}")
    print()


def main():
    named_rows = [
        ("AP_n7", 7, tuple(range(1, 7))),
        ("open_gap_n7_S573", 7, (1, 5, 6, 11, 16, 17)),
        ("AP_n8", 8, tuple(range(1, 8))),
        ("nonunit_hole_n8_A", 8, (1, 2, 3, 4, 5, 7, 18)),
        ("nonunit_hole_n8_B", 8, (1, 3, 4, 5, 7, 13, 18)),
        ("AP_n14", 14, tuple(range(1, 14))),
        ("Vstar_n14", 14, (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)),
        ("doubled_apex_edge_n14", 14, tuple(range(1, 13)) + (26,)),
        ("unit_shift_AP_n14", 14, tuple(range(2, 15))),
    ]
    summaries = [row_summary(name, n, V) for name, n, V in named_rows]

    print("=== C=2n-1 shell strata around the paper frontier ===")
    for n, C, strata in shell_strata_table():
        print(f"n={n:2d} C={C:2d} strata=" + ", ".join(f"gcd {g}:{v}" for g, v in strata.items()))
    print()

    print("=== AP support identity check ===")
    print("For AP_n={1,...,n-1}: speed shells are perfect, pair-sum shell 1 is missing, and D=n is the first unshielded low denominator.")
    failures = [row for row in ap_identity_report() if row[2:] != (True, True, True)]
    print(f"verified n=4..20; failures={failures}")
    print()

    print("=== Named-row support ledgers ===")
    for summary in summaries:
        print_row_summary(summary)

    print("=== Cross-row lessons ===")
    print("1. Perfect speed-shell support and pair-sum shell support are different ledgers.")
    print("   AP rows hit every speed shell once, but their pair-sum shell support always misses shell 1.")
    print("2. Unit-visible speed-shell holes are easy: inverse clocks give the 2/C witness unless a zero residue intervenes.")
    print("3. The hard bounded rows hide in nonunit speed-shell holes or lifted denominators: they pass unit visibility but still need D/U/N, lift, or endpoint labels.")
    print("4. Low denominator shielding is the fold-sieve layer: AP and V* leave D=n unshielded at the floor; unit-shift AP shields D=n by speed n and becomes loose.")
    print()

    print("=== Tournament Analysis over support lenses ===")
    metrics = lens_metrics(summaries)
    names = list(LENS_META)
    for name in sorted(names, key=lambda x: proof_key((x, metrics[x])), reverse=True):
        meta = metrics[name]
        print(
            f"{name}: maturity={meta['maturity']} payload={meta['payload']} "
            f"separates={meta['separates']} cost={meta['cost']}"
        )
    proof_adj, _ = tournament_from_key(names, metrics, proof_key)
    cost_adj, _ = tournament_from_key(names, metrics, cost_key)
    score_hist = Counter(len(proof_adj[name]) for name in names)
    hpath_order = [name for name in sorted(names, key=lambda x: proof_key((x, metrics[x])), reverse=True)]
    print(f"score_hist={format_counter(score_hist)}")
    print(f"directed_3_cycles={directed_3_cycles(names, proof_adj)}")
    print(f"SCC_sizes={scc_sizes(names, proof_adj)}")
    print(f"Hamiltonian_paths={hamiltonian_path_count(names, proof_adj)}")
    print("Hamiltonian_path=" + " > ".join(hpath_order))
    print(f"edge_flips_vs_cost_only={edge_flip_count(names, proof_adj, cost_adj)}")
    print()

    print("=== Assumption challenge ===")
    print("Vertices used here are proof lenses and support fibers, not runners or arcs.")
    print("The quotient preserves whether a row has a unit-visible missing shell, an unshielded low pinch denominator, or a lifted exact witness denominator.")
    print("It destroys exact phase order and unmarked tournament class data; if a fiber mixes floor/open-gap/loose rows, the lens must be lifted with endpoint-owner or D/U/N labels.")


if __name__ == "__main__":
    main()
