#!/usr/bin/env python3
"""Independent exhaustive checks for bounded-black-defect cycle counting.

No external modules, randomness, repository operations, or network access.
Cycles are successor permutations, equivalently words anchored at vertex zero.
All labelled tournaments of orders 1 through 5 are scanned. The target cycle
claims are checked on every nontrivial strong tournament, self-converse or not.
All permutations through order 8 are scanned to check full-cycle root fibers.
"""

from collections import Counter, defaultdict
from itertools import permutations
from math import factorial, gcd
from pathlib import Path
import argparse
import json


HERE = Path(__file__).resolve().parent


def compose(f, g):
    return tuple(f[g[v]] for v in range(len(f)))


def power(f, k):
    ans = tuple(range(len(f)))
    for _ in range(k):
        ans = compose(f, ans)
    return ans


def cycle_type(p):
    seen = set()
    sizes = []
    for v in range(len(p)):
        if v not in seen:
            w, size = v, 0
            while w not in seen:
                seen.add(w)
                size += 1
                w = p[w]
            sizes.append(size)
    return tuple(sorted(sizes))


def strong(adj):
    n = len(adj)
    for start in range(n):
        seen, todo = {start}, [start]
        while todo:
            v = todo.pop()
            for w in range(n):
                if adj[v][w] and w not in seen:
                    seen.add(w)
                    todo.append(w)
        if len(seen) != n:
            return False
    return True


def successor(word):
    r = [0] * len(word)
    for i, v in enumerate(word):
        r[v] = word[(i + 1) % len(word)]
    return tuple(r)


def reversal(word):
    r = [0] * len(word)
    for v, w in zip(word, reversed(word)):
        r[v] = w
    return tuple(r)


def check_tournaments(max_order):
    rows = []
    detailed = []
    totals = Counter()
    for n in range(1, max_order + 1):
        vertices = tuple(range(n))
        perms = list(permutations(vertices))
        pairs = [(i, j) for i in vertices for j in range(i + 1, n)]
        stats = Counter(n=n, labelled=2 ** len(pairs))
        observed = Counter()
        for mask in range(1 << len(pairs)):
            adj = [[False] * n for _ in vertices]
            for bit, (v, w) in enumerate(pairs):
                if (mask >> bit) & 1:
                    adj[v][w] = True
                else:
                    adj[w][v] = True

            is_strong = strong(adj)
            aut = {
                p for p in perms
                if all(adj[p[v]][p[w]] == adj[v][w] for v, w in pairs)
            }
            anti = {
                p for p in perms
                if all(adj[p[v]][p[w]] != adj[v][w] for v, w in pairs)
            }
            stats["strong"] += is_strong
            stats["self_converse"] += bool(anti)
            stats["strong_self_converse"] += is_strong and bool(anti)
            if not is_strong:
                continue
            # The singleton has a Hamiltonian path but no directed cycle.
            if n == 1:
                assert len(aut) == 1
                continue
            stats["tested_nontrivial"] += 1
            assert not anti or len(anti) == len(aut)
            assert len(aut) % 2 == 1
            assert all(all(s % 2 for s in cycle_type(a)) for a in aut)

            paths = [
                p for p in perms
                if all(adj[p[i]][p[i + 1]] for i in range(n - 1))
            ]
            blue_paths = {p for p in paths if reversal(p) in anti}
            black_paths = set(paths) - blue_paths

            def representative(p):
                return min(tuple(a[v] for v in p) for a in aut)

            path_orbits = defaultdict(set)
            for p in paths:
                path_orbits[representative(p)].add(p)
            assert all(len(orbit) == len(aut) for orbit in path_orbits.values())
            assert all(not (orbit & blue_paths) or orbit <= blue_paths
                       for orbit in path_orbits.values())
            blue_orbits = {representative(p) for p in blue_paths}
            h, b = len(path_orbits), len(blue_orbits)
            delta = h - b
            assert len(paths) == h * len(aut)
            assert len(blue_paths) == b * len(aut)
            assert len(black_paths) == delta * len(aut)
            # Converse-reversal acts as an involution only in the SC case.
            if anti:
                assert delta % 2 == 0
            else:
                assert b == 0

            cycles = [
                p for p in paths
                if p[0] == 0 and adj[p[-1]][p[0]]
            ]
            assert cycles, "Camion's theorem check failed"
            black_cycles, all_blue_cycles = [], []
            chosen_black_cut = {}
            q_hist = Counter()
            for cyc in cycles:
                r = successor(cyc)
                powers = {power(r, k) for k in range(n)}
                d = len(aut & powers)
                q = n // d
                assert len(powers) == n and n % d == 0 and d % 2 == 1
                assert aut & powers == {power(r, k * q) for k in range(d)}
                # Equal cyclic-subgroup orbits induce regular tournaments.
                orbits = {
                    frozenset(power(r, k * q)[v] for k in range(d))
                    for v in vertices
                }
                assert len(orbits) == q
                assert all(len(orbit) == d for orbit in orbits)
                assert all(sum(adj[v][w] for w in orbit) == (d - 1) // 2
                           for orbit in orbits for v in orbit)
                minimum_outdegree = min(sum(row) for row in adj)
                assert minimum_outdegree >= (d - 1) // 2
                # Integer-only form of min_out >= (n/(delta+2)-1)/2.
                assert (2 * minimum_outdegree + 1) * (delta + 2) >= n
                cuts = [cyc[i:] + cyc[:i] for i in range(n)]
                cut_reps = {representative(p) for p in cuts}
                blue_reps = {representative(p) for p in cuts if p in blue_paths}
                black_reps = cut_reps - blue_reps
                assert len(cut_reps) == q
                assert len(blue_reps) in (0, gcd(2, q))
                assert len(black_reps) <= delta
                assert q <= delta + 2
                for i in range(n):
                    for j in range(n):
                        assert (representative(cuts[i]) ==
                                representative(cuts[j])) == ((i - j) % q == 0)
                if blue_reps:
                    a = next(i for i, cut in enumerate(cuts) if cut in blue_paths)
                    for i, cut in enumerate(cuts):
                        assert (cut in blue_paths) == (2 * (i - a) % q == 0)
                if black_reps:
                    black_cycles.append(r)
                    pick = min(p for p in cuts if p in black_paths)
                    assert pick not in chosen_black_cut
                    chosen_black_cut[pick] = r
                else:
                    all_blue_cycles.append(r)
                    for i in range(n):
                        assert compose(reversal(cuts[(i + 1) % n]),
                                       reversal(cuts[i])) == compose(r, r)
                    assert compose(r, r) in aut
                q_hist[q] += 1
                stats["cycle_checks"] += 1
                stats["cut_checks"] += n
            assert len(black_cycles) <= len(black_paths) == delta * len(aut)
            assert len(all_blue_cycles) <= n * len(aut)
            sharper_root_bound = 1 if n % 2 else n // 2
            assert len(all_blue_cycles) <= sharper_root_bound * len(aut)
            assert len(cycles) <= (delta + n) * len(aut)
            observed[(len(paths), len(aut), b, delta, len(cycles),
                      len(black_cycles), len(all_blue_cycles))] += 1
            stats["black_cycles"] += len(black_cycles)
            stats["all_blue_cycles"] += len(all_blue_cycles)
            detailed.append({
                "n": n, "mask": mask, "self_converse": bool(anti),
                "H": len(paths), "A": len(aut),
                "h": h, "b": b, "delta": delta, "HC": len(cycles),
                "black_cut_cycles": len(black_cycles),
                "all_blue_cycles": len(all_blue_cycles),
                "q_histogram": dict(sorted(q_hist.items())),
            })

        row = dict(stats)
        row["invariant_profiles"] = [
            {
                "H": key[0], "A": key[1], "b": key[2], "delta": key[3],
                "HC": key[4], "black_cut_cycles": key[5],
                "all_blue_cycles": key[6], "labelled_count": count,
            }
            for key, count in sorted(observed.items())
        ]
        rows.append(row)
        for k, value in stats.items():
            if k != "n":
                totals[k] += value
        print("TOURNAMENT", json.dumps(row, sort_keys=True))
    return {"rows": rows, "totals": dict(totals), "individual_cases": detailed}


def check_full_cycle_roots():
    rows = []
    for n in range(1, 9):
        perms = list(permutations(range(n)))
        cycles = [r for r in perms if cycle_type(r) == (n,)]
        assert len(cycles) == factorial(n - 1)
        fibers = defaultdict(list)
        for r in cycles:
            fibers[compose(r, r)].append(r)
        expected_type = (n,) if n % 2 else (n // 2, n // 2)
        expected_fiber = 1 if n % 2 else n // 2
        checked = 0
        for g in perms:
            actual = len(fibers[g])
            expected = expected_fiber if cycle_type(g) == expected_type else 0
            assert actual == expected, (n, g, actual, expected)
            assert actual <= n
            if n % 2 and actual:
                assert fibers[g][0] == power(g, (n + 1) // 2)
            checked += 1
        nonzero = {g: roots for g, roots in fibers.items() if roots}
        row = {
            "n": n,
            "all_permutations_checked": checked,
            "full_n_cycles": len(cycles),
            "square_image_size": len(nonzero),
            "nonempty_fiber_size": expected_fiber,
            "square_cycle_type": expected_type,
        }
        rows.append(row)
        print("ROOTS", json.dumps(row, sort_keys=True))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=5,
                        help="Exhaustive labelled tournament cutoff (default: 5).")
    args = parser.parse_args()
    assert 1 <= args.max_order <= 6, "This regression test is intended for n <= 6."
    results = {
        "tournaments": check_tournaments(args.max_order),
        "full_cycle_square_roots": check_full_cycle_roots(),
        "status": "PASS",
        "definitions": {
            "cycle_count": "Directed Hamiltonian cycles up to cyclic rotation.",
            "blue": "The path's positional reversal is an anti-automorphism.",
            "delta": "H/|Aut(T)| minus number of blue path orbits.",
        },
    }
    (HERE / "verification_results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n")
    print("TOTALS", json.dumps(results["tournaments"]["totals"], sort_keys=True))
    print("PASS: all assertions succeeded.")


if __name__ == "__main__":
    main()
