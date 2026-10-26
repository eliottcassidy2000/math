#!/usr/bin/env python3
"""Independent bounded tc audit; Python stdlib only, no input files or Git.

Not a replay of any n=10 computation.  No producer code is imported.
All checks use explicit failures and survive Python's -O flag.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def matrix(n):
    return [[0 for _ in range(n)] for _ in range(n)]


def tournament(n, bits):
    a = matrix(n)
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = (bits >> k) & 1
        a[j][i] = 1 - a[i][j]
    return a


def preserves(a, p):
    return all(a[i][j] == a[p[i]][p[j]]
               for i in range(len(a)) for j in range(len(a)))


def brute_aut(a):
    return [p for p in permutations(range(len(a))) if preserves(a, p)]


def graph_aut(a):
    """Enumerate actual vertex maps, without twin compression or codec decoding."""
    n = len(a)
    colors = [0] * n
    while True:
        signatures = [(colors[i], tuple(sorted(Counter(
            (a[i][j], colors[j]) for j in range(n)).items())))
            for i in range(n)]
        names = {s: k for k, s in enumerate(sorted(set(signatures)))}
        refined = [names[s] for s in signatures]
        if len(set(refined)) == len(set(colors)):
            colors = refined
            break
        colors = refined
    domains = [{j for j in range(n) if colors[i] == colors[j]}
               for i in range(n)]

    def visit(p, unused):
        if len(p) == n:
            result = tuple(p[i] for i in range(n))
            require(preserves(a, result), "backtracking produced a non-automorphism")
            yield result
            return
        choices = {
            i: [j for j in sorted(domains[i] & unused)
                if all(a[i][k] == a[j][v] and a[k][i] == a[v][j]
                       for k, v in p.items())]
            for i in range(n) if i not in p}
        i = min(choices, key=lambda k: (len(choices[k]), k))
        for j in choices[i]:
            yield from visit(dict(p, **{}) | {i: j}, unused - {j})

    yield from visit({}, set(range(n)))


def codec(a, kind):
    width = 7 if kind == "joint" else 9
    g = matrix(width * len(a))

    def edge(u, v, c=1):
        g[u][v] = g[v][u] = c

    edges = ([(0, 1), (0, 2), (2, 3), (0, 4), (4, 5), (5, 6)]
             if kind == "joint" else
             [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5),
              (2, 6), (2, 7), (2, 8)])
    for v in range(len(a)):
        for i, j in edges:
            edge(width * v + i, width * v + j)
    for u in range(len(a)):
        for v in range(len(a)):
            if a[u][v]:
                edge(width * u + 1, width * v + (3 if kind == "joint" else 2),
                     2 if kind == "joint" else 1)
    return g


def hamiltonian_words(a):
    return {p for p in permutations(range(len(a)))
            if all(a[x][y] for x, y in zip(p, p[1:]))}


def odd_cycles(a):
    # Minimum vertex first removes rotations only, not reversal or supports.
    out = []
    for length in range(3, len(a) + 1, 2):
        for p in permutations(range(len(a)), length):
            if p[0] == min(p) and all(
                    a[p[i]][p[(i + 1) % length]] for i in range(length)):
                out.append(p)
    return out


def independence_at_two(cycles):
    # Enumerate disjoint cycle collections directly; retain distinct cycle words.
    ans = 0
    for mask in range(1 << len(cycles)):
        used = set()
        count = 0
        for j, cycle in enumerate(cycles):
            if mask & (1 << j):
                if used.intersection(cycle):
                    break
                used.update(cycle)
                count += 1
        else:
            ans += 2 ** count
    return ans


def main():
    # Independently validate the generic graph routine against brute force.
    graph_controls = 0
    for n in range(5):
        for bits in range(1 << (n * (n - 1) // 2)):
            g = matrix(n)
            for k, (i, j) in enumerate(combinations(range(n), 2)):
                g[i][j] = g[j][i] = (bits >> k) & 1
            require(set(graph_aut(g)) == set(brute_aut(g)), "generic graph control")
            graph_controls += 1

    frame = codec(tournament(1, 0), "joint")
    require(len(brute_aut(frame)) == 1, "seven-vertex frame not asymmetric")
    rows = []
    path_orbits = 0
    lifted_maps = 0
    for n in range(5):
        spectrum = set()
        count = 1 << (n * (n - 1) // 2)
        for bits in range(count):
            a = tournament(n, bits)
            aut = brute_aut(a)
            paths = hamiltonian_words(a)
            require(len(paths) == independence_at_two(odd_cycles(a)), "finite OCF check")
            remain = set(paths)
            orbit_count = 0
            while remain:
                p = min(remain)
                orbit = {tuple(s[v] for v in p) for s in aut}
                require(orbit <= paths, "word action not closed")
                require(len(orbit) == len(aut), "non-free word action")
                require(orbit <= remain, "word orbits overlap incorrectly")
                remain.difference_update(orbit)
                orbit_count += 1
            require(len(paths) == orbit_count * len(aut), "orbit formula")
            spectrum.add(orbit_count)
            path_orbits += orbit_count
            joint = codec(a, "joint")
            joint_maps = set(graph_aut(joint))
            expected = {tuple(7 * s[v] + r for v in range(n) for r in range(7))
                        for s in aut}
            require(joint_maps == expected, "joint codec added or lost automorphisms")
            for s in permutations(range(n)):
                lift = tuple(7 * s[v] + r for v in range(n) for r in range(7))
                require(preserves(joint, lift) == preserves(a, s), "joint lift equivariance")
                lifted_maps += 1
        rows.append({"n": n, "labelled_tournaments": count, "tc_spectrum": sorted(spectrum)})

    # Enumerate full Gamma automorphisms, rather than inserting the 12^n factor.
    c3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    kernel_controls = []
    for name, a in [("empty", []), ("singleton", tournament(1, 0)),
                    ("transitive_two", tournament(2, 0)), ("C3", c3)]:
        n = len(a)
        g = codec(a, "gamma")
        maps = list(graph_aut(g))
        projected = Counter()
        for p in maps:
            require(all(p[9 * v] % 9 == 0 for v in range(n)), "Gamma center type")
            projected[tuple(p[9 * v] // 9 for v in range(n))] += 1
        aut = set(brute_aut(a))
        require(set(projected) == aut, "Gamma quotient is not Aut(T)")
        require(set(projected.values()) == {12 ** n}, "Gamma leaf-kernel order")
        H = len(hamiltonian_words(a))
        tc = H // len(aut)
        naive = Fraction(H, len(maps))
        require(Fraction(12 ** n * H, len(maps)) == tc, "kernel-corrected quotient")
        kernel_controls.append({"name": name, "n": n, "H": H,
                                "AutT": len(aut), "AutGamma": len(maps),
                                "kernel": 12 ** n, "naive": str(naive), "tc": tc})
    require(kernel_controls[1]["naive"] == "1/12", "minimal denominator hostile")
    # Directed Hamiltonian cycles modulo rotation are deliberately not word paths.
    require(len(odd_cycles(c3)) == 1 and len(brute_aut(c3)) == 3, "cycle caveat")
    print(json.dumps({"status": "PASS", "scope": "all 76 labelled tournaments, n=0..4",
                      "graph_controls": graph_controls, "frame_brute_order": 1,
                      "small_orders": rows, "explicit_word_orbits": path_orbits,
                      "lifted_permutations": lifted_maps, "gamma_controls": kernel_controls,
                      "cycle_caveat": "C3 has one rotation-quotiented cycle fixed by all 3 automorphisms",
                      "n10": "NOT RUN OR CERTIFIED", "all_order_tc_exclusion": "OPEN"},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
