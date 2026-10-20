"""Independent forest/permutation audit of the path-cover response transform.

All labelled simple graphs on n=1,...,5, with no inherited filters.
q_R(c) counts unordered covers by c oriented paths (singletons one orientation).
F_R(t) sums t^(number of red consecutive pairs) over all vertex permutations.
Checks are explicit, including under python -O. Standard library only.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
from math import comb, factorial
from pathlib import Path
import json


def check(test, message):
    if not test:
        raise RuntimeError(message)


def forest_covers(n, red_edges):
    """Direct edge-subset enumeration, with no Hamilton path routine."""
    q = [0] * (n + 1)
    accepted = 0
    for mask in range(1 << len(red_edges)):
        degree = [0] * n
        parent = list(range(n))

        def root(v):
            while parent[v] != v:
                v = parent[v]
            return v

        edge_count = 0
        valid = True
        for index, (u, v) in enumerate(red_edges):
            if not (mask >> index) & 1:
                continue
            edge_count += 1
            degree[u] += 1
            degree[v] += 1
            a, b = root(u), root(v)
            if degree[u] > 2 or degree[v] > 2 or a == b:
                valid = False
                break
            parent[a] = b
        if valid:
            components = n - edge_count
            nonsingleton_components = len({root(v) for v in range(n) if degree[v]})
            q[components] += 2 ** nonsingleton_components
            accepted += 1
    return q, accepted


def permutation_response(n, edge_bits, red_mask):
    """Literal permutation polynomial; independent from forest enumeration."""
    coefficients = [0] * n
    for word in permutations(range(n)):
        red_count = sum(bool(red_mask & edge_bits[tuple(sorted((u, v)))])
                        for u, v in zip(word, word[1:]))
        coefficients[red_count] += 1
    return coefficients


def response_from_cover(n, q):
    coefficients = [0] * n
    for c in range(1, n + 1):
        power = n - c
        for j in range(power + 1):
            coefficients[j] += factorial(c) * q[c] * comb(power, j) * (-1) ** (power - j)
    return coefficients


def main():
    checked = Counter()
    forest_subsets = accepted_forests = 0
    examples = {}
    for n in range(1, 6):
        edges = list(combinations(range(n), 2))
        edge_bits = {edge: 1 << i for i, edge in enumerate(edges)}
        full_mask = (1 << len(edges)) - 1
        records = []
        for red_mask in range(full_mask + 1):
            red_edges = [edge for i, edge in enumerate(edges) if red_mask >> i & 1]
            q, accepted = forest_covers(n, red_edges)
            polynomial = permutation_response(n, edge_bits, red_mask)
            check(response_from_cover(n, q) == polynomial, f"cover transform n={n} R={red_mask}")
            check(sum(polynomial) == factorial(n), "total permutations")
            check(sum(j * value for j, value in enumerate(polynomial)) ==
                  2 * len(red_edges) * factorial(n - 1), "first derivative at one")
            for c in range(1, n + 1):
                degree = n - c
                inverse = sum(comb(j, degree) * polynomial[j] for j in range(degree, n))
                check(inverse == factorial(c) * q[c], f"inverse coefficient n={n} R={red_mask} c={c}")
            if red_mask == 0:
                check(q[n] == 1 and sum(q) == 1 and polynomial[0] == factorial(n), "empty positive control")
            if red_mask == full_mask:
                check(polynomial[-1] == factorial(n) and sum(polynomial[:-1]) == 0, "complete positive control")
                check(all(factorial(c) * q[c] == factorial(n) * comb(n - 1, c - 1)
                          for c in range(1, n + 1)), "complete cover profile")
            records.append((q, polynomial))
            checked[n] += 1
            forest_subsets += 1 << len(red_edges)
            accepted_forests += accepted
        for red_mask, (q, polynomial) in enumerate(records):
            qbar, fbar = records[full_mask ^ red_mask]
            check(fbar == polynomial[::-1], f"response reciprocity n={n} R={red_mask}")
            check(polynomial[0] == qbar[1] == sum((-1) ** (n - c) * factorial(c) * q[c] for c in range(1, n + 1)),
                  "alternating complement Hamilton count")
            for d in range(1, n + 1):
                transformed = sum((-1) ** (n - c) * comb(c - 1, d - 1) * factorial(c) * q[c]
                                  for c in range(d, n + 1))
                check(transformed == factorial(d) * qbar[d], f"full complement cover transform n={n} R={red_mask} d={d}")
        if n == 4:
            for name, chosen in [("claw", [(0, 1), (0, 2), (0, 3)]),
                                 ("path4", [(0, 1), (1, 2), (2, 3)]),
                                 ("cycle4", [(0, 1), (1, 2), (2, 3), (0, 3)])]:
                red_mask = sum(edge_bits[tuple(sorted(edge))] for edge in chosen)
                q, f = records[red_mask]
                examples[name] = {"edges": chosen, "q_c_for_c_1_to_n": q[1:], "F_coefficients_low_to_high": f,
                                  "blue_H_ordered": f[0], "red_H_ordered": f[-1]}
            check(examples["claw"]["F_coefficients_low_to_high"] == [0, 12, 12, 0], "claw hostile extremes")
    check(dict(checked) == {1: 1, 2: 2, 3: 8, 4: 64, 5: 1024}, "explicit graph universe")
    check(forest_subsets == sum(3 ** comb(n, 2) for n in range(1, 6)), "forest-subset universe accounting")
    receipt = {
        "status": "FINITE-EXACT independent verification of elementary identities",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "universe": "All labelled simple undirected graphs on [n], n=1,...,5; no filters",
        "graphs_by_order": dict(checked), "total_graphs": sum(checked.values()),
        "edge_subsets_examined_for_linear_forests": forest_subsets,
        "accepted_linear_forests_including_repetitions_across_graphs": accepted_forests,
        "checks": [
            "F_R(t)=sum_c c! q_R(c) (t-1)^(n-c)",
            "c! q_R(c)=[x^(n-c)] F_R(1+x)",
            "F_Rbar(t)=t^(n-1) F_R(1/t)",
            "H_ordered(Rbar)=sum_c (-1)^(n-c) c! q_R(c)",
            "d! q_Rbar(d)=sum_(c>=d) (-1)^(n-c) binom(c-1,d-1) c! q_R(c)",
            "F_R(1)=n!; F_R'(1)=2|E(R)|(n-1)!"
        ],
        "independent_paths": "Direct edge-subset/union-find forest enumeration versus literal vertex permutations",
        "examples": examples,
        "scope": "Ordered Hamilton paths and unordered oriented path covers; not tournament counts or knot invariants"
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
