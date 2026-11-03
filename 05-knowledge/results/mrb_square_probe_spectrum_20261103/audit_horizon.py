#!/usr/bin/env python3
"""Independent exhaustive final-graph permutation counts; no producer imports."""
import hashlib
import itertools
import json
from math import isqrt
from pathlib import Path


def check(ok, error):
    if not ok:
        raise ValueError(error)


def graph(source, word):
    n = len(source)
    size = n + len(word)
    a = [[False] * size for _ in range(size)]
    for i in range(n):
        for j in range(n):
            a[i][j] = bool(source[i] & 2**j)
    for step, mask in enumerate(word):
        j = n + step
        for i in range(j):
            a[i][j] = bool(mask & 2**i)
            a[j][i] = not a[i][j]
    return a


def count(a):
    return sum(all(a[x][y] for x, y in zip(p, p[1:]))
               for p in itertools.permutations(range(len(a))))


def main():
    data = json.loads(Path(__file__).with_name("horizon.json").read_text())
    report = {"audit": "same agent; independent Boolean adjacency and full permutations"}
    censuses = []
    for test in data["tensor_tests"]:
        n, k = len(test["source"]), test["horizon"]
        values = [count(graph(test["source"], word))
                  for word in itertools.product(*(range(2**(n+j)) for j in range(k)))]
        digest = hashlib.sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()
        check(digest == test["counts_sha256"], "independent horizon digest")
        censuses.append({"source": test["source"], "horizon": k, "counts": len(values), "sha256": digest})
    report["permutation_censuses"] = censuses
    host = data["common_host"]
    labels = list(map(int, host["labels"]))
    for witness in host["witnesses"]:
        selected = [labels[i] for i in witness["selected_pool_indices"]]
        a = graph(host["source"], witness["word"])
        for i, j in itertools.combinations(range(len(selected)), 2):
            v = selected[i] + selected[j]
            check(a[i][j] == (isqrt(v)**2 == v), "independent host decode")
        check(count(a) == witness["H"], "host count")
    report["arithmetic_host_selected_graphs"] = len(host["witnesses"])
    check([count(graph((2, 0), (s, 1))) for s in range(4)] == [3, 1, 5, 1], "nonterminal pair")
    check([count(graph((0,), (s, 1))) for s in range(2)] == [3, 1], "nonterminal singleton")
    report["terminal_constraint_hostiles"] = "both independently confirmed"
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
