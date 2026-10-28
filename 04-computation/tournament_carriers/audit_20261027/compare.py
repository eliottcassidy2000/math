"""Referee: carrier-by-carrier comparison of the referee census (census_n7.json) with the
author's carrier_census_n7.json (data file only; no author code is imported).
Usage: python compare.py <author carrier_census_n7.json>"""
import itertools
import json
import sys

import numpy as np


def canon(n, arcs, perms):
    a = np.array(arcs, dtype=np.int64)
    codes = ((np.int64(1) << (perms[:, a[:, 0]] * n + perms[:, a[:, 1]])).sum(axis=1))
    return int(codes.min())


mine = json.load(open("census_n7.json"))
theirs = json.load(open(sys.argv[1]))
for n in sorted(mine, key=int):
    nn = int(n)
    perms = np.array(list(itertools.permutations(range(nn))), dtype=np.int64)
    M = {canon(nn, r["arcs"], perms): r for r in mine[n]}
    T = {canon(nn, r["arcs"], perms): r for r in theirs[n]}
    same_set = set(M) == set(T)
    diffs = []
    for c, r in M.items():
        t = T.get(c)
        if t is None:
            diffs.append(("missing", r["arcs"]))
            continue
        for k_m, k_t in (("linear", "linear"), ("weak", "weak"), ("count", "marked_types"),
                         ("A_fixes_vertex", "A_fixes_vertex"), ("order", "order")):
            if k_m in r and k_t in t and t[k_t] is not None and r[k_m] != t[k_t]:
                diffs.append((k_m, r["arcs"], r[k_m], t[k_t]))
    print(f"n={n}: referee {len(M)} classes, author {len(T)} classes, identical class sets: {same_set}, "
          f"flag disagreements: {len(diffs)}")
    for d in diffs[:10]:
        print("   ", d)
