#!/usr/bin/env python3
"""order_plus_path_F8_macmini_20261005.py -- |F_8|: classes of 8-tournaments having a Hamiltonian path whose
complement is acyclic ("order + path"); the envelope of the Syracuse 8-window alphabet (S13: 194 classes).
Reproduce: python3 order_plus_path_F8_macmini_20261005.py   (nauty gentourng; ~tens of minutes)"""
import itertools, subprocess, sys, time
from collections import defaultdict
import networkx as nx
m = 8
r = subprocess.run(["gentourng", "-q", str(m)], capture_output=True, text=True)
hosts = []
for line in r.stdout.split("\n"):
    line = line.strip()
    if len(line) == m * (m - 1) // 2 and set(line) <= {"0", "1"}:
        arcs = set(); k = 0
        for i in range(m):
            for j in range(i + 1, m):
                arcs.add((i, j) if line[k] == "1" else (j, i)); k += 1
        hosts.append(frozenset(arcs))
assert len(hosts) == 6880
def is_order_plus_path(arcs):
    succ = defaultdict(set)
    for a, b in arcs: succ[a].add(b)
    # enumerate Hamiltonian paths by DFS
    def dfs(path, used):
        if len(path) == m:
            pathset = {(path[i], path[i + 1]) for i in range(m - 1)}
            rest = nx.DiGraph([a for a in arcs if a not in pathset]); rest.add_nodes_from(range(m))
            return nx.is_directed_acyclic_graph(rest)
        for v in succ[path[-1]]:
            if v not in used and dfs(path + [v], used | {v}):
                return True
        return False
    return any(dfs([s], {s}) for s in range(m))
t0 = time.time()
F = sum(1 for h in hosts if is_order_plus_path(h))
print(f"|F_8| = {F} of 6880 classes  ({time.time()-t0:.0f}s);  |F_m| for m = 3..8: 2, 4, 11, 37, 143, {F}")
