#!/usr/bin/env python3
# mac-mini-2026-11-08-firstrealizers: per-run summary of sc_walk_runs.out against the exact inventories.
# For each run: the strong self-converse values found; for orders 12 and 13 they are compared with the exact strong
# inventory below 1125 (sc_odd13_1124.out); for order 9 with the census (cls_9.txt in WORK, optional).
# usage: walk_summary.py [WORK]   (run in the packet directory)
import sys, os, collections
WORK = sys.argv[1] if len(sys.argv) > 1 else "."
exact = collections.defaultdict(set)
for line in open("sc_odd13_1124.out"):
    t, h, s = line.split()
    if s == "1": exact[{"E": 12, "O": 13}[t]].add(int(h))
p9 = os.path.join(WORK, "cls_9.txt")
if os.path.exists(p9):
    for line in open(p9):
        hd = line.split("|")[0].split()
        if int(hd[3]) and int(hd[4]): exact[9].add(int(hd[1]))
runs, cur = [], None
for line in open("sc_walk_runs.out"):
    if line.startswith("# run:"):
        cur = (tuple(map(int, line.split()[2:])), set()); runs.append(cur)
    elif line.startswith("W "): cur[1].add(int(line.split()[1]))
targets = [1131, 1145, 1231, 1245, 1263, 1391, 1459, 1475, 1535, 1567, 1583, 1699, 1891]
by_order = collections.defaultdict(set)
for (n, lo, hi, steps, seed), vals in runs:
    by_order[n] |= vals
    print(f"run n={n} H in [{lo},{hi}] steps={steps} seed={seed}: {len(vals)} strong self-converse values")
for n in sorted(by_order):
    vals = by_order[n]
    line = f"order {n}, all runs: {len(vals)} values"
    if n in exact:
        cap = 3357 if n == 9 else 1124
        lo_vals = {h for h in vals if h <= cap}
        line += (f"; at H <= {cap}: {len(lo_vals)} vs exact {len(exact[n])}, identical: {lo_vals == exact[n]}"
                 f"; above {cap}: {len(vals - lo_vals)}")
    line += f"; the 13 unresolved values found: {sorted(set(targets) & vals)}"
    print(line)
