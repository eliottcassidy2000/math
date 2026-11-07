#!/usr/bin/env python3
"""independent_recount.py (klein-2026-11-06): recount, with a pure-Python subset DP written independently of hspec.c,
the Hamiltonian-path count h, the number of cyclic triangles c3, strong connectivity and (orders <= 10) the number of
Hamiltonian cycles of the witnesses recorded in census_<n>.txt (the W lines: the first tournament read with each value,
and the S lines: the first strong one), and of the level-floor minimisers in level_floor_witnesses.txt.
usage: independent_recount.py data_dir nmax [sample_per_order seed]
All witnesses of orders <= 8 are recounted; for larger orders a seeded random sample (default 300) plus the witnesses of
the maximum value and of the least new value of the order."""
import sys, random, gzip, os
from collections import defaultdict
D = sys.argv[1]; NMAX = int(sys.argv[2]); SAMPLE = int(sys.argv[3]) if len(sys.argv) > 3 else 300; SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261106
def open_census(path):
    if os.path.exists(path): return open(path)
    return gzip.open(path + ".gz", "rt")
def parse(n, line):
    """gentourng default format: upper triangle row by row, '1' at (i, j) means i -> j"""
    out = [set() for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if line[idx] == '1': out[i].add(j)
            else: out[j].add(i)
            idx += 1
    return out
def ham_paths(n, out):
    # dp over (visited set, last vertex): number of paths; iterate by increasing subset size
    dp = [dict() for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        if not row: continue
        for v, c in row.items():
            for w in out[v]:
                if not (S >> w) & 1:
                    T = S | (1 << w); dp[T][w] = dp[T].get(w, 0) + c
    return sum(dp[(1 << n) - 1].values())
def ham_cycles(n, out):
    # cycles through vertex 0: paths starting at 0 covering everything, ending at an in-neighbour of 0
    dp = [dict() for _ in range(1 << n)]; dp[1][0] = 1
    for S in range(1, 1 << n, 2):
        row = dp[S]
        if not row: continue
        for v, c in row.items():
            for w in out[v]:
                if not (S >> w) & 1:
                    T = S | (1 << w); dp[T][w] = dp[T].get(w, 0) + c
    full = dp[(1 << n) - 1]
    return sum(c for v, c in full.items() if 0 in out[v])
def c3_of(n, out):
    return sum(1 for i in range(n) for j in range(n) for k in range(n) if i < j and i < k and j != k and j in out[i] and k in out[j] and i in out[k])
def strong(n, out):
    def reach(start, adj):
        seen = {start}; stack = [start]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in seen: seen.add(w); stack.append(w)
        return seen
    inn = [set() for _ in range(n)]
    for v in range(n):
        for w in out[v]: inn[w].add(v)
    return len(reach(0, out)) == n and len(reach(0, inn)) == n
rng = random.Random(SEED); total = 0; bad = []
for n in range(3, NMAX + 1):
    T = defaultdict(lambda: defaultdict(lambda: [0, 0])); C = defaultdict(dict); W = {}; Sw = {}
    with open_census(f"{D}/census_{n}.txt") as fh:
        for line in fh:
            a = line.split()
            if not a: continue
            if a[0] == 'T': T[int(a[1])][int(a[2])][int(a[3])] += int(a[4])
            elif a[0] == 'C': C[int(a[1])][int(a[2])] = int(a[3])
            elif a[0] == 'W': W[int(a[1])] = a[2]
            elif a[0] == 'S': Sw[int(a[1])] = a[2]
    hs = sorted(W)
    if n <= 8: chosen = hs
    else:
        chosen = set(rng.sample(hs, min(SAMPLE, len(hs)))); chosen.add(max(hs))
        prev = set()
        with open_census(f"{D}/census_{n-1}.txt") as fh:
            for line in fh:
                if line.startswith('T '): prev.add(int(line.split()[1]))
        chosen.add(min(h for h in hs if h not in prev)); chosen = sorted(chosen)
    for h in chosen:
        for src, wit in (("W", W.get(h)), ("S", Sw.get(h))):
            if wit is None: continue
            out = parse(n, wit); hh = ham_paths(n, out); k = c3_of(n, out); st = strong(n, out)
            ok = (hh == h) and (k in T[h]) and (T[h][k][1 if st else 0] > 0) and (src == "W" or st)
            if C and st and n <= 10:
                hc = ham_cycles(n, out); ok = ok and (hc in C[h])
            total += 1
            if not ok: bad.append((n, src, h, hh, k, st))
    print(f"order {n}: {len(chosen)} values, witnesses recounted so far {total}, mismatches {len(bad)}")
# level-floor minimisers
for line in open(f"{D}/level_floor_witnesses.txt"):
    a = line.split()
    if not a or a[0] != 'order': continue
    n = int(a[1]); h = int(a[3]); k = int(a[5]); st = int(a[7].rstrip(':')); wit = a[8]
    out = parse(n, wit); ok = ham_paths(n, out) == h and c3_of(n, out) == k and strong(n, out) == bool(st)
    total += 1
    if not ok: bad.append(("floor", n, h, k))
print(f"level-floor minimisers recounted; total witnesses {total}; mismatches: {bad if bad else 'none'}")
print("ALL RECOUNTS AGREE" if not bad else "MISMATCH")
