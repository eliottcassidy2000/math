#!/usr/bin/env python3
"""self_lines_structure_20261020.py -- opus-2026-10-20-S46
Which permutations g realise a self-line?  For the kappa-twist (g . T_t = T_{kappa t}) and the
path-reversal twist, tabulate the consistent g by cycle type and solution dimension, and the
relation of g to the base path (how many path arcs g maps to path arcs).  Also the connectivity
of the class-level line multigraph (all lines; blue lines on SC classes; black lines) from the
census files.  Usage: python3 self_lines_structure_20261020.py n [n ...]
"""
import sys, itertools, os
from collections import Counter

def cycle_type(g):
    n = len(g); seen = [False] * n; ct = []
    for v in range(n):
        if seen[v]: continue
        L = 0; w = v
        while not seen[w]: seen[w] = True; w = g[w]; L += 1
        ct.append(L)
    return tuple(sorted(ct, reverse=True))

def structure(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pid = {p: k for k, p in enumerate(pairs)}; P = len(pairs)
    ispath = [j == i + 1 for (i, j) in pairs]; ONE = P
    def solve(g, twist):
        parent = list(range(P + 1)); par = [0] * (P + 1)
        def find(x):
            if parent[x] == x: return x, 0
            r, p = find(parent[x]); parent[x] = r; par[x] ^= p; return r, par[x]
        def union(x, y, c):
            rx, px = find(x); ry, py = find(y)
            if rx == ry: return (px ^ py) == c
            parent[rx] = ry; par[rx] = px ^ py ^ c; return True
        for k in range(P):
            if ispath[k] and not union(k, ONE, 0): return None
        ginv = [0] * n
        for v in range(n): ginv[g[v]] = v
        for k, (u, v) in enumerate(pairs):
            a, b = ginv[u], ginv[v]; e2 = pid[(a, b) if a < b else (b, a)]
            if not union(k, e2, twist(k, a, b)): return None
        roots = set(find(k)[0] for k in range(P + 1)); rootONE = find(ONE)[0]
        return len(roots) - 1   # dimension
    tw_kappa = lambda k, a, b: ((0 if ispath[k] else 1) + (1 if a > b else 0)) & 1
    tw_rev = lambda k, a, b: ((1 if ispath[k] else 0) + (1 if a > b else 0)) & 1
    hist = {"kappa": Counter(), "rev": Counter()}; total = {"kappa": 0, "rev": 0}
    for g in itertools.permutations(range(n)):
        ct = cycle_type(g)
        pathkeep = sum(1 for i in range(n - 1) if abs(g[i] - g[i + 1]) == 1)   # path arcs mapped onto path pairs
        for name, tw in (("kappa", tw_kappa), ("rev", tw_rev)):
            d = solve(g, tw)
            if d is not None:
                hist[name][(ct, d, pathkeep)] += 1; total[name] += 2 ** d
    for name in ("kappa", "rev"):
        print(f"n={n} {name}-twist: P = {total[name]}; consistent g by (cycle type, dim, #path arcs kept): {dict(sorted(hist[name].items()))}")
    sys.stdout.flush()

def connectivity(n):
    base = os.path.dirname(os.path.abspath(__file__))
    fn = os.path.join(base, f"tiling_line_metagraph_census_20261020_lines_n{n}.tsv")
    cf = os.path.join(base, f"tiling_line_metagraph_census_20261020_classes_n{n}.tsv")
    if not os.path.exists(fn): return
    sc = {}
    with open(cf, encoding="utf-8") as f:
        hdr = f.readline().strip().split("\t")
        for line in f:
            x = line.rstrip("\n").split("\t"); d = dict(zip(hdr, x)); sc[int(d["class"])] = int(d["SC"])
    C = len(sc)
    def comps(edges, nodes):
        parent = {v: v for v in nodes}
        def find(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for a, b in edges: parent[find(a)] = find(b)
        return len(set(find(v) for v in nodes))
    allE = []; blueE = []; blackE = []
    with open(fn, encoding="utf-8") as f:
        f.readline()
        for line in f:
            a, b, bl, bk = map(int, line.split("\t"))
            if a != b:
                allE.append((a, b))
                if bl: blueE.append((a, b))
                if bk: blackE.append((a, b))
    scn = [c for c in sc if sc[c]]
    print(f"n={n} connectivity: all lines -> {comps(allE, list(sc))} component(s) on {C} classes; blue lines on the {len(scn)} SC classes -> {comps(blueE, scn)} component(s); black lines -> {comps(blackE, list(sc))} component(s) on all classes")
    sys.stdout.flush()

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [5, 6, 7]
    for n in ns:
        connectivity(n)
        if n <= 9: structure(n)
