"""Module types of the {-3}-self-dual classes in a selfdual3 record (selfdual3_fewcycle_20261104.c, mode selfdual3).
For every decomposable class of order >= 7: all proper modules of size >= 2 and their types (C3, TT_k, or score list).
Usage: python3 selfdual3_modules_20261104.py selfdual3_fewcycle_20261104_records/sd3_10.out
kind-pasteur-2026-11-04-S10 (scratch name sd3_modules.py)."""
import re, sys, itertools, collections
def parse(path):
    out = []
    for l in open(path):
        if l.startswith("SD3 n="):
            n = int(re.search(r"n=(\d+)", l).group(1))
            d = int(re.search(r"decomposable=(\d)", l).group(1))
            H = int(re.search(r"H=(\d+)", l).group(1)); hc = int(re.search(r"hc=(\d+)", l).group(1))
            sc = int(re.search(r"selfconverse=(\d)", l).group(1))
            rows = [int(x, 16) for x in l.split("rows:")[1].split()]
            out.append((n, d, H, hc, sc, rows))
    return out
def is_module(n, rows, M):
    for x in range(n):
        if M >> x & 1: continue
        o = i = False
        for m in range(n):
            if M >> m & 1:
                if rows[x] >> m & 1: o = True
                else: i = True
        if o and i: return False
    return True
def subtype(n, rows, M):
    vs = [v for v in range(n) if M >> v & 1]
    k = len(vs)
    sc = sorted(sum(1 for w in vs if rows[v] >> w & 1) for v in vs)
    if k == 3 and sc == [1, 1, 1]: return "C3"
    if sc == list(range(k)): return f"TT{k}"
    return f"{k}:{sc}"
for n, d, H, hc, sc, rows in parse(sys.argv[1]):
    if n < 7 or not d: continue
    mods = [M for M in range(1, (1 << n) - 1) if 2 <= bin(M).count("1") <= n - 1 and is_module(n, rows, M)]
    # maximal proper modules
    maxi = [M for M in mods if not any(M != N and (M & N) == M for N in mods)]
    types = collections.Counter(subtype(n, rows, M) for M in mods)
    print(f"n={n} H={H} hc={hc} sc={sc}  #modules={len(mods)}  module types: {dict(types)}")
