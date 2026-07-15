#!/usr/bin/env python3
"""THM-790 n=8 check (opus-S305): classify all 2^21 tilings (certified by
A000568(8) = 6880), then test the parity-law predictions:
  blue tilings = 2^{(21+3)/2} = 4096  (=> 2048 blue lines);
  blue |Dx| in {0,16,32,48} (even n: 0 mod 16), max = 8(n-2) = 48;
  the transitive pipe: x_max = 168 -> 120, landing mixed;
plus the flow tables (type histograms, black flux) for the record."""
import numpy as np
from collections import defaultdict

n = 8
V = list(range(1, n+1))
tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]
m = len(tiles)
ti = {t: i for i, t in enumerate(tiles)}
refl = np.array([ti[(n-y+1, n-x+1)] for (x, y) in tiles])
POP = np.array([bin(i).count('1') for i in range(256)], dtype=np.int64)
FULLV = (1 << n) - 1
N = 1 << m
CH = 1 << 15

def feats(bits):
    """bits: (B,) int64 tilings -> per-vertex packed keys sorted + arc sums."""
    B = bits.shape[0]
    adj = np.zeros((B, n), dtype=np.int64)
    for v in range(2, n+1): adj[:, v-1] |= 1 << (v-2)          # path arcs
    for i, (x, y) in enumerate(tiles):
        b = (bits >> i) & 1
        adj[:, x-1] |= b << (y-1)
        adj[:, y-1] |= (1 - b) << (x-1)
    s = POP[adj]                                                # (B, n) scores
    inv = np.zeros_like(adj)
    for v in range(n):
        inv[:, v] = FULLV ^ (1 << v) ^ adj[:, v]                # in-mask
    c3 = np.zeros((B, n), dtype=np.int64)
    for v in range(n):
        for u in range(n):
            if u == v: continue
            has = (adj[:, v] >> u) & 1
            c3[:, v] += has * POP[adj[:, u] & inv[:, v]]
    so1 = np.zeros((B, n), dtype=np.int64); si1 = np.zeros_like(so1)
    so2 = np.zeros_like(so1); co1 = np.zeros_like(so1); ci1 = np.zeros_like(so1)
    so3 = np.zeros_like(so1)
    for v in range(n):
        for u in range(n):
            if u == v: continue
            hout = (adj[:, v] >> u) & 1
            hin = (adj[:, u] >> v) & 1
            so1[:, v] += hout * s[:, u];  si1[:, v] += hin * s[:, u]
            so2[:, v] += hout * s[:, u]**2
            so3[:, v] += hout * s[:, u]**3
            co1[:, v] += hout * c3[:, u]; ci1[:, v] += hin * c3[:, u]
    # second-round neighbour sums (2-step structure)
    to1 = np.zeros_like(so1); ti1 = np.zeros_like(so1); tc1 = np.zeros_like(so1)
    for v in range(n):
        for u in range(n):
            if u == v: continue
            hout = (adj[:, v] >> u) & 1
            hin = (adj[:, u] >> v) & 1
            to1[:, v] += hout * so1[:, u]
            ti1[:, v] += hin * si1[:, u]
            tc1[:, v] += hout * ci1[:, u]
    key = s.astype(np.int64)
    for extra, base in ((c3, 64), (so1, 64), (si1, 64), (so2, 512), (so3, 4096),
                        (co1, 256), (ci1, 256), (to1, 512), (ti1, 512), (tc1, 2048)):
        key = key * base + extra
    key.sort(axis=1)
    # arc-type sums (invariant scalars)
    a1 = np.zeros(B, dtype=np.int64); a2 = np.zeros(B, dtype=np.int64)
    for u in range(n):
        for v in range(n):
            if u == v: continue
            h = (adj[:, u] >> v) & 1
            common = POP[adj[:, u] & adj[:, v]]
            a1 += h * (s[:, u] * 8 + s[:, v]) * (common + 1)
            a2 += h * common * common
    # vectorized Hamiltonian-path count (level-by-level subset DP)
    adjbit = [[(adj[:, v] >> u) & 1 for u in range(n)] for v in range(n)]
    import itertools as _it
    levels = {1: {}}
    for v in range(n):
        arr = np.zeros(B, dtype=np.int64); arr += 0
        levels[1][(1 << v, v)] = None  # placeholder replaced below
    levels = {1: {(1 << v, v): np.ones(B, dtype=np.int64) for v in range(n)}}
    for size in range(1, n):
        nxt = {}
        for (mask, v), arr in levels[size].items():
            for u in range(n):
                bu = 1 << u
                if mask & bu: continue
                w = adjbit[v][u]
                tgt = (mask | bu, u)
                contrib = arr * w
                if tgt in nxt: nxt[tgt] += contrib
                else: nxt[tgt] = contrib.copy()
        del levels[size]
        levels[size + 1] = nxt
    Hv = np.zeros(B, dtype=np.int64)
    for (mask, v), arr in levels[n].items(): Hv += arr
    # deleted-vertex sub-H multiset (eight 7-vertex DPs), sorted per tiling
    subH = np.zeros((B, n), dtype=np.int64)
    for w in range(n):
        verts = [v for v in range(n) if v != w]
        lv = {(1 << v, v): np.ones(B, dtype=np.int64) for v in verts}
        for size in range(1, n-1):
            nxt = {}
            for (mask, v), arr in lv.items():
                for u in verts:
                    bu = 1 << u
                    if mask & bu: continue
                    contrib = arr * adjbit[v][u]
                    tgt = (mask | bu, u)
                    if tgt in nxt: nxt[tgt] += contrib
                    else: nxt[tgt] = contrib.copy()
            lv = nxt
        col = np.zeros(B, dtype=np.int64)
        for (mask, v), arr in lv.items(): col += arr
        subH[:, w] = col
    subH.sort(axis=1)
    return key, a1, a2, s, Hv, subH

# pass 1: bucket all tilings
buckets = {}
cls_of = np.zeros(N, dtype=np.int32)
rep = []
for lo in range(0, N, CH):
    bits = np.arange(lo, min(lo+CH, N), dtype=np.int64)
    key, a1, a2, s, Hv, subH = feats(bits)
    for j in range(bits.shape[0]):
        k = (key[j].tobytes(), int(a1[j]), int(a2[j]), int(Hv[j]), subH[j].tobytes())
        c = buckets.get(k)
        if c is None:
            c = len(buckets); buckets[k] = c; rep.append(int(bits[j]))
        cls_of[lo+j] = c
C = len(buckets)
print(f"n=8: {C} classes (A000568(8) = 6880): {'OK' if C == 6880 else 'FAIL'}")
assert C == 6880

# class data
def adj_of(b):
    adj = [0]*(n+1)
    for v in range(2, n+1): adj[v] |= 1 << (v-2)
    for i, (x, y) in enumerate(tiles):
        if b >> i & 1: adj[x] |= 1 << (y-1)
        else: adj[y] |= 1 << (x-1)
    return adj
def ham(adj):
    dp = defaultdict(int)
    for v in V: dp[(1 << (v-1), v)] = 1
    for _ in range(n-1):
        nd = defaultdict(int)
        for (mask, v), c in dp.items():
            for u in V:
                b = 1 << (u-1)
                if not mask & b and adj[v] & b: nd[(mask | b, u)] += c
        # merge sizes: keep incremental frontier only
        dp = nd if nd else dp
        if nd: last = nd
    return sum(c for (mask, v), c in last.items() if mask == (1 << n) - 1)

fib = np.bincount(cls_of, minlength=C)
# grid-sym per tiling, vectorized
bits_all = np.arange(N, dtype=np.int64)
gs = np.ones(N, dtype=bool)
for i in range(m):
    j = int(refl[i])
    if j < i: continue
    gs &= (((bits_all >> i) & 1) == ((bits_all >> j) & 1))
n_blue_tilings = int(gs.sum())
gs_per_class = np.zeros(C, dtype=np.int64)
np.add.at(gs_per_class, cls_of[gs], 1)
x_of = np.zeros(C, dtype=np.int64); sc_of = {}
H_of = np.zeros(C, dtype=np.int64)
for c in range(C):
    adj = adj_of(rep[c])
    sc = sorted(bin(adj[v]).count('1') for v in V)
    sc_of[c] = tuple(sc)
    x_of[c] = sum((2*si - (n-1))**2 for si in sc)
    H_of[c] = ham(adj)
assert all(H_of[c] % fib[c] == 0 for c in range(C)), "fiber*Aut=H fails"
typ = np.where(gs_per_class == fib, 2, np.where(gs_per_class == 0, 0, 1))  # 2=pb,1=mx,0=pk
print(f"blue tilings = {n_blue_tilings} (predicted 4096: {'OK' if n_blue_tilings == 4096 else 'FAIL'})")
print(f"pure blue {int((typ==2).sum())}, mixed {int((typ==1).sum())}, pure black {int((typ==0).sum())}")
tnames = {2: 'pureblue', 1: 'mixed', 0: 'pureblack'}
for t in (2, 1, 0):
    hist = defaultdict(int)
    for c in np.where(typ == t)[0]: hist[int(x_of[c])] += 1
    print(f"  {tnames[t]:<9} axis histogram: {dict(sorted(hist.items(), reverse=True))}")

# lines
partner = bits_all ^ (N - 1)
lo_side = bits_all < partner
la, lb = cls_of[bits_all[lo_side]], cls_of[partner[lo_side]]
lblue = gs[bits_all[lo_side]]
dx = np.abs(x_of[la] - x_of[lb])
print(f"blue lines = {int(lblue.sum())} (predicted 2048: {'OK' if int(lblue.sum()) == 2048 else 'FAIL'})")
bh = defaultdict(int); kh = defaultdict(int)
for d, isb in zip(dx, lblue):
    (bh if isb else kh)[int(d)] += 1
print(f"BLUE |dx| histogram: {dict(sorted(bh.items()))}")
print(f"  parity law (all = 0 mod 16): {'OK' if all(k % 16 == 0 for k in bh) else 'FAIL'}")
print(f"  max = {max(bh)} (predicted 48: {'OK' if max(bh) == 48 else 'FAIL'})")
print(f"BLACK |dx| histogram: {dict(sorted(kh.items()))}")
# type-pair line counts
tp = defaultdict(int)
for A, Bc, isb in zip(la, lb, lblue):
    key = ('BLUE' if isb else 'BLACK',) + tuple(sorted([tnames[int(typ[A])], tnames[int(typ[Bc])]]))
    tp[key] += 1
for k in sorted(tp): print(f"  {k}: {tp[k]}")
# blue never touches pure black?
viol = sum(v for k, v in tp.items() if k[0] == 'BLUE' and 'pureblack' in k[1:])
print(f"blue-avoids-pureblack: {'OK' if viol == 0 else 'FAIL'}")
# the transitive pipe
tc = int(np.where(H_of == 1)[0][0])
tl = [(('BLUE' if isb else 'BLACK'), int(x_of[Bc if A == tc else A]))
      for A, Bc, isb in zip(la, lb, lblue) if A == tc or Bc == tc]
print(f"transitive node x = {int(x_of[tc])}; lines: {tl} (predicted [('BLUE', 120)])")
