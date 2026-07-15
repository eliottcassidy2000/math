#!/usr/bin/env python3
"""orbit_decode_n10_burnside_kps_S128c19.py -- kind-pasteur S128 cont.18.
n=9 via BURNSIDE-ENUMERATION: for each of 9! perms, solve the affine GF(2) system
'sigma T(t) = T(kappa t)' (28 vars, 36 eqs, bit-packed) and ENUMERATE its solutions.
Union (with multiplicity = |Aut|) = the quasi-fixed set X(9) -- never touching 2^28.
Then: dedupe/validate, Klein orbits, SC-type each rep. Also corrected NS counts."""
import sys
from math import comb
from itertools import permutations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)
n = 10
tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
tidx = {t: i for i, t in enumerate(tiles)}
m = len(tiles)
gmap = [tidx[(n - y + 1, n - x + 1)] for (x, y) in tiles]
pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
def slot(a, b):
    if b - a == 1:
        return None
    return tidx[(b + 1, a + 1)]
print("m=%d; sweeping %d perms..." % (m, 3628800))
sols = Counter()
np_ = 0
for sig in permutations(range(n)):
    np_ += 1
    rows = []   # each row: (mask, rhs)
    ok = True
    for (u, v) in pairs:
        p, q = sig[u], sig[v]
        uu, vv = (u, v) if p < q else (v, u)
        if p > q:
            p, q = q, p
        a, b = (uu, vv) if uu < vv else (vv, uu)
        e = slot(a, b)
        if e is None:
            f1c, f1m = 0, 0
        else:
            f1c, f1m = 1, 1 << e
        if uu > vv:
            f1c ^= 1
        e2 = slot(p, q)
        if e2 is None:
            f2c, f2m = 0, 0
        else:
            f2c, f2m = 0, 1 << e2
        mask = f1m ^ f2m
        rhs = f1c ^ f2c
        if mask:
            rows.append((mask, rhs))
        elif rhs:
            ok = False
            break
    if not ok:
        continue
    # gauss (bit-packed)
    piv = {}
    okc = True
    for mask, rhs in rows:
        while mask:
            hb = mask.bit_length() - 1
            if hb in piv:
                pm, pr = piv[hb]
                mask ^= pm
                rhs ^= pr
            else:
                piv[hb] = (mask, rhs)
                break
        else:
            if rhs:
                okc = False
                break
    if not okc:
        continue
    rank = len(piv)
    dim = m - rank
    if dim > 18:
        print("  dim %d too big at perm %d -- abort" % (dim, np_)); sys.exit(1)
    # back-substitute: reduce pivots to rref
    keys = sorted(piv, reverse=True)
    rref = []
    for k in keys:
        mask, rhs = piv[k]
        for k2, (m2, r2) in list(piv.items()):
            pass
    # simpler: iterate to clear (small ranks fine)
    changed = True
    pv = dict(piv)
    for k in sorted(pv, reverse=True):
        mask, rhs = pv[k]
        for k2 in sorted(pv, reverse=True):
            if k2 < k and (mask >> k2) & 1:
                m2, r2 = pv[k2]
                mask ^= m2
                rhs ^= r2
        pv[k] = (mask, rhs)
    freecols = [c for c in range(m) if c not in pv]
    # particular solution: free = 0 -> x[k] = rhs for pivot cols
    base = 0
    for k, (mask, rhs) in pv.items():
        if rhs:
            base |= 1 << k
    # nullspace basis: for each free col f: x_f = 1, pivots adjust
    basis = []
    for f in freecols:
        vvec = 1 << f
        for k, (mask, rhs) in pv.items():
            if (mask >> f) & 1:
                vvec |= 1 << k
        basis.append(vvec)
    for sel in range(1 << dim):
        x = base
        s = sel
        i = 0
        while s:
            if s & 1:
                x ^= basis[i]
            s >>= 1
            i += 1
        sols[x] += 1
    if np_ % 400000 == 0:
        print("  ...%d perms, distinct sols %d, mass %d" % (np_, len(sols), sum(sols.values())), flush=True)
W = sum(sols.values())
X = list(sols)
mult = Counter(sols.values())
print("W(9) = %d ; |X|(9) = %d ; multiplicity histogram (|Aut| profile) = %s" % (W, len(X), dict(mult)))
gs = sum(1 for t in X if all(((t >> i) & 1) == ((t >> gmap[i]) & 1) for i in range(m)))
print("gridsym-qf at n=10: %d (even n: expect > 0)" % gs)
# Klein orbits + SC typing
def beats(tv):
    B = [[False] * n for _ in range(n)]
    for k in range(2, n + 1):
        B[k - 1][k - 2] = True
    for (x, y), i in tidx.items():
        if tv[i]:
            B[x - 1][y - 1] = True
        else:
            B[y - 1][x - 1] = True
    return B
def iso(B1, B2):
    d1 = [sum(B1[u]) for u in range(n)]
    d2 = [sum(B2[u]) for u in range(n)]
    if sorted(d1) != sorted(d2):
        return False
    buckets = defaultdict(list)
    for v in range(n):
        buckets[d2[v]].append(v)
    order = sorted(range(n), key=lambda u: (d1[u], u))
    assign = [-1] * n
    used = [False] * n
    def bt(i):
        if i == n:
            return True
        u = order[i]
        for v in buckets[d1[u]]:
            if used[v]:
                continue
            good = True
            for j in range(i):
                w = order[j]
                if B1[u][w] != B2[v][assign[w]] or B1[w][u] != B2[assign[w]][v]:
                    good = False
                    break
            if good:
                assign[u] = v; used[v] = True
                if bt(i + 1):
                    return True
                used[v] = False; assign[u] = -1
        return False
    return bt(0)
seen = set(); oSC = oNS = 0
Xs = set(X)
Xngs = [t for t in X if not all(((t >> i) & 1) == ((t >> gmap[i]) & 1) for i in range(m))]
print("|X_ngs|(10) = %d ; selfK(10) = %d ; orbits target = %d" % (len(Xngs), len(Xngs)//2, len(Xngs)//4))
for t in Xngs:
    if t in seen:
        continue
    tv = [(t >> i) & 1 for i in range(m)]
    gt = sum(tv[gmap[i]] << i for i in range(m))
    kt = t ^ ((1 << m) - 1); gkt = gt ^ ((1 << m) - 1)
    seen |= {t, gt, kt, gkt}
    B1 = beats(tv)
    BT = [[B1[u][v] for u in range(n)] for v in range(n)]
    if iso(B1, BT):
        oSC += 1
    else:
        oNS += 1
print("n=10: orbits=%d = orbitsSC %d + orbitsNS %d" % (oSC + oNS, oSC, oNS))
print("SEQUENCES: orbitsSC = 2,0,9,3,69,%d ; orbitsNS = 0,3,13,98,783,%d" % (oSC, oNS))
NS10 = (9733056 - 8784) // 2
print("NS rates incl n=10: 0.136, 0.071, 0.029, 0.0083, %d/%d=%.5f"
      % (oNS, NS10, oNS / NS10))
print("DONE")
