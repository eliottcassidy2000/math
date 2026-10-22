#!/usr/bin/env python3
"""(b) Paper 2: census of formal ADE multisets by rank and exhaustive 4-moment collision audit.
Signature s = (c1, c2, 3*c3, c4) (integers) from the printed formulas (equivalently the
orthogonal-frame counts O_1..O_4). Unbounded-knapsack DP with numpy arrays, per rank n."""
import sys, time
import numpy as np
from fractions import Fraction as Fr
t0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 43

def cA(k): return (k*(k+1), -k*(k+1)*(2*k-1), 4*k*(k+1)*(5*k**2-7*k+3), -2*k*(k+1)*(14*k**3-37*k**2+39*k-15))
def cD(k): return (2*k*(k-1), -2*k*(k-1)*(4*k-7), 8*k*(k-1)*(20*k**2-74*k+69), -4*k*(k-1)*(112*k**3-672*k**2+1378*k-963))
cE = {6: (72, -1512, 3*50976, -2011824), 7: (126, -4158, 3*223272, -14196924), 8: (240, -13680, 3*1281600, -143445600)}
# 3*c3 for A: 4k(k+1)(...) ; for D: 8k(k-1)(...)  (already multiplied by 3 above)

def alphabet(nmax, which='ADE'):
    T = [(f"A{k}", k, cA(k)) for k in range(1, nmax+1)]
    if which == 'ADE':
        T += [(f"D{k}", k, cD(k)) for k in range(4, nmax+1)]
        T += [(f"E{k}", k, cE[k]) for k in (6, 7, 8) if k <= nmax]
    return T

# ---- census from the generating function
def census(nmax):
    p = [1] + [0]*nmax
    parts = list(range(1, nmax+1)) + list(range(4, nmax+1)) + [6, 7, 8]
    for r in parts:
        for n in range(r, nmax+1): p[n] += p[n-r]
    return p
p = census(NMAX)
print("census p_n (n=1..%d):" % NMAX, p[1:])
print("sum_{n<=42} p_n =", sum(p[1:43]) if NMAX >= 42 else None, " p_43 =", p[43] if NMAX >= 43 else None,
      " sum_{n<=43} =", sum(p[1:44]) if NMAX >= 43 else None)
printed_table = {1:1,2:2,3:3,4:6,5:9,6:16,7:24,8:39,9:57,10:88,11:128,12:193,13:276,14:403,15:570,16:815,17:1137,18:1599,
 19:2207,20:3058,21:4177,22:5713,23:7728,24:10458,25:14021,26:18786,27:24985,28:33186,29:43801,30:57721,31:75655,32:98981,
 33:128897,34:167524,35:216842,36:280117,37:360537,38:463108,39:592930,40:757618,41:965192,42:1227224,43:1556185}
mism = [n for n in printed_table if n <= NMAX and printed_table[n] != p[n]]
print("printed Table 3.2 counts agree with generating function:", not mism, mism)

def audit(nmax, which='ADE', M=4, report=True):
    T = alphabet(nmax, which)
    A = [np.zeros((0, M), dtype=np.int64) for _ in range(nmax+1)]
    A[0] = np.zeros((1, M), dtype=np.int64)
    for name, r, v in sorted(T, key=lambda x: -x[1]):
        vv = np.array(v[:M], dtype=np.int64)
        for n in range(r, nmax+1):
            if len(A[n-r]):
                A[n] = np.concatenate([A[n], A[n-r] + vv])
    first = None
    res = {}
    for n in range(1, nmax+1):
        u, cnt = np.unique(A[n], axis=0, return_counts=True)
        dup = u[cnt > 1]
        res[n] = (len(A[n]), len(u), int((cnt > 1).sum()), int(cnt.max()), dup)
        if first is None and len(dup): first = n
    return res, first

if __name__ == '__main__':
    res, first = audit(NMAX, 'ADE', 4)
    tot42 = 0
    for n in range(1, NMAX+1):
        size, nu, ndup, mx, dup = res[n]
        assert size == p[n], (n, size, p[n])
        if n <= 42: tot42 += size
        if ndup or n >= 40:
            print(f"  rank {n}: {size} multisets, {nu} distinct signatures, {ndup} colliding signatures, max multiplicity {mx}")
    print("total multisets ranks 1..42 enumerated:", tot42)
    print("first rank with a 4-moment collision:", first)
    if NMAX >= 43:
        print("rank-43 colliding signature(s) (c1,c2,3c3,c4):", res[43][4].tolist())
    print(f"(ADE audit time {time.time()-t0:.1f}s)")
    # smaller-M thresholds and other alphabets
    for which, Mm, nm in [('ADE', 1, 12), ('ADE', 2, 20), ('ADE', 3, 30), ('A', 1, 12), ('A', 2, 20), ('A', 3, 40), ('A', 4, 60)]:
        r2, f2 = audit(nm, which, Mm)
        print(f"  alphabet {which}, M={Mm}: minimal collision rank = {f2}; colliding signatures there = {r2[f2][4].tolist() if f2 else None}")
    print(f"(total time {time.time()-t0:.1f}s)")
