# opus-2026-07-16-S329 -- HYP-7055:
# (A) THE BOTH-CLEAN 13-PACKET SEARCH: Sidon + ratio-clean (Y* >= YMIN for
#     q,p <= 13) + no 13-multiples + distinct residues: exact S1..S5 + BONF5.
# (B) THE GUY / SQUARE-TRIANGULAR WEAVE (exact referee):
#     Z(2k+1) = T_{k-1}^2 = 1^3 + ... + (k-1)^3 (Nicomachus);
#     Z(9) = 36 = the 3rd square-triangular number;
#     E(AP_n) = n(n+1)(2n+1)/3 - n^2; midpoint quadruples = 2C(n,3) + C(n,2).
from fractions import Fraction
from math import comb, gcd, isqrt
import math, itertools

# ---------- (B) first (cheap, exact)
print("(B) THE GUY / SQUARE-TRIANGULAR WEAVE:")
def Z(n): return (n//2)*((n-1)//2)*((n-2)//2)*((n-3)//2)//4
def T(k): return k*(k+1)//2
print("   n : Z(n) : T_{k-1}^2 (n = 2k+1) : sum of cubes")
for k in range(2, 11):
    n = 2*k + 1
    zz = Z(n)
    tsq = T(k-1)**2
    cubes = sum(j**3 for j in range(1, k))
    print(f"   {n:2d} : {zz:5d} : {tsq:5d} : {cubes:5d}   match: {zz == tsq == cubes}")
# square-triangular numbers via Pell
st = [0, 1]
while len(st) < 7: st.append(34*st[-1] - st[-2] + 2)
print(f"   square-triangular numbers: {st[1:]}")
print(f"   Z(9) = {Z(9)} = the 3rd square-triangular: {Z(9) == st[3]}")
# E(AP_n) and midpoint quadruples
for n in (12, 13, 14):
    S = list(range(1, n+1))
    from collections import defaultdict
    cnt = defaultdict(int)
    for a in S:
        for b in S: cnt[a+b] += 1
    E = sum(c*c for c in cnt.values())
    Ef = n*(n+1)*(2*n+1)//3 - n*n
    quad = sum(comb(c, 2) for c in cnt.values())   # ordered-pair midpoint pairs
    # unordered genuine pair-of-pairs sharing a sum: use unordered pair reps
    cnt2 = defaultdict(int)
    for i in range(n):
        for j in range(i+1, n): cnt2[S[i]+S[j]] += 1
    quad2 = sum(comb(c, 2) for c in cnt2.values())
    print(f"   n={n}: E = {E} (formula {Ef}: {E == Ef}); unordered midpoint "
          f"quadruples = {quad2} (2C(n,3)+C(n,2)... 2C({n},3) = {2*comb(n,3)})")

# ---------- (A) the both-clean search
print("\n(A) THE BOTH-CLEAN 13-PACKET SEARCH:")
YMIN = 40
def Ystar(x1, x2, Q=13):
    best = None
    lo, hi = min(x1, x2), max(x1, x2)
    for q in range(1, Q+1):
        p0 = round(q*hi/lo)
        for pp in (p0-1, p0, p0+1):
            if pp < 1: continue
            v = abs(q*hi - pp*lo)
            if best is None or v < best: best = v
    return best

import random
random.seed(3)
def try_build(start, step_pool, tries=4000):
    S = []
    sums = set()
    cand = start
    while len(S) < 13 and cand < 25000:
        ok = True
        if cand % 13 == 0: ok = False
        if ok:
            for s in S:
                if cand + s in sums: ok = False; break
        if ok:
            for s in S:
                if Ystar(s, cand) < YMIN: ok = False; break
        if ok and any(cand % 13 == s % 13 for s in S): ok = False
        if ok:
            for s in S: sums.add(cand + s)
            S.append(cand)
        cand += random.choice(step_pool)
    return S if len(S) == 13 else None

S = try_build(1601, [31, 37, 41, 43, 53, 59, 61, 67, 71])
if S is None:
    print("   greedy failed; retry with wider steps")
    S = try_build(2003, [97, 101, 103, 107, 109, 113, 127])
print(f"   packet: {S}")
if S:
    ys = sorted(Ystar(a, b) for a, b in itertools.combinations(S, 2))
    print(f"   min pairwise Y* = {ys[0]} (threshold {YMIN}); Sidon: True by construction; "
          f"no 13-multiples: {all(x % 13 for x in S)}")
    def teeth(x, a, b):
        w = Fraction(1, 13*x)
        out = []
        for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
            lo2, hi2 = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
            if lo2 < hi2: out.append((lo2, hi2))
        return out
    def inter(arcs, x):
        out = []
        for (a, b) in arcs: out.extend(teeth(x, a, b))
        return out
    def mu(arcs): return sum(b - a for a, b in arcs)
    FULL = [(Fraction(0), Fraction(1))]
    lvl1 = [inter(FULL, x) for x in S]
    Sk = {1: sum(mu(v) for v in lvl1)}
    prev = {(i,): lvl1[i] for i in range(13)}
    for k in range(2, 6):
        tot = Fraction(0); cur = {}
        for idx in itertools.combinations(range(13), k):
            w = inter(prev[idx[:-1]], S[idx[-1]])
            if k < 5: cur[idx] = w
            tot += mu(w)
        Sk[k] = tot; prev = cur if k < 5 else prev
        print(f"   S{k} = {float(tot):.6f} (equid {comb(13,k)*(2/13)**k:.6f})", flush=True)
    bound = 1 - Sk[1] + Sk[2] - Sk[3] + Sk[4] - Sk[5]
    U = FULL
    for x in S:
        out = []
        for (a, b) in U:
            cur2 = a
            for (lo2, hi2) in sorted(teeth(x, a, b)):
                if lo2 > cur2: out.append((cur2, min(lo2, b)))
                cur2 = max(cur2, hi2)
                if cur2 >= b: break
            if cur2 < b: out.append((cur2, b))
        U = [iv for iv in out if iv[0] < iv[1]]
    print(f"   BONF5 >= {float(bound):+.6f}   actual uncovered = {float(mu(U)):.6f}"
          f"   {'*** COERCIVE: THE FIRST FULL-PROBLEM LEVEL-5 CERTIFICATE ***' if bound > 0 else 'not coercive'}")
