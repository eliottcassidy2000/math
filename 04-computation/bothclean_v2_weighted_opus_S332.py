# opus-2026-07-16-S332 -- HYP-7135 / THM-928 (C): THE CORRECTED BOTH-CLEAN
# DFS (post-MISTAKE-151). Filter fixes vs S329-S331:
#   (F1) WEIGHTED CAPPED pair cleanliness (THM-864 shape):
#        min over 1<=q,p<=13 of |q*b - p*a| * (q+p-1)  >=  THETA.
#        Both coordinates capped (no beat-Dirichlet vacuity); the weight
#        (q+p-1) matches THM-864's error 13*kappa*rho/(y*(p+q-1)); the
#        gcd-quantization loophole dies because q=p=1 demands |b-a| >= THETA.
#   (F2) Sidon (distinct pairwise sums) -- additive quadruples out.
#   (F3) no-3AP -- the q=2,(r,s)=(1,1) exact triple resonance out.
#   (F4) line-clean at height 3: no exact triple relation q*c = r*a + s*b
#        with 1<=q<=3, |r|,|s|<=3 (mac-mini THM-926: all-nonzero resonance
#        lines decay like 1/|k1 k2 k3| -- only small heights matter).
#   (F5) no 13-multiples, distinct residues mod 13 not required (too strong).
# On any 13-packet found: exact BONF5 (S1..S5 + true uncovered measure).
from fractions import Fraction
from math import comb
import math, itertools, sys

THETA = 100
LO, HI = 300, 40000

def weighted_clean(a, b):
    lo, hi = (a, b) if a < b else (b, a)
    for q in range(1, 14):
        for p in range(1, 14):
            if abs(q * hi - p * lo) * (q + p - 1) < THETA:
                return False
    return True

def line_clean(S, x):
    for a, b in itertools.combinations(S, 2):
        for qq in range(1, 4):
            for r in range(-3, 4):
                for s in range(-3, 4):
                    if r == 0 and s == 0: continue
                    if qq * x == r * a + s * b: return False
                    # x may also play the a/b role vs pairs incl itself later;
                    # symmetric roles: r*x + s*a = qq*b etc.
                    if qq * b == r * a + s * x: return False
                    if qq * a == r * b + s * x: return False
    return True

def no3ap(S, x):
    ss = set(S)
    for s in S:
        if 2 * s - x in ss or 2 * x - s in ss: return False
        if (x + s) % 2 == 0 and (x + s) // 2 in ss: return False
    return True

def dfs():
    sys.setrecursionlimit(10000)
    best = [0]
    def rec(S, sums, start):
        if len(S) == 13: return list(S)
        if len(S) > best[0]:
            best[0] = len(S)
            print(f"   depth {len(S)}: {S}", flush=True)
        x = start
        while x < HI:
            if x % 13 != 0:
                ok = all(x + s not in sums for s in S)
                if ok: ok = no3ap(S, x)
                if ok: ok = all(weighted_clean(s, x) for s in S)
                if ok: ok = line_clean(S, x)
                if ok:
                    r = rec(S + [x], sums | {x + s for s in S}, x + 1)
                    if r: return r
            x += 1
        return None
    return rec([], set(), LO)

print(f"(1) CORRECTED both-clean DFS: THETA = {THETA} (weighted capped), "
      f"window [{LO}, {HI}]", flush=True)
S = dfs()
if not S:
    print("   NO weighted-clean 13-packet in the window (exhaustive DFS)")
else:
    print(f"   FOUND: {S}")
    def teeth_iv(x, a, b):
        w = Fraction(1, 13*x)
        out = []
        for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
            lo2, hi2 = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
            if lo2 < hi2: out.append((lo2, hi2))
        return out
    def inter(arcs, x):
        out = []
        for (a, b) in arcs: out.extend(teeth_iv(x, a, b))
        return out
    def muv(arcs): return sum(b - a for a, b in arcs)
    FULL = [(Fraction(0), Fraction(1))]
    lvl1 = [inter(FULL, x) for x in S]
    Sk = {1: sum(muv(v) for v in lvl1)}
    prev = {(i,): lvl1[i] for i in range(13)}
    for k in range(2, 6):
        tot = Fraction(0); cur = {}
        for idx in itertools.combinations(range(13), k):
            w = inter(prev[idx[:-1]], S[idx[-1]])
            if k < 5: cur[idx] = w
            tot += muv(w)
        Sk[k] = tot; prev = cur if k < 5 else prev
        print(f"   S{k} = {float(tot):.6f} (equid {comb(13,k)*(2/13)**k:.6f})",
              flush=True)
    bound = 1 - Sk[1] + Sk[2] - Sk[3] + Sk[4] - Sk[5]
    U = FULL
    for x in S:
        out = []
        for (a, b) in U:
            cur2 = a
            for (lo2, hi2) in sorted(teeth_iv(x, a, b)):
                if lo2 > cur2: out.append((cur2, min(lo2, b)))
                cur2 = max(cur2, hi2)
                if cur2 >= b: break
            if cur2 < b: out.append((cur2, b))
        U = [iv for iv in out if iv[0] < iv[1]]
    verdict = ('*** COERCIVE: THE FIRST ACCESSIBLE-SCALE FULL-PROBLEM '
               'LEVEL-5 CERTIFICATE ***' if bound > 0 else
               'not coercive (weighted cleanliness at THETA=100 still insufficient)')
    print(f"   BONF5 >= {float(bound):+.6f}   actual uncovered = "
          f"{float(muv(U)):.6f}   {verdict}")
