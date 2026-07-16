# opus-2026-07-16-S330 -- HYP-7065: (1) TRUE DFS search for both-clean
# 13-packets at small scale (Sidon + Y* >= Y0 for q,p <= 13 + no 13-mult);
# (2) if found: exact BONF5 -- possibly the first accessible-scale
# full-problem level-5 certificate.
from fractions import Fraction
from math import comb
import math, itertools, sys

Y0 = 30
QMAX = 20
LO, HI = 300, 40000   # S331 relaunch: wider window

def Ystar_ok(a, b):
    lo, hi = (a, b) if a < b else (b, a)
    for q in range(1, QMAX+1):
        p0 = round(q*hi/lo)
        for pp in (p0-1, p0, p0+1):
            if pp < 1: continue
            if abs(q*hi - pp*lo) < Y0: return False
    return True

def no3ap_ok(S, x):
    for s in S:
        if 2*s - x in S or 2*x - s in S or (x + s) % 2 == 0 and (x + s)//2 in S:
            return False
    return True

def linform_ok(S, x, cmax=3, thresh=25):
    # no small |c1*x + c2*a + c3*b| with coefficients in [-cmax, cmax]
    for a in S:
        for b in S:
            if b < a: continue
            for c1 in range(1, cmax+1):
                for c2 in range(-cmax, cmax+1):
                    for c3 in range(-cmax, cmax+1):
                        v = c1*x + c2*a + c3*b
                        if v != 0 and abs(v) < thresh: return False
                        if v == 0 and not (c2 == 0 and c3 == 0): return False
    return True

def dfs():
    sols = []
    cand_all = [x for x in range(LO, HI) if x % 13 != 0]
    sys.setrecursionlimit(10000)
    best_len = [0]
    def rec(S, sums, start):
        if len(S) == 13:
            sols.append(list(S)); return True
        if len(S) > best_len[0]:
            best_len[0] = len(S)
            print(f"   depth {len(S)}: {S}", flush=True)
        for x in range(start, HI):
            if x % 13 == 0: continue
            ok = True
            for s in S:
                if x + s in sums: ok = False; break
            if ok:
                for s in S:
                    if not Ystar_ok(s, x): ok = False; break
            if ok and not no3ap_ok(S, x): ok = False
            if ok and not linform_ok(S, x): ok = False
            if ok:
                ns = sums | {x + s for s in S}
                if rec(S + [x], ns, x + Y0):
                    return True
        return False
    rec([], set(), LO)
    return sols

print(f"(1) DFS both-clean search: Y0 = {Y0}, window [{LO}, {HI}]")
sols = dfs()
if not sols:
    print("   NO both-clean 13-packet in the window (exhaustive DFS)")
else:
    S = sols[0]
    print(f"   FOUND: {S}")
    # (2) exact BONF5
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
    print(f"   BONF5 >= {float(bound):+.6f}   actual = {float(mu(U)):.6f}   "
          f"{'*** THE FIRST ACCESSIBLE-SCALE FULL-PROBLEM CERTIFICATE ***' if bound > 0 else 'not coercive (q,p <= 13 cleanliness insufficient)'}")
