"""
LRC(14) -- the rho* density route (THM-527 / OPEN-Q-108): VERIFICATION + REFUTATION
===================================================================================
Self-contained, exact-rational. Independent of the repo's scripts.

rho*(P,E) = meas{ x in [0,1) : ||p x|| >= 1/14 for all p in P
                            AND maxgap{ frac(e x) : e in E } > 2/7 }
THM-527: rho*(P,E) > 0  ==>  M(S) >= 1/14  (SUFFICIENT) for S = P u {Vmax - e}.

Route's goal (THM-527 part G item 1 / OPEN-Q-108):  inf rho* >= c0 > 0  ==> LRC(14).

RESULTS proved/verified below:
 [1] Reproduce exactly the canon values mu_k and the consecutive floor 1/84.
 [2] Unconditional union bound: rho* >= mu(E) - |P|/7  ==>  rho* > 0 for k >= 12.
 [3] REFUTATION: inf rho* = 0. Explicit ADMISSIBLE primitive covering 13-sets
     (multiple of every q in 2..14) with rho* = 0.
 [4] Those sets are still LONELY (M ~ 1/9 >> 1/14): rho*=0 is a BLIND SPOT of the
     criterion, NOT a counterexample to LRC(14). Witness denominators are
     unrelated to Vmax (consistent with THM-566: no uniform bounded denominator).
Conclusion: the "lower-bound the density rho*" strategy cannot prove LRC(14).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce

# ----------------------------- core exact machinery -----------------------------
def frac(x): return x - (x.numerator // x.denominator)

def circ_maxgap(ph):
    pts = sorted(set(ph))
    if len(pts) == 1: return F(1)
    g = [pts[i+1]-pts[i] for i in range(len(pts)-1)]; g.append(F(1)-pts[-1]+pts[0])
    return max(g)

def breakpoints(E, P):
    bps = {F(0), F(1)}; Es = list(E)
    for a in range(len(Es)):
        for b in range(a+1, len(Es)):
            d = abs(Es[a]-Es[b])
            if d == 0: continue
            for c in range(d+1):
                bps.add(F(c, d)); bps.add(F(7*c+2, 7*d)); bps.add(F(7*c-2, 7*d))
    for e in Es:
        if e:
            for c in range(e+1): bps.add(F(c, e))
    for p in P:
        for a in range(p+1):
            bps.add(F(14*a+1, 14*p)); bps.add(F(14*a-1, 14*p))
    return sorted(b for b in bps if 0 <= b <= 1)

def rho_star(P, E):
    bps = breakpoints(E, P); tot = F(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi <= lo: continue
        mid = (lo+hi)/2; ok = True
        for p in P:
            f = frac(p*mid)
            if min(f, 1-f) < F(1, 14): ok = False; break
        if ok and circ_maxgap([frac(e*mid) for e in E]) > F(2, 7): tot += hi-lo
    return tot

def norm(x):
    f = x - (x.numerator // x.denominator); return min(f, 1-f)

def M_exact(S):
    cands = set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j], abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): cands.add(F(m, d))
        for m in range(2*S[i]+1): cands.add(F(m, 2*S[i]))
    best = F(0); arg = F(0)
    for t in cands:
        if 0 < t < 1:
            v = min(norm(s*t) for s in S)
            if v > best: best, arg = v, t
    return best, arg

def is_covering(S): return all(any(s % q == 0 for s in S) for q in range(2, 15))
def primitive(S):   return reduce(gcd, S) == 1

# ----------------------------------- [1] canon -----------------------------------
print("[1] EXACT canon check: mu_k = rho*(empty, {0..k-1})")
claim = {3:F(1), 4:F(19,21), 5:F(9,14), 13:F(829,4620)}
ok1 = True
for k in sorted(claim):
    mu = rho_star(set(), list(range(k)))
    m = (mu == claim[k]); ok1 &= m
    print(f"    mu_{k:<2} = {str(mu):>10}  claim {claim[k]}  {'OK' if m else 'MISMATCH'}")
floor_consec = min((rho_star(set(P), list(range(k))), k, P)
                   for k in range(3,14) for P in combinations(range(1,14), 13-k))
print(f"    consecutive floor = {floor_consec[0]} at k={floor_consec[1]}, P={floor_consec[2]}  "
      f"(claim 1/84: {'OK' if floor_consec[0]==F(1,84) else 'MISMATCH'})")

# ------------------------------- [2] union bound --------------------------------
print("\n[2] UNCONDITIONAL union bound rho* >= mu(E) - |P|/7  (proves rho*>0 for k>=12):")
for k in (11, 12, 13):
    mu = rho_star(set(), list(range(k)))
    lb = mu - F(13-k, 7)
    print(f"    k={k}: mu_consec={float(mu):.4f} - {13-k}/7 = {lb} = {float(lb):.4f}"
          f"  -> {'PROVES rho*>0' if lb>0 else 'inconclusive'}")

# --------------------------- [3]+[4] refutation ----------------------------------
print("\n[3]+[4] REFUTATION: admissible primitive covering 13-sets with rho*=0, yet lonely")
shapes = [((1,2,3,12),[0,2,3,4,5,6,7,8,10]),
          ((1,2,3,6), [0,1,2,4,5,6,7,8,10]),
          ((1,2,3),   [0,1,2,3,4,5,6,7,8,10])]
allgood = True
for P, E in shapes:
    r = rho_star(set(P), E)
    S = None
    for Vmax in range(20, 5000):
        cl = [Vmax-e for e in E]
        cand = sorted(set(P) | set(cl))
        if min(cl) > 13 and len(cand) == 13 and sorted(s for s in cand if s<=13)==sorted(P) \
           and primitive(cand) and is_covering(cand):
            S = cand; break
    M, t = M_exact(S)
    lonely = M >= F(1,14)
    allgood &= (r == 0 and S is not None and lonely)
    print(f"    P={P} E={E}")
    print(f"       rho*={r}  |  covering set S={S}")
    print(f"       primitive={primitive(S)} covering={is_covering(S)} | "
          f"M(S)={M}={float(M):.4f} at t={t} (>=1/14: {lonely})")
print(f"\n==> inf rho* = 0 (route's c0>0 conjecture REFUTED); sets remain lonely. "
      f"All checks pass: {ok1 and floor_consec[0]==F(1,84) and allgood}")
