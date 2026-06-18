"""
Part 3: resolve the a3/a4 count discrepancy (claim a3:127, a4:4116; I got 34/2074
with V<220). Re-run with the claim's stated V<=400 scope, and confirm min M / 0 viol.
Also rigorously settle the SCOPE point: {1..12,m} is an S1 family (one speed>13),
so the AP theorem re-closes the already-proved S1 case, not the open S3 residual.
Finally: verify the two headline witness FLOORS are honest lower bounds on M
(witness gives min_v ||v tau|| which is <= M(S) by definition of M as a max over tau).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2): C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x * t) for x in S)
        if v > b: b = v
    return b

def is_cov(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def gcd_all(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g

print("="*70)
print("F. a3/a4 with wider scope V<=400 (claim a3:127, a4:4116)")
print("="*70)
# a3 scaling {t,2t,...,12t,V}
a3_count = 0; a3_min = None; a3_viol = 0
for t in range(2, 20):
    for V in range(14, 401):
        S = sorted(set([t*i for i in range(1,13)] + [V]))
        if len(S) != 13: continue
        if gcd_all(S) != 1: continue
        if not is_cov(S): continue
        above13 = [v for v in S if v > 13]
        if len(above13) < 2: continue
        Vmin = min(S); Vmax = max(S)
        if Vmax < 13*Vmin: continue
        a3_count += 1
        Mv = Mval(S)
        if a3_min is None or Mv < a3_min[0]: a3_min = (Mv, t, V)
        if Mv < F(1,14): a3_viol += 1
print(f"  a3 (t<20, V<=400): count={a3_count}, minM={a3_min[0] if a3_min else None}, viol={a3_viol}")
print(f"    (claim count 127 — scope-dependent; min M=1/13 & 0 viol are the load-bearing claims)")

a4_count = 0; a4_min = None; a4_viol = 0
for a in range(1, 14):
    for d in range(1, 14):
        ap = [a + d*i for i in range(12)]
        if min(ap) < 1: continue
        for V in range(14, 401):
            S = sorted(set(ap + [V]))
            if len(S) != 13: continue
            if gcd_all(S) != 1: continue
            if not is_cov(S): continue
            above13 = [v for v in S if v > 13]
            if len(above13) < 2: continue
            Vmin = min(S); Vmax = max(S)
            if Vmax < 13*Vmin: continue
            a4_count += 1
            Mv = Mval(S)
            if a4_min is None or Mv < a4_min[0]: a4_min = (Mv, a, d, V)
            if Mv < F(1,14): a4_viol += 1
print(f"  a4 (a,d<14, V<=400): count={a4_count}, minM={a4_min[0] if a4_min else None}, viol={a4_viol}")

print()
print("="*70)
print("G. SCOPE: case parameter k=#{v>13} for {1..12,m}")
print("="*70)
# Definitive: every {1..12,m} has exactly one element > 13 (m>=14). So k_param=1.
# The project's case split: S1 (k<=1) PROVED. So {1..12,m} in S1.
for m in [14, 26, 182, 364, 1000]:
    S = list(range(1,13)) + [m]
    above = [v for v in S if v > 13]
    print(f"  m={m}: speeds>13 = {above}  =>  k_param={len(above)} => S{'1' if len(above)<=1 else '3'}")
print("  CONCLUSION: AP theorem's {1..12,m} family is ENTIRELY within S1 (k<=1),")
print("  which was already PROVED. It does NOT touch the open S3 residual.")
print("  (The 'k=3 slice' VERIFICATION on board [14,30] DOES touch S3 but is finite/VERIFIED, not PROVED.)")

print()
print("="*70)
print("H. Are witness levels honest lower bounds on M? (definitional check)")
print("="*70)
# M(S) = max_tau min_v ||v tau||. For any fixed tau0, min_v ||v tau0|| <= M(S).
# So witness gives a valid lower bound. Confirm numerically: witness level <= Mval.
ok = True
for k in range(1, 40):
    m = 182*k
    S = list(range(1,13)) + [m]
    w1 = min(nrm(F(v)*F(2,27)) for v in S) if (m%27) not in (0,13,14) else None
    w2 = min(nrm(F(v)*F(14*k,182*k+1)) for v in S)
    Mv = Mval(S)
    best = w2 if w1 is None else max(w1, w2)
    if best > Mv:
        ok = False
        print(f"  !!! k={k}: witness {best} EXCEEDS Mval {Mv} — witness invalid!")
print(f"  All witness levels <= Mval (valid lower bounds) for k<40: {ok}")

print()
print("="*70)
print("I. Headline: is ANY covering {1..12,m} an LRC(14) counterexample? (M<1/14)")
print("="*70)
cex = []
for k in range(1, 60):
    m = 182*k
    S = list(range(1,13)) + [m]
    Mv = Mval(S)
    if Mv < F(1,14): cex.append((m, Mv))
print(f"  covering {{1..12,m}}, m=182k, k<60, with M<1/14: {len(cex)}")
print(f"  => headline 'no {{1..12,m}} is a counterexample' holds on tested range: {len(cex)==0}")
print(f"  Note {{1,...,13}} (m=13) has M=1/14 exactly but is NON-covering (no mult of 14): "
      f"M({{1..13}})={Mval(list(range(1,14)))}, covering={is_cov(list(range(1,14)))}")
