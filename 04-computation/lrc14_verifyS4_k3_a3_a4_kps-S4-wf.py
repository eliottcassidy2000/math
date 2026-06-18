"""
Part 2: stress the FRAGILE claims of the closure:
 - 'no wraparound' inequality 12*(m/13) < m+1 (load-bearing for witness 2).
 - a3 scaling family {t,2t,...,12t,V}: claim min M = 1/13, 0 violations.
 - a4 general AP {a,a+d,...,a+11d,V}: claim min M = 1/13, 0 violations.
 - k=3 board [14,30]: claim 11632 covering primitive S3 sets, min M = 2/23.
 - The reported k=3 worst set [1,2,3,4,5,7,8,11,12,13,18,20,28] M=2/23.
 - Check whether the AP-theorem scope claim "closes every set {1..12,m}" is
   really airtight, and whether any covering {1..12,m} is actually a 13-set
   counterexample (the headline).
Also: re-examine the SCOPE. Does {1,...,12,m} being closed actually mean
anything for LRC(14)? It must be a primitive covering 13-set in S3. Check.
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

base = list(range(1, 13))

print("="*70)
print("A. 'No wraparound' inequality 12*(m/13) < m+1 for witness 2 (m=182k)")
print("="*70)
allok = True
for k in range(1, 5000):
    m = 182 * k
    q = m // 13  # 14k
    if not (12 * q < m + 1):
        allok = False
        print("  FAIL at k=", k)
# algebraic: 12q = 12*14k = 168k ; m+1 = 182k+1 ; 168k < 182k+1 always for k>=0
print(f"  12*(m/13) < m+1 for all k in [1,5000): {allok}")
print(f"  algebraic margin (m+1) - 12q = (182k+1) - 168k = 14k+1 > 0 always.")

print()
print("="*70)
print("B. Is {1,...,12,m} (m=182k) a primitive covering S3 set? (scope check)")
print("="*70)
# Primitive: gcd of all = 1 (1 is in the set, so always primitive).
# S3: k>=2 (>=2 speeds > 13) AND Vmax >= 13*Vmin.
# But {1..12,m} has only ONE speed > 13 (namely m). So #{v>13} = 1 => k=1 => this is S1!
for k in [1, 2, 3]:
    m = 182 * k
    S = base + [m]
    above13 = [v for v in S if v > 13]
    print(f"  m={m}: #{{v>13}}={len(above13)} -> k_caseparam={len(above13)}; "
          f"this is S{'1' if len(above13)<=1 else '3'}, gcd={gcd_all(S)}")
print("  NOTE: {1..12,m} has exactly ONE speed >13 (m). So k_param=1 => it is in S1,")
print("        which was ALREADY PROVED. The AP theorem closes an S1 family, not S3.")

print()
print("="*70)
print("C. a3 scaling family {t,2t,...,12t,V}: claim min M=1/13, 0 violations")
print("="*70)
# Enumerate t in some range, V coprime-ish, covering primitive S3, exact M.
a3_min = None; a3_count = 0; a3_viol = []
for t in range(2, 12):
    for V in range(14, 220):
        S = sorted(set([t*i for i in range(1,13)] + [V]))
        if len(S) != 13: continue
        if gcd_all(S) != 1: continue
        if not is_cov(S): continue
        above13 = [v for v in S if v > 13]
        if len(above13) < 2: continue  # need S3 (k>=2)
        Vmin = min(S); Vmax = max(S)
        if Vmax < 13*Vmin: continue     # S3 residual condition
        a3_count += 1
        Mv = Mval(S)
        if a3_min is None or Mv < a3_min[0]:
            a3_min = (Mv, t, V, S)
        if Mv < F(1,14):
            a3_viol.append((t, V, Mv, S))
print(f"  a3 covering primitive S3 sets found: {a3_count}")
if a3_min:
    print(f"  min M = {a3_min[0]} (={float(a3_min[0]):.5f}) at t={a3_min[1]}, V={a3_min[2]}")
print(f"  violations (M<1/14): {len(a3_viol)}")
if a3_viol: print("   ", a3_viol[:5])

print()
print("="*70)
print("D. a4 general AP {a,a+d,...,a+11d,V}: claim min M=1/13, 0 violations")
print("="*70)
a4_min = None; a4_count = 0; a4_viol = []
for a in range(1, 12):
    for d in range(1, 12):
        ap = [a + d*i for i in range(12)]
        for V in range(14, 220):
            S = sorted(set(ap + [V]))
            if len(S) != 13: continue
            if min(S) < 1: continue
            if gcd_all(S) != 1: continue
            if not is_cov(S): continue
            above13 = [v for v in S if v > 13]
            if len(above13) < 2: continue
            Vmin = min(S); Vmax = max(S)
            if Vmax < 13*Vmin: continue
            a4_count += 1
            Mv = Mval(S)
            if a4_min is None or Mv < a4_min[0]:
                a4_min = (Mv, a, d, V, S)
            if Mv < F(1,14):
                a4_viol.append((a, d, V, Mv, S))
print(f"  a4 covering primitive S3 sets found: {a4_count}")
if a4_min:
    print(f"  min M = {a4_min[0]} (={float(a4_min[0]):.5f}) at a={a4_min[1]}, d={a4_min[2]}, V={a4_min[3]}")
print(f"  violations (M<1/14): {len(a4_viol)}")
if a4_viol: print("   ", a4_viol[:5])

print()
print("="*70)
print("E. k=3 board [14,30]: enumerate covering primitive S3 k=3 sets, min M")
print("="*70)
# k=3 means EXACTLY 3 speeds > 13. Small part subset of [1,13]; 3 large in [14,30].
# Need |S|=13 => 10 small from [1,13] + 3 large from [14,30].
# Covering, primitive, S3 (Vmax>=13 Vmin, k>=2 auto).
small_universe = list(range(1, 14))
k3_min = None; k3_count = 0; k3_viol = []
checked = 0
for smalls in combinations(small_universe, 10):
    for larges in combinations(range(14, 31), 3):
        S = sorted(smalls + larges)
        if len(S) != 13: continue
        if gcd_all(S) != 1: continue
        if not is_cov(S): continue
        Vmin = min(S); Vmax = max(S)
        if Vmax < 13*Vmin: continue  # S3
        # k>=2 auto (3 larges)
        checked += 1
        Mv = Mval(S)
        k3_count += 1
        if k3_min is None or Mv < k3_min[0]:
            k3_min = (Mv, S)
        if Mv < F(1,14):
            k3_viol.append((Mv, S))
print(f"  k=3 covering primitive S3 sets on board [14,30]: {k3_count}")
if k3_min:
    print(f"  min M = {k3_min[0]} (={float(k3_min[0]):.5f})  at S={k3_min[1]}")
print(f"  claim: 11632 sets, min M = 2/23 = {F(2,23)} ({float(F(2,23)):.5f})")
print(f"  count match: {k3_count==11632};  minM match 2/23: {k3_min and k3_min[0]==F(2,23)}")
print(f"  violations (M<1/14): {len(k3_viol)}")
if k3_viol: print("   first few:", k3_viol[:3])

# Spot-check the reported worst set
worst_set = [1,2,3,4,5,7,8,11,12,13,18,20,28]
print()
print(f"  Reported worst set {worst_set}:")
print(f"    covering={is_cov(worst_set)}, primitive(gcd={gcd_all(worst_set)}), "
      f"Vmin={min(worst_set)}, Vmax={max(worst_set)}, 13*Vmin={13*min(worst_set)}")
print(f"    M = {Mval(worst_set)} (claim 2/23 = {F(2,23)})")
