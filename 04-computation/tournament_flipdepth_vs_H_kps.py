#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
BRIDGE the minimal-flip DEPTH (min feedback arc set) to the project's core invariant H (Ham paths / OCF).

kind-pasteur-2026-07-01-S10. The metagraph geodesic from a class to TRANSITIVE = its MIN FEEDBACK ARC SET
(MFAS) = the minimal-flip depth (HYP-3803 frame).  Transitive: MFAS=0, H=1.  R_7: MFAS=6=phi(14), H=175.
Paley: MFAS=7, H=189 (max).  Does MFAS track H (both measure 'distance from transitive')?  Map n=7 landmarks
+ a sample; find the MAX MFAS; and honestly bound the 6=phi(14) alignment (n=7-specific? check R_5 vs phi(10)).
"""
import sys, itertools, random
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def MFAS(A,n):
    best=10**9
    for order in itertools.permutations(range(n)):
        pos={v:k for k,v in enumerate(order)}
        b=sum(1 for i in range(n) for j in range(n) if A[i][j] and pos[i]>pos[j])
        if b<best: best=b
        if best==0: break
    return best
def Hpaths(A,n):
    c=0
    for p in itertools.permutations(range(n)):
        if all(A[p[t]][p[t+1]] for t in range(n-1)): c+=1
    return c
def c3(A,n):
    c=0
    for i,j,k in itertools.combinations(range(n),3):
        if (A[i][j]and A[j][k]and A[k][i]) or (A[i][k]and A[k][j]and A[j][i]): c+=1
    return c
def circ(conn,n):
    return [[1 if (i!=j and (j-i)%n in conn) else 0 for j in range(n)] for i in range(n)]
def rand_tour(n,rng):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<0.5: A[i][j]=1
            else: A[j][i]=1
    return A

n=7
transitive=[[1 if i<j else 0 for j in range(n)] for i in range(n)]
R7=circ({1,2,3},n); Paley=circ({1,2,4},n)
print("="*88); print(" n=7 LANDMARKS: minimal-flip depth (MFAS) vs H (Ham paths) vs #3-cycles"); print("="*88)
print(f"  {'tournament':>22} {'MFAS(depth)':>11} {'H':>6} {'c3':>5}")
for name,A in [("transitive",transitive),("R_7 rot {1,2,3}",R7),("Paley QR {1,2,4}",Paley)]:
    print(f"  {name:>22} {MFAS(A,n):>11} {Hpaths(A,n):>6} {c3(A,n):>5}")

print("\n=== sample 500 random 7-tournaments: MFAS vs H correlation, and MAX MFAS ===")
rng=random.Random(11); data=[]; maxmf=(0,None)
for _ in range(500):
    A=rand_tour(n,rng); mf=MFAS(A,n); H=Hpaths(A,n); data.append((mf,H))
    if mf>maxmf[0]: maxmf=(mf,A)
# Spearman-ish: sort by MFAS, check H monotonic trend via bucket means
from collections import defaultdict
buck=defaultdict(list)
for mf,H in data: buck[mf].append(H)
print(f"  {'MFAS':>5} {'#tours':>7} {'meanH':>9} {'minH':>6} {'maxH':>6}")
for mf in sorted(buck):
    hs=buck[mf]; print(f"  {mf:>5} {len(hs):>7} {sum(hs)/len(hs):>9.1f} {min(hs):>6} {max(hs):>6}")
# Pearson correlation MFAS vs H
import statistics as st
mfs=[d[0] for d in data]; Hs=[d[1] for d in data]
mmf=st.mean(mfs); mH=st.mean(Hs)
cov=sum((a-mmf)*(b-mH) for a,b in data)/len(data)
r=cov/(st.pstdev(mfs)*st.pstdev(Hs))
print(f"  Pearson r(MFAS, H) = {r:.3f}  (both measure distance-from-transitive => expect positive)")
print(f"  MAX MFAS in sample = {maxmf[0]}  (Paley MFAS=7, transitive=0)")

print("\n=== honesty: is MFAS=phi(2n) a general law, or n=7-specific? R_n rotational, n=5,7 ===")
for nn in (5,7):
    Rn=circ(set(range(1,(nn-1)//2+1)),nn)  # 'beat next (n-1)/2'
    import math
    phi=sum(1 for a in range(1,2*nn) if math.gcd(a,2*nn)==1)
    print(f"  n={nn}: MFAS(R_{nn})={MFAS(Rn,nn)}, phi({2*nn})={phi}  -> {'EQUAL' if MFAS(Rn,nn)==phi else 'NOT equal'}")
print("  => 6=phi(14) at n=7 is a SPECIFIC alignment (the LRC-relevant case), not a general identity.")
print("DONE.")
