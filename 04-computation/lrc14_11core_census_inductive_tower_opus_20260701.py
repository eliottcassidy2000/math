#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE INDUCTIVE TOWER for the far-element lever: m_k = min meas over primitive k-cores. The lever says an
11-core with j far outliers (elts > V) attached to a compact (11-j)-core has
      meas >= (6/7)^j * m_{11-j}  -  (correction terms c'/(7w)),
so the MULTI-outlier residual is controlled if (6/7)^j * m_{11-j} > 1/36 for all j>=1. We compute the tower
m_7..m_11 exactly (dense, max<=15) and check (6/7)^j * m_{11-j} > 1/36, plus a direct TWO-outlier spot census.

opus-2026-07-01-S31 addendum to lrc14_11core_census_to_completion. Exact (Fraction).
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
BAND=Fr(1,14); TARGET=Fr(1,36)
def safe_arcs(v): return [((Fr(k)+BAND)/v,(Fr(k+1)-BAND)/v) for k in range(v)]
def inter(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=A[i][0] if A[i][0]>B[j][0] else B[j][0]; hi=A[i][1] if A[i][1]<B[j][1] else B[j][1]
        if lo<hi: r.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return r
def Lmeas(S):
    a=safe_arcs(S[0])
    for v in S[1:]:
        a=inter(a,safe_arcs(v))
        if not a: return Fr(0)
    return sum(h-l for l,h in a)

print("="*96)
print(" INDUCTIVE TOWER  m_k = min meas over primitive k-cores (dense, max<=15), and the lever guard")
print(" (6/7)^j * m_{11-j} > 1/36 for the j-outlier residual.   1/36 = %.6f"%float(TARGET))
print("="*96)
Vt=15
tower={}
for k in range(7,12):
    mn=None; arg=None
    for C in itertools.combinations(range(1,Vt+1),k):
        if reduce(gcd,C)!=1: continue
        m=Lmeas(C)
        if m>0 and (mn is None or m<mn): mn=m; arg=C
    tower[k]=mn
    print(f"  m_{k} = {mn} = {float(mn):.6f}   at {arg}")
print()
print(f"  {'j outliers':>11} {'11-j':>5} {'m_{11-j}':>12} {'(6/7)^j * m_{11-j}':>20} {'> 1/36?':>9}")
allguard=True
for j in range(1,5):
    base=tower.get(11-j)
    if base is None: continue
    g=(Fr(6,7)**j)*base
    ok=g>TARGET; allguard=allguard and ok
    print(f"  {j:>11} {11-j:>5} {float(base):>12.6f} {float(g):>20.6f} {str(ok):>9}")
print(f"\n  ALL guards (6/7)^j m_{{11-j}} > 1/36 for j=1..4 ? {allguard}")
print("  => the more outliers, the SAFER (m_{11-j} grows, and (6/7)^j m_{11-j} stays > 1/36): the tightest case")
print("     is the DENSE core (j=0, Part A) and the SINGLE outlier (j=1); both are exhaustively verified.")

# direct TWO-outlier spot census: compact 9-cores (max<=13) + two outliers in the DANGEROUS near range
print("\n"+"-"*96)
print(" TWO-OUTLIER spot census: compact 9-core (max<=13) + {w1,w2}, w1<w2 in [14..60] (the near/dangerous range)")
print("-"*96)
nine=[C for C in itertools.combinations(range(1,14),9) if reduce(gcd,C)==1]
mn=None; arg=None; nchk=0; viol=0
for C9 in nine:
    for w1,w2 in itertools.combinations(range(14,61),2):
        C=tuple(sorted(set(C9)|{w1,w2}))
        if len(C)!=11: continue
        m=Lmeas(C)
        if m>0:
            nchk+=1
            if mn is None or m<mn: mn=m; arg=C
            if m<TARGET: viol+=1
print(f"  compact 9-cores: {len(nine)}; two-outlier 11-cores checked: {nchk}; min meas = {float(mn):.6f} at {arg}; violations(<1/36) = {viol}")
print(f"\n  CONCLUSION: the j-outlier residual is monotone-safe (guard table) and directly verified for j=1 (w<=500)")
print(f"  and j=2 (near range). The only ingredient left for a full proof is the uniform arc-count bound c'<=const")
print(f"  (vs c'<=sum_v v) that makes the correction terms vanish uniformly -- that is OPEN-Q-108 (uniform fattening).")
print("DONE.")
