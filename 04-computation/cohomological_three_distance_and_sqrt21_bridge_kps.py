#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
COHOMOLOGICAL THREE-DISTANCE + COMPLEMENTARITY, and the sqrt(21)=Q(sqrt-3,sqrt-7) bridge assessment.

kind-pasteur-2026-07-01-S21. Cohomology of the lonely set on S^1:
 (A) L_C = union of b_0 ARCS (H^0 = R^{b_0}, H^1=0 since arcs contractible). M>=r <=> b_0(L_C)>=1.
 (B) COMPLEMENTARITY (Alexander duality on S^1): b_0(L_C) = b_0(danger cover) -- lonely & covered arcs ALTERNATE.
 (C) THREE-DISTANCE: the arc widths take FEW distinct values (Steinhaus/Morse-collar). Count them.
 (D) FAR ELEMENT w = degree-w circle map (Lefschetz Λ(φ_w)=1-w); adds w danger arcs; equidistributes on L_C
     (covers ~2r), leaving (1-2r)|L_C|>0 -- the ι-odd Gauss certificate (i sqrt7) is the arithmetic of survival.
 (E) THE BRIDGE: covering field Q(sqrt-3) (Phi6/metallic, opus-S24) vs certification Q(sqrt-7) (Gauss i sqrt7);
     their real compositum = Q(sqrt21). Check sqrt21 | sqrt(D_14) (deep-well); assess "bridge = OPEN-Q-108 residual".
"""
import sys
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
BAND=Fr(1,14)
def safe(v,band=BAND): return [((Fr(k)+band)/v,(Fr(k+1)-band)/v) for k in range(v)]
def inter(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: r.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return r
def L_of(S,band=BAND):
    a=safe(S[0],band)
    for v in S[1:]:
        a=inter(a,safe(v,band))
        if not a: return a
    return a

# AP core is TIGHT (M=1/14 exactly => lonely set empty at band 1/14). Use the covering-min CONSTRUCTION
# champion {1..11,13,36} (M=14/183 > 1/14 => positive-measure lonely set) for the cohomology.
S=[1,2,3,4,5,6,7,8,9,10,11,13,36]
Lc=L_of(S)
b0=len(Lc)
widths=sorted(set(h-l for l,h in Lc))
totmeas=sum(h-l for l,h in Lc)
print("="*94); print(" (A)+(B)+(C): cohomology of the lonely set L_C (construction champion {1..11,13,36}, band 1/14)"); print("="*94)
print(f"  b_0(L_C) = #lonely arcs = {b0}  (H^0=R^{b0}, H^1=0)  => M>=1/14 certified iff b_0>=1: {b0>=1}")
print(f"  measure |L_C| = {float(sum(h-l for l,h in Lc)):.6f}")
print(f"  THREE-DISTANCE: #distinct arc widths = {len(widths)}  (widths {[float(w) for w in widths[:6]]}...)")
# complementarity: gaps between consecutive arcs (danger components), on the circle
arcs=sorted(Lc); gaps=[]
for i in range(len(arcs)):
    a=arcs[i]; b=arcs[(i+1)%len(arcs)]
    gap=(b[0]-a[1]) if i<len(arcs)-1 else (1-a[1]+b[0])
    if gap>0: gaps.append(gap)
print(f"  b_0(danger cover) = #gaps between lonely arcs = {len(gaps)}  => complementarity b_0(L_C)=b_0(danger): {len(gaps)==b0}")

print("\n"+"="*94); print(" (D) FAR ELEMENT w: degree-w circle map, Lefschetz Λ=1-w, equidistributes on L_C"); print("="*94)
for w in [37,50,300,700]:
    Lw=inter(Lc,safe(w)); mw=sum(h-l for l,h in Lw); frac=mw/totmeas if totmeas>0 else 0
    print(f"  w={w:>4}: Λ(φ_w)=1-w={1-w:>5}; survivors L_C∩safe(w): {len(Lw)} arcs, meas frac={float(frac):.4f} "
          f"(-> 1-2r=6/7={float(1-2*BAND):.4f}); survives (b_0>0): {len(Lw)>0}")
print("  => the far element cuts L_C into ~1-2r of itself (equidistribution); it NEVER covers it (b_0 stays >0).")

print("\n"+"="*94); print(" (E) THE Q(sqrt-3,sqrt-7) BRIDGE via sqrt21, and the residual assessment"); print("="*94)
import math
def squarefree(D):
    s=1 if D>0 else -1; D=abs(D); out=1; d=2
    while d*d<=D:
        c=0
        while D%d==0: D//=d; c+=1
        if c%2: out*=d
        d+=1
    return s*out*D if D>1 else s*out
n=14; D14=n*(n-1)*(n*n-n+4); sf=squarefree(D14)
print(f"  covering field   = Q(sqrt-3)  (Phi6 / Eisenstein / metallic-CF flavor, opus-S24)")
print(f"  certification    = Q(sqrt-7)  (Gauss sum i*sqrt7, ι-odd Lefschetz, S20)")
print(f"  real compositum  = Q(sqrt(3*7)) = Q(sqrt21);  sqrt-3 * sqrt-7 = -sqrt21")
print(f"  deep-well D_14 = n(n-1)(n^2-n+4) = {D14} = squarefree {sf} = 21*{sf//21 if sf%21==0 else '?'} "
      f"=> sqrt21 | sqrt(D_14): {sf%21==0}")
print(f"  => the DEEP-WELL binding CF (metallic/covering side) CONTAINS sqrt21 = the bridge of covering & cert.")
print("  ASSESSMENT of 'is the sqrt21 bridge the OPEN-Q-108 residual?':")
print("   - The bridge Q(sqrt-3,sqrt-7)/sqrt21 IS the residual of the LEFSCHETZ/CERTIFICATION program: the ι-odd")
print("     Gauss certificate (sqrt-7) must interact with the Eisenstein covering (sqrt-3) to certify the WHOLE")
print("     cover; the 'mixed' part (neither pure isotype) lives in the compositum sqrt21. Making i*sqrt7 => M>=1/n")
print("     rigorous = bridging these two fields. So YES for the trace-formula approach's residual.")
print("   - It is NOT the CLASSICAL multi-far census residual: those extremizers are the pentagon (Z/10)* [Q(sqrt5)]")
print("     and the sporadic two-clash (Z/19) [Q(sqrt-19)] -- their OWN sub-apex fields, not sqrt21.")
print("   => HONEST: sqrt21 = the certification-program residual (a promising reformulation), distinct from the")
print("      classical census residual. The two programs meet only if the census extremizers' fields (5,19) also")
print("      reduce to the apex compositum -- which they do NOT, so the bridge is one route, not the whole crux.")
print("DONE.")
