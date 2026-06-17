#!/usr/bin/env python3
"""
lrc14_perturbation_gap_recover_kps  (kind-pasteur, PROVE side)

EXACT decomposition of the opened lonely measure under single perturbation e->w:

    L(C\{e} u {w}) = meas( G_e \ D_w ) = meas(G_e) - meas(G_e ∩ D_w)

where G_e = e-exclusive gap of tight config C, D_w = danger set of new speed w.

For each tight config and each dropped e, we sweep w and report:
  - best (largest) recovery meas(G_e ∩ D_w)  -> minimal opened L for that e
  - the global minimum over (e,w)
This pins down WHICH drop and WHICH w realize inf, with EXACT fractions, and
shows the gap-vs-recovery tradeoff explicitly.

Also: for the champion (e=12,w=36 in AP) print the residual lonely arcs explicitly
to see geometric structure (why exactly 1/1260).
"""
from fractions import Fraction as F

def danger(v):
    out=[]; w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(a,b) for a,b in out if b>a]
def union(arcs):
    arcs=sorted((a,b) for a,b in arcs if b>a)
    if not arcs: return []
    res=[]; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=max(ch,hi)
        else: res.append((cl,ch)); cl,ch=lo,hi
    res.append((cl,ch)); return res
def total(arcs): return sum(b-a for a,b in union(arcs))
def complement(arcs):
    u=union(arcs); res=[]; prev=F(0)
    for lo,hi in u:
        if lo>prev: res.append((prev,lo))
        prev=max(prev,hi)
    if prev<1: res.append((prev,F(1)))
    return res
def intersect(A,B):
    A=union(A); B=union(B); res=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: res.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return res
def L_exact(S):
    arcs=[]
    for v in set(S): arcs+=danger(v)
    return F(1)-total(arcs)
def gap_of_drop(C,e):
    others=[v for v in C if v!=e]; cov=[]
    for v in others: cov+=danger(v)
    return complement(cov)

TIGHT = {
    "AP": list(range(1,14)),
    "SPORADIC": list(range(1,12))+[13,24],
}
WMAX=400

for name,C in TIGHT.items():
    print("="*78)
    print(f"TIGHT {name} = {C}   (L=0)")
    print("="*78)
    per_e_best=[]
    for e in C:
        Ge=gap_of_drop(C,e); mg=total(Ge)
        # sweep w: minimal STRICTLY-POSITIVE opened L (exclude tight-preserving w)
        best=(F(2),None)   # (opened, w)
        for w in range(2,WMAX+1):
            if w==e: continue
            if w in C and w!=e: continue   # must keep 13 distinct
            opened = mg-total(intersect(Ge,danger(w)))
            if opened>0 and opened<best[0]: best=(opened,w)
        opened,w=best
        per_e_best.append((opened,e,w,mg))
    per_e_best.sort()
    print(f"  {'e':>3s} {'best_w':>6s} {'meas(G_e)':>12s} {'min opened L>0':>16s}  float")
    for opened,e,w,mg in per_e_best:
        print(f"  {e:>3d} {w:>6d} {str(mg):>12s} {str(opened):>16s}  {float(opened):.3e}")
    glob=per_e_best[0]
    print(f"\n  >>> minimal POSITIVE single-perturbation L for {name}: {glob[0]} = {float(glob[0]):.6e}")
    print(f"      via drop e={glob[1]} -> add w={glob[2]}")
    S=[x for x in C if x!=glob[1]]+[glob[2]]
    print(f"      direct L_exact check = {L_exact(S)}")

# ---- residual lonely arcs for the AP champion 12->36 ----
print("\n"+"="*78)
print("CHAMPION residual lonely set: AP, drop 12, add 36  (expect L=1/1260)")
print("="*78)
C=list(range(1,14)); S=[x for x in C if x!=12]+[36]
arcs=[]
for v in S: arcs+=danger(v)
lonely=complement(arcs)
print(f"  L = {L_exact(S)}")
print(f"  residual lonely arcs ({len(lonely)}):")
for lo,hi in lonely:
    print(f"     [{lo}, {hi}]  len={hi-lo}={float(hi-lo):.6e}  center~{float((lo+hi)/2):.6f}")
# what speed-12 danger arcs were NOT recovered by 36?
G12=gap_of_drop(C,12)
print(f"\n  G_12 (12-exclusive gap) meas={total(G12)} = {float(total(G12)):.6e}, arcs:")
for lo,hi in G12: print(f"     [{lo},{hi}] center~{float((lo+hi)/2):.6f}")
print(f"  meas(G_12 cap D_36) recovered = {total(intersect(G12,danger(36)))}")
print(f"  meas(G_12 cap D_24) recovered = {total(intersect(G12,danger(24)))} (the tight w)")
