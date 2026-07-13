#!/usr/bin/env python3
"""
lrc14_thin_thick_klein_S282.py
==============================
klein-2026-07-13-S282 (owner: prove the width-weighted 2nd moment bound).

Test whether the width-weighting CLOSES Q_s=O(M) or reveals a recursive residual.
S281 idea: clustered endpoints = thin arcs (width < 1/w), weight-suppressed. Split R_s arcs into
THICK (w_i >= 1/w) and THIN (< 1/w). Compute Sum_{ell!=0}|f_hat(ell w)|^2 for full, thick-only, thin-only.
 - If thin negligible AND thick gives sqrt-cancellation (~#thick/w^2 not (#thick)^2/w^2): width-weighting works.
 - If thick still ~(#thick)^2/w^2 (no cancellation) OR #thick~r (no reduction): residual persists.
Also: min separation of THICK-arc endpoints under xw (do coarse endpoints still cluster?).
"""
import math
NG=1500000
LMAX=600
def sec(e,x): return int((e*x%1.0)*7.0)%7
def occ(E,x):
    o=0
    for e in E: o|=1<<sec(e,x)
    return o
def Rs_arcs(E,s):
    arcs=[]; inR=False; start=0.0
    for k in range(1,NG):
        x=k/NG; o=occ(E,x); cur=(7-bin(o).count("1")==1) and not((o>>s)&1)
        if cur and not inR: start=(k-0.5)/NG; inR=True
        if (not cur) and inR: arcs.append((start,(k-0.5)/NG)); inR=False
    if inR: arcs.append((start,1.0))
    return arcs
def fhat2_sum(arcs,w):
    # Sum_{ell=1}^{LMAX} 2*|f_hat(ell w)|^2  (real, use both signs => factor 2 for ell>0)
    tot=0.0
    for ell in range(1,LMAX+1):
        N=ell*w; re=im=0.0
        for (a,b) in arcs:
            # integral_a^b e(-2pi N x) dx = (e(-2pi N a)-e(-2pi N b))/(2pi i N)
            ca=math.cos(-2*math.pi*N*a); sa=math.sin(-2*math.pi*N*a)
            cb=math.cos(-2*math.pi*N*b); sb=math.sin(-2*math.pi*N*b)
            # (e(-Na)-e(-Nb))/(2pi i N): real=( (sa-sb) )/(2pi N)?? do complex
            dr=ca-cb; di=sa-sb   # numerator e(-Na)-e(-Nb)
            # divide by (2 pi i N): 1/(i)= -i, so (dr+i di)/(2pi i N)= (dr+idi)(-i)/(2piN)=( di - i dr)/(2piN)
            re+= di/(2*math.pi*N); im+= -dr/(2*math.pi*N)
        tot+= 2*(re*re+im*im)
    return tot
def min_sep_xw(arcs,w):
    pts=sorted(set([a for a,b in arcs]+[b for a,b in arcs]))
    d=1.0
    for i in range(len(pts)):
        for j in range(i+1,len(pts)):
            f=(w*(pts[j]-pts[i]))%1.0; v=min(f,1-f)
            if 0<v<d: d=v
    return d

print("thin/thick decomposition of Sum|f_hat(ell w)|^2 (w=997); does removing thin arcs close it?")
print("="*78)
w=997; thr=1.0/w
print("  {:24s} {:>4} {:>7} {:>7} {:>10} {:>10} {:>10}".format("E'","diam","#arcs","#thick","full*w^2","thick*w^2","thin*w^2"))
for E in [[0,1,2,3,4,5,6],[0,3,7,15,30,55,90],[0,10,27,55,99,150,199]]:
    for s in [1]:  # one representative sector
        arcs=Rs_arcs(E,s)
        thick=[(a,b) for a,b in arcs if b-a>=thr]
        thin=[(a,b) for a,b in arcs if b-a<thr]
        full=fhat2_sum(arcs,w)*w*w
        thk=fhat2_sum(thick,w)*w*w if thick else 0.0
        thn=fhat2_sum(thin,w)*w*w if thin else 0.0
        print("  {:24s} {:>4} {:>7} {:>7} {:>10.2f} {:>10.2f} {:>10.4f}".format(str(E),max(E),len(arcs),len(thick),full,thk,thn))
        # is thick sqrt (~#thick) or square (~#thick^2)?  (full*w^2 ~ Q_s/(2pi)^2 ~ 0.04*Q_s)
print("-"*78)
print("  full*w^2 = Sum|f_hat(ell w)|^2 * w^2 = Q_s/(2pi)^2.  compare thick*w^2 to #thick (sqrt) vs #thick^2/... ")
print("  #thick vs #arcs: if #thick ~ #arcs (no reduction) => width-weighting does NOT reduce the count.")
print()
print("  min-sep of THICK endpoints under xw (do coarse endpoints still cluster?):")
for E in [[0,3,7,15,30,55,90],[0,10,27,55,99,150,199]]:
    arcs=Rs_arcs(E,1); thick=[(a,b) for a,b in arcs if b-a>=thr]
    print("    E'={:26s} #thick={:3d}  min-sep(xw)={:.5f}".format(str(E),len(thick),min_sep_xw(thick,w)))
print("\ndone.")
