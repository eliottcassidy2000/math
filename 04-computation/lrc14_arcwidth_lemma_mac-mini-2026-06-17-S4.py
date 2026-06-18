#!/usr/bin/env python3
"""
lrc14_arcwidth_lemma — mac-mini-2026-06-17-S4

A NON-TAUTOLOGICAL sufficient condition (the genuine missing estimate, large-m side).

For a covering core A with M(A)>1/14, the level-1/14 SAFE SET
    G_A = {tau : ||v tau|| >= 1/14 for all v in A}
is a finite union of closed arcs (positive measure since M(A)>1/14). The parked
runner w=14m has danger TEETH {||14m tau||<1/14}: 14m intervals, each of full width
t(m)=1/(98m), centered at k/(14m). A tooth has width t(m); the safe gaps between
teeth have width 6 t(m).

LEMMA. If the widest arc of G_A has width W(A) > t(m) = 1/(98m), then that arc cannot
lie inside a single tooth, so it contains a point with ||14m tau||>=1/14; at that
point min over ALL of S=A u {14m} is >= 1/14. Hence  M(S) >= 1/14.
=> For all m >= M0(A) := floor(1/(98 W(A))) + 1, looseness is AUTOMATIC; only m<M0(A)
   needs a finite check. This rigorously discharges the infinite m-direction.

This is NOT tautological: W(A) depends only on the CORE A (computed without reference
to S), and the conclusion M(S)>=1/14 is derived, not assumed. Verifies the lemma and
computes M0(A) for the covering cores.
"""
from fractions import Fraction as F
from math import gcd, floor

C=F(1,14)  # the loneliness level

def danger_arcs_level(v, c=C):
    """open arcs where ||v tau|| < c : (k/v - c/v, k/v + c/v), k=0..v-1, on [0,1)."""
    hw=F(c,v); out=[]
    for k in range(v):
        lo=F(k,v)-hw; hi=F(k,v)+hw
        out.append((lo,hi))
    return out

def wrap(intervals):
    out=[]
    for lo,hi in intervals:
        a=lo-(lo - (lo % 1)); # normalize lo into [0,1)
        shift=lo-(lo%1); a=lo-shift; b=hi-shift
        if b<=1: out.append((a,b))
        else: out.append((a,F(1))); out.append((F(0),b-1))
    return out

def union(intervals):
    iv=sorted(intervals); res=[]
    cl,ch=iv[0]
    for lo,hi in iv[1:]:
        if lo<=ch: ch=ch if ch>hi else hi
        else: res.append((cl,ch)); cl,ch=lo,hi
    res.append((cl,ch)); return res

def safe_arcs(A, c=C):
    """closed safe arcs G_A = [0,1) minus union of open danger arcs."""
    dz=[]
    for v in set(A): dz.extend(danger_arcs_level(v,c))
    dz=union(wrap(dz))
    # complement on [0,1) (circle): gaps between consecutive danger arcs
    arcs=[]
    for i in range(len(dz)):
        hi=dz[i][1]; lo=dz[(i+1)%len(dz)][0]
        if i==len(dz)-1: lo=dz[0][0]+1
        w=lo-hi
        if w>0: arcs.append((hi%1, w))
    return arcs

def widest_arc(A,c=C):
    arcs=safe_arcs(A,c)
    if not arcs: return F(0)
    return max(w for _,w in arcs)

# exact M tool
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); Cc=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): Cc.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): Cc.add(F(k,d)); k+=1
    Cc.add(F(1,2)); return Cc
def M(S): return max(g(S,t) for t in cand(S))

print("="*78)
print("ARC-WIDTH LEMMA: W(A) > 1/(98m) => M(A u {14m}) >= 1/14   (large-m, rigorous)")
print("="*78)
cores={
 'drop1 {2..13}'   : [v for v in range(1,14) if v!=1],
 'drop6 {1..5,7..13}': [v for v in range(1,14) if v!=6],
 'drop12 {1..11,13}' : [v for v in range(1,14) if v!=12],
}
for name,A in cores.items():
    MA=M(A); W=widest_arc(A)
    M0=floor(F(1,98)/W)+1 if W>0 else None
    print(f"\n  {name}: M(A)={MA}={float(MA):.5f} (>1/14? {MA>C})")
    print(f"    widest safe arc W(A)={W}={float(W):.6f};  tooth t(m)=1/(98m)")
    print(f"    => M0(A)=floor(1/(98 W))+1 = {M0}  (m>=M0 : W>t(m) => M(S)>=1/14 automatic)")
    # verify lemma: for m=1..M0+5 check the predicate vs actual
    print(f"    verify (m : W>t(m)? predicted-safe | actual M(A u 14m)>=1/14?):")
    allok=True
    for m in range(1,(M0 or 1)+6):
        S=A+[14*m]
        if len(set(S))!=len(A)+1: continue
        pred = W>F(1,98*m)
        act = M(S)>=C
        if pred and not act: allok=False; print("       LEMMA VIOLATED at m=",m)
        if m<=3 or m in (M0-1,M0,M0+1):
            print(f"       m={m:2d}: W>t? {str(pred):5s} (t={float(F(1,98*m)):.5f}) | M(S)={float(M(S)):.5f}>=1/14? {act}")
    print(f"    lemma holds (predicted-safe => actually safe) for all tested m: {allok}")

print("\n"+"="*78)
print("CONSEQUENCE: the infinite m-direction is discharged RIGOROUSLY for each core:")
print("  m>=M0(A): automatic by the arc-width lemma (W(A) computed from A alone).")
print("  m<M0(A):  a FINITE check. Residual LRC(14) = (finitely many cores)x(m<M0).")
print("  This is NON-tautological: W(A) is a property of the CORE, M(S)>=1/14 is derived.")
print("  The remaining gap is only the bounded small-m window + bounding the cores.")
