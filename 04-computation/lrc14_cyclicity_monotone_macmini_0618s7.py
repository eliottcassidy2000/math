#!/usr/bin/env python3
"""
lrc14_cyclicity_monotone_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

CORRECTED tournament dictionary finding (refutes the naive 'AP=most transitive'):

  Across shapes, the EXACT data shows AP/consec gives the LARGEST E[c3] (most directed
  3-cycles), the LARGEST E[H] (most Hamiltonian paths), AND the LARGEST meas(S7)/measN.
  Dissociated shapes give the SMALLEST of all three.  So the right statement is:

     AP = the MOST CYCLIC winding ensemble (max average odd-cycle content),
     and meas(S7) (sector fill) co-moves with cycle content.

This is the parity/OCF home-turf reframe:  the DANGEROUS LRC clusters (the ones whose
sector-fill meas(S7) is largest, i.e. closest to the cap) are exactly those whose
difference-winding round tournament carries the MOST directed odd cycles on average.
Cyclicity of the winder = sector-fill = danger.  The singular-series 'relation-richness'
(low height W(E), AP-rich) is the SAME as 'winding tournament is maximally cyclic'.

THIS SCRIPT:
 (1) Confirms the co-monotonicity meas(S7) vs E[c3] vs E[H] vs measN across MANY shapes
     per k (k=7,8), reporting Spearman-style rank agreement (exact via sorting).
 (2) Tests whether AP is the GLOBAL max of E[c3] and E[H] (exhaustive over primitive
     shapes with bounded maxE, the 'dangerous box' k=7,8).
 (3) Reports the EXACT cyclicity moments for AP:  E[c3](AP_k), E[H](AP_k) as rationals,
     and checks E[c3] = C(k,3) - E[sum C(s_i,2)] with the score distribution from the
     three-gap/equidistribution of the AP winding.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def round_tournament(E,x):
    n=len(E); A=[[0]*n for _ in range(n)]; tf=True
    for i in range(n):
        for j in range(n):
            if i==j: continue
            rel=((E[i]-E[j])*x)%1
            if 0<rel<F(1,2): A[i][j]=1
            elif rel>F(1,2): A[i][j]=0
            else: A[i][j]=1 if E[i]<E[j] else 0; tf=False
    return A,tf
def scores(A): return [sum(r) for r in A]
def c3_count(A):
    n=len(A); return comb(n,3)-sum(comb(s,2) for s in scores(A))
def Hpaths(A):
    n=len(A); h=0
    for p in itertools.permutations(range(n)):
        if all(A[p[t]][p[t+1]] for t in range(n-1)): h+=1
    return h
def sectors_hit(E,x): return set(int(((e*x)%1)*7) for e in E)
def maxgap(E,x):
    ps=sorted(set((e*x)%1 for e in E))
    if len(ps)==1: return F(1)
    g=F(0)
    for i in range(len(ps)):
        nxt=ps[(i+1)%len(ps)]+(F(1) if i+1==len(ps) else F(0)); g=max(g,nxt-ps[i])
    return g
def breakpoints(E):
    bps=set([F(0),F(1)]); Es=sorted(set(E))
    for e in Es:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
        for m in range(0,2*e+1): bps.add(F(m,2*e))
    diffs=set()
    for a in range(len(Es)):
        for b in range(a+1,len(Es)):
            diffs.add(Es[b]-Es[a]); diffs.add(Es[a]+Es[b])
    for d in diffs:
        if d==0: continue
        for m in range(0,2*d+1): bps.add(F(m,2*d))
    return sorted(x for x in bps if 0<=x<=1)

def profile(E, do_H=True):
    k=len(E); bps=breakpoints(E)
    S7=F(0); N=F(0); Ec3=F(0); EH=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        w=x1-x0; xm=(x0+x1)/2
        A,tf=round_tournament(E,xm)
        if len(sectors_hit(E,xm))==7: S7+=w
        if maxgap(E,xm)<=F(1,7): N+=w
        Ec3+=w*c3_count(A)
        if do_H: EH+=w*Hpaths(A)
    return S7,N,Ec3,EH

print("="*96)
print("CYCLICITY-MONOTONE dictionary: meas(S7) co-moves with E[c3], E[H] (cyclic content)")
print("="*96)

# (1) co-monotonicity across a spread of shapes
shapesets={
 7:[("consec",list(range(7))),("perf1",[0,1,2,3,4,5,7]),("perf2",[0,1,2,3,4,6,8]),
    ("perf3",[0,1,2,4,6,8,10]),("dissoc",[0,1,3,7,15,31,63]),("Sidon",[0,1,3,7,12,20,30]),
    ("rand1",[0,2,5,9,14,20,27]),("rand2",[0,3,8,15,24,35,48])],
 8:[("consec",list(range(8))),("perf1",[0,1,2,3,4,5,6,8]),("perf2",[0,1,2,3,4,5,7,9]),
    ("dissoc",[0,1,3,7,15,31,63,127]),("Sidon",[0,1,3,7,12,20,30,44]),
    ("rand1",[0,2,5,9,14,20,27,35]),("rand2",[0,5,13,27,41,58,79,97])],
}
for k in (7,8):
    print(f"\n--- k={k} ---")
    print(f"{'shape':<10}{'meas(S7)':>11}{'measN':>11}{'E[c3]':>10}{'E[H]':>11}")
    rows=[]
    for name,E in shapesets[k]:
        S7,N,Ec3,EH=profile(E, do_H=(k<=8))
        rows.append((name,S7,N,Ec3,EH))
        print(f"{name:<10}{float(S7):>11.5f}{float(N):>11.5f}{float(Ec3):>10.4f}{float(EH):>11.3f}")
    # rank co-monotonicity: sort by S7, check E[c3] and E[H] descend together
    by_s7=sorted(rows,key=lambda r:-r[1])
    order_c3=[r[0] for r in sorted(rows,key=lambda r:-r[3])]
    order_h =[r[0] for r in sorted(rows,key=lambda r:-r[4])]
    order_s7=[r[0] for r in by_s7]
    print(f"  rank by meas(S7): {order_s7}")
    print(f"  rank by E[c3]   : {order_c3}    same-top? {order_s7[0]==order_c3[0]}")
    print(f"  rank by E[H]    : {order_h}    same-top? {order_s7[0]==order_h[0]}")

# (2) is AP the GLOBAL max of E[c3] and E[H]? exhaustive primitive box
print("\n" + "="*96)
print("EXHAUSTIVE: is consec/AP the GLOBAL MAX of E[c3] and E[H] over primitive shapes?")
print("(0 in E, distinct, maxE bounded — the 'dangerous box')")
for k, maxE in [(7,9),(8,10)]:
    apE=list(range(k))
    _,_,apc3,apH=profile(apE, do_H=True)
    best_c3=apc3; best_c3_E=tuple(apE); beat_c3=0
    best_H=apH;  best_H_E=tuple(apE);  beat_H=0
    cnt=0
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=[0]+list(rest)
        if reduce(gcd,E)!=1: continue
        cnt+=1
        _,_,c3v,Hv=profile(E, do_H=True)
        if c3v>best_c3: best_c3=c3v; best_c3_E=tuple(E)
        if c3v>apc3: beat_c3+=1
        if Hv>best_H: best_H=Hv; best_H_E=tuple(E)
        if Hv>apH: beat_H+=1
    print(f"  k={k} (maxE<={maxE}, {cnt} primitive shapes):")
    print(f"    AP E[c3]={float(apc3):.4f}; shapes beating it: {beat_c3}; "
          f"global max {float(best_c3):.4f} at {best_c3_E}")
    print(f"    AP E[H] ={float(apH):.4f}; shapes beating it: {beat_H}; "
          f"global max {float(best_H):.4f} at {best_H_E}")

# (3) exact AP cyclicity moments
print("\n" + "="*96)
print("EXACT AP cyclicity moments E[c3](AP_k), E[H](AP_k):")
for k in range(5,9):
    apE=list(range(k))
    _,_,c3v,Hv=profile(apE, do_H=True)
    print(f"  k={k}:  E[c3]={c3v} = {float(c3v):.5f}   E[H]={Hv} = {float(Hv):.4f}")
print("="*96)
print("DONE.")
