#!/usr/bin/env python3
"""
covering_min_open_edge_klein.py  --  klein-2026-06-30-S53

The OPEN EDGE of the covering-min conjecture: does the construction {1..n-2, n(n-1)} = n/Phi6(n)
stay the covering-min for all large n?  Past thin-edge patterns:
  - HYP-3551: tightest covering = densest coverable core {1..n-2} (M=1/(n-1)) + minimal killer n(n-1).
  - HYP-3726: floor margin 1/(n(2n-1)) = 1/hexagonal; HYP-3701: construction beaten n=7,8,9, margin ~1/n^2.
  - HYP-3763 (mine): single-drop covering escapes exceed n/Phi6 for n=10,12,14 with margin shrinking ~1/n^2.

This maps, exactly, the competition across n between:
  (a) CONSTRUCTION  {1..n-2, n(n-1)}                       M = n/Phi6(n)
  (b) DROP-2 mediant {1..n-1}\{2} + killer (mult of n; +even)   the small-n beater
  (c) DROP-(n-2) double-killer  {1..n-3,n-1} + lcm(n-2,n)       my 7/89 family
  (d) densest-core lower reference 1/(n-1)
and reports, per n, the winner, the runner-up, and the runner-up margin over the construction.
"""
from fractions import Fraction as F
from math import gcd

def Phi6(n): return n*n - n + 1
def lcm(a,b): return a*b//gcd(a,b)
def dist0(x,D):
    x%=D; return min(x,D-x)

def M_and_arg(S):
    Dmax = 2*max(S)+2; best=F(0); arg=None
    for D in range(2,Dmax+1):
        for a in range(1,D):
            m=D
            for s in S:
                d=dist0(s*a,D)
                if d<m: m=d
                if m==0: break
            if F(m,D)>best: best=F(m,D); arg=(D,a)
    return best,arg

def is_covering(S,n): return all(any(s%q==0 for s in S) for q in range(2,n+1))
def prim(S):
    g=0
    for s in S: g=gcd(g,s)
    return g==1

def best_drop2(n, mmax=8):
    """{1..n-1}\{2} + one killer (mult of n, made even so it covers 2). min M over the killer."""
    base=[j for j in range(1,n) if j!=2]
    best=None
    for m in range(1,mmax+1):
        for kill in ({n*m} | {2*n*m}):
            S=base+[kill]
            if len(set(S))!=n-1 or not prim(S) or not is_covering(S,n): continue
            Mv,arg=M_and_arg(S)
            if best is None or Mv<best[0]: best=(Mv,kill,arg)
    return best

def best_double_killer(n, mmax=6):
    """{1..n-3, n-1} + killer covering (n-2) and n, i.e. mult of lcm(n-2,n). min M."""
    base=[j for j in range(1,n-2)]+[n-1]
    L=lcm(n-2,n); best=None
    for m in range(1,mmax+1):
        S=base+[L*m]
        if len(set(S))!=n-1 or not prim(S) or not is_covering(S,n): continue
        Mv,arg=M_and_arg(S)
        if best is None or Mv<best[0]: best=(Mv,L*m,arg)
    return best

if __name__=="__main__":
    print(f"{'n':>3} {'construction':>13} {'drop2-mediant':>16} {'double-killer':>16} {'1/(n-1)':>9}  winner / runner-up margin over constr")
    print("-"*112)
    for n in range(7,19):
        thr=F(n,Phi6(n)); constr=list(range(1,n-1))+[n*(n-1)]
        Mc,_=M_and_arg(constr)
        d2=best_drop2(n); dk=best_double_killer(n)
        cands={'constr':Mc}
        if d2: cands['drop2']=d2[0]
        if dk: cands['dblkill']=dk[0]
        win=min(cands,key=lambda k:cands[k]); winM=cands[win]
        # runner-up among non-winners
        others={k:v for k,v in cands.items() if k!=win}
        ru=min(others,key=lambda k:others[k]) if others else None
        rum = (cands[ru]-winM) if ru else None
        d2s=f"{str(d2[0]):>10}={float(d2[0]):.4f}" if d2 else "  --"
        dks=f"{str(dk[0]):>10}={float(dk[0]):.4f}" if dk else "  --"
        flag = "  <== CONSTR NOT min!" if win!='constr' else ""
        print(f"{n:>3} {str(Mc):>7}={float(Mc):.4f} {d2s:>16} {dks:>16} {float(F(1,n-1)):>9.4f}  "
              f"win={win}{flag}; runner-up={ru} margin={float(rum):+.5f}" )
    print()
    print("constr = n/Phi6 (densest core {1..n-2}+killer n(n-1)); drop2 = mediant beater; dblkill = my 7/89 family")
