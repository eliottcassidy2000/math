#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_routeB_tournament_layercake_opus_0620s7.py   (opus-2026-06-20-S7)

FINAL Route-B synthesis: the SURVIVAL/LAYER-CAKE certificate + the TOURNAMENT bridge,
and a search for the SINGLE coupling that unifies the three cuts.

THE REDUCTION (PROVED algebra, THM-556 + layer-cake):
    U4(E) = p0 + p5 + 5 p6
          = 1 - G_1 + G_5 + 4 G_6              (S3, exact, by Abel summation)
    where G_b(E) = P_x(N_E(x) >= b),  N_E(x) = # inner sectors {1..6} missed by {frac(e_i x)}.
  Layer-cake: any phi with phi(0)=phi[0], phi(t)=phi[0]+sum_{b>=1}(phi[b]-phi[b-1])1{t>=b}.
  For U4, phi=(1,0,0,0,0,1,5), so the survival coefficients are
    c_1=-1, c_2=c_3=c_4=0, c_5=+1, c_6=+4.   (CONVEX phi => but signs NOT monotone.)

THE CUT CERTIFICATE (VERIFIED on the binding rows, S4/S5/S6, refuted as standalone at k=12):
    (I)   G_1(E) >= G_1(consec)   [consec MAX the covered/coherent measure p0]
    (II)  G_5(E) <= G_5(consec)   [consec MAX P(N>=5)]
    (III) G_6(E) <= G_6(consec)   [consec MAX P(N>=6)]
  => U4(E) <= U4(consec) on the binding rows k=8,9,10.

THE TOURNAMENT BRIDGE (HYP-2605/R4):  at each x, T(x) is a round tournament on the
k runners; N_E(x) = 6 - #(inner sectors hit).  N large <=> orbit concentrated <=>
T(x) HIERARCHICAL (transitive-like, big score spread, low c3).  N=6 <=> orbit in one
1/7-sector <=> T(x) MAXIMALLY transitive.  So:
  G_b(E) = P_x( T(x) is at least b-deep hierarchical ).
  The survival ladder G_1>=G_2>=...>=G_6 = the FILTRATION of x by tournament hierarchy depth.
  Cut certificate = "consec's tournament-hierarchy filtration is extremal at BOTH ends:
  it spends the LEAST x-measure being merely-hierarchical (G_1 min) yet the MOST x-measure
  being DEEPLY hierarchical (G_5,G_6 max)."  That two-ended squeeze is the non-separability.

THIS SCRIPT:
  (A) re-verify the layer-cake identity U4 = 1-G1+G5+4G6 exactly on a sample.
  (B) test the UNIFYING coupling hope: is the N_E distribution of consec a
      MONOTONE REARRANGEMENT extremum?  Concretely test the two SHARP one-sided orders
      that the cut signs demand SIMULTANEOUSLY:
        - upper tail dominance:   G_b(consec) >= G_b(E) for ALL b>=2   (deep tail)
        - lower step recession:   G_1(consec) <= G_1(E)                (shallow step)
      i.e. consec is the UNIQUE shape that is tail-heavy above b=2 AND step-light at b=1.
      Count, over the bank, shapes satisfying BOTH (should be only consec).
  (C) connect to c3: compute E_x[c3(T(x))] per shape; does consec MINIMIZE mean c3
      (most hierarchical on average)?  (links to THM-555 cycle moments.)
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def dist_and_scores(E):
    """return (p-vector, mean_c3) exactly. mean_c3 = E_x[c3(T(x))], c3=C(k,3)-sum C(s_i,2).
       At cell midpoint, T(x): i->j iff frac((e_i-e_j)x) in (0,1/2). score s_i=outdeg."""
    E=sorted(set(E)); k=len(E)
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    # for c3 we also need 1/2-crossings of differences; but mean_c3 only needs the
    # tournament at generic x. Add difference half-crossings for exact c3 measure:
    diffs=set(abs(E[i]-E[j]) for i in range(k) for j in range(i+1,k))
    for d in diffs:
        if d==0: continue
        for a in range(0,2*d+1): bps.add(F(a,2*d))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[F(0)]*7; mc3=F(0); C3=comb(k,3)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; w=hi-lo
        # missed-sector count
        hit=set()
        for e in E:
            v=(e*mid)%1; hit.add((v.numerator*7)//v.denominator)
        t=sum(1 for j in range(1,7) if j not in hit); p[t]+=w
        # tournament scores
        s=[0]*k
        for a in range(k):
            for b in range(k):
                if a==b: continue
                fr=((E[a]-E[b])*mid)%1
                if F(0)<fr<F(1,2): s[a]+=1
        c3=C3-sum(comb(si,2) for si in s)
        mc3+=w*c3
    return p,mc3

def survival(p): return [sum(p[t] for t in range(a,7)) for a in range(7)]
def U4(p): return p[0]+p[5]+5*p[6]
def consec(k): return list(range(k))
def primitive(E):
    g=reduce(gcd,[e for e in E if e>0]); return tuple(e//g for e in E)

if __name__=="__main__":
    print("="*78); print("(A) layer-cake identity U4 = 1 - G1 + G5 + 4 G6 (exact spot check)"); print("="*78)
    for E in [list(range(8)), [0,2,3,5,7,8,11,13], list(range(9))]:
        p,_=dist_and_scores(E); G=survival(p)
        lhs=U4(p); rhs=1-G[1]+G[5]+4*G[6]
        print(f"  E={E}: U4={lhs}  1-G1+G5+4G6={rhs}  EQUAL={lhs==rhs}")

    print("\n"+"="*78)
    print("(B) UNIFYING two-ended squeeze + (C) mean c3 (hierarchy) extremality")
    print("="*78)
    for k in (8,9,10):
        C=consec(k); pc,mc3c=dist_and_scores(C); Gc=survival(pc); Uc=U4(pc)
        span = 13
        seen=set(); n=0
        both_squeeze=0   # shapes with G1<=consec AND G_b>=consec for all b>=2 (besides consec)
        beats_U4=0
        c3_smaller=0     # shapes with mean c3 < consec (more hierarchical than consec)
        worst_c3=None; minc3=mc3c
        EPS=F(1,10**15)
        for rest in itertools.combinations(range(1,span+1),k-1):
            E=[0]+list(rest); pr=primitive(E)
            if pr in seen: continue
            seen.add(pr); n+=1
            p,mc3=dist_and_scores(E); G=survival(p); U=U4(p)
            if U>Uc+EPS: beats_U4+=1
            # two-ended squeeze
            cond = (G[1]<=Gc[1]+EPS) and all(G[b]>=Gc[b]-EPS for b in range(2,7))
            if cond and pr!=tuple(C): both_squeeze+=1
            if mc3<mc3c-EPS:
                c3_smaller+=1
                if mc3<minc3: minc3=mc3; worst_c3=E
        print(f"  k={k}: bank={n} primitive (span<= {span})")
        print(f"     consec mean_c3 = {float(mc3c):.5f}  (c3=0 is fully transitive/hierarchical)")
        print(f"     U4 beats consec: {beats_U4}")
        print(f"     (B) non-consec shapes with the TWO-ENDED squeeze (G1<=consec & all deep G_b>=consec): {both_squeeze}")
        print(f"         => {'consec is the UNIQUE two-ended extremum (clean!)' if both_squeeze==0 else 'NOT unique; squeeze is necessary-joint'}")
        print(f"     (C) shapes with mean_c3 < consec (more hierarchical on avg): {c3_smaller}"
              + (f"  most-hier at {worst_c3} c3={float(minc3):.4f}" if worst_c3 else "  => consec MINIMIZES mean c3 (most hierarchical)!"))
