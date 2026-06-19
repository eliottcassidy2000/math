#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Part 2: stochastic orders, S_1, palindrome, the (d) counterexamples."""
import itertools
from fractions import Fraction
from functools import reduce
from math import gcd, comb

def dist_p(E):
    E=sorted(set(E))
    bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        N=sum(1 for j in range(1,7) if j not in hit)
        p[N]+=(hi-lo)
    return p

def moments(p): return [sum(p[t]*comb(t,r) for t in range(7)) for r in range(7)]
def g_poly(k):
    g=[]
    for t in range(7):
        if k==8: g.append(Fraction((t-1)*(t-2)*(t-4)*(t-5),40))
        elif k in(9,10): g.append(Fraction(-(t-2)*(t-3)*(t-6),36))
        else: g.append(Fraction((t-3)*(t-4),12))
    return g
def L_y(E,k):
    p=dist_p(E); g=g_poly(k); return sum(p[t]*g[t] for t in range(7)),p
def consec(k): return list(range(k))
def gen_sets(k,mx):
    for rest in itertools.combinations(range(1,mx+1),k-1): yield [0]+list(rest)

# per-sector J({j}) = meas{x: sector j avoided by all e}
def J_single(E,j):
    """measure of x where inner sector j is missed by orbit."""
    E=sorted(set(E))
    bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        if j not in hit: tot+=(hi-lo)
    return tot

if __name__=="__main__":
    import sys
    sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

    # claim(6): palindrome J[j]=J[6-j] for consec
    print("=== palindrome of per-sector avoid (consec) ===")
    for k in [8,9]:
        C=consec(k); Js=[J_single(C,j) for j in range(7)]
        print(f"k={k}: J[1..6]={[str(Js[j]) for j in range(1,7)]}")
        print(f"   J[1]==J[5]? {Js[1]==Js[5]}  J[2]==J[4]? {Js[2]==Js[4]}  (palindrome j<->6-j)")

    # claim(d): specific counterexamples to monotone descent
    print("\n=== (d) decrement counterexamples ===")
    for E,Edec in [([0,1,2,3,4,5,7,13],[0,1,2,3,4,5,7,12]),
                   ([0,1,2,3,4,5,6,12],[0,1,2,3,4,5,6,11])]:
        p1=dist_p(E); p2=dist_p(Edec)
        print(f"  {E}: p0={float(p1[0]):.4f}   ->  {Edec}: p0={float(p2[0]):.4f}   descent decreases p0? {p2[0]<p1[0]}")

    # claims(4)(5): stochastic dominance & S_1 over k=8 maxval<=14
    print("\n=== stochastic orders & S_1 (k=8, maxval<=14, primitive) ===")
    k=8;mx=14
    C=consec(k); pc=dist_p(C); Sc=moments(pc); S1c=Sc[1]
    # cdf for usual stoch dom: N_E <= N_C means P(N_E>=t) <= P(N_C>=t) all t
    def tail(p,t): return sum(p[s] for s in range(t,7))
    usd=0; icx=0; s1_smaller=0; nsets=0; S2_beat=0
    for E in gen_sets(k,mx):
        if reduce(gcd,E)!=1: continue
        nsets+=1
        p=dist_p(E); S=moments(p)
        if all(tail(p,t)<=tail(pc,t) for t in range(7)): usd+=1
        if S[1]<S1c: s1_smaller+=1
        if S[2]>Sc[2]: S2_beat+=1
    print(f"  {nsets} sets; usual-stoch-dom (N_E<=N_consec) holds for {usd}")
    print(f"  S_1 strictly smaller than consec for {s1_smaller} sets (consec NOT min if >0)")
    print(f"  S_2 (y_2>0 weight) beaten by {S2_beat} competitors (consec not S_2-extremal if >0)")
    print(f"  Sc=[{','.join(str(x) for x in Sc)}]")
    # sign of y_r: derive y from g via binomial inversion, check alternation
    g=g_poly(8)
    # y_r = sum_t (-1)^(t-r) C(t,r)? No: g(t)=sum_r y_r C(t,r); invert
    y=[]
    gg=list(g)
    # finite difference inversion: y_r = sum_{t} (-1)^{r-t} C(r,t) g(t)... use Newton forward diff at 0
    for r in range(7):
        yr=sum((-1)**(r-t)*comb(r,t)*g[t] for t in range(r+1))
        y.append(yr)
    print(f"  derived y (k=8) = {[str(v) for v in y]}")
    print(f"  signs: {[ (1 if v>0 else (-1 if v<0 else 0)) for v in y]}")
