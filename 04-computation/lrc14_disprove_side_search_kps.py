#!/usr/bin/env python3
"""
lrc14_disprove_side_search_kps  (kind-pasteur, 2026-06-16, DISPROVE side)

Adversarial stress test of the inf_S L(S) > 0 program for LRC(14).
S = 13 distinct positive integers.  L(S) = lonely measure (exact rational).

Goals:
  (A) Try to drive L -> 0 with primitive 13-sets of growing lcm (resonance/escape).
      Beat 1/1260 = 0.00079365 if possible.
  (B) Find a genuine LRC(14) counterexample: primitive 13-set with
      max-min := max_tau min_v ||v tau|| < 1/14.   (LRC holds iff max-min >= 1/14.)
  (C) Stress compactness: do large-lcm configs sneak below the bounded-lcm extremizers?

All exact (Fraction). max-min computed exactly via the arc-endpoint candidate set.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd
from functools import reduce
import sys

# ---------------- exact lonely measure (verbatim arc-sweep) ----------------
def danger_arcs(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A

def L_exact(S):
    A=[]
    for v in set(S): A.extend(danger_arcs(v))
    A=sorted((a,b) for a,b in A if b>a)
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch: ch=max(ch,b)
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

# ---------------- fast float lonely measure (screen) ----------------
def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=(14*k-1)*inv; hi=(14*k+1)*inv
            if lo<0.0: arcs.append((0.0,hi)); arcs.append((1.0+lo,1.0))
            elif hi>1.0: arcs.append((lo,1.0)); arcs.append((0.0,hi-1.0))
            else: arcs.append((lo,hi))
    arcs=[(a,b) for a,b in arcs if b>a]; arcs.sort()
    tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch:
            if hi>ch: ch=hi
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot

# ---------------- exact max-min  = max_tau min_v ||v tau|| ----------------
# min_v ||v tau|| is piecewise linear in tau; its breakpoints are where some
# ||v tau|| has a kink (tau = j/v) or where two of them cross.
# The maximum of the lower envelope of {||v tau||} occurs at a vertex:
#  either a kink point j/v, or a crossing  ||a tau||=||b tau|| with opposite slopes.
# We collect ALL candidate tau (rationals), evaluate min_v ||v tau|| exactly, take max.
def frac_norm(x):  # ||x|| exact, x a Fraction
    f = x - (x.numerator//x.denominator)   # fractional part in [0,1)
    return min(f, 1-f)

def maxmin_exact(S):
    S=sorted(set(S))
    cands=set()
    # kink points j/v
    for v in S:
        for j in range(v+1):
            cands.add(Fr(j,v))
    # crossings of ||a tau|| and ||b tau||:  a*tau - p = +-(b*tau - q)
    # case +: (a-b)tau = p-q -> tau=(p-q)/(a-b);  case -: (a+b)tau = p+q
    for i in range(len(S)):
        a=S[i]
        for jx in range(i+1,len(S)):
            b=S[jx]
            # ||a tau|| = |a tau - round|, lines a tau - p for integer p over relevant range
            for p in range(a+1):
                for q in range(b+1):
                    if a!=b:
                        t=Fr(p-q,a-b)
                        if 0<=t<=1: cands.add(t)
                    t2=Fr(p+q,a+b)
                    if 0<=t2<=1: cands.add(t2)
    best=Fr(-1); bestt=None
    for t in cands:
        if t<0 or t>1: continue
        m=min(frac_norm(v*t) for v in S)
        if m>best: best=m; bestt=t
    return best,bestt

def lcm(S):
    return reduce(lambda a,b:a*b//gcd(a,b), S, 1)
def is_primitive(S):
    return reduce(gcd,S)==1

THRESH=Fr(1,14)
TARGET=Fr(1,1260)
BASE=list(range(1,14))

def report(tag,S):
    Sf=tuple(sorted(set(S)))
    L=L_exact(S); mm,t=maxmin_exact(S)
    print(f"  [{tag}] S={Sf}")
    print(f"        lcm={lcm(Sf)} prim={is_primitive(Sf)} L={float(L):.9g}={L}  maxmin={float(mm):.9g}={mm} (1/14={float(THRESH):.6g}) at tau={t}")
    return L,mm

if __name__=="__main__":
    print("="*78)
    print("SANITY: known configs")
    print("="*78)
    report("AP 1..13", BASE)
    report("sporadic 1..11,13,24", [1,2,3,4,5,6,7,8,9,10,11,13,24])
    report("inf-cand 12->36", [1,2,3,4,5,6,7,8,9,10,11,13,36])
