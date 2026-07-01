#!/usr/bin/env python3
"""
flip_rank_bipartite_shape_klein.py  --  klein-2026-07-01-S71 (part 2)

THE SHAPE of the minimal realizing subcube = FIX A BALANCED CUT, FLIP THE TWO SIDES.
Discovery (part 1): the flip-rank achievers for n=4,5 have FIXED arcs = a complete bipartite K_{a,b} (the
cut of a balanced vertex bipartition A|B, a=ceil(n/2), b=floor(n/2)) and FREE arcs = the WITHIN-part edges
(the two sub-tournaments on A and on B). Free count = f(n) = C(a,2)+C(b,2).

STRIKING: f(n) = C(ceil(n/2),2)+C(floor(n/2),2) EQUALS ceil(log2 |G_n|) for n=3..7, then f(n) < LB for
n>=8 (so the bipartite config becomes information-theoretically IMPOSSIBLE at n>=8 -- the shape MUST change).

This script: (A) tabulate f(n) vs LB; (B) TEST whether the balanced-bipartite config realizes all iso
classes (settling rho(6)=6 vs 7): for each cross-edge orientation, do the 2^f completions cover G_n?
"""
import itertools, math

def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def perm_maps(n,E):
    idx={e:k for k,e in enumerate(E)}; maps=[]
    for pi in itertools.permutations(range(n)):
        m=[]
        for (i,j) in E:
            a,b=pi[i],pi[j]
            m.append((idx[(a,b)],0) if a<b else (idx[(b,a)],1))
        maps.append(m)
    return maps
def canon(t,maps,ne):
    best=None
    for m in maps:
        v=0
        for k in range(ne):
            bit=(t>>k)&1; tk,fl=m[k]
            if fl: bit^=1
            v|=bit<<tk
        if best is None or v<best: best=v
    return best
def classes_of(n):
    E=edges(n); ne=len(E); maps=perm_maps(n,E); cid={}; cls=[0]*(1<<ne); nc=0
    for t in range(1<<ne):
        c=canon(t,maps,ne)
        if c not in cid: cid[c]=nc; nc+=1
        cls[t]=cid[c]
    return E,ne,cls,nc

def bipartite_realizes(n, verbose=False):
    a=(n+1)//2; b=n//2; A=list(range(a)); B=list(range(a,n))
    E,ne,cls,nc=classes_of(n)
    idx={e:k for k,e in enumerate(E)}
    cross=[idx[(i,j)] for i in A for j in B]          # a*b fixed cut edges
    within=[idx[e] for e in E if e[0] in A and e[1] in A] + [idx[e] for e in E if e[0] in B and e[1] in B]
    f=len(within); lb=math.ceil(math.log2(nc))
    ok=None
    for fa in range(1<<len(cross)):                    # every cut orientation
        fb=0
        for bit,e in enumerate(cross):
            if (fa>>bit)&1: fb|=1<<e
        seen=set()
        for wa in range(1<<f):
            t=fb
            for bit,e in enumerate(within):
                if (wa>>bit)&1: t|=1<<e
            seen.add(cls[t])
            if len(seen)==nc: break
        if len(seen)==nc:
            ok=fa; break
    print(f"n={n}: |G|={nc} split {a}+{b}  free f(n)=C({a},2)+C({b},2)={f}  LB=ceil(log2)={lb}  ->  balanced-bipartite REALIZES all classes: {ok is not None}"
          + (f" (cut orientation #{ok})" if ok is not None else " (NO cut orientation works)"))
    return f, lb, ok is not None

if __name__=="__main__":
    print("(A) f(n) = C(ceil(n/2),2)+C(floor(n/2),2)  vs  LB = ceil(log2 |G_n|):")
    A=[1,1,1,2,4,12,56,456,6880,191536,9733056]
    for n in range(3,11):
        a=(n+1)//2; b=n//2; f=math.comb(a,2)+math.comb(b,2); lb=math.ceil(math.log2(A[n]))
        flag = "=" if f==lb else ("< (bipartite IMPOSSIBLE)" if f<lb else ">")
        print(f"   n={n}: f={f:>2}  {flag}  LB={lb:>2}   |G_n|={A[n]}")
    print("\n(B) does the balanced-bipartite config actually realize all iso classes?")
    for n in [4,5,6]:
        bipartite_realizes(n)
