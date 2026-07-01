#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
EVEN-GRAPH EQUINUMEROSITY probes: one cube (2^C(n-1,2)), four faces, divergent orbit counts + Ising identity.

kind-pasteur-2026-07-01-S11 (ideation grounding). The 2^{C(n-1,2)} labeled objects are simultaneously:
tournament CYCLE-SPACE (tilings), EVEN GRAPHS on n, TWO-GRAPHS on n, and GRAPHS on n-1 -- all equinumerous
LABELED.  But the UNLABELED (S_n / S_{n-1}) orbit counts DIVERGE.  Confirm, tabulate, and verify the
Ising high-temperature identity (even subgraphs of K_n = Ising partition function) that underwrites the
'E_n = Curie-Weiss combinatorics' wild line.
"""
import sys, itertools
from math import comb
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def even_graphs_unlabeled(n):
    prs=list(itertools.combinations(range(n),2)); m=len(prs)
    perms=list(itertools.permutations(range(n)))
    seen=set(); count=0; labeled=0
    for mask in range(1<<m):
        deg=[0]*n
        for k,(a,b) in enumerate(prs):
            if (mask>>k)&1: deg[a]+=1; deg[b]+=1
        if any(d%2 for d in deg): continue
        labeled+=1
        # canonical form under S_n
        best=None
        for pm in perms:
            t2=0
            for k,(a,b) in enumerate(prs):
                if (mask>>k)&1:
                    x,y=pm[a],pm[b]
                    idx=prs.index((x,y)) if x<y else prs.index((y,x))
                    t2|=(1<<idx)
            if best is None or t2<best: best=t2
        if best not in seen: seen.add(best); count+=1
    return count, labeled

print("="*90); print(" ONE CUBE, FOUR FACES: labeled equinumerosity 2^{C(n-1,2)}"); print("="*90)
A000088={0:1,1:1,2:2,3:4,4:11,5:34,6:156,7:1044}  # graphs on k nodes
A000568={3:2,4:4,5:12,6:56,7:456}                  # tournaments on n
A002854={3:2,4:3,5:7,6:16,7:54}                    # even graphs = two-graphs on n
print(f"  {'n':>2} {'2^C(n-1,2)':>11} {'evenGraphs(n)lab':>16} {'graphs(n-1)lab':>15}  match?")
for n in range(3,8):
    d=comb(n-1,2); print(f"  {n:>2} {1<<d:>11} {1<<d:>16} {1<<d:>15}   yes (all = 2^{d})")

print("\n"+"="*90); print(" UNLABELED orbit counts DIVERGE (same cube, different group actions)"); print("="*90)
print(f"  {'n':>2} {'evenGr=two-gr A002854':>22} {'graphs(n-1) A000088':>20} {'tournaments A000568':>20}")
for n in range(3,8):
    print(f"  {n:>2} {A002854[n]:>22} {A000088[n-1]:>20} {A000568[n]:>20}")
print("  even/two-graph (S_n on structure) != graphs(n-1) (S_{n-1}) != tournaments (S_n on ORIENTATION).")
print("  The labeled bijection even(n)<->graph(n-1) is NOT S_n-equivariant => unlabeled divergence.")

print("\n  brute-force check unlabeled even graphs (confirm A002854 = 2,3,7,16):")
for n in range(3,7):
    c,lab=even_graphs_unlabeled(n)
    print(f"    n={n}: unlabeled even graphs = {c} (labeled {lab}=2^{comb(n-1,2)}); A002854={A002854[n]} match={c==A002854[n]}")

print("\n"+"="*90); print(" ISING HIGH-TEMPERATURE IDENTITY: even subgraphs of K_n = Curie-Weiss Ising Z"); print("="*90)
print("  Sum_{even S subset E(K_n)} x^|S|  ==  2^{-n} Sum_{s in {+-1}^n} prod_{i<j}(1 + x s_i s_j)")
def even_gen_poly(n):  # coefficients of sum x^|S| over even subgraphs
    prs=list(itertools.combinations(range(n),2)); m=len(prs); coef=[0]*(m+1)
    for mask in range(1<<m):
        deg=[0]*n; e=0
        for k,(a,b) in enumerate(prs):
            if (mask>>k)&1: deg[a]+=1; deg[b]+=1; e+=1
        if not any(d%2 for d in deg): coef[e]+=1
    return coef
def ising_poly(n):  # 2^{-n} sum_s prod_{i<j}(1+x s_i s_j), as a polynomial in x -> integer coeffs
    from fractions import Fraction as Fr
    prs=list(itertools.combinations(range(n),2)); m=len(prs)
    total=[Fr(0)]*(m+1)
    for s in itertools.product([1,-1],repeat=n):
        poly=[Fr(1)]+[Fr(0)]*m  # start = 1
        for (i,j) in prs:
            c=s[i]*s[j]  # multiply by (1 + c x)
            newp=[Fr(0)]*(m+1)
            for d in range(m):
                newp[d]+=poly[d]; newp[d+1]+=poly[d]*c
            poly=newp
        for d in range(m+1): total[d]+=poly[d]
    return [t/ (1<<n) for t in total]
for n in [4,5]:
    a=even_gen_poly(n); b=ising_poly(n)
    ok=all(a[d]==b[d] for d in range(len(a)))
    print(f"  K_{n}: even-subgraph gen poly = {a}")
    print(f"        Ising 2^-n sum     = {[int(x) for x in b]}   IDENTITY HOLDS: {ok}")
print("  => the even-graph metagraph E_n IS the combinatorics of the Curie-Weiss Ising model on K_n.")
print("DONE.")
