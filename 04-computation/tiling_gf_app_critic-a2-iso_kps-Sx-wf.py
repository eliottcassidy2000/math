#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Critic part 2: (i) E[alpha2] full at n=8 via per-subset decomposition;
(ii) iso-class counts A000568, self-converse SC, V_merged=(A000568+SC)/2 n=3..6."""
import sys
from itertools import combinations, permutations, product
from math import comb
from fractions import Fraction as Fr
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# ============ PART (i): E[alpha2] via per-subset linearity ============
# alpha2 = # unordered pairs of vertex-disjoint odd directed cycles.
# At n<=8 the relevant odd-cycle pairs are 3-3 and 3-5 (and 5-5 needs >=10 vtx).
# By linearity of expectation over disjoint vertex subsets:
#   E[alpha2] = sum over {disjoint subset pairs (A,B)} Pr[A is dir odd cycle]*Pr[B is dir odd cycle]
# where "A is a dir cycle" sums over the (|A|-1)!/2 cyclic orders of A.
# Probability a specified directed cycle exists on integer-labeled vertices:
#  each arc (u,v): if |u-v|==1 it's a forced base-path arc (prob 1 if it matches base dir, else 0);
#  else free fair coin prob 1/2.
# Base path: edge between consecutive integers k,k-1 is oriented k->k-1.

def cycle_prob(cyc):
    """Pr that the directed closed walk cyc (tuple of vertices, last->first) is all present."""
    p = Fr(1)
    k = len(cyc)
    for i in range(k):
        u = cyc[i]; v = cyc[(i+1)%k]
        if abs(u-v) == 1:
            # base-path arc oriented high->low
            need_high_to_low = (u > v)  # arc u->v exists in base iff u=v+1
            present = (u == v+1)
            p *= Fr(1) if present else Fr(0)
            if not present:
                return Fr(0)
        else:
            p *= Fr(1,2)
    return p

def prob_subset_is_dircycle(subset):
    """sum over directed Hamiltonian cycles of the subset of the existence probability."""
    s = list(subset)
    first = s[0]; rest = s[1:]
    total = Fr(0)
    seen = set()
    for perm in permutations(rest):
        cyc = (first,)+perm
        # avoid double counting reversal: canonicalize by smaller second element vs last
        rev = (first,)+tuple(reversed(perm))
        key = min(cyc, rev)
        if key in seen: continue
        seen.add(key)
        total += cycle_prob(cyc)
    return total

def E_alpha2(n, sizes_pairs):
    """sizes_pairs: list of (k1,k2) odd cycle-length pairs to include."""
    total = Fr(0)
    verts = list(range(1, n+1))
    contributions = {}
    for (k1,k2) in sizes_pairs:
        sub = Fr(0)
        if k1 == k2:
            for A in combinations(verts, k1):
                Aset=set(A)
                restv=[v for v in verts if v not in Aset]
                for B in combinations(restv, k2):
                    if min(B) < min(A):  # unordered, avoid double count for equal sizes
                        continue
                    sub += prob_subset_is_dircycle(A)*prob_subset_is_dircycle(B)
        else:
            for A in combinations(verts, k1):
                Aset=set(A)
                restv=[v for v in verts if v not in Aset]
                for B in combinations(restv, k2):
                    sub += prob_subset_is_dircycle(A)*prob_subset_is_dircycle(B)
        contributions[(k1,k2)] = sub
        total += sub
    return total, contributions

print("E[alpha2] per-subset decomposition:")
for n in range(3, 9):
    pairs = [(3,3)]
    if n >= 8: pairs += [(3,5)]
    tot, contr = E_alpha2(n, pairs)
    print(f" n={n}: total={tot}  parts={ {k:str(v) for k,v in contr.items()} }")
print("  claim: n=8 full E[alpha2] = 173/8 (3-3) + 447/32 (3-5) = 1139/32")
n=8
tot,contr=E_alpha2(8,[(3,3),(3,5)])
print(f"  computed n=8: 3-3={contr[(3,3)]} 3-5={contr[(3,5)]} total={tot}  ==1139/32? {tot==Fr(1139,32)}")

# ============ PART (ii): iso-class counts ============
def canon_tournament(adj, n):
    best=None
    verts=list(range(n))
    for perm in permutations(verts):
        bits=0
        for i in range(n):
            for j in range(n):
                if i!=j and adj[perm[i]][perm[j]]:
                    bits |= 1 << (i*n+j)
        if best is None or bits<best: best=bits
    return best

def all_tournaments(n):
    pairs=list(combinations(range(n),2))
    for bv in product((0,1),repeat=len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for (i,j),b in zip(pairs,bv):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

print()
print("iso-class counts (canonical form under S_n):")
for n in range(3,7):
    classes=set(); sc=set()
    for adj in all_tournaments(n):
        c=canon_tournament(adj,n); classes.add(c)
        # complement (transpose) tournament
        radj=[[adj[j][i] for j in range(n)] for i in range(n)]
        rc=canon_tournament(radj,n)
        if rc==c: sc.add(c)
    A=len(classes); SC=len(sc); Vm=(A+SC)//2
    print(f" n={n}: A000568={A}  SC(self-converse)={SC}  V_merged=(A+SC)/2={Vm}")
print("  claim: A000568=2,4,12,56 ; SC=2,2,8,12 ; V_merged=2,3,10,34 (n=3..6)")
