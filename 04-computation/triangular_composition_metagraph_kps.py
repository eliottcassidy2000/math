#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The triangular-number composition algebra vs the n=4 tournament metagraph.
kind-pasteur-2026-06-14-S3.  (User observation: T(x)=x*b(a(x)), a=+1, b=/2; the 4
subsets of {a,b} = the 4 tournament iso classes on 4 vertices.)

VERIFY the correspondence as a GRADED POSET:
  Boolean square 2^{a,b}:  rank |S|:  emptyset (0) < {a},{b} (1) < {a,b} (2)
  vs tournament classes on n=4 graded by H (= #directed Hamiltonian paths = OCF).
Also: score sequences, complementation (T->T^op) orbit structure (self-comp vs
complement-pair), to test the match (extremes = self-complementary, middle pair =
complement-orbit) and the involution (mis)match. Then n=5 (12 classes) to test
whether the Boolean-lattice match extends or is special to n=4.
"""
import sys, itertools
from math import comb
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def tournaments(n):
    """all labeled tournaments on n vertices as adjacency (i->j) bool matrices,
    one bit per pair (i<j): bit=1 means i->j else j->i."""
    pairs=list(itertools.combinations(range(n),2)); m=len(pairs)
    for bits in range(1<<m):
        adj=[[False]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: adj[i][j]=True
            else: adj[j][i]=True
        yield adj

def scores(adj,n): return tuple(sorted(sum(adj[i]) for i in range(n)))

def ham_paths(adj,n):
    """count directed Hamiltonian paths."""
    cnt=0
    for perm in itertools.permutations(range(n)):
        if all(adj[perm[k]][perm[k+1]] for k in range(n-1)): cnt+=1
    return cnt

def canon(adj,n):
    """canonical form under vertex relabeling (iso invariant)."""
    best=None
    for perm in itertools.permutations(range(n)):
        bits=0; k=0
        for i in range(n):
            for j in range(i+1,n):
                # edge between perm-positions
                pi,pj=perm[i],perm[j]
                if adj[pi][pj]: bits|=1<<k
                k+=1
        if best is None or bits<best: best=bits
    return best

def complement(adj,n):
    return [[adj[j][i] for j in range(n)] for i in range(n)]

def classify(n):
    seen={}
    for adj in tournaments(n):
        c=canon(adj,n)
        if c not in seen:
            seen[c]={'adj':adj,'H':ham_paths(adj,n),'score':scores(adj,n)}
    # complement pairing
    for c,info in seen.items():
        comp=canon(complement(info['adj'],n),n)
        info['comp']=comp
        info['selfcomp']=(comp==c)
    return seen

def main():
    for n in (4,5):
        cls=classify(n)
        print(f"\n=== n={n}: {len(cls)} iso classes (A000568); arc count C(n,2)=T(n-1)={comb(n,2)} ===",flush=True)
        # group by H
        byH={}
        for c,info in cls.items():
            byH.setdefault(info['H'],[]).append(info)
        print(f"   H-graded poset (rank = H):",flush=True)
        for H in sorted(byH):
            grp=byH[H]
            scs=[i['score'] for i in grp]
            sc_self=[(i['score'],'self-comp' if i['selfcomp'] else 'comp-pair') for i in grp]
            print(f"     H={H:2d}: {len(grp)} class(es)  scores={scs}  "
                  f"[{', '.join('SC' if i['selfcomp'] else 'pair' for i in grp)}]",flush=True)
        nsc=sum(1 for i in cls.values() if i['selfcomp'])
        print(f"   self-complementary: {nsc}; complement-pairs: {(len(cls)-nsc)//2}; "
              f"merged classes (A000568+SC)/2 = {(len(cls)+nsc)//2}",flush=True)

    print("\n=== the n=4 correspondence (user's observation) ===",flush=True)
    print("   Boolean square 2^{a,b}, rank=|S|:  ∅(0) < {a},{b}(1) < {a,b}(2)",flush=True)
    print("   tournament G_4, rank=H:            H=1 < H=3,3 < H=5",flush=True)
    print("   MATCH: rank-0 ∅ <-> transitive (H=1, self-comp);",flush=True)
    print("          rank-2 {a,b} <-> strong (H=5, self-comp);",flush=True)
    print("          rank-1 {a},{b} <-> the complement-pair (H=3 each).",flush=True)
    print("   => extremes (∅,{a,b}) = the 2 self-comp classes; middle = the comp-pair.",flush=True)

if __name__=="__main__":
    main()
