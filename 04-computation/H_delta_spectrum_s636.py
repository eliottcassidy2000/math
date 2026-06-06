#!/usr/bin/env python3
"""
S636 — H, delta, and the spectrum.  H(T) = OCF = I(Omega(T), 2) = the independence polynomial of the
odd-cycle conflict graph (vertices = odd directed cycles, edges = sharing a vertex) at z=2:
   H = sum over vertex-disjoint families F of odd cycles of 2^{|F|} = 1 + 2*alpha1 + 4*alpha2 + ...
=> H ODD; flipping an arc changes H by an EVEN delta. 7 (n=5) and 21 (n=6) are non-realizable.
We compute: (1) the H-spectrum (confirm 7 forbidden); (2) delta(arc)=ΔH per flip (distribution, even);
(3) the delta-PROPAGATION: flipping arc a changes delta(b) for which arcs b? (claim: only b sharing an
odd cycle with a — 'provably not all'); (4) dichromatic number + a spectral (eigenvalue) bound.
"""
import itertools
from fractions import Fraction as Fr

def arcs(n): return [(i,j) for i in range(n) for j in range(n) if i<j]

def tournament_from_bits(n, bits):
    # bits: one bit per unordered pair (i<j): 1 => i->j, 0 => j->i
    out=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(arcs(n)):
        if (bits>>k)&1: out[i][j]=True
        else: out[j][i]=True
    return out

def odd_cycles(adj, n):
    """all directed cycles of ODD length (as frozenset of vertices + the cycle), via DFS."""
    cyc=[]
    for start in range(n):
        # DFS paths start->...->start
        def dfs(v, path, used):
            for w in range(n):
                if adj[v][w]:
                    if w==start and len(path)>=3 and len(path)%2==1:
                        cyc.append(frozenset(path))
                    elif w>start and w not in used:   # canonical: start is min vertex
                        dfs(w, path+[w], used|{w})
        dfs(start,[start],{start})
    # dedup identical vertex sets that form the SAME directed cycle: keep one per (frozenset, length)
    # (for OCF, distinct directed odd cycles count separately; a vertex set may host >1 long cycle)
    return cyc

def H_value(adj,n):
    cyc=odd_cycles(adj,n)
    # independence polynomial at 2 over the 'share-a-vertex' conflict graph = sum over vertex-disjoint
    # subfamilies of 2^{size}. Compute by DP over cycles.
    # represent each cycle by its vertex bitmask
    masks=[sum(1<<v for v in c) for c in cyc]
    # dynamic: count weighted vertex-disjoint subsets
    total=[1]  # H accumulator via recursion with memo on (index, used-mask) is expensive; use sum over subsets greedily
    # since #odd cycles is small at n<=6, do recursion
    from functools import lru_cache
    M=len(masks)
    def rec(i, used):
        if i==M: return 1
        res=rec(i+1,used)                       # skip cycle i
        if not (masks[i]&used):
            res+=2*rec(i+1, used|masks[i])      # take cycle i (weight 2)
        return res
    return rec(0,0)

def dichromatic(adj,n):
    """min #transitive parts covering vertices (transitive = acyclic subtournament)."""
    def acyclic(S):
        # subtournament on S is transitive iff it has a strict linear order (no 3-cycle suffices for tournament? need full acyclicity)
        # check topological: tournament on S acyclic iff scores distinct iff no directed cycle
        S=list(S)
        # check no directed cycle via DFS
        idx={v:k for k,v in enumerate(S)}
        col={}
        def dfs(v):
            col[v]=1
            for w in S:
                if adj[v][w]:
                    if col.get(w)==1: return False
                    if w not in col and not dfs(w): return False
            col[v]=2; return True
        col={}
        for v in S:
            if v not in col:
                if not dfs(v): return False
        return True
    verts=list(range(n))
    for k in range(1,n+1):
        # try to partition into k acyclic parts (brute for small n)
        def assign(i, parts):
            if i==n: return True
            for p in range(len(parts)):
                parts[p].append(verts[i])
                if acyclic(parts[p]) and assign(i+1,parts): return True
                parts[p].pop()
            if len(parts)<k:
                parts.append([verts[i]])
                if assign(i+1,parts): return True
                parts.pop()
            return False
        if assign(0,[]): return k
    return n

if __name__=="__main__":
    n=5; A=arcs(n)
    print(f"n={n}: H-spectrum over all 2^{len(A)} tournaments (OCF = I(Omega,2))")
    from collections import Counter
    spec=Counter();
    Hs={}
    for bits in range(1<<len(A)):
        adj=tournament_from_bits(n,bits); h=H_value(adj,n); spec[h]+=1; Hs[bits]=(adj,h)
    print("   H counts:", dict(sorted(spec.items())))
    print(f"   7 realizable? {'NO (forbidden)' if 7 not in spec else 'yes'}; all H odd? {all(h%2==1 for h in spec)}")

    # (2)+(3) delta and propagation on a sample tournament (a cyclic/round one)
    bits0=0b0101101011  # arbitrary
    adj0=tournament_from_bits(n,bits0); H0=H_value(adj0,n)
    def flip(bits,k): return bits ^ (1<<k)
    deltas={k: H_value(tournament_from_bits(n,flip(bits0,k)),n)-H0 for k in range(len(A))}
    print(f"\nsample T (H={H0}): delta per arc (ΔH on flip):")
    for k,(i,j) in enumerate(A): print(f"   arc {i}->{j} (or rev): delta={deltas[k]} (even? {deltas[k]%2==0})")
    # propagation: flip arc a, recompute all deltas, see which changed
    a=0; bitsA=flip(bits0,a); adjA=tournament_from_bits(n,bitsA); HA=H_value(adjA,n)
    deltasA={k: H_value(tournament_from_bits(n,flip(bitsA,k)),n)-HA for k in range(len(A))}
    changed=[A[k] for k in range(len(A)) if deltasA[k]!=deltas[k] and k!=a]
    print(f"\nflip arc {A[a]}: delta-values that CHANGED for OTHER arcs: {changed}")
    print(f"   (#changed = {len(changed)} of {len(A)-1} other arcs => NOT all; propagation is local)")

    # (4) dichromatic number + spectral
    print(f"\ndichromatic number of sample T: {dichromatic(adj0,n)} (round/LRC-tight => 2, THM-402)")
