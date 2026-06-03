#!/usr/bin/env python3
"""Validity check: LRC-tight regular = INTERVAL circulant (additive), never Paley/QR
(multiplicative). interval vs Paley coincide at m=3, split at m>=7; Paley not round for
m>=7; the unique VT round tournament is the interval circulant; dichromatic number
(chi) of round tournaments. opus-2026-06-03-S591."""
from itertools import permutations, combinations
def circ(m,S):  # Cayley tournament on Z/m with connection set S (i->j iff j-i in S)
    S=set(x%m for x in S)
    return [[1 if (j-i)%m in S else 0 for j in range(m)] for i in range(m)]
def interval(m): return circ(m, range(1,(m-1)//2+1))
def QR(m): return sorted({(x*x)%m for x in range(1,m)})
def paley(m): return circ(m, QR(m))
def is_tournament(A,m): return all(A[i][j]+A[j][i]==1 for i in range(m) for j in range(m) if i!=j)
def is_round(A,m):
    # round/locally-transitive: each out-neighborhood is a contiguous arc under SOME rotation labeling
    # test: exists cyclic order where out-set of every vertex is an interval. Use: round <=> locally transitive
    # locally transitive: out-nbhd and in-nbhd each induce a transitive subtournament
    for v in range(m):
        out=[u for u in range(m) if A[v][u]]
        # transitive subtournament on out?
        if not transitive(A,out): return False
    return True
def transitive(A,verts):
    for a,b,c in permutations(verts,3):
        if A[a][b] and A[b][c] and A[c][a]: return False
    return True
def canon(A,m):
    best=None
    for p in permutations(range(m)):
        k=tuple(A[p[i]][p[j]] for i in range(m) for j in range(m))
        if best is None or k<best: best=k
    return best
def dichromatic(A,m):  # min # acyclic(transitive) parts covering V
    from itertools import product
    for k in range(1,m+1):
        for col in product(range(k),repeat=m):
            if set(col)!=set(range(k)): continue
            ok=True
            for c in range(k):
                part=[v for v in range(m) if col[v]==c]
                if not transitive(A,part): ok=False; break
            if ok: return k
    return m
def main():
    print("interval vs Paley (m=1 mod... Paley needs m=3 mod 4 prime):")
    for m in [3,7,11]:
        I=interval(m); P=paley(m)
        iso = canon(I,m)==canon(P,m)
        print(f"  m={m}: interval round={is_round(I,m)}; Paley round={is_round(P,m)}; "
              f"interval~=Paley(iso)={iso}; QR(m)={QR(m)} vs interval-S={list(range(1,(m-1)//2+1))}")
    print()
    print("dichromatic number chi of round tournaments (interval circulant) vs Paley:")
    for m in [3,5,7,9,11]:
        I=interval(m)
        print(f"  m={m}: chi(interval R_{m})={dichromatic(I,m)}"+ (f"; chi(Paley)={dichromatic(paley(m),m)}" if m in (3,7,11) else ""))
if __name__=='__main__': main()
