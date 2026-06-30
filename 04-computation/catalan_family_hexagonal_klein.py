#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE THREE NESTED CATALAN-FAMILY SEQUENCES + the hexagonal/A2 & Fibonacci connections (klein-S25).

seq3 = Catalan(n)          [A000108]  1,2,5,14,42
seq1 = C(2n,n-1)=n*Cat(n)  [A001791]  1,4,15,56,210
seq2 = C(2n+1,n-1)         [A002054]  1,5,21,84,330  = (n/2)*Cat(n+1)
Relations, the Catalan-triangle nesting, the Pascal/Fibonacci analogy, and a test of whether the Catalan
family counts walks in the A2 (hexagonal/triangular) Weyl chamber -- tying to the hexagonal covering bridge.
"""
from math import comb, isqrt, pi, sqrt
def cat(n): return comb(2*n,n)//(n+1)

print("="*84); print(" THE THREE SEQUENCES = the Catalan family modulated by linear factors"); print("="*84)
for n in range(1,7):
    s3=cat(n); s1=n*cat(n); s2=(n*cat(n+1))//2
    print(f"  n={n}: Catalan={s3:>4}  n*Catalan(seq1)={s1:>5}  (n/2)Cat(n+1)(seq2)={s2:>5}  central C(2n,n)={comb(2*n,n):>5} = Cat+seq1={s3+s1}")
print("  => seq3=Cat, seq1=n*Cat, seq2=(n/2)Cat(n+1); Cat+seq1 = central binomial C(2n,n).")

print("\n"+"="*84); print(" NESTING: the three are adjacent near-central Pascal columns C(2n, n-j), j=0,1,2"); print("="*84)
print("  (Catalan = C(2n,n)-C(2n,n-1) = the j=0 normalized 'reflected' count)")
for n in range(2,7):
    row=[comb(2*n, n-j) for j in range(0,3)]
    print(f"  n={n}: C(2n,n)={row[0]}, C(2n,n-1)={row[1]}(=seq1), C(2n,n-2)={row[2]}; seq2=C(2n+1,n-1)={comb(2*n+1,n-1)}=C(2n,n-1)+C(2n,n-2)")

print("\n"+"="*84); print(" FIBONACCI analogy: both families are read off PASCAL's triangle"); print("="*84)
def fib(n):
    a,b=0,1
    for _ in range(n): a,b=b,a+b
    return a
print(f"  Fibonacci = SHALLOW DIAGONAL SUMS of Pascal: F(n)=sum_k C(n-1-k,k); {[fib(n) for n in range(1,11)]}; rate phi=1.618")
print(f"  Catalan family = NEAR-CENTRAL ENTRIES of Pascal (the 'fat diagonal'); rate 4.")
print(f"  Both nested in Pascal; Fibonacci is the THIN (shallow) reading, Catalan the CENTRAL reading.")

print("\n"+"="*84); print(" A2 (hexagonal/triangular) Weyl-chamber walk test -- does the Catalan family count them?"); print("="*84)
# sl3 standard-rep weights in fundamental-weight coords: (1,0),(-1,1),(0,-1). Count walks 0->0 of length n
# staying DOMINANT (a>=0,b>=0) = multiplicity of trivial in V^{ox n}. Compare to Catalan family.
from collections import defaultdict
steps=[(1,0),(-1,1),(0,-1)]
def a2_walks(N):
    cur={(0,0):1}; out=[]
    for step in range(1,N+1):
        nxt=defaultdict(int)
        for (a,b),c in cur.items():
            for da,db in steps:
                na,nb=a+da,b+db
                if na>=0 and nb>=0: nxt[(na,nb)]+=c
        cur=nxt; out.append(cur.get((0,0),0))
    return out
w=a2_walks(12)
print(f"  sl3 A2 dominant walks 0->0 by length: {w}")
print(f"  (nonzero only at multiples of 3; the A2 'Catalan' = sub-sequence) ; Catalan = {[cat(n) for n in range(6)]}")
a2sub=[w[k] for k in range(2,12,3)]
print(f"  A2 walk sub-sequence (len 3,6,9,12): {a2sub}  -- compare A2-Catalan / our family")

print("\n"+"="*84); print(" HEXAGONAL COVERING OPTIMALITY (Kershner 1939) -- the continuous side of the bridge"); print("="*84)
theta=2*pi/sqrt(27)
print(f"  thinnest plane covering by congruent disks = the HEXAGONAL lattice, density theta = 2pi/sqrt(27) = {theta:.6f}")
print(f"  (Kershner 1939; Fejes Toth). This is the OPTIMAL lattice covering, wallpaper p6m -- a THEOREM.")
print(f"  So the continuous-side optimality is settled; the open bridge is: is the LRC covering the hexagonal one?")
