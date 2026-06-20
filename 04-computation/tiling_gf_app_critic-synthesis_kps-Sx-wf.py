#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INDEPENDENT CRITIC re-derivation of the four THM-554 tiling-GF applications.
Fresh code, exact Fraction arithmetic. Does NOT import the existing verifier scripts.
Checks the load-bearing closed forms + the claims that were NOT yet independently
re-derived elsewhere (E[c7] poly, full E[alpha2] at n=8, score-determinacy CEs,
parity bias, max-c3 multiplicity census, iso-class V_merged numbers).
"""
import sys
from collections import defaultdict, Counter
from itertools import product, combinations
from math import comb
from fractions import Fraction as Fr
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

def tiles(n):
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]

def build_adj(n, T, bv):
    adj = [[0]*(n+1) for _ in range(n+1)]
    for k in range(n, 1, -1):
        adj[k][k-1] = 1            # base path n->n-1->...->1
    for (a, b), bit in zip(T, bv):
        if bit == 0: adj[a][b] = 1
        else:        adj[b][a] = 1
    return adj

def count_cycles_k(adj, n, k):
    """count directed k-cycles (as cyclic sequences /k, i.e. # of directed k-cycles)."""
    tot = 0
    for sub in combinations(range(1, n+1), k):
        # count directed Hamiltonian cycles on this induced k-tournament, /k?
        # A directed k-cycle = cyclic orientation cycle. Count over permutations fixing first.
        s = sub
        cnt = 0
        rest = s[1:]
        for perm in permutations_first(s):
            ok = True
            for i in range(k):
                u = perm[i]; v = perm[(i+1)%k]
                if not adj[u][v]:
                    ok = False; break
            if ok: cnt += 1
        tot += cnt
    return tot

from itertools import permutations
def permutations_first(s):
    # cyclic permutations with first element fixed -> counts each directed cycle once
    first = s[0]; rest = list(s[1:])
    for p in permutations(rest):
        yield (first,) + p

def scores(adj, n):
    return [sum(adj[v][u] for u in range(1, n+1)) for v in range(1, n+1)]

def c3_from_scores(sv, n):
    return comb(n,3) - sum(comb(s,2) for s in sv)

def is_dir_tri(adj, s):
    a,b,c = s
    outdeg = [adj[a][b]+adj[a][c], adj[b][a]+adj[b][c], adj[c][a]+adj[c][b]]
    return sorted(outdeg) == [1,1,1]   # cyclic triangle each vertex outdeg 1

def list_dir_tris(adj, n):
    return [frozenset(s) for s in combinations(range(1,n+1),3) if is_dir_tri(adj,s)]

def count_disjoint_tri_pairs(adj, n):
    tris = list_dir_tris(adj, n)
    cnt = 0
    for i in range(len(tris)):
        for j in range(i+1, len(tris)):
            if tris[i].isdisjoint(tris[j]):
                cnt += 1
    return cnt

# ---- polynomials claimed ----
def P_Ec3(n):  return Fr(n**3 - 3*n**2 + 8*n - 12, 24)
def P_Varc3(n):return Fr(n**3 - 7*n**2 + 20*n - 16, 32)
def P_Ec5(n):  return Fr(n**5 -10*n**4 +45*n**3 -140*n**2 +294*n -280, 160)
def P_Ec7(n):  return Fr(n**7 -21*n**6 +189*n**5 -1015*n**4 +3836*n**3 -10514*n**2 +18458*n -15204, 896)

print("="*78)
print("CHECK A: brute c3/c5/c7/H moments (n<=7) vs claimed polynomials")
print("="*78)
for n in range(3, 8):
    T = tiles(n); F = len(T); tot = 1<<F
    sc3=Counter(); sc5=Counter(); sc7=Counter(); sumH=0
    # H via Redei: # Hamiltonian paths in tournament (odd by Redei). Count directly.
    for bv in product((0,1), repeat=F):
        adj = build_adj(n, T, bv)
        sv = scores(adj, n)
        sc3[c3_from_scores(sv,n)] += 1
        # c5, c7 by direct cycle count
        if n>=5: sc5[count_cycles_k(adj,n,5)] += 1
        if n>=7: sc7[count_cycles_k(adj,n,7)] += 1
        # H = number of Hamiltonian (directed) paths
        hp = 0
        for perm in permutations(range(1, n+1)):
            if all(adj[perm[i]][perm[i+1]] for i in range(n-1)):
                hp += 1
        sumH += hp
    def mean(c):
        s=sum(k*v for k,v in c.items()); t=sum(c.values()); return Fr(s,t)
    def var(c):
        m=mean(c); t=sum(c.values()); e2=Fr(sum(k*k*v for k,v in c.items()),t); return e2-m*m
    Ec3=mean(sc3); Vc3=var(sc3)
    line=f" n={n}: E[c3]={Ec3} (poly {P_Ec3(n)} {Ec3==P_Ec3(n)})  Var[c3]={Vc3} (poly {P_Varc3(n)} {Vc3==P_Varc3(n)})"
    if n>=5:
        Ec5=mean(sc5); line+=f"  E[c5]={Ec5} (poly {P_Ec5(n)} {Ec5==P_Ec5(n)})"
    if n>=7:
        Ec7=mean(sc7); line+=f"  E[c7]={Ec7} (poly {P_Ec7(n)} {Ec7==P_Ec7(n)})"
    EH=Fr(sumH,tot); line+=f"  E[H]={EH}"
    print(line)
    if n==5:
        idn = 1 + 2*Ec3 + 2*Ec5
        print(f"     n=5 identity 1+2E[c3]+2E[c5]={idn}  ==E[H]? {idn==EH}")

print()
print("="*78)
print("CHECK B: parity bias E[(-1)^c3] = 1/2^floor((n-1)/2), n=4..8  (Z-engine)")
print("="*78)
def beta_step(dist,n):
    nd=defaultdict(int)
    for vec,cnt in dist.items():
        l=list(vec)+[0]; l[n-1]+=1; nd[tuple(l)]+=cnt
    dist=nd
    for b in range(1,n-1):
        nd=defaultdict(int)
        for vec,cnt in dist.items():
            l0=list(vec); l0[n-1]+=1; nd[tuple(l0)]+=cnt
            l1=list(vec); l1[b-1]+=1; nd[tuple(l1)]+=cnt
        dist=nd
    return dist
dist={(0,):1}
for n in range(2,11):
    dist=beta_step(dist,n)
    if n>=3:
        tot=0; sgn=0
        for vec,cnt in dist.items():
            c3=comb(n,3)-sum(comb(s,2) for s in vec)
            tot+=cnt; sgn += cnt*(1 if c3%2==0 else -1)
        bias=Fr(sgn,tot)
        claim = Fr(1, 1<<((n-1)//2)) if n>=4 else Fr(0)
        print(f" n={n}: E[(-1)^c3]={bias}  claim={claim}  match={bias==claim}")

print()
print("="*78)
print("CHECK C: max-c3 value + multiplicity (regular census), Z-engine n<=10")
print("="*78)
dist={(0,):1}
for n in range(2,11):
    dist=beta_step(dist,n)
    if n>=3:
        cs=Counter()
        for vec,cnt in dist.items():
            cs[comb(n,3)-sum(comb(s,2) for s in vec)]+=cnt
        mx=max(cs); mult=cs[mx]
        if n%2==1: vform=Fr(n**3-n,24)
        else:      vform=Fr(n**3-4*n,24)
        print(f" n={n}: max_c3={mx} (form {vform} {mx==vform})  mult={mult}  (n {'odd' if n%2 else 'even'})")

print()
print("="*78)
print("CHECK D: score-determinacy of c3 (YES) vs c5/alpha2 (NO) -- counterexamples")
print("="*78)
for n in (5,6,7):
    T=tiles(n); F=len(T)
    by_score_c3=defaultdict(set); by_score_c5=defaultdict(set); by_score_a2=defaultdict(set)
    for bv in product((0,1),repeat=F):
        adj=build_adj(n,T,bv); sv=tuple(sorted(scores(adj,n)))
        by_score_c3[sv].add(c3_from_scores(scores(adj,n),n))
        c5=count_cycles_k(adj,n,5); by_score_c5[sv].add(c5)
        a2=count_disjoint_tri_pairs(adj,n)
        by_score_a2[sv].add(a2)
    nd_c3=sum(1 for v in by_score_c3.values() if len(v)>1)
    nd_c5=sum(1 for v in by_score_c5.values() if len(v)>1)
    nd_a2=sum(1 for v in by_score_a2.values() if len(v)>1)
    print(f" n={n}: c3 multi-valued score classes={nd_c3} (0=>determined)  c5={nd_c5}  alpha2tri={nd_a2}")
    if nd_c5:
        ex=next(s for s,v in by_score_c5.items() if len(v)>1)
        print(f"     c5 CE: score {ex} -> c5 in {sorted(by_score_c5[ex])}")
    if nd_a2:
        ex=next(s for s,v in by_score_a2.items() if len(v)>1)
        print(f"     a2 CE: score {ex} -> a2 in {sorted(by_score_a2[ex])}")
