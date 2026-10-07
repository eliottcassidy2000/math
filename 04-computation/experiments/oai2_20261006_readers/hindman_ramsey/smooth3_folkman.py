#!/usr/bin/env python3
"""Folkman / Hindman FS u FP with generators in S = 3-smooth numbers, 2 colours, SAT."""
import itertools, sys
from pysat.solvers import Solver
def smooth3(N):
    out=[]; a=1
    while a<=N:
        b=a
        while b<=N:
            out.append(b); b*=3
        a*=2
    return sorted(out)
def Omega(n):
    c=0
    while n%2==0: n//=2; c+=1
    while n%3==0: n//=3; c+=1
    return c
def lam(n): return 1-2*(Omega(n)%2)

# A. witness 2-colouring of [1..12] for distinct-summand Schur
def schur_witness(N):
    S=smooth3(N)
    s=Solver(name='cadical195'); trip=[]
    for i,a in enumerate(S):
        for b in S[i+1:]:
            if a+b>N: break
            trip.append((a,b,a+b))
            s.add_clause([a,b,a+b]); s.add_clause([-a,-b,-(a+b)])
    for v in range(1,N+1): s.add_clause([v,-v])
    ok=s.solve(); m=s.get_model() if ok else None
    # enumerate all solutions restricted to involved vertices
    sols=[]
    if ok:
        inv=sorted({v for t in trip for v in t})
        while s.solve():
            m=s.get_model(); sol=tuple(1 if m[v-1]>0 else 0 for v in inv); sols.append(sol)
            s.add_clause([-m[v-1] for v in inv])
        return inv, sols
    return None, []
inv,sols=schur_witness(12)
print("A. distinct-summand Schur, N=12: involved", inv, "#solutions", len(sols))
for so in sols[:8]: print("   ", dict(zip(inv,so)))

# B. hybrid colouring check: chi = lambda on S, -1 off S; no mono {a,b,a+b,ab}, a,b in S (a==b allowed)
N=10**15; S=smooth3(N); Sset=set(S)
def chi(n): return lam(n) if n in Sset else -1
bad=0; tot=0
for i,a in enumerate(S):
    for b in S[i:]:
        if a*b>N: break
        tot+=1
        if chi(a)==chi(b)==chi(a+b)==chi(a*b): bad+=1
print(f"B. hybrid colouring: mono {{a,b,a+b,ab}} with a<=b in S, ab<=1e15: {bad} of {tot}")
# 3-colouring: lambda on S, 0 off S; Schur triples a,b in S
bad3=0;tot3=0
for i,a in enumerate(S):
    for b in S[i:]:
        if a+b>N: break
        tot3+=1
        c=a+b
        col=lambda n:(lam(n) if n in Sset else 0)
        if col(a)==col(b)==col(c): bad3+=1
print(f"   3-colouring lambda|S, 0 off S: mono Schur triples a<=b in S, a+b<=1e15: {bad3} of {tot3}")

# C. Folkman m=3 with generators in S, two colours
def folkman3(N, allow_dup_sums=True):
    S=smooth3(N)
    s=Solver(name='cadical195'); nc=0
    for a,b,c in itertools.combinations(S,3):
        if a+b+c>N: continue
        fs=sorted({a,b,c,a+b,a+c,b+c,a+b+c})
        s.add_clause(fs); s.add_clause([-x for x in fs]); nc+=1
    ok=s.solve(); m=s.get_model() if ok else None
    s.delete(); return ok,nc,m
for e in [6,8,10,12,14,16,18,20]:
    N=2**e
    ok,nc,m=folkman3(N)
    print(f"C. Folkman m=3, generators 3-smooth, N=2^{e}: {'SAT (2-colourable)' if ok else 'UNSAT'}  configs={nc}")
    if not ok: break
