#!/usr/bin/env python3
r"""
lrc_chic_defect_sweep_kps_S77.py  (kind-pasteur-2026-07-07-S77)

Solidify the defect-half conjecture (mu>M => chi_c < 1/M) by a broad sweep, and MAP which
mu>M instances are the nontrivial HARD case (chi = 1/M) vs the trivial EASY case (chi < 1/M).

For each S: M, mu, integer chi, and a sub-1/M circular coloring (confirming chi_c < 1/M):
  - EASY (chi < 1/M): the integer chi-coloring already gives chi_c <= chi < 1/M.
  - HARD (chi = 1/M = n): test the (2n-1,2)-coloring => chi_c <= n-1/2 < n.
Report any counterexample (mu>M but no sub-1/M coloring found).
"""
import io, itertools
from math import gcd
from fractions import Fraction as F
from contextlib import redirect_stdout
from pysat.solvers import Cadical153 as Sat
import networkx as nx

def cd(a,b,p): d=abs(a-b)%p; return min(d,p-d)
def M_of(S,Nmax=140):
    best=F(0)
    for N in range(2,Nmax+1):
        for a in range(1,N):
            if gcd(a,N)!=1: continue
            m=min(min((a*s)%N,N-(a*s)%N) for s in S)
            if F(m,N)>best: best=F(m,N)
    return best
def mu_of(S,Nmax=46):
    mx=max(S); best=F(0)
    for N in range(2*mx+1,Nmax+1):
        G=nx.Graph(); G.add_nodes_from(range(N))
        for x in range(N):
            for s in S: G.add_edge(x,(x+s)%N)
        a=max(len(c) for c in nx.find_cliques(nx.complement(G)))
        if F(a,N)>best: best=F(a,N)
    return best
def has_coloring(S,p,q,Tmax=7):
    def V(r,a): return r*p+a+1
    for T in range(1,Tmax+1):
        for D in range(0,p):
            cls=[]
            for r in range(T):
                cls.append([V(r,a) for a in range(p)])
                for a in range(p):
                    for b in range(a+1,p): cls.append([-V(r,a),-V(r,b)])
            for r in range(T):
                for s in S:
                    j=r+s; carry=(j//T)*D; jm=j%T
                    for a in range(p):
                        for b in range(p):
                            if cd((b+carry)%p,a,p)<q: cls.append([-V(r,a),-V(jm,b)])
            with redirect_stdout(io.StringIO()):
                sv=Sat(bootstrap_with=cls); ok=sv.solve()
                col=[next(a for a in range(p) if V(r,a) in set(v for v in sv.get_model() if v>0)) for r in range(T)] if ok else None
                sv.delete()
            if col is not None:
                def C(x,col=col,T=T,D=D): return (col[x%T]+(x//T)*D)%p
                if all(cd(C(x+s),C(x),p)>=q for x in range(-50,350) for s in S):
                    return (T,D,col)
    return None
def integer_chi(S,invM):
    import math
    for n in range(2, math.ceil(invM)+2):
        if has_coloring(S,n,1): return n
    return None

# enumerate primitive sets (gcd 1, contains small elements), |S|=3,4, max<=8
sets=[]
for k in (3,4):
    for c in itertools.combinations(range(1,9),k):
        g=0
        for x in c: g=gcd(g,x)
        if g==1: sets.append(list(c))
print(f"sweeping {len(sets)} primitive sets |S|=3,4 (max<=8)")
print(f"{'S':16s}{'M':>7s}{'1/M':>6s}{'mu':>7s}{'chi':>4s}{'case':>6s}{'defect witness':>22s}")
hard=[]; counterex=[]; nmu=0
for S in sets:
    M=M_of(S); invM=1/float(M); mu=mu_of(S)
    if float(mu)<=float(M)+1e-12: continue   # only mu>M
    nmu+=1
    chi=integer_chi(S,invM)
    if chi is None: continue
    if chi < invM-1e-9:
        case='easy'; wit=f'chi={chi}<{invM:.2f}'
    else:
        # hard: chi = 1/M (integer). test (2chi-1,2)
        n=chi; res=has_coloring(S,2*n-1,2)
        case='HARD'
        if res: wit=f'({2*n-1},2) SAT: chi_c<={2*n-1}/2'
        else:
            wit='*** NO sub-1/M found ***'; counterex.append((tuple(S),chi,mu,M))
        hard.append((tuple(S),n,res is not None))
    print(f"{str(S):16s}{str(M):>7s}{invM:6.2f}{str(mu):>7s}{str(chi):>4s}{case:>6s}{wit:>22s}")
print(f"\nmu>M sets: {nmu}; HARD (chi=1/M): {len(hard)}; hard closed by (2n-1,2): {sum(1 for h in hard if h[2])}")
print(f"COUNTEREXAMPLES (mu>M, no sub-1/M coloring): {len(counterex)} {counterex if counterex else '(none -- conjecture holds)'}")
