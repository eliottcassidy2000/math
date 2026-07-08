#!/usr/bin/env python3
r"""
lrc_chic_linearization_locus_kps_S76.py  (kind-pasteur-2026-07-07-S76)

THE LINEARIZATION LOCUS: chi_c(G_S) = 1/M(S)  <=>  mu(S) = M(S).

Framework (proved): for any finite S (distance graph G_S = Cay(Z,+-S)),
    1/mu(S)  =  chi_f(G_S)  <=  chi_c(G_S)  <=  1/M(S),
the LEFT eq because G_S is vertex-transitive (chi_f = 1/independence-ratio = 1/mu), and the
RIGHT because the LINEAR coloring c(x) = a*x mod N (t=a/N the loneliness witness,
M = m/N, min_s ||a s/N|| = m/N) is an (N,m)-coloring, p/q = N/m = 1/M.  Since mu >= M
always, the sandwich is nonempty; and

  * mu = M  =>  1/mu = 1/M  =>  chi_c = 1/M  (SQUEEZE, proved -- the circular rung is
    FAITHFUL exactly when Motzkin density meets loneliness);
  * mu > M  =>  CONJECTURE chi_c < 1/M (the linearization DEFECT lives exactly on the
    Haralambis mu>M separation locus).  Verified: GW (THM-658, chi_c<=13.5<14),
    Lucas {1,3,4,7} (chi_c<=25/6<5), {1,3,4,5} (chi_c<=4<4.5).

This file tabulates M, mu, and a quasi-periodic chi_c upper bound for a bank of small S,
testing the equivalence.
"""
import io, itertools
from math import gcd, floor
from fractions import Fraction as F
from contextlib import redirect_stdout
from pysat.solvers import Cadical153 as Sat
import networkx as nx

def cd(a,b,p): d=abs(a-b)%p; return min(d,p-d)

def M_of(S, Nmax=240):
    best=F(0)
    for N in range(2,Nmax+1):
        for a in range(1,N):
            if gcd(a,N)!=1: continue
            m=min(min((a*s)%N,N-(a*s)%N) for s in S)
            if F(m,N)>best: best=F(m,N)
    return best

def mu_of(S, Nmax=56):
    """max independent-set density = max_{N>2max S} alpha(C_N,S)/N (exact via clique of complement)."""
    mx=max(S); best=F(0)
    for N in range(2*mx+1, Nmax+1):
        G=nx.Graph(); G.add_nodes_from(range(N))
        for x in range(N):
            for s in S: G.add_edge(x,(x+s)%N)
        alpha=max(len(c) for c in nx.find_cliques(nx.complement(G)))
        if F(alpha,N)>best: best=F(alpha,N)
    return best

def quasi(S,p,q,T,Delta):
    def V(r,a): return r*p+a+1
    cls=[]
    for r in range(T):
        cls.append([V(r,a) for a in range(p)])
        for a in range(p):
            for b in range(a+1,p): cls.append([-V(r,a),-V(r,b)])
    for r in range(T):
        for s in S:
            j=r+s; carry=(j//T)*Delta; jm=j%T
            for a in range(p):
                for b in range(p):
                    if cd((b+carry)%p,a,p)<q: cls.append([-V(r,a),-V(jm,b)])
    with redirect_stdout(io.StringIO()):
        sv=Sat(bootstrap_with=cls); ok=sv.solve()
        col=[next(a for a in range(p) if V(r,a) in set(v for v in sv.get_model() if v>0)) for r in range(T)] if ok else None
        sv.delete()
    return col

def chic_ub(S, invM, Tmax=11):
    lo=int(floor(invM)); best=None
    for q in range(1,7):
        for p in range((lo-1)*q+1, int(invM*q)+1):
            if p/q>=invM-1e-12 or gcd(p,q)!=1 or p<=0: continue
            for T in range(1,Tmax):
                hit=False
                for Delta in range(0,p):
                    col=quasi(S,p,q,T,Delta)
                    if col is not None:
                        def C(x): return (col[x%T]+(x//T)*Delta)%p
                        if all(cd(C(x+s),C(x),p)>=q for x in range(-100,600) for s in S):
                            if best is None or p/q<best[1]: best=((p,q),p/q); hit=True; break
                if hit: break
    return best

INSTS = {
    '{1,2,3}':[1,2,3], '{1,2,5}':[1,2,5], '{2,3,4}':[2,3,4], '{1,3,5,7}':[1,3,5,7],
    '{1,4,6}':[1,4,6], 'Lucas{1,3,4,7}':[1,3,4,7], '{1,3,4,5}':[1,3,4,5],
    '{1,3,4,6}':[1,3,4,6], '{1,4,5,6}':[1,4,5,6], '{2,3,4,9}':[2,3,4,9],
}
print("THE LINEARIZATION LOCUS: chi_c(G_S) = 1/M(S) <=> mu(S) = M(S)?")
print(f"{'S':16s}{'M':>7s}{'1/M':>7s}{'mu':>7s}{'chi_f=1/mu':>11s}{'mu>M':>6s}{'chi_c UB':>12s}{'defect':>7s}")
n_ok = 0; n_tot = 0
for name,S in INSTS.items():
    M=M_of(S); invM=1/float(M); mu=mu_of(S); chif=1/float(mu)
    ub=chic_ub(S,invM)
    ubs=f'{ub[0][0]}/{ub[0][1]}={ub[1]:.2f}' if ub else f'>={invM:.2f}'
    defect = bool(ub) and ub[1] < invM-1e-9
    muGtM = float(mu) > float(M)
    n_tot += 1
    if defect == muGtM: n_ok += 1
    print(f"{name:16s}{str(M):>7s}{invM:7.2f}{str(mu):>7s}{chif:11.2f}{('YES' if muGtM else 'no(=)'):>6s}{ubs:>12s}{('YES' if defect else 'no'):>7s}")
print(f"\n  equivalence (defect <=> mu>M) holds on {n_ok}/{n_tot} instances")
print("  => CONJECTURE: chi_c(G_S)=1/M(S) iff mu(S)=M(S); the linearization defect is EXACTLY")
print("     the Haralambis mu>M separation locus. Easy half (mu=M => chi_c=1/M) PROVED by the")
print("     sandwich 1/mu <= chi_c <= 1/M. GW = the |S|=13 flagship (THM-658).")
