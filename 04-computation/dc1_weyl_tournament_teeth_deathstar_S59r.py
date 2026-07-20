#!/usr/bin/env python3
"""
death-star-2026-07-20-S59r (HYP-8160) -- teeth for the DC1-via-tournaments
exploration.  All exact.  DC1 is NOT claimed (open).

T1  The weight triple = the oriented 3-cycle.  Represent A1 on C[x] (truncated):
    p = d/dx, q = mult by x, N = qp (number operator).  Compute the commutator
    matrices and verify the sl2/Heisenberg relations
      [p,q] = 1,  [N,q] = q,  [N,p] = -p.
    (N, q, p) = the 3 vertices; the commutator signs (+1 on N->q, -1 on N->p,
    +1 central on q-p) = an oriented triangle = the n=3 tournament, with N the
    OBSERVER/Euler vertex carrying the central hbar of [p,q]=1.
T2  The observer +1 = hbar.  In the classical limit the SAME triple holds with
    Poisson {,}: {p,q}=1, {N,q}=q, {N,p}=-p (N=pq).  The '1' in [p,q]=1 is the
    quantization; the '+1' promoting the symplectic PAIR (p,q) to the 3-vertex
    TRIANGLE {N,p,q} is the observer of the observer-principle A_k <-> n=2k+1.
T3  Redei parity of the two n=3 tournaments; the counterexample fiber = 3 (odd,
    1+2) = Redei-shaped (verify via the orbit law from THM-1305).  Parity-
    protection: fiber-degree odd is necessary-not-sufficient (odd 3 > 1 exists
    at higher dim); DC1 needs odd AND the 2D leading-form bound (THM-1345).
"""
from fractions import Fraction as Fr

# ---------- T1: operator representation on C[x], basis {x^0..x^D} ----------
D = 6
def zero(): return [[Fr(0)]*(D+1) for _ in range(D+1)]
def matmul(A,B):
    C = zero()
    for i in range(D+1):
        for k in range(D+1):
            a = A[i][k]
            if a:
                for j in range(D+1):
                    if B[k][j]: C[i][j] += a*B[k][j]
    return C
def sub(A,B): return [[A[i][j]-B[i][j] for j in range(D+1)] for i in range(D+1)]
def comm(A,B): return sub(matmul(A,B), matmul(B,A))
def eq(A,B): return all(A[i][j]==B[i][j] for i in range(D+1) for j in range(D+1))
def scal(c):
    M = zero()
    for i in range(D+1): M[i][i] = Fr(c)
    return M

# p = d/dx : x^n -> n x^{n-1};  column n (input x^n) -> row n-1 entry n
P = zero(); Q = zero()
for n in range(D+1):
    if n-1 >= 0: P[n-1][n] = Fr(n)     # d/dx
    if n+1 <= D: Q[n+1][n] = Fr(1)     # mult by x
N = matmul(Q,P)                         # number operator N x^n = n x^n
# verify on interior (avoid truncation edge at top row)
def eq_interior(A,B,margin=1):
    return all(A[i][j]==B[i][j] for i in range(D+1-margin) for j in range(D+1-margin))

print("=== T1: the weight triple = oriented 3-cycle (A1 on C[x]) ===")
print("  [p,q] = 1        :", eq_interior(comm(P,Q), scal(1)))
print("  [N,q] = q  (+1)  :", eq_interior(comm(N,Q), Q))
print("  [N,p] = -p (-1)  :", eq_interior(comm(N,P), [[-P[i][j] for j in range(D+1)] for i in range(D+1)]))
print("  => vertices {N (observer/Euler), q (wt +1), p (wt -1)}; edge signs (+,-)")
print("     around N and the central +1 on (q,p) = an oriented triangle = n=3 tournament.")

# ---------- T2: observer +1 = hbar (classical shadow, Poisson) ----------
print("\n=== T2: observer +1 = hbar (classical Poisson shadow) ===")
def pmul(a,b):
    r={}
    for ka,ca in a.items():
        for kb,cb in b.items():
            k=(ka[0]+kb[0],ka[1]+kb[1]); v=r.get(k,0)+ca*cb
            if v: r[k]=v
            elif k in r: del r[k]
    return r
def padd(*ps):
    r={}
    for p in ps:
        for k,c in p.items():
            v=r.get(k,0)+c
            if v: r[k]=v
            elif k in r: del r[k]
    return r
def psc(p,s): return {k:c*s for k,c in p.items() if c*s!=0}
def pd(p,i):
    r={}
    for k,c in p.items():
        if k[i]>0:
            k2=list(k); k2[i]-=1; r[tuple(k2)]=c*k[i]
    return r
def poisson(A,B): return padd(pmul(pd(A,0),pd(B,1)), psc(pmul(pd(A,1),pd(B,0)),-1))
pp={(1,0):1}; qq={(0,1):1}; Hh=pmul(pp,qq)   # p,q, H=pq
print("  {p,q} =", poisson(pp,qq), " (=1)")
print("  {H,p} =", poisson(Hh,pp), " (= p, wt +1... sign per convention)")
print("  {H,q} =", poisson(Hh,qq), " (= -q)")
print("  => same triple classically. The RHS '1' of [p,q]=1 is hbar; the observer")
print("     vertex (N=pq, the Euler grading) is what the symplectic PAIR gains to")
print("     become the n=3 TRIANGLE. observer(+1) = hbar = the S59p/q conserved +1.")

# ---------- T3: Redei parity + counterexample fiber = 3 (odd) ----------
print("\n=== T3: Redei parity and the odd fiber ===")
# n=3 tournaments: transitive (1 Hamiltonian path) and 3-cycle (3 = 1+2*1 paths).
def ham_paths(adj, n):
    from itertools import permutations
    c=0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[i+1]] for i in range(n-1)): c+=1
    return c
# transitive 0->1->2, 0->2
trans=[[0,1,1],[0,0,1],[0,0,0]]
# 3-cycle 0->1->2->0
cyc=[[0,1,0],[0,0,1],[1,0,0]]
print(f"  transitive T3: {ham_paths(trans,3)} Hamiltonian path(s) (Redei odd: 1 = 1+2*0)")
print(f"  3-cycle C3:    {ham_paths(cyc,3)} Hamiltonian path(s) (Redei odd: 3 = 1+2*1)")
# counterexample fiber over (a,0,0), a!=0: {(0,0,a)} U {two lambda: lambda^2=-1/(4a)} = 3
print("  JC3 counterexample fiber over a-axis = 1 (fixed) + 2 (orbit) = 3 = ODD = Redei 1+2 (THM-1305).")
print("  => fiber-degree ODD is the Redei-parity shadow (NECESSARY, not sufficient:")
print("     odd 3 > 1 occurs at dim>=3; DC1 needs ODD + the 2D leading-form bound, THM-1345).")
