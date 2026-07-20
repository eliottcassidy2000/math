#!/usr/bin/env python3
"""
kind-pasteur-2026-07-20-S128c99 -- HYP-8130: THE REFINED k=4/5 HUNT (backlog c98(i)).

Generalized design ansatz (relaxations over the c98 solver, which was EMPTY at
k=4,5): A = (u^k, kappa*u^{k-1}*x, -x^k) with kappa FREE (was kappa=k);
B1 = u*W;  B2 = S + kappa*x*W with S a FREE low-degree polynomial (was S=y);
B3 = (R - x^k*u*W)/u^k with R in a WIDER box span{x^m u^n} including the
x^2u-stratification (base-multiplicity tuning: Bezout 2k - mult must equal the
field degree).  Then impose u^k | (R - x^k uW) [linear] and det J == nonzero
constant [quadratic].  Also a two-term-middle variant A2 = c1 u^3 x + c2 u^2 x^2
at k=4 (B1 = u^2 W design).  k=3 rerun with kappa free = extended control.
Empties here + fleet empties = evidence for the ORDER-{1,3} line-congruence
conjecture: a z-affine Keller map of C^3 has field degree 1 or 3.
"""
import sympy as sp
from sympy import symbols, expand, Poly, cancel, together, fraction
import time

x, y, z, v = symbols('x y z v')
u = 1 + x*y

def hunt(k, Wdeg, Sdeg, Rmx, Rnx, label, two_term=False):
    t0 = time.time()
    print(f"\n==== {label}: k={k}, Wdeg<={Wdeg}, Sdeg<={Sdeg}, R: m<={Rmx}, n<={Rnx}, two_term={two_term} ====", flush=True)
    kap = sp.symbols('kap')
    c1, c2 = sp.symbols('c1 c2')
    W = sp.Integer(0); syms = []
    for i in range(Wdeg+1):
        for j in range(Wdeg+1-i):
            cs = sp.symbols(f'w{i}_{j}'); syms.append(cs); W += cs*x**i*y**j
    S = sp.Integer(0)
    for i in range(Sdeg+1):
        for j in range(Sdeg+1-i):
            cs = sp.symbols(f's{i}_{j}'); syms.append(cs); S += cs*x**i*y**j
    R = sp.Integer(0)
    for m_ in range(Rmx+1):
        for n_ in range(Rnx+1):
            cs = sp.symbols(f'r{m_}_{n_}'); syms.append(cs); R += cs*x**m_*u**n_
    if two_term:
        upow = 2
        A = [u**k, c1*u**3*x + c2*u**2*x**2, -x**k]
        B1 = u**upow*W
        B2 = S + x*W*(c1*u + c2*x)
        extra = [c1, c2]
    else:
        upow = 1
        A = [u**k, kap*u**(k-1)*x, -x**k]
        B1 = u*W
        B2 = S + kap*x*W
        extra = [kap]
    # divisibility u^k | (R - x^k * B1)
    P = sp.expand(R - x**k*B1)
    Pv = sp.together(P.subs(y, (v-1)/x))
    num, den = sp.fraction(sp.cancel(Pv))
    num = sp.expand(num)
    pn = sp.Poly(num, v)
    lin = []
    for j in range(min(k, pn.degree()+1)):
        cj = pn.nth(j)
        for c_ in sp.Poly(cj, x).all_coeffs():
            if c_ != 0: lin.append(c_)
    lin = list(dict.fromkeys(lin))
    sol1 = sp.solve(lin, syms, dict=True)
    if not sol1:
        if lin:
            print(f"  divisibility unsolvable -- EMPTY", flush=True); return
        sub1 = {}
    else:
        sub1 = sol1[0]
    W1, S1, R1, B11, B21 = [t.subs(sub1) for t in (W, S, R, B1, B2)]
    B31 = sp.cancel((R1 - x**k*B11)/u**k)
    b3n, b3d = sp.fraction(sp.together(B31))
    if b3d.free_symbols & {x, y}:
        print(f"  B3 not polynomial (den {b3d}) -- abort", flush=True); return
    F1 = A[0]*z + B11; F2 = A[1]*z + B21; F3 = A[2]*z + B31
    rem = sorted(set().union(*[e.free_symbols for e in (F1, F2, F3)]) - {x, y, z}, key=str)
    print(f"  free params after divisibility: {len(rem)}", flush=True)
    J = sp.Matrix([F1, F2, F3]).jacobian([x, y, z])
    d = sp.expand(J.det())
    dc = d - d.subs({x:0, y:0, z:0})
    eqs = [c_ for c_ in sp.Poly(dc, x, y, z).coeffs() if c_ != 0]
    eqs = list(dict.fromkeys(eqs))
    print(f"  Keller system: {len(eqs)} eqs, {len(rem)} unknowns [{time.time()-t0:.0f}s]", flush=True)
    try:
        sols = sp.solve(eqs, rem, dict=True)
    except Exception as e:
        print(f"  solve error: {e}; groebner fallback...", flush=True)
        try:
            G = sp.groebner(eqs, *rem)
            sols = sp.solve(list(G.exprs), rem, dict=True)
        except Exception as e2:
            print(f"  groebner failed: {e2} -- UNRESOLVED", flush=True); return
    print(f"  branches: {len(sols)}", flush=True)
    nz = 0
    for so in sols:
        dv = sp.simplify(d.subs(so))
        if dv != 0:
            nz += 1
            print(f"   *** NONZERO-DET BRANCH: det = {dv}", flush=True)
            print(f"       W = {sp.expand(W1.subs(so))}", flush=True)
            print(f"       S = {sp.expand(S1.subs(so))}", flush=True)
            print(f"       R = {sp.factor(R1.subs(so))}", flush=True)
            print(f"       extra = {[ (str(e), so.get(e)) for e in extra if e in so ]}", flush=True)
    if nz == 0:
        print(f"  ALL branches det==0 or none -- EMPTY at this box [{time.time()-t0:.0f}s]", flush=True)

# extended control: k=3 with kappa free and S free -- must re-find kappa=3, S=y, owner W
hunt(3, 4, 2, 3, 3, "k=3 EXTENDED CONTROL (kappa, S free)")
# the refined k=4 hunts
hunt(4, 5, 2, 5, 4, "k=4 REFINED (kappa free, wide R)")
hunt(4, 5, 2, 4, 4, "k=4 TWO-TERM MIDDLE", two_term=True)
# k=5 with kappa free
hunt(5, 6, 2, 5, 4, "k=5 REFINED (kappa free)")
print("\nDONE.", flush=True)
