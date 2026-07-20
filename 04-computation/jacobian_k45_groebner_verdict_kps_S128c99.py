#!/usr/bin/env python3
"""
kind-pasteur-2026-07-20-S128c99 -- HYP-8130: k=4/5 refined-box VERDICT via
Rabinowitsch + Groebner mod p = 32003.  GB == [1] for (Keller eqs + t*det - 1)
<=> NO solution with det != 0 over GF(32003)-bar => (good-prime specialization)
strong evidence of emptiness over C in the box.  Control: k=3 must be
CONSISTENT (GB != [1]) -- the owner map lives there.
"""
import sympy as sp
from sympy import symbols, expand, Poly, together, fraction, cancel
import time

x, y, z, v, t_ = symbols('x y z v t_')
u = 1 + x*y
PRIME = 32003

def verdict(k, Wdeg, Sdeg, Rmx, Rnx, label, two_term=False):
    t0 = time.time()
    print(f"\n==== {label}: k={k} ====", flush=True)
    kap = sp.symbols('kap'); c1, c2 = sp.symbols('c1 c2')
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
        A = [u**k, c1*u**3*x + c2*u**2*x**2, -x**k]; B1 = u**2*W; B2 = S + x*W*(c1*u + c2*x); extra = [c1, c2]
    else:
        A = [u**k, kap*u**(k-1)*x, -x**k]; B1 = u*W; B2 = S + kap*x*W; extra = [kap]
    P = sp.expand(R - x**k*B1)
    Pv = sp.together(P.subs(y, (v-1)/x)); num, den = sp.fraction(sp.cancel(Pv)); num = sp.expand(num)
    pn = sp.Poly(num, v)
    lin = []
    for j in range(min(k, pn.degree()+1)):
        cj = pn.nth(j)
        for c_ in sp.Poly(cj, x).all_coeffs():
            if c_ != 0: lin.append(c_)
    lin = list(dict.fromkeys(lin))
    sol1 = sp.solve(lin, syms, dict=True)
    sub1 = sol1[0] if sol1 else {}
    W1, S1, R1, B11, B21 = [e.subs(sub1) for e in (W, S, R, B1, B2)]
    B31 = sp.cancel((R1 - x**k*B11)/u**k)
    F1 = A[0]*z + B11; F2 = A[1]*z + B21; F3 = A[2]*z + B31
    rem = sorted(set().union(*[e.free_symbols for e in (F1, F2, F3)]) - {x, y, z}, key=str)
    J = sp.Matrix([F1, F2, F3]).jacobian([x, y, z])
    d = sp.expand(J.det())
    d0 = d.subs({x:0, y:0, z:0})
    dc = d - d0
    eqs = [c_ for c_ in sp.Poly(dc, x, y, z).coeffs() if c_ != 0]
    eqs = list(dict.fromkeys(eqs))
    print(f"  system: {len(eqs)} eqs, {len(rem)} unknowns; det(0) = d0 (Rabinowitsch on d0) [{time.time()-t0:.0f}s]", flush=True)
    gens = rem + [t_]
    def clearden(e):
        _, pe = sp.Poly(sp.expand(e), *gens, domain='QQ').clear_denoms()
        return pe.as_expr()
    try:
        G = sp.groebner([clearden(e) for e in eqs] + [clearden(d0*t_ - 1)], *gens, modulus=PRIME, order='grevlex')
        exprs = list(G.exprs)
        if exprs == [sp.Integer(1)] or exprs == [1]:
            print(f"  GB = [1]  =>  EMPTY over GF({PRIME})-bar (no nonzero-det solution in box) [{time.time()-t0:.0f}s]", flush=True)
        else:
            print(f"  GB size {len(exprs)} != [1] => CONSISTENT (solutions exist mod {PRIME}) [{time.time()-t0:.0f}s]", flush=True)
    except Exception as e:
        print(f"  groebner failed: {e}", flush=True)

verdict(3, 4, 2, 3, 3, "k=3 CONTROL (must be CONSISTENT)")
verdict(4, 5, 2, 5, 4, "k=4 REFINED (kappa free, wide R)")
verdict(4, 5, 2, 4, 4, "k=4 TWO-TERM MIDDLE", two_term=True)
verdict(5, 6, 2, 5, 4, "k=5 REFINED (kappa free)")
print("\nDONE.", flush=True)
