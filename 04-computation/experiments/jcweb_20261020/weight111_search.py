#!/usr/bin/env python3
"""HYP-9343 (c)(3): the weight type (1,-1,-1) on A^3 (W = 2, exponent 1). Invariants v = xy, s = xz; an equivariant map is
F = (x l(v,s), A(v,s)/x, B(v,s)/x) with A, B in (v,s) (integrality); quotient map (lA, lB); Keller iff
   l {A,B} + A {l,B} + B {A,l} = c  (c != 0).
Search: l of degree <= 1, A and B of degree <= D (D = 2, then 3) with unknown coefficients; Groebner basis of the
coefficient ideal of (equation - c) with c a free symbol; report whether solutions with nonconstant l or nonlinear (A,B) exist,
and whether any is non-injective (a collision over {l = 0}). Also the same for weights (1,-1,-2) in the box deg l <= 1,
deg A <= 2, deg B <= 3 as a control (the Alpoge map needs deg A = 3, deg B = 4, so the box should contain only trivial solutions)."""
import sympy as sp, itertools, sys, time
v, s, c = sp.symbols('v s c')
def monos(D, lowest=0):
    return [v**i * s**j for i in range(D + 1) for j in range(D + 1 - i) if i + j >= lowest]
def run(weights, D_A, D_B, B_ideal_lowest):
    t0 = time.time()
    lm = monos(1); am = monos(D_A, 1); bm = monos(D_B, B_ideal_lowest)
    if weights == (1, -1, -2):
        # B in (v^2, t): drop the monomial v (degree-1 in v only); s plays the role of t
        bm = [m for m in bm if m != v]
    la = sp.symbols(f'l0:{len(lm)}'); aa = sp.symbols(f'a0:{len(am)}'); ba = sp.symbols(f'b0:{len(bm)}')
    l = sum(x * m for x, m in zip(la, lm)); A = sum(x * m for x, m in zip(aa, am)); B = sum(x * m for x, m in zip(ba, bm))
    def br(f, g): return sp.diff(f, v) * sp.diff(g, s) - sp.diff(f, s) * sp.diff(g, v)
    if weights == (1, -1, -1):
        Om = sp.expand(l * br(A, B) + A * br(l, B) + B * br(A, l))
    else:
        Om = sp.expand(l * br(A, B) + A * br(l, B) + 2 * B * br(A, l))
    eqs = sp.Poly(Om - c, v, s).coeffs()
    gens = list(la) + list(aa) + list(ba) + [c]
    G = sp.groebner(eqs, *gens, order='grevlex')
    basis = list(G.exprs)
    # is c forced to zero?  is l forced constant?  is (A,B) forced linear?
    c_forced = any(sp.simplify(b - c) == 0 for b in basis) or (c in basis)
    info = {"basis size": len(basis), "c forced 0": c_forced, "time": round(time.time() - t0, 1)}
    # sample solutions: solve for the variety dimension via sympy solve on a reduced system if small
    try:
        sols = sp.solve(eqs, gens[:-1], dict=True)
        info["solutions"] = len(sols)
        nontriv = []
        for so in sols:
            lv = sp.expand(l.subs(so)); Av = sp.expand(A.subs(so)); Bv = sp.expand(B.subs(so))
            if sp.Poly(lv, v, s).total_degree() >= 1 or sp.Poly(Av, v, s).total_degree() >= 2 or sp.Poly(Bv, v, s).total_degree() >= 2:
                nontriv.append((lv, Av, Bv))
        info["nontrivial (nonconstant l or nonlinear A/B)"] = len(nontriv)
        for nt in nontriv[:6]: info.setdefault("examples", []).append([str(e) for e in nt])
    except Exception as e:
        info["solve"] = f"skipped ({type(e).__name__})"
    return info
print("weights (1,-1,-1), deg l <= 1, deg A, B <= 2:", run((1, -1, -1), 2, 2, 1), flush=True)
print("weights (1,-1,-2) control, deg l <= 1, deg A <= 2, deg B <= 3:", run((1, -1, -2), 2, 3, 2), flush=True)
print("weights (1,-1,-1), deg l <= 1, deg A, B <= 3:", run((1, -1, -1), 3, 3, 1), flush=True)
