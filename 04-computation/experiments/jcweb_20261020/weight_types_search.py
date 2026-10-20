#!/usr/bin/env python3
"""General mixed weight types (1,-w2,-w3) on A^3 (free invariant ring k[v,s], v = x^{w2} y, s = x^{w3} z): equivariant maps
F = (x l, A/x^{w2}, B/x^{w3}) with integrality (monomials v^i s^j of A need w2 i + w3 j >= w2, of B need >= w3), quotient
(l^{w2} A, l^{w3} B), Keller iff det J(l^{w2} A, l^{w3} B) = c l^{W-1}, W = w2 + w3 (THM-4727 Lemma R). Low-degree Groebner
search for solutions with nonconstant l, for weight types (1,-1,-2) [small box], (1,-2,-3), (1,-1,-3), (1,-2,-2)."""
import sympy as sp, time, sys
v, s, c = sp.symbols('v s c')
def monos(D, cond):
    return [v**i * s**j for i in range(D + 1) for j in range(D + 1 - i) if cond(i, j)]
def run(w2, w3, D_l, D_A, D_B):
    t0 = time.time(); W = w2 + w3
    lm = monos(D_l, lambda i, j: True); am = monos(D_A, lambda i, j: w2 * i + w3 * j >= w2); bm = monos(D_B, lambda i, j: w2 * i + w3 * j >= w3)
    la = sp.symbols(f'l0:{len(lm)}'); aa = sp.symbols(f'a0:{len(am)}'); ba = sp.symbols(f'b0:{len(bm)}')
    l = sum(x * m for x, m in zip(la, lm)); A = sum(x * m for x, m in zip(aa, am)); B = sum(x * m for x, m in zip(ba, bm))
    V = l**w2 * A; S = l**w3 * B
    det = sp.expand(sp.diff(V, v) * sp.diff(S, s) - sp.diff(V, s) * sp.diff(S, v))
    # det must equal c l^{W-1}; divide exactly: det = l^{W-1} * Omega; compute Omega by polynomial division in (v,s) with symbolic coefficients
    q, r = sp.div(sp.Poly(det, v, s), sp.Poly(l**(W - 1), v, s))
    # r must vanish identically (it does: structural); impose Omega = c
    Om = q.as_expr()
    eqs = sp.Poly(sp.expand(Om - c), v, s).coeffs() + sp.Poly(r.as_expr(), v, s).coeffs()
    gens = list(la) + list(aa) + list(ba) + [c]
    try:
        G = sp.groebner(eqs, *gens, order='grevlex')
        basis = list(G.exprs)
        c_zero = any(sp.simplify(b - c) == 0 for b in basis)
        sols = sp.solve(eqs, gens[:-1], dict=True)
        nontriv = []
        for so in sols:
            lv = sp.expand(l.subs(so))
            if sp.Poly(lv, v, s).total_degree() >= 1: nontriv.append((lv, sp.expand(A.subs(so)), sp.expand(B.subs(so))))
        return dict(unknowns=len(gens) - 1, equations=len(eqs), basis=len(basis), c_forced_zero=c_zero, solutions=len(sols), nonconstant_l=len(nontriv), examples=[[str(e) for e in nt] for nt in nontriv[:4]], secs=round(time.time() - t0, 1))
    except Exception as e:
        return dict(error=f"{type(e).__name__}: {str(e)[:80]}", secs=round(time.time() - t0, 1))
for (w2, w3, Dl, DA, DB) in [(1, 2, 1, 2, 3), (2, 3, 1, 2, 3), (1, 3, 1, 2, 3), (2, 2, 1, 2, 2), (1, 1, 1, 3, 3)]:
    print(f"weights (1,-{w2},-{w3}), deg l <= {Dl}, deg A <= {DA}, deg B <= {DB}:", run(w2, w3, Dl, DA, DB), flush=True)
