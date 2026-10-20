#!/usr/bin/env python3
"""THM-4727 Lemma R (ramification exponent): for a polynomial map F: C^n -> C^n equivariant for weights (1,-w_2,...,-w_n) on
the source and (1,-w'_2,...,-w'_n) on the target (orbit coordinates x_1, X_1 of weight 1), with invariants v_i = x_1^{w_i} x_i,
V_i = X_1^{w'_i} X_i, the quotient map Fbar (V = Fbar(v)) and the orbit scaling l = X_1/x_1 (an invariant) satisfy
   det JF = x_1^{W - W'} * l^{1-W'} * det J Fbar,   W = sum w_i, W' = sum w'_i.
Checked symbolically on random equivariant maps built from random polynomials (l, Fbar) with the integrality conditions."""
import sympy as sp, random
random.seed(7)
def check(ws, wps, deg=2, trials=3):
    n = len(ws) + 1
    xs = sp.symbols(f'x1:{n+1}'); vs = sp.symbols(f'v2:{n+1}')
    x1 = xs[0]
    ok = True
    for _ in range(trials):
        # random polynomial l(v) and Fbar(v) of degree <= deg with small integer coefficients
        mons = [sp.Mul(*[v**e for v, e in zip(vs, es)]) for es in sp.utilities.iterables.iproduct(*[range(deg + 1)] * (n - 1)) if sum(es) <= deg]
        rnd = lambda: sum(random.randint(-2, 2) * m for m in mons) + random.randint(1, 3)
        l = rnd(); Fb = [rnd() for _ in range(n - 1)]
        # build F: X1 = x1 l(v), X_i = V_i / X1^{w'_i}; this is polynomial only if divisibility holds -- so we test the
        # identity as rational functions (it is an identity of rational functions regardless)
        sub = {vs[i]: x1**ws[i] * xs[i + 1] for i in range(n - 1)}
        X1 = x1 * l.subs(sub)
        F = [X1] + [sp.cancel(Fb[i].subs(sub) / X1**wps[i]) for i in range(n - 1)]
        JF = sp.Matrix(F).jacobian(list(xs)).det()
        JFb = sp.Matrix(Fb).jacobian(list(vs)).det().subs(sub)
        W, Wp = sum(ws), sum(wps)
        lhs = sp.cancel(JF); rhs = sp.cancel(x1**(W - Wp) * l.subs(sub)**(1 - Wp) * JFb)
        ok &= sp.simplify(lhs - rhs) == 0
    return ok
for ws, wps in [((1, 2), (1, 2)), ((1, 2), (2, 1)), ((1, 1), (1, 1)), ((2, 3), (2, 3)), ((1, 1, 1), (1, 1, 1)), ((1,), (1,)), ((2,), (2,)), ((1, 2), (1, 1))]:
    print(f"source weights (1,{','.join('-'+str(w) for w in ws)}), target weights (1,{','.join('-'+str(w) for w in wps)}): identity holds = {check(ws, wps)}")
