#!/usr/bin/env python3
"""Weight type (1,-1,-1): sweep-type ansatz A = p(w)/l + alpha, B = q(w)/l + beta, w = l*u, u = 1 + v, with p, q in (w) of
degree <= 3 and fixed linear l; Keller iff Omega = l{A,B} + A{l,B} + B{A,l} = c != 0 (THM-4727 R with W = 2).
For each l the system in the 6 + 2 unknowns is bilinear; solve with sympy (per-case timeout via alarm)."""
import sympy as sp, signal, itertools
v, s, c = sp.symbols('v s c')
def handler(signum, frame): raise TimeoutError
signal.signal(signal.SIGALRM, handler)
p1, p2, p3, q1, q2, q3, al, be = sp.symbols('p1 p2 p3 q1 q2 q3 alpha beta')
def br(f, g): return sp.diff(f, v) * sp.diff(g, s) - sp.diff(f, s) * sp.diff(g, v)
for lname, l in [("2 - 3v - s", 2 - 3*v - s), ("1 - v", 1 - v), ("1 - s", 1 - s), ("1 + v - s", 1 + v - s), ("v + s", v + s), ("s", s)]:
    u = 1 + v; w = l * u
    p = p1*w + p2*w**2 + p3*w**3; q = q1*w + q2*w**2 + q3*w**3
    A = sp.cancel(p / l) + al; B = sp.cancel(q / l) + be
    assert sp.denom(sp.together(A)) == 1 and sp.denom(sp.together(B)) == 1
    Om = sp.expand(l * br(A, B) + A * br(l, B) + B * br(A, l))
    eqs = sp.Poly(Om - c, v, s).coeffs()
    signal.alarm(150)
    try:
        sols = sp.solve(eqs, [p1, p2, p3, q1, q2, q3, al, be], dict=True)
        signal.alarm(0)
        good = [so for so in sols if sp.simplify(c.subs(so)) != 0 or c not in so]
        nz = []
        for so in sols:
            # c is determined? substitute and see whether Omega is a nonzero constant
            Omv = sp.expand(Om.subs(so))
            if Omv.free_symbols <= {c} | set([p1,p2,p3,q1,q2,q3,al,be]) and sp.Poly(Omv, v, s).total_degree() == 0 and sp.simplify(Omv) != 0:
                nz.append(so)
        print(f"l = {lname}: solutions {len(sols)}, with Omega a nonzero constant: {len(nz)}", flush=True)
        for so in nz[:4]:
            print("   ", {k: v_ for k, v_ in so.items()}, " A =", sp.factor(A.subs(so)), " B =", sp.factor(B.subs(so)), " Omega =", sp.simplify(Om.subs(so)))
    except TimeoutError:
        print(f"l = {lname}: solve timed out (150 s)", flush=True)
