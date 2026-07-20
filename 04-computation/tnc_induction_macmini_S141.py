#!/usr/bin/env python3
"""
An induction for the Toral Nullcone Conjecture, on the coefficient ladder
                                                        (mac-mini-S141)
================================================================================
Owner: "work induction to prove GMC(2) and TNC."

THE RESTATEMENT THAT MAKES INDUCTION POSSIBLE.  For Lambda = u^{-M} R(u),
    CT(Lambda^m) = CT( u^{-Mm} R(u)^m ) = [u^{Mm}] R(u)^m.
So the whole infinite family is a COEFFICIENT condition on powers of one polynomial:

    TNC(M,N):   [u^{Mm}] R(u)^m = 0  for every m >= 1,   deg R = d = M+N, R(0) != 0
                ==>  R is degenerate.

Normalise R(0) = 1 (scaling R rescales t, cf. THM-1550).  Write R = 1 + r_1 u + ... + r_d u^d.

THE LADDER.  m = 1 gives [u^M]R = r_M = 0 outright.  The question this file asks is whether
the higher m keep peeling coefficients in a TRIANGULAR way -- each new m killing one new
r_j given the previous ones -- which is exactly what an induction needs.
"""
from fractions import Fraction as F
import itertools

def polymul(a, b, cap):
    out = [0]*(min(len(a)+len(b)-1, cap+1))
    for i, x in enumerate(a):
        if not x: continue
        for j, y in enumerate(b):
            if i+j > cap: break
            out[i+j] += x*y
    return out

def coeff_MM(rs, M, m):
    """[u^{Mm}] R(u)^m  for R = 1 + sum r_i u^i (rs is [r_1,...,r_d])."""
    cap = M*m
    R = [1] + list(rs)
    P = [1]
    for _ in range(m):
        P = polymul(P, R, cap)
    return P[cap] if cap < len(P) else 0

print("=" * 78)
print("PART A -- the restatement, and the m = 1 seed")
print("=" * 78)
print("  CT(Lambda^m) = [u^{Mm}] R(u)^m.   m = 1 gives  [u^M] R = r_M = 0.")
print("  So the FIRST coefficient is killed for free.  Checking the ladder numerically.")
print()

try:
    import sympy as sp
    HAVE = True
except ImportError:
    HAVE = False
    print("  sympy unavailable -- symbolic part skipped")

if HAVE:
    for (M, N) in ((1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2)):
        d = M + N
        r = sp.symbols(f'r1:{d+1}')
        eqs = []
        MM = d + 3
        for m in range(1, MM+1):
            eqs.append(sp.expand(coeff_MM(list(r), M, m)))
        # peel: at each stage, which single r_j does the next equation isolate?
        known = {}
        ladder = []
        for m, e in enumerate(eqs, start=1):
            e2 = sp.expand(e.subs(known))
            if e2 == 0:
                ladder.append((m, "vacuous")); continue
            free = sorted(e2.free_symbols, key=lambda s: int(str(s)[1:]))
            # is it LINEAR in its highest-index unknown?
            top = free[-1]
            deg = sp.Poly(e2, top).degree()
            if deg == 1:
                sol = sp.solve(e2, top)[0]
                known[top] = sp.simplify(sol)
                ladder.append((m, f"solves {top} (linear)"))
            else:
                ladder.append((m, f"deg {deg} in {top} -- LADDER BREAKS"))
        solved = [str(k) for k in known]
        broke = [x for x in ladder if "BREAKS" in x[1]]
        print(f"  TNC({M},{N})  d={d}:  peeled {len(known)}/{d} coefficients "
              f"{solved}")
        print(f"      ladder: {[f'm={m}:{w}' for m, w in ladder[:6]]}")
        if not broke and len(known) == d:
            print(f"      => FULL TRIANGULAR PEEL: every r_j forced. "
                  f"Residual: {[sp.simplify(v) for v in known.values()]}")
        elif broke:
            print(f"      => breaks first at {broke[0][0]}")
        print()

print("=" * 78)
print("PART B -- what the peel actually forces: brute search over small integer R")
print("=" * 78)
print(f"{'M':>3} {'N':>3} {'d':>3} {'R with [u^Mm]R^m=0 for m=1..8, coeffs in [-3,3]':>52}")
for (M, N) in ((1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2)):
    d = M + N
    hits = []
    for rs in itertools.product(range(-3, 4), repeat=d):
        if all(coeff_MM(list(rs), M, m) == 0 for m in range(1, 9)):
            hits.append(rs)
    shown = [h for h in hits if any(h)][:4]
    print(f"{M:>3} {N:>3} {d:>3} {str(shown) if shown else 'ONLY R = 1 (all r_i = 0)':>52}")
print("  A nonzero hit here would be a TNC COUNTEREXAMPLE at that bidegree.")

print()
print("=" * 78)
print("PART C -- the Laplace layer (GMC(2) side): a Borel/Watson reduction")
print("=" * 78)
print("  HYP-8350 asks: L(p^m) = 0 for all m >= 1  =>  p = 0,  L(f) = int_0^inf f e^{-v} dv.")
print("  Exponential generating function:")
print("      F(t) := sum_m L(p^m) t^m / m!  'is'  int_0^inf e^{t p(v)} e^{-v} dv.")
print("  If all L(p^m) = 0 for m >= 1 then F is asymptotically the CONSTANT 1 at t = 0.")
print()
print("  CONVERGENCE.  On the ray v = rho e^{i phi} (|phi| < pi/2, so Re v > 0 keeps e^{-v}")
print("  decaying, and L(v^k) = k! is unchanged by contour rotation since the integrand is")
print("  entire), the integral converges when Re( t a_D e^{i D phi} ) < 0, an OPEN HALF-PLANE")
print("  in t.  Sweeping phi over (-pi/2, pi/2) rotates that half-plane through D*pi, so F")
print("  continues analytically to a sector of opening")
print("      pi + D*pi = pi(1 + D).")
print()
print(f"{'D':>3} {'sector opening':>16} {'> pi (Watson threshold)?':>26}")
import math
for D in range(1, 7):
    op = math.pi*(1+D)
    print(f"{D:>3} {f'{op:.4f} = pi*{1+D}':>16} {str(op > math.pi):>26}")
print()
print("  So for EVERY D >= 1 the sector exceeds Watson's threshold pi, which is the")
print("  hypothesis of the Watson-Nevanlinna theorem: a function analytic in a sector of")
print("  opening > pi with Gevrey-1 asymptotics is BOREL-DETERMINED by its asymptotic")
print("  series.  Here that series is identically 1, so F == 1, so the pushforward of")
print("  e^{-v}dv under p is delta_0, so p == 0.")
print()
print("  WHAT IS AND IS NOT ESTABLISHED.  The sector computation is elementary and is done")
print("  above.  The Gevrey-1 bound |L(p^m)| <= C A^m (m!)^{1} needed to APPLY")
print("  Watson-Nevanlinna is NOT verified here -- and it is exactly the estimate HYP-8350")
print("  has been missing all along, now in a named classical form rather than an ad hoc")
print("  one.  Checking the growth numerically:")
print(f"{'p':>14} {'m':>3} {'|L(p^m)|':>16} {'|L(p^m)|/(Dm)!':>18}")
from math import factorial
def L(c): return sum(a*factorial(k) for k, a in enumerate(c))
def ppow(a, m):
    r = [F(1)]
    for _ in range(m):
        out = [F(0)]*(len(r)+len(a)-1)
        for i, x in enumerate(r):
            for j, y in enumerate(a): out[i+j] += x*y
        r = out
    return r
for co, nm in (([F(-1), F(1)], "v-1"), ([F(1), F(0), F(1)], "v^2+1"),
               ([F(0), F(2), F(0), F(1)], "v^3+2v")):
    D = len(co)-1
    for m in (2, 5, 8):
        val = L(ppow(co, m))
        print(f"{nm:>14} {m:>3} {abs(val):>16} {float(abs(val))/factorial(D*m):>18.6f}")
print("  |L(p^m)| stays bounded by (Dm)! times an O(1) factor -- consistent with Gevrey-D,")
print("  i.e. Gevrey-1 in the variable t^{1/D}.  That is the right normalisation and it")
print("  matches the sector opening pi(1+D) being measured in the SAME variable.")
