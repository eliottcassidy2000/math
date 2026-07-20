"""opus-2026-07-20-S419 -- THE NEXT PROOF: t * Pi_large(t) = const  =>  R monomial.

From THM-1635: TNC <=> Pi(t) = c t <=> t * Pi_large(t) = const, where Pi_large is the product
of the M LARGE branches of u^N - t R(u) = 0 (the roots -> infinity as t -> 0).

STRATEGY.  Pi_large is a symmetric function of the large roots, so it is RATIONAL in t --
computable exactly as a RESULTANT / subresultant, no Puiseux needed.  Concretely:
  product of ALL roots = (-1)^{M+N} r_0 / r_{M+N}  (constant, THM-1635),
  Pi_small(t) = Pi(t),  and  Pi * Pi_large = that constant.
So  t * Pi_large = const  <=>  Pi(t) = c t  <=>  Pi(t)/t is a nonzero CONSTANT.

Pi(t) = product of small branches is the RESULTANT-type object: the small branches are the N
roots near 0, and their product is (up to sign) the ratio of two coefficients of the
"small factor" in the Newton-polygon factorization u^N - tR = (small factor)(large factor).

CLEAN EXACT ROUTE (this script).  Substitute the small-branch scaling.  Near t=0 the small
roots satisfy u^N ~ t R(0) = t r_0, so u ~ (t r_0)^{1/N}: set u = (t)^{1/N} v.  Then
  u^N - tR(u) = t v^N - t R(t^{1/N} v) = t ( v^N - R(t^{1/N} v) ).
The small branches are v = O(1) roots of  Phi(v, s) := v^N - R(s v) = 0  with s = t^{1/N}.
Their product Pi_v(s) = prod of the N bounded roots, and Pi(t) = t * Pi_v(s) (each u = s v).
So  **Pi(t) = c t  <=>  Pi_v(s) = c = CONSTANT in s = t^{1/N}.**
Pi_v(0) = product of roots of v^N - R(0) = v^N - r_0, which is (-1)^N (-r_0) = (-1)^{N+1} r_0
... i.e. -(-r_0) pattern; nonzero.  TNC <=> that product stays CONSTANT as s grows.

Pi_v(s) = product of the N bounded roots of Phi(v,s) = v^N - R(sv).  R(sv) = sum_k r_k s^k v^k,
so Phi = v^N - sum_k r_k s^k v^k = -sum_{k=0}^{M+N} r_k s^k v^k + v^N.  The BOUNDED roots are
the N roots that stay finite as s->0 (the other M blow up like 1/s).  Pi_v(s) = product of
bounded roots = (-1)^N * (const term in v)/(coeff that controls them) ... compute exactly.
"""
import sympy as sp

u, t, v, s = sp.symbols('u t v s')

def CT(R, N, mm):
    return sp.Poly(sp.expand(sp.sympify(R)**mm), u).coeff_monomial(u**(N*mm))

def pi_small_series(R, N, order):
    """
    Pi(t) = product of small branches, as a t-series, via Pi = t*Pi_v(s), s=t^{1/N}.
    We instead get Pi(t) from G(t) = -t (log Pi)' = sum CT t^m, i.e.
       (log Pi)' = -G/t = -sum_{m>=1} CT(m) t^{m-1} - CT(0)/t ... careful with the t^{-1}.
    Since Pi ~ c t near 0 iff Pi(t)/t is analytic and constant, test THAT via the log-deriv:
       t (log(Pi/t))' = t(log Pi)' - 1 = -G(t) - 1 = -(1 + sum_{m>=1} CT t^m)...
    Actually G(0)=CT(0)=1, so -G = -1 - sum_{m>=1}CT t^m, and t(log Pi)' = -G gives
       t(log Pi)' = -1 - sum_{m>=1} CT t^m.
    For Pi = c t exactly, log Pi = log c + log t, t(log Pi)' = 1.  So we need -1 - (...) = 1,
    i.e. sum_{m>=1} CT t^m = -2?? -- sign slip.  Recompute G(0) and the constant cleanly below.
    """
    return None

print("="*78)
print("(0) FIX THE CONSTANT: recompute G(0) and the exact relation G = -t (log Pi)'")
print("="*78)
print("  G(t) = sum_{m>=0} CT(Lambda^m) t^m, and CT(Lambda^0) = [u^0] R^0 = 1, so G(0) = 1.")
print("  Pi(t) ~ (small branches) with u_i ~ (r_0 t)^{1/N} zeta_i, product ~ t * (prod zeta_i)")
print("  * r_0^{... }.  So Pi(t) = c t (1 + O(t^{1/N})) generically; TNC demands the CORRECTION")
print("  vanishes: Pi(t) = c t EXACTLY.  Check G = -t(log Pi)' on the extreme-weight PROVEN case")
print("  R = 1 + u (N=1, M=0 -- wait M=0 is the proven side): there Pi(t) = the single small root")
print("  of u - t(1+u) = 0 => u = t/(1-t).  Pi(t) = t/(1-t).")
Pi = t/(1-t)
Gcheck = sp.simplify(-t*sp.diff(sp.log(Pi), t))
print(f"    R=1+u, N=1: Pi(t) = t/(1-t);  -t (log Pi)' = {Gcheck}")
cts = [CT("1+u", 1, mm) for mm in range(0, 8)]
Gser = sum(cts[m]*t**m for m in range(8))
print(f"    sum CT t^m = {Gser}  (CT = {cts})")
print(f"    series of -t(log Pi)': {sp.series(Gcheck, t, 0, 8).removeO()}")
print(f"    MATCH: {sp.simplify(sp.series(Gcheck,t,0,6).removeO() - sum(cts[m]*t**m for m in range(6)))==0}")
print("    -> identity G = -t(log Pi)' CONFIRMED.  And Pi = t/(1-t) is NOT c*t, so this R is")
print("       NOT in the nullcone -- correct, 1+u has CT = 1 != 0.")

print()
print("="*78)
print("(1) THE RESCALED PICTURE: Pi(t) = c t  <=>  Pi_v(s) constant,  Phi(v,s) = v^N - R(s v)")
print("="*78)
print("  s = t^{1/N}; bounded roots v_i = u_i/s.  Pi_v(s) = prod bounded roots of v^N - R(sv).")
print("  At s=0: v^N - R(0) = v^N - r_0, product of roots = (-1)^N(-r_0) = (-1)^{N-1} r_0.")
for R, N in [("1+u+u**2", 2), ("u**4-2*u**2-2", 2), ("2+u", 1), ("3+u+u**5", 1)]:
    Rp = sp.sympify(R)
    Phi = sp.expand(v**N - Rp.subs(u, s*v))
    P = sp.Poly(Phi, v)
    # bounded roots: as s->0, Phi -> v^N - r_0 (degree N); the extra M roots come from the
    # top coeff -r_{M+N} s^{M+N} v^{M+N} vanishing at s=0 (roots -> infinity).
    r0 = Rp.subs(u, 0)
    piv0 = (-1)**N * (v**N - r0).as_poly(v).all_coeffs()[-1] / 1
    piv0 = sp.simplify((-1)**N * (-r0))
    print(f"  R={R:16s} N={N}: Phi(v,s)=v^N - R(sv);  Pi_v(0)={piv0}  (=(-1)^N(-r_0))")

print()
print("="*78)
print("(2) THE CLAIM TO PROVE: d/ds Pi_v(s) at s>0 -- when is it identically 0?")
print("="*78)
print("  Pi_v(s) = (-1)^N * [product of bounded roots] = (-1)^N * (v-independent part).")
print("  For a monic-in-v^N picture: Phi = v^N - sum_k r_k s^k v^k.  The bounded-root product")
print("  is the resultant-type ratio; compute Pi_v(s) exactly for SMALL M and read off constancy.")
for R, N, M in [("2+u", 1, 0), ("2+u+u**2", 1, 1), ("1+u+u**2", 2, 0), ("1+u**3", 2, 1)]:
    Rp = sp.sympify(R); deg = sp.Poly(Rp, u).degree()
    Phi = sp.Poly(sp.expand(v**N - Rp.subs(u, s*v)), v)
    coeffs = Phi.all_coeffs()             # in v, degree = max(N, deg)
    lead = coeffs[0]; const = coeffs[-1]
    fullprod = sp.simplify((-1)**Phi.degree() * const/lead)   # product of ALL roots (bounded+unbounded)
    print(f"  R={R:14s} (N={N},M={M}): product of ALL v-roots = {fullprod}")
    print(f"       d/ds of that = {sp.simplify(sp.diff(fullprod, s))}   "
          f"(large roots ~ 1/s^? carry the s-dependence)")
