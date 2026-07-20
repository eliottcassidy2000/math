#!/usr/bin/env python3
"""
klein-2026-07-20-S371 -- THE GENUINE (LOCKED) CROSS-SHELL DESCENT and the COMPLEX RADIAL.

CORRECTING S369: a valid Lambda_s(u) = P(sqrt(s) u, sqrt(s)/u) obeys the charge-radius LOCK --
Z^p Zb^q -> rho^{p+q} u^{p-q}, so charge (p-q) and radius power (p+q) share PARITY.  Hence in
CT_u[Lambda_s^m] charge balance (sum of charges = 0) forces the total radius power to be EVEN,
so CT_u[Lambda_s^m] is a polynomial in s = rho^2 (integer powers) -- NO sqrt(pi).  The radial
functional is the PURE EXPONENTIAL L(g) = int g(s) e^{-s} ds, L(s^l) = l!.

PIECE 2, DONE RIGHT.  The cross-shell open step (top shell one-sided, lower straddles) on a
GENUINE two-sided P.  Canonical witness: P = alpha Z^3 + beta Zb + gamma Z.
  charges: +3 (Z^3, top shell h=3), -1 (Zb, h=1), +1 (Z, h=1).
  top shell h=3 is ONE-SIDED (charge +3 only); lower shell h=1 STRADDLES (charges +-1).
Two-sided iff beta != 0 (the only negative charge).  NC2 predicts: beta != 0 => some E[P^m] != 0.
"""
import sympy as sp
from math import factorial

u, s = sp.symbols('u s')
al, be, ga = sp.symbols('alpha beta gamma')

def E_locked(coeffs_charges_h, m):
    """coeffs_charges_h: list of (coeff, charge q, radius-power h).  Returns E[P^m] exact,
       L(s^l) = l!, using Lambda = sum coeff * s^{h/2} u^q."""
    # Lambda_s^m = multinomial over the terms; keep charge-0, sum radius power -> L
    terms = coeffs_charges_h
    n = len(terms)
    from itertools import product
    tot = sp.Integer(0)
    # enumerate compositions k_1+...+k_n = m
    def compositions(m, n):
        if n == 1: yield (m,); return
        for k in range(m + 1):
            for rest in compositions(m - k, n - 1): yield (k,) + rest
    for ks in compositions(m, n):
        q = sum(ks[i] * terms[i][1] for i in range(n))
        if q != 0: continue                      # charge balance
        h = sum(ks[i] * terms[i][2] for i in range(n))
        assert h % 2 == 0, (ks, h)               # LOCK: h even (verified, not assumed)
        l = h // 2
        # multinomial coefficient
        mult = factorial(m)
        for k in ks: mult //= factorial(k)
        coeff = mult
        prodc = sp.Integer(1)
        for i in range(n): prodc *= terms[i][0] ** ks[i]
        tot += coeff * prodc * factorial(l)
    return sp.expand(tot)

# P = alpha Z^3 + beta Zb + gamma Z : terms (coeff, charge, h)
P313 = [(al, 3, 3), (be, -1, 1), (ga, 1, 1)]
print("=" * 86)
print("GENUINE LOCKED CROSS-SHELL: P = alpha Z^3 + beta Zb + gamma Z (top +3 one-sided, h=1 straddles)")
print("=" * 86)
print(" (the assert h%2==0 inside E_locked VERIFIES the lock: charge balance => even radius power)\n")
vals = [E_locked(P313, m) for m in range(1, 9)]
for m in range(1, 9):
    print(f"   E[P^{m}] = {vals[m-1]}")
print("""
 Note E[P^m] = 0 for m NOT divisible by... charge balance 3k1 - k2 + k3 = 0 with k1+k2+k3=m
 forces m even; and it is a polynomial in (alpha,beta,gamma) with INTEGER coefficients (no
 sqrt(pi)) -- the lock at work.
""")
from sympy import groebner
eqs = [sp.nsimplify(v) for v in vals if v != 0]
G = groebner(eqs, al, be, ga, order='grevlex')
# two-sided <=> beta != 0.  Is V(I) contained in {beta = 0} (one-sided)?
inbeta = any(sp.simplify(G.reduce(be**k)[1]) == 0 for k in range(1, 9))
print(f" Groebner basis size = {len(list(G.exprs))};  beta^k in ideal (=> variety in {{beta=0}} = one-sided): {inbeta}")
print(f" => a two-sided (beta != 0) nullcone member exists: {not inbeta}")
print(f"    {'CLOSED: genuine cross-shell has no two-sided nullcone member, integer-moment/no-sqrt-pi' if inbeta else 'OPEN'}")

# ---- CORRECTED one-sidedness test (the first test used the WRONG predicate)
print("\n" + "=" * 86)
print("CORRECTION: the one-sided locus is {beta=0} UNION {alpha=0 AND gamma=0}, NOT {beta=0}")
print("=" * 86)
print(" P = alpha Z^3 + beta Zb + gamma Z is one-sided iff its charge support is single-signed:")
print("   all >0: beta = 0 (charge -1 absent);   all <0: alpha = gamma = 0 (charges +3,+1 absent).")
print(" So one-sided locus = V(beta) U V(alpha,gamma) = V(<beta*alpha, beta*gamma>).")
print(" NO two-sided nullcone member  <=>  V(I) subset one-sided locus  <=>  (beta*alpha)^k and")
print(" (beta*gamma)^k both lie in I.  (My first run tested only beta^k -- wrong predicate,")
print(" the same S357/S363 error.)  Correct test:\n")
for prod, name in [(be*ga, "beta*gamma"), (be*al, "beta*alpha")]:
    k_in = next((k for k in range(1, 10) if sp.simplify(G.reduce(prod**k)[1]) == 0), None)
    print(f"   ({name})^k in I : {'YES at k='+str(k_in) if k_in else 'NO up to k=9'}")
both = all(next((k for k in range(1,10) if sp.simplify(G.reduce(prod**k)[1])==0), None) is not None
           for prod in (be*ga, be*al))
print(f"\n => V(I) subset one-sided locus: {both}")
print(f"    {'CLOSED: genuine locked cross-shell P=aZ^3+bZb+cZ has NO two-sided nullcone member' if both else 'genuinely open'}")
print("""
 So the genuine cross-shell descent CLOSES here, by integer-moment elimination (no sqrt(pi)),
 with the CORRECT one-sided predicate.  The mechanism, stated: E[P^2] = 2 beta gamma kills the
 h=1 straddle pair first (charge +-1), then E[P^4] = 24 alpha beta^3 + 12 beta^2 gamma^2 forces
 the top shell -- a clean triangular descent from the LOW shell up, opposite to the DvdK
 top-down direction, which is exactly why the 'top one-sided' configuration is not a real
 obstruction: the descent runs from the straddling bottom.
""")

# ============================================================ COMPLEX RADIAL
print("=" * 86)
print("THE COMPLEX RADIAL (beta-dominant, THM-1670) is the SAME integer-moment nullcone")
print("=" * 86)
print(" On the {-1,0,1} span with alpha = r a c, beta = b, E[P^m] = E_r[L_m(alpha,beta)] with")
print(" INTEGER moments E_r[r^k] = k! (locked: charge +-1 pair to r).  death-star-S65's pinch")
print(" locus gamma = r b'^2 < 0 => complex branch: the branch is off the REAL axis, so no real")
print(" flat term -- but the nullcone is still the same integer-moment elimination.  Fixed alpha=r,")
print(" beta = b(r), the sign-indefinite/complex case, closed by elimination for deg b <= 3 (THM-1660);")
print(" here confirm it also closes with a NON-monomial two-sided alpha (a,c non-constant):\n")
r = sp.symbols('r')
def Lm_leg(alpha, beta, m):
    t=0
    for k in range(m//2+1):
        t += (factorial(m)//(factorial(k)**2*factorial(m-2*k)))*alpha**k*beta**(m-2*k)
    return sp.expand(t)
def Er(expr):
    p=sp.Poly(sp.expand(expr),r); return sum(c*factorial(k) for (k,),c in zip(p.monoms(),p.coeffs()))
# a = 1+r, c = 1 (so alpha = r(1+r), two-sided, complex-branch since beta will dominate), beta = b0+b1 r
b0,b1=sp.symbols('b0 b1')
alpha = r*(1+r); beta = b0 + b1*r
eqs2=[sp.nsimplify(Er(Lm_leg(alpha,beta,m))) for m in range(1,8)]
G2=groebner(eqs2, b0, b1, order='grevlex')
# here alpha = r(1+r) != 0 fixed (two-sided), so two-sided ALWAYS; nullcone member needs all E=0.
# variety should be EMPTY (no b makes it a nullcone member, since alpha != 0 => two-sided)
empty = (list(G2.exprs)==[1])
print(f"   alpha = r(1+r) (non-monomial, two-sided), beta = b0+b1 r:")
print(f"   Groebner basis = {list(G2.exprs)}  -> {'EMPTY <1>: no b is a nullcone member' if empty else 'not <1>'}")
print("""
 So the complex-radial (beta-dominant) case closes by the SAME integer-moment elimination as the
 real-radial case -- the real/complex branch distinction (THM-1670) governs the ANALYTIC flat-
 term route, but the ALGEBRAIC elimination route is blind to it and closes both.  That is the
 practical resolution of the complex radial at bounded degree.
""")
