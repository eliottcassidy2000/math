#!/usr/bin/env python3
"""
The two GMC(2) residual pieces: the radial Liouville step, and the cross-shell coupling
                                                        (mac-mini-S146)
================================================================================
Owner: "work the 2 GMC(2) residual pieces."

Reduction recap.  For 2 real Gaussians, write Z = r w, W = Zbar = r w^{-1} (w = e^{i theta},
s = r^2).  The Gaussian measure is (1/2pi) e^{-s} ds d theta.  Any P in C[Z,W] is a Laurent
polynomial in w with s-dependent coefficients:
    Lambda_s(w) = sum_k w^k s^{|k|/2} lambda_k(s),   lambda_k a polynomial in s,
and E[P^m] = L( CT_w(Lambda_s^m) ),  L(f) = int_0^inf f e^{-s} ds,  L(s^j) = j!.
The two residuals left after THM-1665:

PIECE 1 (HYP-8350, the radial Liouville step).  THM-1665 reduced the charge-0/radial layer
to  Psi(t) = int_0^inf e^{-v}[1/(1-tp(v)) - 1] dv == 0  =>  p == 0.  Here we CLOSE it for
REAL p by the discontinuity (jump) argument -- DvdK's Theorem-2 mechanism transplanted from
CT to L.

PIECE 2 (HYP-8470, cross-shell coupling).  CT_w(Lambda_s^m) = sum over charge-balanced
m-tuples (k_1..k_m), sum k_i = 0, of s^{(sum|k_i|)/2} prod lambda_{k_i}(s).  A PARITY LEMMA
(sum|k| == sum k == 0 mod 2) shows the s-power is ALWAYS an integer -- no half-shells in
E[P^m].  We establish the exact shell structure, verify span-2 decouples (THM-1600), and
locate where genuine cross-shell coupling first appears.
"""
from fractions import Fraction as F
from math import factorial, exp
import itertools, numpy as np

def L(coeffs):
    return sum(c*factorial(k) for k, c in enumerate(coeffs))
def ppow(a, m):
    r = [F(1)]
    for _ in range(m):
        out = [F(0)]*(len(r)+len(a)-1)
        for i, x in enumerate(r):
            for j, y in enumerate(a): out[i+j] += x*y
        r = out
    return r

# ================================================================= PIECE 1
print("=" * 78)
print("PIECE 1 -- the radial Liouville step: Psi == 0 => p == 0 for REAL p")
print("=" * 78)
print("  h(t) := int_0^inf e^{-v}/(1 - t p(v)) dv.  Psi == 0 (THM-1665) means h == 1.")
print("  For REAL p of degree D>=1 the boundary-value JUMP across real t>0 is")
print("     h(t+i0) - h(t-i0) = (2 pi i / t) * sum_{v: p(v)=1/t} e^{-v}/|p'(v)|,")
print("  each term POSITIVE.  h == 1 (analytic, single-valued -- the sector opening")
print("  pi(1+D) > pi wraps both boundary values) forces the jump = 0, hence NO v with")
print("  p(v) = 1/t and p'(v) != 0, for all large 1/t.  A nonconstant polynomial is")
print("  unbounded, contradiction.  So p is constant, and L(p)=p=0 => p==0.  QED (real p).")
print()
print("  Numerical witness of the positive jump density  D(x) = sum_{p(v)=x} e^{-v}/|p'(v)|:")
print(f"{'p(v)':>16} {'x=value':>9} {'roots v>=0 of p=x':>22} {'jump density D(x)':>18}")
def jump_density(coeffs, x):
    # real roots v>=0 of p(v) - x
    c = [float(a) for a in coeffs]; c[0] -= x
    roots = np.roots(list(reversed(c)))
    tot = 0.0; vs = []
    dp = np.polyder(list(reversed([float(a) for a in coeffs])))
    for rt in roots:
        if abs(rt.imag) < 1e-8 and rt.real >= -1e-9:
            v = rt.real; vs.append(round(v, 3))
            tot += exp(-v)/abs(np.polyval(dp, v))
    return vs, tot
for coeffs, nm in (([F(-1), F(1)], "v-1"), ([F(1), F(-3), F(1)], "v^2-3v+1"),
                   ([F(0), F(0), F(1)], "v^2")):
    for x in (2.0, 5.0):
        vs, d = jump_density(coeffs, x)
        print(f"{nm:>16} {x:>9} {str(vs):>22} {d:>18.6f}")
print("  D(x) > 0 wherever p hits x with p' != 0  =>  h != const  =>  no real nullcone p.")
print()
print("  Exhaustive check: NO real p (deg<=4, coeffs in [-3,3]) with L(p^m)=0 for m=1..12,")
print("  other than p == 0:")
found = 0; tested = 0
for D in (1, 2, 3, 4):
    for co in itertools.product(range(-3, 4), repeat=D+1):
        if co[-1] == 0: continue           # genuine degree D
        tested += 1
        if all(L(ppow([F(x) for x in co], m)) == 0 for m in range(1, 13)):
            found += 1
print(f"    tested {tested} nonzero real polynomials; L-nullcone members found: {found}")

# ================================================================= PIECE 2
print()
print("=" * 78)
print("PIECE 2 -- cross-shell coupling: the parity lemma and the shell structure")
print("=" * 78)
print("  CT_w(Lambda_s^m) = sum over balanced m-tuples of s^{(sum|k_i|)/2} prod lambda_{k_i}.")
print("  PARITY LEMMA: sum|k_i| == sum k_i (mod 2) = 0, so the s-power is ALWAYS an integer.")
print("  => E[P^m] never sees a half-integer moment; the 'half-shell' worry is a red herring.")
print()
# verify the parity lemma by brute enumeration of balanced tuples
def check_parity(charges, m):
    bad = 0
    for tup in itertools.product(charges, repeat=m):
        if sum(tup) == 0 and sum(abs(k) for k in tup) % 2 != 0: bad += 1
    return bad
for charges in ([-1, 0, 1], [-2, -1, 0, 1, 2], [-3, 0, 3]):
    bad = sum(check_parity(charges, m) for m in range(1, 6))
    print(f"  charges {str(charges):>16}: balanced tuples with odd sum|k| (m<=5): {bad}"
          f"  {'(parity lemma holds)' if bad == 0 else '*** FAILS ***'}")

print()
print("  THE SHELL DECOMPOSITION.  E[P^m] = sum_j L( s^j * c_j(m) ) where j = (sum|k|)/2 is")
print("  the CHARGE-DEGREE and c_j(m) collects prod lambda_{k_i} over balanced tuples with")
print("  that charge-degree.  For SPAN 2 (charges {-1,0,1}) only j-pairs occur:")
print("      E[P^m] = sum_j  m!/(j!^2 (m-2j)!)  L( s^j (lam1 lam_{-1})^j lam0^{m-2j} ).")
print("  With CONSTANT lambdas this DECOUPLES (THM-1600: m=1 => lam0=0, then rs=0). Check:")
def EPm_span2(l1, l0, lm1, m):
    tot = F(0); rs = [x for x in (l1[i]*lm1[j] for i in range(len(l1)) for j in range(len(lm1)))]
    # do it properly with polynomial multiply in s
    def pm(a, b):
        o = [F(0)]*(len(a)+len(b)-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): o[i+j] += x*y
        return o
    def pw(a, k):
        r = [F(1)]
        for _ in range(k): r = pm(r, a)
        return r
    prod = pm(l1, lm1)
    for j in range(m//2+1):
        mult = F(factorial(m), factorial(j)**2 * factorial(m-2*j))
        term = pm([F(0)]*j + [F(1)], pm(pw(prod, j), pw(l0, m-2*j)))
        tot += mult * L(term)
    return tot
print(f"{'(lam1,lam0,lam-1)':>20} {'E[P^m], m=1..6':>34} {'nullcone?':>10}")
for (a, b, c) in (([1], [0], [1]), ([1], [1], [1]), ([1], [0], [0]), ([2], [-1], [3])):
    vals = [EPm_span2([F(x) for x in a], [F(x) for x in b], [F(x) for x in c], m)
            for m in range(1, 7)]
    nc = all(v == 0 for v in vals)
    print(f"{str((a,b,c)):>20} {str(vals):>34} {str(nc):>10}")

print()
print("  WHERE CROSS-SHELL COUPLING FIRST BITES.  It is NOT parity (always integer shells).")
print("  It is that at each m, E[P^m] = sum_j L(s^j c_j) is a SINGLE number mixing all")
print("  charge-degrees j -- no per-shell separation at fixed m.  Separation must come from")
print("  varying m.  The TOP shell (all charges +-k_max) has the fastest-growing moment")
print("  (m k_max/2)!, so it dominates for large m: this is the shell-descent handle.")
print(f"{'span':>6} {'charges':>16} {'top-shell moment growth in m':>30}")
for span, charges in ((2, [-1, 0, 1]), (3, [-2, -1, 0, 1, 2]), (4, [-3, 0, 3])):
    kmax = max(charges)
    print(f"{span:>6} {str(charges):>16} {'(m*kmax/2)! = ('+str(kmax)+'m/2)!':>30}")
print()
print("SUMMARY")
print("  PIECE 1 (HYP-8350): CLOSED for real p by the jump argument -- no nonconstant real p")
print("  is in the L-nullcone, so the charge-0/radial layer of GMC(2) closes for Hermitian P.")
print("  PIECE 2 (HYP-8470): the shells are ALL INTEGER (parity lemma), span-2 decouples")
print("  (THM-1600), and the genuine coupling is the fixed-m mixing sum_j L(s^j c_j), with")
print("  the top shell (m kmax/2)! the descent handle. Framed, not closed.")
