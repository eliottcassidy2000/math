#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S64
=======================
B2 / E2 / DEDEKIND objects on the genuinely hard residual of LRC14.

CENTERPIECE (verified n=3..15): the LRC14 covering-min MARGIN is the observer's DEDEKIND SUM.
  margin(n) := n/Phi6(n) - 1/n  =  (n-1)/(n*Phi6)  =  -12 * s(n, Phi6(n)) / n^2
where Phi6(n)=n^2-n+1 and s(h,k) is the Dedekind sum.  The closed form
  12 * s(n, Phi6) = -(Phi6-1)/Phi6            [because n^3 = -1 (mod Phi6): n is the order-6 element]
is the substance; it makes the margin (an LRC quantity) equal to a modular/Dedekind quantity.

B2 vs A2/E2 DICHOTOMY:
  - order-4 element h (h^2 = -1 mod k) = the SQUARE lattice / B2 / imaginary unit i (90 deg): s(h,k)=0.
  - order-6 element n (n^2-n+1 = 0 mod Phi6, n^3=-1) = the HEXAGONAL lattice / A2 / Eisenstein omega
    (60 deg): s(n,Phi6) = -(Phi6-1)/(12 Phi6) != 0.
  The LRC14 margin is POSITIVE exactly because the observer lives on the HEXAGONAL side (nonzero
  Dedekind anomaly).  A square/B2 observer would give ZERO margin (covering-min = floor, degenerate).

E2 (weight-2 Eisenstein):  F_d(tau) = E2(tau) - d*E2(d*tau) is a genuine weight-2 modular form on
  Gamma_0(d); for d=14, F_14 has constant term 1-14 = -(n-1).  The three forms F_2, F_7, F_14 span the
  3-dim Eisenstein subspace of M_2(Gamma_0(14)) (4 cusps => dim c-1 = 3); the 1-dim cusp form f_14
  (genus(X_0(14))=1) is the OBSTRUCTION = the hard residual (control f_14 at the apex cusp d=7).

RESIDUAL ATTACK:  the margin = a Dedekind sum; eroding the margin to 0 (a patch beating the covering-min)
  would require the observer Dedekind anomaly to vanish -- i.e. restore the order-4/square structure.
  A single patch integer has CRT-linked residues; by Dedekind-Rademacher reciprocity (the repo's
  far-coherence tool, HYP-2808) its Dedekind sums at different moduli cannot all be tuned to the
  anomaly-free value at once.  So the hexagonal anomaly (=the margin) persists -- the mechanism behind
  "the hole moves but never vanishes."
"""
from fractions import Fraction as F
from math import gcd

def dsum(h, k):
    """Exact Dedekind sum s(h,k) = sum_{i=1}^{k-1} ((i/k))((hi/k))."""
    def saw(x):
        x = x - int(x)
        return F(0) if x == 0 else x - F(1, 2)
    return sum((saw(F(i, k))*saw(F((h*i) % k, k)) for i in range(1, k)), F(0))

def sigma1(k): return sum(d for d in range(1, k+1) if k % d == 0)

print("="*78)
print("(1) CLOSED FORM  12 s(n,Phi6) = -(Phi6-1)/Phi6   (from n^3 = -1 mod Phi6, order-6)")
print("="*78)
for n in range(3, 16):
    P = n*n-n+1
    if gcd(n, P) != 1: continue
    s = dsum(n % P, P)
    print(f"  n={n:2d} Phi6={P:3d}: 12 s(n,Phi6)={str(12*s):8} vs -(Phi6-1)/Phi6={str(F(-(P-1),P)):8}"
          f"  match={12*s==F(-(P-1),P)}  n^3 mod Phi6={pow(n,3,P)} (=-1)")

print()
print("="*78)
print("(2) THE MARGIN IDENTITY  margin = n/Phi6 - 1/n = -12 s(n,Phi6)/n^2")
print("="*78)
for n in range(3, 16):
    P = n*n-n+1
    if gcd(n, P) != 1: continue
    margin = F(n, P)-F(1, n)
    ded = -12*dsum(n % P, P)/(n*n)
    print(f"  n={n:2d}: margin={str(margin):9}={float(margin):.6f}  -12 s/n^2={str(ded):9}  match={margin==ded}")
print("  => LRC14 (n=14) margin 13/2562 IS the observer Dedekind sum s(14,183)=-91/1098.")

print()
print("="*78)
print("(3) B2 (square, order-4, s=0) vs A2/E2 (hexagonal, order-6, s!=0): the margin's source")
print("="*78)
print("  SQUARE/B2: order-4 h (h^2=-1 mod k) -> s(h,k)=0 -> would give ZERO margin (degenerate):")
for k in [5, 13, 17, 25, 29]:
    hs = [h for h in range(2, k) if (h*h) % k == k-1]
    if hs: print(f"    k={k:2d}: h={hs[0]} (h^2=-1): s={dsum(hs[0],k)}  [B2, margin 0]")
print("  HEXAGONAL/A2: order-6 n (n^3=-1 mod Phi6) -> s=-(Phi6-1)/12Phi6 != 0 -> POSITIVE margin:")
for n in [4, 7, 8, 14]:
    P = n*n-n+1
    print(f"    Phi6({n})={P:3d}: s(n,Phi6)={dsum(n%P,P)}  [A2, margin={F(n,P)-F(1,n)}]")

print()
print("="*78)
print("(4) E2 weight-2 Eisenstein forms on Gamma_0(14): F_d = E2(tau)-d E2(d tau), d|14")
print("="*78)
for d in [2, 7, 14]:
    coeffs = [1-d] + [-24*(sigma1(k) - (d*sigma1(k//d) if k % d == 0 else 0)) for k in range(1, 9)]
    print(f"  F_{d:2d} = E2-{d}E2({d}t): const={1-d:4d}  q-exp[0..8]={coeffs}")
print("  const(F_14) = 1-14 = -(n-1) = -13.  F_2,F_7,F_14 span the 3-dim Eisenstein subspace")
print("  (X_0(14) has 4 cusps => dim = c-1 = 3); the 1-dim cusp form f_14 (genus 1) = the residual.")
print("  The apex cusp d=7 (Atkin-Lehner W_7, width 2) is where f_14 must be controlled = the hard core.")
print("\nDONE.")
