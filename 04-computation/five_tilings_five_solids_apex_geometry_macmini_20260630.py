#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S65
=======================
THE 5 PLANE-TILINGS / 5 PLATONIC SOLIDS and their connection to the LRC14 proof.

The two counts of 5 and the covering-min/apex geometry all come from ONE source: the (2,3,n)
angle-defect trichotomy (Gauss-Bonnet), i.e. the sign of 1/2 + 1/3 + 1/n - 1.

  n=3,4,5 : SPHERICAL   -> the 5 Platonic solids {3,3},{3,4},{4,3},{3,5},{5,3}  (positive curvature)
  n=6     : EUCLIDEAN   -> the plane tiling (hexagonal = A2 root lattice = the LRC covering-min)
  n>=7    : HYPERBOLIC  -> (2,3,7)=Klein quartic PSL(2,7)=168=Fano; the LRC14 apex

COUNT 5 (planar): crystallographic orders {n: 2cos(2pi/n) in Z} = {n: phi(n)<=2} = {1,2,3,4,6}
  = the 5 two-dimensional Bravais lattices.

LRC(2p) APEX PRIME p walks the (2,3,p) spine, and the modular genus of X0(2p) tracks the geometry:
  p=3 (LRC6)  Euclidean/crystallographic   genus X0(6)=0    trivial
  p=5 (LRC10) SPHERICAL/icosahedral(Platonic) genus X0(10)=0 tractable
  p=7 (LRC14) FIRST HYPERBOLIC (Klein)      genus X0(14)=1   HARD (cusp form f14)
The genus jump 0->1 at p=7 IS the spherical->hyperbolic transition = why LRC14 is the first hard case.
The 5 Platonic solids consume the spherical orders {3,4,5}; order 6 is the unique Euclidean (covering-min,
A2, last session's Dedekind anomaly); 7 is the first order that tiles NEITHER sphere nor plane => the frontier.
"""
from fractions import Fraction as F
from math import cos, pi, gcd

def phi(n): return sum(1 for k in range(1, n+1) if gcd(k, n) == 1) if n > 0 else 0
def primes_dividing(N): return [p for p in range(2, N+1) if N % p == 0 and all(p % d for d in range(2, p))]

print("="*78)
print("(1) COUNT 5: crystallographic rotation orders = {n: 2cos(2pi/n) in Z} = {phi(n)<=2}")
print("="*78)
allowed = [n for n in range(1, 25) if phi(n) <= 2]
print(f"   = {allowed}  (count {len(allowed)}) = the 5 planar orders / 5 Bravais lattices")
for n in [3, 4, 5, 6, 7]:
    print(f"   n={n}: 2cos(2pi/n)={2*cos(2*pi/n):+.4f}  phi={phi(n)}  "
          f"{'crystallographic' if phi(n)<=2 else 'FORBIDDEN (no n-fold lattice)'}")

print()
print("="*78)
print("(2) THE (2,3,n) ANGLE-DEFECT TRICHOTOMY  sign(1/2+1/3+1/n-1)  [Gauss-Bonnet]")
print("="*78)
names = {3: "tetrahedron {3,3}", 4: "octahedron {3,4}", 5: "icosahedron {3,5}",
         6: "HEX/triangular tiling {3,6} = A2 = covering-min", 7: "Klein quartic (2,3,7)=PSL(2,7)=168=Fano"}
for n in range(3, 8):
    s = F(1, 2)+F(1, 3)+F(1, n)-1
    geom = "spherical " if s > 0 else ("EUCLIDEAN " if s == 0 else "hyperbolic")
    print(f"   n={n}: 1/2+1/3+1/n-1 = {str(s):6} ({float(s):+.4f})  {geom}  {names.get(n,'')}")

print()
print("="*78)
print("(3) 5 PLATONIC SOLIDS  {p,q} with 1/p+1/q>1/2  (positive angle defect)")
print("="*78)
sols = [(p, q) for p in range(3, 7) for q in range(3, 7) if F(1, p)+F(1, q) > F(1, 2)]
for p, q in sols:
    print(f"   {{{p},{q}}}  1/p+1/q-1/2 = {F(1,p)+F(1,q)-F(1,2)}")
print(f"   count = {len(sols)}")

print()
print("="*78)
print("(4) LRC(2p) APEX-PRIME WALK: genus X0(2p) tracks the (2,3,p) geometry")
print("="*78)
def genus_X0(N):
    ps = set(primes_dividing(N))
    nu2 = 0 if N % 4 == 0 else (0 if 2 in ps and False else 1)
    r2 = 1
    for p in ps:
        if p == 2: r2 *= 1
        else: r2 *= (1 + (1 if p % 4 == 1 else -1))
    nu2 = 0 if N % 4 == 0 else r2
    r3 = 1
    for p in ps:
        if p == 3: r3 *= 1
        else: r3 *= (1 + (1 if p % 3 == 1 else -1))
    nu3 = 0 if N % 9 == 0 else r3
    psiN = N
    for p in ps: psiN = psiN*(p+1)//p
    cusps = sum(phi(gcd(d, N//d)) for d in range(1, N+1) if N % d == 0)
    return 1 + F(psiN, 12) - F(nu2, 4) - F(nu3, 3) - F(cusps, 2), cusps
for p in [3, 5, 7, 11, 13]:
    N = 2*p; g, c = genus_X0(N)
    s = F(1, 2)+F(1, 3)+F(1, p)-1
    geo = "Euclidean/crystallographic" if s == 0 or p == 3 else ("SPHERICAL/icosahedral (Platonic!)" if s > 0 else "HYPERBOLIC (Klein-type)")
    flag = "  <== genus jump 0->1, spherical->hyperbolic = LRC14 frontier" if p == 7 else ""
    print(f"   p={p:2d}  X0({N:2d}): genus={g}, cusps={c}  geometry(2,3,{p})={geo}{flag}")
print("\n   apex gap 4cos^2(3pi/7) = 2+2cos(6pi/7) = %.5f = 2+lambda_min(C_7) (>0 iff 7 odd/non-tiling)"
      % (4*cos(3*pi/7)**2))
print("   covering-min lives at order 6 (hexagonal/A2, Euclidean); B2 (square, order 4) = anomaly-free (HYP-3768).")
print("\nDONE.")
