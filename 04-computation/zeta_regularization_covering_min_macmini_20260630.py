#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S67
=======================
THE COVERING-MIN IS A ZETA-REGULARIZATION CARRIER (owner's seed, verified + synthesized).

T = 1+2+...+(n-1) = n(n-1)/2  (the finite speed-sum).
  Phi6 = n^2-n+1 = 2T+1                 (sum-of-naturals identity)
  killer = n(n-1) = 2T = -1 (mod Phi6)  (the reflection anchor)
  s(n,Phi6) = -T/(12T+6) = -1/(12+6/T)  (Mobius in T) -> -1/12 = zeta(-1)   (Ramanujan 1+2+3+...)
  margin = n/Phi6 - 1/n = -12 s/n^2 ,  n^2*margin -> -12 zeta(-1) = 1.
  zeta(-1) = -1/12 = -B2/2 ; 12 = -1/zeta(-1) ; 6 = 1/B2 ; psi(14)=24=2*12 = eta exponent (Delta=eta^24).
  The -1/12 anomaly is HEXAGONAL (order 6, n^3=-1); the SQUARE (order 4, h^2=-1) gives s=0 (no anomaly)
  => the LRC margin EXISTS because the covering-min is hexagonal.
"""
from fractions import Fraction as F
from math import gcd

def dsum(h, k):
    def saw(x):
        x = x - int(x)
        return F(0) if x == 0 else x - F(1, 2)
    return sum((saw(F(i, k))*saw(F((h*i) % k, k)) for i in range(1, k)), F(0))

def psi(N):
    r = N
    for p in set(f for f in range(2, N+1) if N % f == 0 and all(f % d for d in range(2, f))):
        r = r*(p+1)//p
    return r

print("(A) SEED IDENTITIES (exact):")
for n in range(3, 15):
    T = n*(n-1)//2; P = n*n-n+1; s = dsum(n % P, P)
    ok = (2*T+1 == P) and ((n*(n-1)) % P == P-1) and (s == F(-T, 12*T+6))
    print(f"  n={n:2d}: T={T:3d} Phi6=2T+1={P:3d} killer=2T=-1 s=-T/(12T+6)={str(s):9}  all-exact={ok}")

print("\n(B) MARGIN carries the regularization: n^2*margin -> -12 zeta(-1) = 1")
for n in [10, 50, 200, 1000, 5000]:
    P = n*n-n+1; m = F(n, P)-F(1, n)
    print(f"  n={n:5d}: n^2*margin = {float(n*n*m):.6f}")

print("\n(C) the regularization constants: zeta(-1)=-1/12=-B2/2, 12=-1/zeta(-1), 6=1/B2, psi(14)=%d=2*12 (eta exponent)" % psi(14))

print("\n(D) hexagonal (order 6) -> s -> -1/12 = zeta(-1);  square (order 4, h^2=-1) -> s = 0 (no anomaly)")
for n in [4, 7, 14]:
    P = n*n-n+1; print(f"  hex order-6 n={n:2d}: s={dsum(n % P, P)}  -> -1/12")
for k in [5, 13, 17]:
    h = [x for x in range(2, k) if (x*x) % k == k-1][0]
    print(f"  square order-4 k={k:2d}: h={h} (h^2=-1): s={dsum(h, k)}  (no regularization => zero margin)")
print("\nDONE.")
