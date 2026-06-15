#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The slope-s diagonal family of Pascal's triangle, its growth-constant ladder, and
the hard-core-gas / independence-polynomial reading.  kind-pasteur-2026-06-15.
(User dispatch: extend the family 2^n / Fibonacci / construction-3 and find the
recursive patterns.)

a_s(n) = sum_k C(n - s*k, k)  = diagonal sum of Pascal's triangle at slope s.
Recurrence (Pascal identity):  a_s(n) = a_s(n-1) + a_s(n-s-1).
  s=0 -> 2^n;  s=1 -> Fibonacci;  s=2 -> Narayana's cows (A000930);
  s=3 -> A003269;  s=4 -> A003520; ...
Combinatorial: a_s(n) counts placements with the recurrence "position n is empty
(a(n-1)) or occupied, forcing the previous s positions empty (a(n-s-1))" = the
1-D hard-core lattice gas with exclusion radius s at fugacity 1 = independent sets
of the s-th power of a path = binary strings with all 1s >= s+1 apart (up to a
boundary/index convention).
Growth constant beta_s = dominant root of x^{s+1} = x^s + 1:
  2 (s=0), phi (s=1, golden), 1.46557 (s=2, supergolden), 1.38028 (s=3),
  1.32472 (s=4, PLASTIC, since x^5-x^4-1 = (x^3-x-1)(x^2-x+1)), ... -> 1.
At general fugacity x the same recurrence is the path-power independence polynomial
I(P_n^s, x) = I(P_{n-1}^s,x) + x*I(P_{n-1-s}^s,x); s=1 is the repo path I.P.
(THM-485 two-temperatures: x=1 Fibonacci/phi, x=2 H/Jacobsthal); the OCF H=I(Omega,2)
is the same object-type at x=2.  HYP-614's phi is the s=1 rung of beta_s.
"""
import sys
from math import comb
import numpy as np
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def a_s(s,n): return sum(comb(n-s*k,k) for k in range(0,n//(s+1)+2) if n-s*k>=k>=0)

def beta(s):
    if s==0: return 2.0
    c=[0]*(s+2); c[0]=1; c[1]=-1; c[-1]=-1   # x^{s+1} - x^s - 1
    return max(r.real for r in np.roots(c) if abs(r.imag)<1e-9 and r.real>0.5)

def main():
    print("=== a_s(n) (n=0..12) ===")
    for s in range(6):
        print(f"  s={s}: {[a_s(s,n) for n in range(13)]}")
    print("\n=== recurrence a_s(n)=a_s(n-1)+a_s(n-s-1) ===")
    for s in range(6):
        print(f"  s={s}: {all(a_s(s,n)==a_s(s,n-1)+a_s(s,n-s-1) for n in range(s+1,40))}")
    print("\n=== growth constant beta_s = dominant root of x^(s+1)=x^s+1 ===")
    for s in range(7):
        b=beta(s); seq=[a_s(s,n) for n in range(50,53)]
        print(f"  s={s}: beta={b:.6f}  measured ratio={seq[-1]/seq[-2]:.6f}")
    print("\n  plastic check: rho^3=rho+1 => rho^5=rho^4+1 ?",
          "x^5-x^4-1 == (x^3-x-1)(x^2-x+1):",
          np.allclose((np.poly1d([1,0,-1,-1])*np.poly1d([1,-1,1])).c, [1,-1,0,0,0,-1]))
    print("\n=== auxiliary (s=0): central binomials C(n,floor(n/2)) = A001405 ===")
    print("  ", [comb(n,n//2) for n in range(12)])
    print("  (these carry the sub-exponential ~beta^n / sqrt(n) correction; the 'width'/dominant column)")

if __name__=="__main__":
    main()
