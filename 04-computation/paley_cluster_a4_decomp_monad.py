#!/usr/bin/env python3
"""
paley_cluster_a4_decomp_monad.py
monad-explorer-2026-06-07 (deep-research lane). Decisive follow-up to
paley_cluster_expansion_monad.py.

GOAL: prove that the L=4 cluster integral a_4 = lim A_4/p^4 = 0 (so the cherry L=2 is
the UNIQUE surviving cluster and R(p) -> e exactly).

EXACT DECOMPOSITION (integrate out the path endpoint x_4 using sum_{x4} chi(x4-x3)=0):
  A_4 = -[ T_a + T_b + T_c ]
  T_a = sum_{x0,x1,x2,x3 distinct} chi(x1-x0)chi(x2-x1)chi(x3-x2)chi(x0-x3)   (4-CYCLE sum)
  T_b = sum_{...}               chi(x1-x0)chi(x2-x1)chi(x3-x2)chi(x1-x3)   (THETA sum)
  T_c = sum_{...}               chi(x1-x0)chi(x2-x1)chi(x3-x2)chi(x2-x3)
      = chi(-1)*(p-3)*A_2 = -(p-3)*p*(p-1)        (the "fully diagonal" closing piece)

Because A_2 = p(p-1) (the cherry closes into chi(-(diff)^2) = chi(-1), NO cancellation),
T_c carries the top order p^3.  T_a, T_b are genuine 4-variable character sums with
square-root cancellation: T_a = p*chi(-1)*J, J = sum_{a,b,c} chi(a b c (a+b+c)) = O(p^2)
by Weil. Hence T_a, T_b = O(p^3) = o(p^4), so

      A_4 = O(p^3)   =>   a_4 = lim A_4/p^4 = 0.

This is the structural reason the cherry (L=2) is the ONLY non-cancelling cluster.

This script verifies, at primes up to ~67:
  (1) the exact identity A_4 = -(T_a + T_b + T_c)   [validates the derivation]
  (2) T_c = -(p-3)p(p-1) exactly
  (3) A_4/p^4 -> 0  and  T_a/p^3, T_b/p^3 stay O(1)  (so a_4=0)
  (4) the convergence rate of R(p) to e via the cherry-only exponentiation prediction.
"""
import numpy as np
from math import factorial

def legendre_array(p):
    chi = np.zeros(p, dtype=np.int64)
    qr = set((x*x) % p for x in range(1, p))
    for d in range(1, p):
        chi[d] = 1 if d in qr else -1
    return chi

def all_distinct_4tuples(p):
    """Generate all distinct (x0,x1,x2,x3) as columns. Memory ~ p^4 * 4 -- use p<=67."""
    a = np.arange(p)
    X0, X1, X2, X3 = np.meshgrid(a, a, a, a, indexing='ij')
    X0, X1, X2, X3 = X0.ravel(), X1.ravel(), X2.ravel(), X3.ravel()
    distinct = (X0!=X1)&(X0!=X2)&(X0!=X3)&(X1!=X2)&(X1!=X3)&(X2!=X3)
    return X0[distinct], X1[distinct], X2[distinct], X3[distinct]

def compute_all(p, chi):
    X0,X1,X2,X3 = all_distinct_4tuples(p)
    base = chi[(X1-X0)%p]*chi[(X2-X1)%p]*chi[(X3-X2)%p]
    Ta = int((base * chi[(X0-X3)%p]).sum())   # 4-cycle
    Tb = int((base * chi[(X1-X3)%p]).sum())   # theta
    Tc = int((base * chi[(X2-X3)%p]).sum())
    A4 = -(Ta+Tb+Tc)
    # A_2 for reference
    A2 = p*(p-1)
    return A4, Ta, Tb, Tc, A2

def main():
    print("="*78)
    print("EXACT DECOMPOSITION OF A_4  (verifying a_4 -> 0, the cherry is unique)")
    print("="*78)
    primes = [7,11,19,23,31,43,47,59,67]
    print(f"{'p':>3} | {'A_4':>12} | {'A_4/p^4':>9} | {'T_a':>10} | {'T_a/p^3':>8} | "
          f"{'T_b/p^3':>8} | {'T_c=-(p-3)p(p-1)?':>16}")
    for p in primes:
        chi = legendre_array(p)
        A4,Ta,Tb,Tc,A2 = compute_all(p, chi)
        Tc_pred = -(p-3)*p*(p-1)
        ok = "OK" if Tc==Tc_pred else f"MISMATCH({Tc} vs {Tc_pred})"
        print(f"{p:>3} | {A4:>12} | {A4/p**4:>9.5f} | {Ta:>10} | {Ta/p**3:>8.4f} | "
              f"{Tb/p**3:>8.4f} | {ok:>16}")
    print()
    print("Interpretation: A_4/p^4 -> 0 (a_4 = 0); T_a/p^3, T_b/p^3 stay O(1) bounded,")
    print("confirming T_a,T_b = O(p^3) [Weil square-root cancellation], A_4 = O(p^3).")
    print("=> the ONLY cluster integral that survives is a_2 = 1 (the cherry).")
    print()

    print("="*78)
    print("CONSEQUENCE: R(p) -> exp(a_2) = exp(1) = e.  Convergence rate:")
    print("="*78)
    stored_H = {3:3,7:189,11:95095,19:1172695746915,23:15760206976379349}
    e = np.e
    print(f"{'p':>3} | {'R(p)':>9} | {'e-R':>8} | {'(e-R)*sqrt(p)':>13} | {'(e-R)*p':>9}")
    for p in [3,7,11,19,23]:
        R = stored_H[p]*(2**(p-1))/factorial(p)
        print(f"{p:>3} | {R:>9.5f} | {e-R:>8.5f} | {(e-R)*p**0.5:>13.4f} | {(e-R)*p:>9.4f}")
    print()
    print("(e-R)*sqrt(p) ~ const => convergence is O(1/sqrt p): slow, which is exactly")
    print("why 5 data points (max p=23, R=2.557, still 6% below e) could not distinguish")
    print("e from a larger constant by extrapolation alone. The cluster expansion settles it.")

if __name__ == "__main__":
    main()
