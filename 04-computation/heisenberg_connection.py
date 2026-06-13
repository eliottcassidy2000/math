#!/usr/bin/env python3
"""
MAJOR DISCOVERY: β_m for Paley P_{2m+1} = b_2 of Heisenberg Lie algebra h_{2n+1}

THM-130: β_m(P_p) = m(m-3)/2 where m = (p-1)/2

Santharoubane's theorem (1983): b_2(h_{2n+1}) = C(2n,2) - 1 = n(2n-1) - 1

Under re-indexing n = (m-1)/2:
  n(2n-1) - 1 = ((m-1)/2)(m-2) - 1 = (m²-3m+2-2)/2 = m(m-3)/2

EXACT MATCH.

OEIS: A014106 = 0, 5, 14, 27, 44, 65, ...

Reference: L.J. Santharoubane, "Cohomology of Heisenberg Lie algebras",
           Proc. Amer. Math. Soc. 87 (1983), 23-28.

This connection appears to be NEW — no paper links Paley tournament
path homology to Heisenberg Lie algebra cohomology.

opus-2026-03-13-S71b
"""

print("="*70)
print("HEISENBERG LIE ALGEBRA ↔ PALEY TOURNAMENT PATH HOMOLOGY")
print("="*70)

print("""
The Heisenberg Lie algebra h_{2n+1} is the (2n+1)-dimensional Lie algebra
with generators x_1,...,x_n, y_1,...,y_n, z satisfying:
  [x_i, y_j] = δ_{ij} z,   all other brackets = 0

Its Betti numbers (Santharoubane 1983):
  b_k(h_{2n+1}) = C(2n,k) - C(2n,k-2)

In particular:
  b_0 = 1
  b_1 = 2n
  b_2 = C(2n,2) - 1 = n(2n-1) - 1
""")

# Verify the match
print(f"{'m':>4} {'p=2m+1':>8} {'β_m=m(m-3)/2':>14} {'n=(m-1)/2':>10} {'b_2(h_{2n+1})':>14} {'Match':>6}")
print("-" * 60)

for m in range(3, 22, 2):
    p = 2 * m + 1
    beta_m = m * (m - 3) // 2
    n = (m - 1) // 2
    b2 = n * (2 * n - 1) - 1
    match = "✓" if beta_m == b2 else "✗"
    print(f"{m:4d} {p:8d} {beta_m:14d} {n:10d} {b2:14d} {match:>6}")

print("""
KEY OBSERVATIONS:

1. The match is EXACT for all m (algebraic identity, not just numerical).

2. The Heisenberg algebra h_{2n+1} has dimension 2n+1 = m.
   The Paley tournament P_p has m = (p-1)/2 quadratic residues.
   So the "Heisenberg dimension" equals the QR count.

3. The Heisenberg algebra is the SIMPLEST non-abelian nilpotent Lie algebra.
   Its structure (central extension of R^{2n} by R via a symplectic form)
   mirrors the QR structure: the Legendre symbol (a/p) is a ±1 form on Z_p*.

4. The FULL Santharoubane formula gives all Betti numbers of h_{2n+1}:
     b_k = C(2n,k) - C(2n,k-2)
   For THM-130:
     β_m = b_2 = C(2n,2) - 1 = C(m-1,2) - 1
     β_{m+1} = C(m+1,2) = C(m+1,2)

   Is β_{m+1} also a Heisenberg Betti number?
""")

# Check: does β_{m+1} = C(m+1,2) match b_3(h_{2n+1})?
print("Checking β_{m+1} against Heisenberg b_k:")
print(f"{'m':>4} {'β_{m+1}=C(m+1,2)':>18} {'b_3(h_{m})':>12} {'b_2(h_{m+2})':>14}")
print("-" * 60)

from math import comb

for m in range(3, 22, 2):
    beta_mp1 = m * (m + 1) // 2
    n = (m - 1) // 2
    # b_3(h_{2n+1}) = C(2n,3) - C(2n,1)
    b3 = comb(2*n, 3) - 2*n if 2*n >= 3 else 0
    # b_2(h_{m+2}) = b_2(h_{2n+3}) where n' = n+1
    n2 = n + 1
    b2_next = n2 * (2*n2 - 1) - 1
    print(f"{m:4d} {beta_mp1:18d} {b3:12d} {b2_next:14d}")

print("""
β_{m+1} = C(m+1,2) does NOT match b_3(h_{2n+1}).
It DOES match b_2(h_{2n+3}) = b_2(h_{m+2})!

So: β_m = b_2(h_m),  β_{m+1} = b_2(h_{m+2})

This means the SHIFT β_{m+1} - β_m = p-1 = 2m corresponds to:
  b_2(h_{m+2}) - b_2(h_m) = [n'(2n'-1)-1] - [n(2n-1)-1]
  where n' = n+1
  = (n+1)(2n+1) - n(2n-1) = 2n²+3n+1 - 2n²+n = 4n+1 ...
  Wait: 4n+1 ≠ 2m = 2(2n+1).

Actually (n+1)(2(n+1)-1)-1 - (n(2n-1)-1) = (n+1)(2n+1) - n(2n-1)
= 2n²+3n+1 - 2n²+n = 4n+1. And 2m = 2(2n+1) = 4n+2.
So the shift is 4n+1, not 4n+2 = 2m. Close but NOT exact.

The connection works perfectly for β_m but the β_{m+1} match is off by 1.
This suggests β_m is the "core" Heisenberg connection, and β_{m+1}
receives additional contributions from the eigenspace structure (the p-1
nonzero eigenspaces each contributing 1).
""")

# Deeper: the Heisenberg Poincaré polynomial
print("="*70)
print("HEISENBERG POINCARÉ POLYNOMIAL")
print("="*70)

for m in [3, 5, 9]:
    n = (m - 1) // 2
    print(f"\nm={m} (n={n}): h_{2*n+1} = h_{m}")
    poincare = []
    for k in range(2*n + 2):
        bk = comb(2*n, k) - (comb(2*n, k-2) if k >= 2 else 0)
        poincare.append(bk)
    print(f"  Betti: {poincare}")
    print(f"  χ(h_m) = {sum((-1)**k * poincare[k] for k in range(len(poincare)))}")

print("""
CONCLUSION:

The formula β_m = m(m-3)/2 for Paley P_p is PRECISELY b_2(h_m),
the second Betti number of the Heisenberg Lie algebra of dimension m.

This suggests a deep structural connection: the k=0 eigenspace of
the Paley path complex at degree m might be (co)homologically
related to the Heisenberg Lie algebra h_m.

The Heisenberg algebra h_m has a natural SYMPLECTIC structure
(the bracket is defined by a symplectic form on R^{m-1}).
The Legendre symbol on Z_p* also has a form-like structure.

OPEN QUESTION: Is there a functorial connection between:
  (a) The k=0 eigenspace chain complex of P_p at degrees m, m+1
  (b) The Chevalley-Eilenberg complex of h_m at degrees 2, 3

If so, the FULL Santharoubane formula b_k = C(m-1,k) - C(m-1,k-2)
might give Betti numbers for MORE degrees of the Paley path complex.
""")
