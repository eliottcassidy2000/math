#!/usr/bin/env python3
"""
mod2_cumulant_connection.py - What does THM-094 (F mod 2) imply for cumulants?

THM-094 (kind-pasteur): F(T,x) = (1+x)^{n-1} mod 2 for ALL tournaments.
THM-095 (opus-S46e): coeff(t_{2k+1} in κ_{2k}) = 2/C(n,2k).

QUESTION: Does the mod-2 universality of F imply anything about cumulants mod 2?

The moments μ_r = E[fwd^r] are determined by F_k = #{perms with fwd=k}.
If F_k ≡ C(n-1,k) mod 2, then:
  n! · E[fwd^r] = Σ_k k^r · F_k ≡ Σ_k k^r · C(n-1,k) mod 2

The RHS is the r-th moment of the binomial(n-1, 1/2) distribution times 2^{n-1}.

Let me compute this and see what cumulant identities emerge mod 2.

Author: opus-2026-03-07-S46e
"""
from math import comb, factorial
from fractions import Fraction

def F_mod2(n, k):
    """F_k(T) mod 2 = C(n-1,k) mod 2 (by THM-094)"""
    return comb(n-1, k) % 2

def moment_mod2(n, r):
    """n! · E[fwd^r] mod 2 = Σ_k k^r · F_k mod 2"""
    total = 0
    for k in range(n):
        total += pow(k, r) * comb(n-1, k)
    return total % 2

print("=" * 60)
print("MOMENTS OF F mod 2")
print("=" * 60)

for n in range(3, 10):
    print(f"\nn={n}:")
    for r in range(1, 8):
        # Exact moment sum
        exact = sum(pow(k, r) * comb(n-1, k) for k in range(n))
        print(f"  Σ k^{r} C({n-1},k) = {exact}, mod 2 = {exact % 2}")

# Now compute the ACTUAL cumulant implications
print("\n" + "=" * 60)
print("CUMULANT CONSEQUENCES mod 2")
print("=" * 60)

for n in range(3, 10):
    # Mean = Σ k·C(n-1,k) / 2^{n-1} ... but this is for binomial
    # For tournaments: mean = (n-1)/2 (exactly)
    # variance = (n+1)/12 + 4t₃/(n(n-1))

    # The variance mod 2:
    # 2/C(n,2) · t₃ + Bernoulli term
    # At mod 2: 2/C(n,2) ≡ ? mod 2

    cn2 = comb(n, 2)
    # 2/C(n,2) is a fraction. In the integer moment:
    # n! · variance = n! · [(n+1)/12 + 4t₃/(n(n-1))]
    # = n!·(n+1)/12 + n!·4t₃/(n(n-1))
    # = n!·(n+1)/12 + 4·(n-2)!·t₃

    nfact_var_const = factorial(n) * (n+1) // 12
    nfact_var_t3_coeff = 4 * factorial(n-2)

    print(f"\nn={n}: n!·Var constant part = {nfact_var_const} (mod 2 = {nfact_var_const % 2})")
    print(f"  n!·Var t₃ coeff = {nfact_var_t3_coeff} (mod 2 = {nfact_var_t3_coeff % 2})")

    # If both are even, then n!·Var ≡ 0 mod 2 regardless of t₃!
    if nfact_var_const % 2 == 0 and nfact_var_t3_coeff % 2 == 0:
        print(f"  → n!·Var ≡ 0 mod 2 for ALL tournaments (consistent with THM-094)")

# More generally: compute n!·μ_r mod 2 from THM-094
print("\n" + "=" * 60)
print("n!·E[fwd^r] mod 2 from THM-094")
print("=" * 60)

for n in range(3, 10):
    print(f"\nn={n}:")
    # From THM-094: F_k ≡ C(n-1,k) mod 2
    # So n!·E[fwd^r] = Σ_k k^r F_k
    # This is tournament-independent mod 2!
    for r in range(1, min(n+1, 8)):
        val = sum(pow(k, r) * comb(n-1, k) for k in range(n))
        # The actual value for a specific tournament would be Σ k^r F_k(T)
        # which ≡ val mod 2 by THM-094
        print(f"  n!·E[fwd^{r}] ≡ {val % 2} mod 2 (universal)")

# KEY: Does n!·E[fwd^r] ≡ 0 mod 2 for all r, or can it be 1?
print("\n" + "=" * 60)
print("PATTERN: Is n!·E[fwd^r] always even?")
print("=" * 60)

for n in range(2, 12):
    vals = []
    for r in range(1, 12):
        val = sum(pow(k, r) * comb(n-1, k) for k in range(n))
        vals.append(val % 2)
    print(f"n={n:2d}: r=1..11: {vals}")
    # Is there any odd value?
    if any(v == 1 for v in vals):
        print(f"  ⚠ NOT always even!")

# Now: connection to centered moments and cumulants
print("\n" + "=" * 60)
print("CENTERED MOMENTS mod 2")
print("=" * 60)

for n in range(3, 9):
    mean_num = n - 1  # mean = (n-1)/2, so 2·mean = n-1
    # centered fwd = fwd - (n-1)/2
    # n! · E[(fwd - mean)^r] = Σ_k (k - (n-1)/2)^r · F_k
    # This involves half-integers when n is even...
    # Better: 2^r · n! · E[(fwd-mean)^r] = Σ_k (2k-(n-1))^r · F_k
    # is always an integer

    print(f"\nn={n}:")
    for r in range(2, min(n+1, 7)):
        val = sum(pow(2*k - (n-1), r) * comb(n-1, k) for k in range(n))
        # This is 2^r · n! · μ_r(centered)
        # which by THM-094 is universal mod 2
        print(f"  2^{r}·n!·μ_{r} = {val}, mod 2 = {val % 2}")

print("\n" + "=" * 60)
print("INTERPRETATION")
print("=" * 60)
print("""
THM-094 (F mod 2 universal) implies that ALL integer moment sums
  Σ_k k^r · F_k(T) ≡ Σ_k k^r · C(n-1,k) mod 2
are tournament-independent mod 2.

Combined with THM-095 (2/C(n,2k) coefficient), this means:
- The cumulant hierarchy's t₃ dependence enters at order 2 in 2-adic valuation
  (the coefficient 2/C(n,2k) contributes exactly one factor of 2)
- The mod-2 universality "zeroes out" all tournament-dependent terms

This is the mod-2 shadow of the OCF: since H = I(Ω,2) and we're evaluating
at x=2 ≡ 0 mod 2, the independence polynomial collapses!
  H mod 2 = 1 + 0 + 0 + ... = 1 mod 2 (Rédei!)

So: OCF at x=2 gives H. OCF at x=2 mod 2 gives H mod 2 = 1 = Rédei.
The entire cumulant hierarchy collapses mod 2 because OCF uses x=2.
""")
