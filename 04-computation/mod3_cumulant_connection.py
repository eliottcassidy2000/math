#!/usr/bin/env python3
"""
mod3_cumulant_connection.py - Connect THM-086 (mod 3 Taylor zeros) to cumulant hierarchy.

THM-086: c_j(T) = 0 mod 3 for j < val(n), where val(n) = 2*floor((n-1)/2).
For n odd: F(T,x) mod 3 = alpha*(x-1)^{n-1} — single free parameter.
For n even: F(T,x) mod 3 = alpha*(x-1)^{n-2} + beta*(x-1)^{n-1} — two free parameters.

The cumulant hierarchy (THM-093):
  kappa_2 = (n+1)/12 + 4*t3/(n(n-1))
  kappa_4 = -(n+1)/120 + 2*(t5+2*alpha2)/C(n,4) - 48*t3^2/(n(n-1))^2
  kappa_6 = (n+1)/252 + 2*t7/C(n,6) + cross terms

QUESTION: What does the mod 3 structure imply for cumulants?

The Taylor coefficients c_j = [x^j] F(T, x+1) (shifting x -> x+1).
The moments E[fwd^r] = sum_k k^r * F_k / n!
The cumulants are functions of the moments.

If c_j = 0 mod 3 for j < val(n), this constrains the moments mod 3,
which constrains the cumulants mod 3.

Let me work out the exact implications.

Author: opus-2026-03-07-S46e
"""
from fractions import Fraction
from math import factorial, comb
from itertools import permutations, combinations
from collections import defaultdict

# === Part 1: Taylor coefficients and moments ===
# F(T, x) = sum_k F_k * x^k
# c_j(T) = sum_k F_k * C(k, j) (Taylor coefficients around x=1)
# Actually c_j = [x^j] F(T, x+1) = sum_k F_k * C(k, j)

# Moments: mu_r = E[fwd^r] = (1/n!) * sum_k k^r * F_k
# Express k^r in terms of falling factorials: k^r = sum_j S(r,j) * C(k,j) * j!
# where S(r,j) = Stirling numbers of the second kind.
# So mu_r = (1/n!) * sum_j S(r,j) * j! * c_j

# Therefore: mu_r = (1/n!) * sum_{j=0}^{r} S(r,j) * j! * c_j

# If c_j = 0 mod 3 for j < val(n), then:
# mu_r = (1/n!) * (sum_{j=val(n)}^{r} S(r,j)*j!*c_j) mod (3/n!)
# But we need to be careful with the rational arithmetic.

def stirling2(n, k):
    """Stirling numbers of the second kind."""
    if n == 0 and k == 0:
        return 1
    if n == 0 or k == 0:
        return 0
    if k > n:
        return 0
    # Use recurrence: S(n,k) = k*S(n-1,k) + S(n-1,k-1)
    S = [[0]*(k+1) for _ in range(n+1)]
    S[0][0] = 1
    for i in range(1, n+1):
        for j in range(1, min(i, k)+1):
            S[i][j] = j * S[i-1][j] + S[i-1][j-1]
    return S[n][k]

print("Stirling numbers S(r,j) for r=1..6:")
for r in range(1, 7):
    row = [stirling2(r, j) for j in range(r+1)]
    print(f"  r={r}: {row}")

# === Part 2: Work out implications at n=5 ===
print("\n" + "=" * 60)
print("n=5: val(5) = 4. c_0=c_1=c_2=c_3 = 0 mod 3")
print("=" * 60)

# mu_r = (1/120) * sum_j S(r,j)*j!*c_j
# For r <= 3: contributions only from c_0,...,c_3 which are all 0 mod 3.
# So 120*mu_r = 0 mod 3 for r <= 3.
# Since 120 = 8*15 = 8*3*5, and gcd(120, 3) = 3:
# mu_r = (120*mu_r) / 120. If 120*mu_r = 0 mod 3, then 120*mu_r / 120 may not be 0 mod 3
# because 120 has factor 3.

# Let's be more precise. Let N! = n! and let M_r = n! * mu_r = sum_k k^r * F_k.
# M_r is always an integer.
# M_r = sum_{j=0}^{r} S(r,j) * j! * c_j
# If c_j = 0 mod 3 for j < 4, then:
# M_1 = S(1,0)*0!*c_0 + S(1,1)*1!*c_1 = 0 + c_1 (mod 3) = 0 mod 3
# M_2 = S(2,0)*0!*c_0 + S(2,1)*1!*c_1 + S(2,2)*2!*c_2 = 0 + c_1 + 2*c_2 (mod 3) = 0 mod 3
# M_3 = S(3,0)*0!*c_0 + S(3,1)*1!*c_1 + S(3,2)*2!*c_2 + S(3,3)*3!*c_3
#      = 0 + c_1 + 6*c_2 + 6*c_3 = c_1 + 0*c_2 + 0*c_3 (mod 3)
#      = c_1 (mod 3) = 0 mod 3.

print("\nM_r = n!*mu_r mod 3, given c_j=0 mod 3 for j<4:")
for r in range(1, 7):
    terms = []
    total_mod3 = 0
    for j in range(r+1):
        coeff = stirling2(r, j) * factorial(j)
        terms.append(f"  j={j}: S({r},{j})*{j}! = {coeff} (mod 3: {coeff%3})")
        if j < 4:
            total_mod3 += 0  # c_j = 0 mod 3
        # Only c_4 contributes at r >= 4
    print(f"\nr={r}:")
    for t in terms:
        print(t)
    # At which j does the first non-vanishing contribution come?
    first_nonvanishing = None
    for j in range(4, r+1):
        coeff = stirling2(r, j) * factorial(j)
        if coeff % 3 != 0:
            first_nonvanishing = j
            break
    if first_nonvanishing:
        print(f"  First non-vanishing (mod 3) from c_{first_nonvanishing}")
    else:
        print(f"  M_{r} = 0 mod 3 regardless of c_j values!")

# === Part 3: Verify numerically at n=5 ===
print("\n" + "=" * 60)
print("n=5: NUMERICAL VERIFICATION")
print("=" * 60)

n = 5
m_edges = n*(n-1)//2
sample_data = []

for bits in range(1 << m_edges):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    # Compute F_k distribution
    F = [0]*n
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        F[fwd] += 1

    # Taylor coefficients c_j = sum_k F_k * C(k, j)
    c = [0]*n
    for j in range(n):
        c[j] = sum(F[k] * comb(k, j) for k in range(n))

    # Moments M_r = sum_k k^r * F_k (integer moments, = n! * mu_r)
    M = [0]*7
    for r in range(7):
        M[r] = sum(k**r * F[k] for k in range(n))

    sample_data.append({'c': c[:], 'M': M[:], 'F': F[:]})

# Check c_j mod 3 for j < 4
print("c_j mod 3 distribution for j=0,1,2,3:")
for j in range(4):
    mods = [d['c'][j] % 3 for d in sample_data]
    from collections import Counter
    print(f"  c_{j} mod 3: {Counter(mods)}")

print("\nc_4 mod 3 distribution:")
mods4 = [d['c'][4] % 3 for d in sample_data]
print(f"  c_4 mod 3: {Counter(mods4)}")

# Check M_r mod 3
print("\nM_r mod 3 (integer moments):")
for r in range(1, 7):
    mods = [d['M'][r] % 3 for d in sample_data]
    print(f"  M_{r} mod 3: {Counter(mods)}")

# Check M_r mod 9
print("\nM_r mod 9:")
for r in range(1, 7):
    mods = [d['M'][r] % 9 for d in sample_data]
    print(f"  M_{r} mod 9: {Counter(mods)}")

# === Part 4: What do the cumulants look like mod 3? ===
print("\n" + "=" * 60)
print("n=5: CUMULANTS mod 3")
print("=" * 60)

# kappa_2 = (n+1)/12 + 4*t3/(n*(n-1))
# At n=5: kappa_2 = 1/2 + t3/5
# So 120*kappa_2 = 60 + 24*t3
# mod 3: 120*kappa_2 = 0 + 0*t3 = 0 mod 3

# kappa_4: more complex. Let me compute directly.
print("kappa_2 * n! for each tournament mod 3:")
k2_mod3 = Counter()
k4_mod3 = Counter()
for d in sample_data:
    F = d['F']
    total = factorial(n)
    mu = Fraction(n-1, 2)
    moments = {}
    for r in range(1, 7):
        moments[r] = sum(Fraction((k - mu)**r * F[k], total) for k in range(n))
    k2 = moments[2]
    k4 = moments[4] - 3 * moments[2]**2
    # Multiply by n!^2 to get integers? No, kappa_2 is rational.
    # Let's look at n! * kappa_2
    nk2 = k2 * total  # integer
    nk4 = k4 * total  # rational (kappa_4 * n!)
    k2_mod3[int(nk2) % 3] += 1
    # k4 is rational. Let me find LCD.
    # kappa_4 * n!: may not be integer.

# Actually, let me work with centered moments instead
print("\nCentered moments n!*mu_r mod 3:")
for r in [2, 3, 4, 5, 6]:
    vals = Counter()
    for d in sample_data:
        F = d['F']
        total = factorial(n)
        mu_val = Fraction(n-1, 2)
        M_centered = sum(Fraction((k - mu_val)**r * F[k]) for k in range(n))
        # M_centered = n! * mu_r (centered)
        # Is this always an integer? mu = (n-1)/2, so (k-mu)^r has denominator 2^r
        # n! * mu_r = sum (k - (n-1)/2)^r * F_k, which has denominator 2^r
        # Multiply by 2^r to get integer
        M_int = M_centered * 2**r
        vals[int(M_int) % 3] += 1
    print(f"  2^{r} * n! * mu_{r} mod 3: {vals}")

# === Part 5: The mod 3 structure at n=7 ===
print("\n" + "=" * 60)
print("n=7: val(7)=6. CUMULANT IMPLICATIONS")
print("=" * 60)

# At n=7: c_0=...=c_5 = 0 mod 3.
# The moments M_r = sum S(r,j)*j!*c_j
# For r<=5: all c_j in the sum are 0 mod 3.
# But S(r,j)*j! may have factors of 3...

print("Which M_r are forced to 0 mod 3 at n=7 (val=6)?")
for r in range(1, 10):
    # M_r = sum_{j=0}^{r} S(r,j)*j!*c_j
    # c_j = 0 mod 3 for j < 6
    # So M_r mod 3 = sum_{j=6}^{r} S(r,j)*j!*c_j mod 3
    contributions = []
    for j in range(min(6, r+1)):
        coeff = stirling2(r, j) * factorial(j)
        if coeff % 3 != 0:
            contributions.append(j)
    forced_zero = (len(contributions) == 0 or all(j < 6 for j in contributions))
    # Actually if all coefficients with j < 6 have S(r,j)*j! divisible by 3,
    # then M_r is NOT forced to 0 mod 3 (could be nonzero from c_6+).
    # We want: for j < 6, S(r,j)*j! contributions are 0 mod 3 (because c_j = 0 mod 3).
    # For j >= 6, c_j may be anything.
    # So M_r is forced to be sum_{j>=6} S(r,j)*j!*c_j mod 3.
    # This is 0 mod 3 iff all S(r,j)*j! for j>=6 are 0 mod 3, or iff r < 6.
    if r < 6:
        print(f"  M_{r}: FORCED 0 mod 3 (all contributing c_j vanish)")
    else:
        nonzero_high = []
        for j in range(6, r+1):
            coeff = stirling2(r, j) * factorial(j)
            if coeff % 3 != 0:
                nonzero_high.append((j, coeff % 3))
        if nonzero_high:
            print(f"  M_{r}: free (nonzero coeffs from c_{nonzero_high[0][0]}+)")
        else:
            print(f"  M_{r}: FORCED 0 mod 3 (all high coefficients divisible by 3!)")

# === Part 6: Key insight ===
print("\n" + "=" * 60)
print("KEY INSIGHT")
print("=" * 60)
print("""
THM-086 says: c_j = 0 mod 3 for j < val(n).

This implies: M_r = n!*E[fwd^r] = 0 mod 3 for r < val(n).
(Because M_r = sum_{j<=r} S(r,j)*j!*c_j, and all c_j in range are 0 mod 3.)

Since val(n) = 2*floor((n-1)/2) = n-1 for n odd, n-2 for n even:

For ODD n: M_1 through M_{n-2} are ALL 0 mod 3.
  -> All moments up to order n-2 are divisible by 3 (in integer form).
  -> All centered moments of order <= n-2 are similarly constrained.
  -> ALL EVEN CUMULANTS kappa_2, kappa_4, ..., kappa_{n-3} are constrained mod 3.

For EVEN n: M_1 through M_{n-3} are ALL 0 mod 3.
  -> Similarly for moments/cumulants up to order n-3.

The cumulant hierarchy kappa_{2k} encodes cycle invariants t_{2k+1}.
The mod 3 vanishing constrains these cumulants.
The TWO structures are COMPATIBLE because:
  - Cumulants are polynomial in cycle counts
  - Taylor coefficients are polynomial in the same data
  - The mod 3 constraint propagates through both representations

QUESTION: Does the mod 3 constraint on cumulants give a NEW proof
of THM-086, or vice versa?
""")
