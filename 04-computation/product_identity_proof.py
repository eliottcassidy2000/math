#!/usr/bin/env python3
"""
product_identity_proof.py — opus-2026-03-13-S70

Investigate the identity ∏_{k=1}^m (1 + Q_k) = F_p for the Interval tournament.

Where:
  Q_k = |S_hat(k)|^2 for S = {1,...,m}, m = (p-1)/2
  F_p = Fibonacci number (the p-th one)

The Q_k for Interval are: Q_k = (sin(m*pi*k/p)/sin(pi*k/p))^2
  = m^2 * F_m(2*pi*k/p) where F_m is the Fejér kernel.

Also: e_1(Q) = sum(Q_k) = m (from the known identity)
      e_m(Q) = prod(Q_k) = 1 (disc(Q) = p^{m-1} implies prod = ... actually verify)

Key identity to prove: prod_{k=1}^m (1 + Q_k) = F_p

Connection to Chebyshev / Morgan-Voyce:
  Q_k = |sum_{s=1}^m omega^{sk}|^2 where omega = e^{2pi*i/p}

  Let z = omega^k. Then S_hat(k) = z + z^2 + ... + z^m = z(1-z^m)/(1-z).
  Q_k = |z(1-z^m)/(1-z)|^2 = |1-z^m|^2/|1-z|^2
       = (2-2cos(2pi*mk/p))/(2-2cos(2pi*k/p))
       = sin^2(m*pi*k/p)/sin^2(pi*k/p)

Product: prod_{k=1}^m (1 + Q_k) = prod (1 + sin^2(mθ_k)/sin^2(θ_k))
  where θ_k = pi*k/p.

  = prod (sin^2(θ_k) + sin^2(mθ_k)) / sin^2(θ_k)

Since p = 2m+1:
  sin^2(θ_k) + sin^2(mθ_k) = sin^2(pi*k/p) + sin^2(m*pi*k/p)

Note: mθ_k = m*pi*k/(2m+1). And (m+1)θ_k = (m+1)*pi*k/(2m+1).
  m*pi*k/(2m+1) + pi*k/(2m+1) = (m+1)*pi*k/(2m+1)
  m*pi*k/(2m+1) = pi/2 * k * 2m/(2m+1)

Actually, let's just verify computationally and look for algebraic structure.
"""

import numpy as np
from math import comb, factorial

def fib(n):
    """Fibonacci number F_n (F_1=1, F_2=1, F_3=2, ...)"""
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a

def morgan_voyce_B(m, x):
    """Morgan-Voyce B polynomial: B_m(x) = sum_{j=0}^{m} C(m+j,2j) x^j"""
    return sum(comb(m+j, 2*j) * x**j for j in range(m+1))

def interval_Q(p):
    """Q_k for Interval tournament on Z_p."""
    m = (p-1)//2
    return [
        (np.sin(m*np.pi*k/p)/np.sin(np.pi*k/p))**2
        for k in range(1, m+1)
    ]

print("="*70)
print("PRODUCT IDENTITY: prod(1+Q_k) = F_p for Interval")
print("="*70)

for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    m = (p-1)//2
    Q = interval_Q(p)
    prod_val = 1.0
    for q in Q:
        prod_val *= (1 + q)
    fp = fib(p)
    ratio = prod_val / fp
    print(f"  p={p:2d}, m={m:2d}: prod(1+Q)={prod_val:20.6f}, F_p={fp:15d}, ratio={ratio:.12f}")

# Now: can we prove this algebraically?
# prod_{k=1}^m (1 + sin^2(mθ)/sin^2(θ)) where θ = πk/p
# = prod (sin^2(θ) + sin^2(mθ)) / prod sin^2(θ)
#
# The denominator: prod_{k=1}^m sin^2(πk/p) = (p/2^{p-1})^2 ??? No.
# Actually, prod_{k=1}^{p-1} sin(πk/p) = p/2^{p-1}.
# So prod_{k=1}^m sin(πk/p) * prod_{k=m+1}^{p-1} sin(πk/p) = p/2^{p-1}
# By symmetry sin(πk/p) = sin(π(p-k)/p), so prod_{k=1}^m = prod_{k=m+1}^{p-1}
# => [prod_{k=1}^m sin(πk/p)]^2 = p/2^{p-1}
# => prod_{k=1}^m sin^2(πk/p) = p/2^{p-1}

print("\n" + "="*70)
print("DENOMINATOR: prod sin^2(πk/p)")
print("="*70)

for p in [5, 7, 11, 13, 17, 19, 23, 29]:
    m = (p-1)//2
    denom = 1.0
    for k in range(1, m+1):
        denom *= np.sin(np.pi*k/p)**2
    expected = p / 2**(p-1)
    print(f"  p={p:2d}: prod sin^2 = {denom:.12e}, p/2^(p-1) = {expected:.12e}, ratio = {denom/expected:.10f}")

# So the product identity becomes:
# prod(sin^2(θ) + sin^2(mθ)) / (p/2^{p-1}) = F_p
# => prod(sin^2(θ) + sin^2(mθ)) = F_p * p / 2^{p-1}

print("\n" + "="*70)
print("NUMERATOR: prod (sin^2(θ) + sin^2(mθ))")
print("="*70)

for p in [5, 7, 11, 13, 17, 19, 23, 29]:
    m = (p-1)//2
    numer = 1.0
    for k in range(1, m+1):
        th = np.pi * k / p
        numer *= np.sin(th)**2 + np.sin(m*th)**2
    expected = fib(p) * p / 2**(p-1)
    print(f"  p={p:2d}: prod(s² + s²m) = {numer:.10e}, F_p*p/2^(p-1) = {expected:.10e}, ratio = {numer/expected:.10f}")

# Alternative: use the Chebyshev representation.
# sin(mθ)/sin(θ) = U_{m-1}(cos θ) where U is Chebyshev of 2nd kind.
# So Q_k = U_{m-1}(cos θ_k)^2 where θ_k = πk/p.
# And 1 + Q_k = 1 + U_{m-1}(cos θ_k)^2.
# = 1 + U_{m-1}(x_k)^2 where x_k = cos(πk/p), k=1,...,m.
#
# The x_k are roots of U_{p-1}(x)/(some polynomial).
# Actually, cos(πk/p) for k=1,...,p-1 are roots of U_{p-1}(x)/2 = 0?
# No: U_{n-1}(cos θ) = sin(nθ)/sin(θ). So U_{p-1}(cos(πk/p)) = sin(πk)/sin(πk/p) = 0 for k=1,...,p-1.
# So x_k = cos(πk/p) for k=1,...,p-1 are ALL roots of U_{p-1}(x).
# Since U_{p-1} has degree p-1, and we have p-1 roots, those are all the roots.
#
# For k=1,...,m, x_k = cos(πk/p) > 0 (since πk/p < π/2 for k < p/2).
# These are the POSITIVE roots of U_{p-1}.

print("\n" + "="*70)
print("CHEBYSHEV CONNECTION: Q_k = U_{m-1}(x_k)^2")
print("="*70)

from numpy.polynomial.chebyshev import chebval

for p in [7, 11, 13]:
    m = (p-1)//2
    print(f"\n  p={p}, m={m}:")
    for k in range(1, m+1):
        xk = np.cos(np.pi*k/p)
        # U_{m-1}(x) computed via recurrence
        U = [1, 2*xk]  # U_0, U_1
        for j in range(2, m):
            U.append(2*xk*U[-1] - U[-2])
        Um1 = U[m-1]
        Q_direct = (np.sin(m*np.pi*k/p)/np.sin(np.pi*k/p))**2
        print(f"    k={k}: x_k={xk:.6f}, U_{m-1}(x)={Um1:.6f}, U^2={Um1**2:.6f}, Q_k={Q_direct:.6f}")

# So we want: prod_{k=1}^m (1 + U_{m-1}(x_k)^2) = F_p
# where x_k are the m positive roots of U_{p-1}(x).
#
# This is a beautiful identity! Let's see if we can relate it to
# the RESULTANT or some polynomial identity.
#
# Define P(x) = 1 + U_{m-1}(x)^2. Then prod P(x_k) = F_p.
# Since x_k are roots of U_{p-1}(x), this is related to the resultant
# of P and U_{p-1}.

# Actually, let me try a different approach.
# Using the identity: 1 + U_{m-1}(cos θ)^2 = 1 + sin^2(mθ)/sin^2(θ)
# = (sin^2θ + sin^2(mθ))/sin^2θ
#
# Now p = 2m+1, so mθ = m*πk/p = πk*(p-1)/(2p).
# And (m+1)θ = πk*(p+1)/(2p).
# sin(mθ) = sin(πk/2 - πk/(2p)) = sin(πk/2)cos(πk/(2p)) - cos(πk/2)sin(πk/(2p))
# Hmm, this is getting messy. Let me try cos(2mθ) = cos((p-1)θ) approach.
# cos((p-1)θ) = cos(pθ - θ) = cos(pθ)cos(θ) + sin(pθ)sin(θ)
# At θ = πk/p: cos(pθ) = cos(πk) = (-1)^k, sin(pθ) = 0.
# So cos((p-1)θ_k) = (-1)^k * cos(θ_k).
# Similarly, sin((p-1)θ_k) = sin(pθ_k - θ_k) = -sin(pθ_k)cos(θ_k) + cos(pθ_k)sin(θ_k)
# = (-1)^k sin(θ_k).
# So sin^2((p-1)θ_k) = sin^2(θ_k). -> sin((p-1)θ) = ±sin(θ). Makes sense since (p-1)θ = πk - θ.
# And sin(πk - θ) = sin(θ) when k is odd, -sin(θ) when k is even (wrong, sin(πk-θ) = (-1)^{k+1} sin(-θ)... )
# sin(πk - θ) = sin(πk)cos(θ) - cos(πk)sin(θ) = -(-1)^k sin(θ) = (-1)^{k+1} sin(θ)
# So sin^2((p-1)θ) = sin^2(θ). This confirms Q_{p-k} = Q_k.

# Key identity to prove:
# prod_{k=1}^m (sin^2(θ_k) + sin^2(mθ_k)) = F_p * p / 2^{p-1}
#
# Let's try using complex exponentials.
# sin^2(θ) = (1-cos(2θ))/2
# sin^2(mθ) = (1-cos(2mθ))/2
# sum = 1 - (cos(2θ)+cos(2mθ))/2 = 1 - cos((m+1)θ)*cos((m-1)θ)
# using product-to-sum: cos A + cos B = 2cos((A+B)/2)cos((A-B)/2)
# cos(2θ)+cos(2mθ) = 2cos((m+1)θ)cos((m-1)θ)
# So sin^2(θ) + sin^2(mθ) = 1 - cos((m+1)θ)cos((m-1)θ)

print("\n" + "="*70)
print("KEY FACTORIZATION: sin²θ + sin²(mθ) = 1 - cos((m+1)θ)cos((m-1)θ)")
print("="*70)

for p in [7, 11, 13]:
    m = (p-1)//2
    print(f"\n  p={p}, m={m}:")
    for k in range(1, m+1):
        th = np.pi*k/p
        lhs = np.sin(th)**2 + np.sin(m*th)**2
        rhs = 1 - np.cos((m+1)*th)*np.cos((m-1)*th)
        print(f"    k={k}: sin²θ+sin²mθ = {lhs:.10f}, 1-cos(m+1)θ·cos(m-1)θ = {rhs:.10f}, match: {abs(lhs-rhs)<1e-12}")

# So prod (1+Q_k) = prod (1 - cos((m+1)θ_k)*cos((m-1)θ_k)) / prod sin^2(θ_k)
#
# Now (m+1)θ_k = (m+1)πk/p = (p+1)/2 * πk/p, and p = 2m+1
# So (m+1)πk/p = πk/2 + πk/(2p)
# And (m-1)θ_k = (m-1)πk/p = πk/2 - 3πk/(2p) ... hmm these aren't clean.
#
# Let me try: at θ_k = πk/p,
# cos((m+1)θ_k) = cos((m+1)πk/p)
# cos((m-1)θ_k) = cos((m-1)πk/p)
#
# Note (m+1) = (p+1)/2 and (m-1) = (p-3)/2.
# cos((p+1)πk/(2p)) * cos((p-3)πk/(2p))
# = cos(πk/2 + πk/(2p)) * cos(πk/2 - 3πk/(2p))
# Using cos(A+B)cos(A-B) = cos²A - sin²B:
# Wait no, cos(A+B)cos(A-C) ≠ simple.
#
# Let me try a substitution. Let α = (m+1)πk/p, β = (m-1)πk/p.
# α + β = 2mπk/p = (p-1)πk/p = πk - πk/p
# α - β = 2πk/p
# cos(α)cos(β) = [cos(α-β) + cos(α+β)]/2
# = [cos(2πk/p) + cos(πk - πk/p)]/2
# = [cos(2πk/p) + (-1)^k cos(πk/p)]/2  (since cos(πk - x) = (-1)^k cos(x))
# = [cos(2θ_k) + (-1)^k cos(θ_k)]/2
# = [(2cos²θ_k - 1) + (-1)^k cos(θ_k)]/2
# = cos²θ_k - 1/2 + (-1)^k cos(θ_k)/2
#
# So 1 - cos(α)cos(β) = 1 - cos²θ_k + 1/2 - (-1)^k cos(θ_k)/2
# = sin²θ_k + 1/2 - (-1)^k cos(θ_k)/2
# = sin²θ_k + (1 - (-1)^k cos(θ_k))/2
# Hmm, but we said 1 - cos(α)cos(β) = sin²θ + sin²mθ. Let me verify...
# sin²θ + sin²mθ = (1-cos2θ)/2 + (1-cos2mθ)/2 = 1 - (cos2θ + cos2mθ)/2
# = 1 - cos((m+1)θ)cos((m-1)θ) ... need to check.
# cos(2θ) + cos(2mθ) = 2cos((m+1)θ)cos((m-1)θ)
# So sin²θ + sin²mθ = 1 - cos((m+1)θ)cos((m-1)θ). ✓

# Substituting the expression for cos(α)cos(β):
# sin²θ + sin²mθ = sin²θ + (1 - (-1)^k cosθ)/2
# So sin²mθ = (1 - (-1)^k cosθ)/2. Let me check:
# sin²mθ at θ=πk/p: sin(mπk/p)² = sin((p-1)πk/(2p))² = sin(πk/2 - πk/(2p))²
# Hmm this should not simplify to (1-(-1)^k cosθ)/2 in general...

# Actually let me recheck. cos2θ + cos2mθ = 2cos((m+1)θ)cos((m-1)θ)?
# cos A + cos B = 2cos((A+B)/2)cos((A-B)/2). A=2θ, B=2mθ.
# (A+B)/2 = (m+1)θ, (A-B)/2 = (1-m)θ = -(m-1)θ.
# So cos2θ + cos2mθ = 2cos((m+1)θ)cos((m-1)θ). ✓

print("\n" + "="*70)
print("PRODUCT OVER (1 - cos(α)cos(β)) for various p")
print("="*70)

for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
    m = (p-1)//2
    prod_num = 1.0
    for k in range(1, m+1):
        th = np.pi*k/p
        val = 1 - np.cos((m+1)*th)*np.cos((m-1)*th)
        prod_num *= val

    # Compare with F_p * p / 2^{p-1}
    fp = fib(p)
    target = fp * p / 2**(p-1)
    ratio = prod_num / target
    print(f"  p={p:2d}: prod = {prod_num:.10e}, F_p*p/2^(p-1) = {target:.10e}, ratio = {ratio:.10f}")

# Now let's try to find a POLYNOMIAL whose roots give us this product.
# The product prod_{k=1}^m (1 - cos(α_k)cos(β_k)) where α_k = (m+1)πk/p, β_k = (m-1)πk/p.
#
# Define f(x) = 1 - cos((m+1)arccos(x)) * cos((m-1)arccos(x))
# At x_k = cos(πk/p), f(x_k) = 1 - cos(α_k)cos(β_k).
# cos((m+1)arccos(x)) = T_{m+1}(x) (Chebyshev)
# cos((m-1)arccos(x)) = T_{m-1}(x)
# So f(x) = 1 - T_{m+1}(x)*T_{m-1}(x).
#
# Now prod_{k=1}^m f(x_k) where x_k are specific roots of U_{p-1}...
# Actually x_k = cos(πk/p) for k=1,...,m are NOT roots of U_{p-1} (those would be cos(πk/(p+1))).
# cos(πk/p) are roots of sin(pθ)/sin(θ) when θ=πk/p, k=1,...,p-1.
# sin(pπk/p) = sin(πk) = 0, and sin(πk/p) ≠ 0. So U_{p-1}(cos(πk/p)) = sin(pπk/p)/sin(πk/p) = 0. ✓
# So x_k = cos(πk/p) for k=1,...,p-1 are the p-1 roots of U_{p-1}(x).
# We only take k=1,...,m (the positive roots).

print("\n" + "="*70)
print("CHEBYSHEV POLYNOMIAL APPROACH")
print("="*70)
print("f(x) = 1 - T_{m+1}(x) * T_{m-1}(x)")
print("Product over positive roots of U_{p-1}")

# Identity: T_{m+1}(x) * T_{m-1}(x) = (T_{2m}(x) + T_2(x)) / 2  (incorrect?)
# Product formula: T_a * T_b = (T_{a+b} + T_{|a-b|})/2
# T_{m+1} * T_{m-1} = (T_{2m} + T_2)/2
# So f(x) = 1 - (T_{2m}(x) + T_2(x))/2
# = 1 - T_{2m}(x)/2 - T_2(x)/2
# = 1 - T_{2m}(x)/2 - (2x² - 1)/2
# = 1 - T_{2m}(x)/2 - x² + 1/2
# = 3/2 - x² - T_{2m}(x)/2
#
# Now 2m = p-1, so T_{2m} = T_{p-1}.
# At x = cos(πk/p): T_{p-1}(cos(πk/p)) = cos((p-1)πk/p) = cos(πk - πk/p) = (-1)^k cos(πk/p) = (-1)^k x_k
# So f(x_k) = 3/2 - x_k² - (-1)^k x_k / 2

for p in [7, 11, 13]:
    m = (p-1)//2
    print(f"\n  p={p}, m={m}:")
    for k in range(1, m+1):
        xk = np.cos(np.pi*k/p)
        f_direct = 1 - np.cos((m+1)*np.pi*k/p)*np.cos((m-1)*np.pi*k/p)
        f_cheb = 3/2 - xk**2 - (-1)**k * xk / 2
        print(f"    k={k}: f_direct={f_direct:.10f}, 3/2-x²-(-1)^k*x/2 = {f_cheb:.10f}, match: {abs(f_direct-f_cheb)<1e-10}")

# BEAUTIFUL! So we need:
# prod_{k=1}^m (3/2 - x_k² - (-1)^k x_k/2) = F_p * p / 2^{p-1}
# where x_k = cos(πk/p) for k=1,...,m.
#
# The (-1)^k factor mixes even and odd k, making this a product over two interleaved sets.
# But maybe we can split into even-k and odd-k contributions.

print("\n" + "="*70)
print("SPLIT BY k PARITY")
print("="*70)

for p in [7, 11, 13, 17, 19, 23, 29]:
    m = (p-1)//2
    prod_even = 1.0
    prod_odd = 1.0
    for k in range(1, m+1):
        xk = np.cos(np.pi*k/p)
        if k % 2 == 0:
            f = 3/2 - xk**2 - xk/2  # (-1)^k = 1
            prod_even *= f
        else:
            f = 3/2 - xk**2 + xk/2  # (-1)^k = -1
            prod_odd *= f
    print(f"  p={p:2d}: prod_odd = {prod_odd:.10f}, prod_even = {prod_even:.10f}, total = {prod_odd*prod_even:.10e}")

# For odd k: f = 3/2 - x² + x/2 = (-2x² + x + 3)/2 = -(2x²-x-3)/2 = -(2x-3)(x+1)/2
# Wait: 2x²-x-3 = (2x-3)(x+1). So 3/2 - x² + x/2 = -(2x-3)(x+1)/(-2) = (2x-3)(x+1)/(-2).
# Hmm let me redo: 3/2 - x² + x/2 = 1/2 * (3 - 2x² + x) = 1/2 * (-2x² + x + 3)
# -2x² + x + 3 = -(2x² - x - 3) = -(2x-3)(x+1) ✓ (since 2*3/2-3 = 0 and 2*(-1)+3 = 1, hmm)
# Actually: 2x²-x-3 at x=3/2: 2*9/4-3/2-3 = 9/2-3/2-3 = 3-3=0. And at x=-1: 2+1-3=0. ✓
# So for odd k: f = -(2x_k-3)(x_k+1)/2
# For even k: f = 3/2 - x² - x/2 = 1/2*(3-2x²-x) = -1/2*(2x²+x-3) = -1/2*(2x+3)(x-1)
# Check: 2x²+x-3 at x=1: 2+1-3=0. At x=-3/2: 2*9/4-3/2-3 = 9/2-3/2-3 = 3-3=0. ✓
# So for even k: f = -(2x_k+3)(x_k-1)/2

print("\n" + "="*70)
print("FACTORED FORM VERIFICATION")
print("="*70)

for p in [7, 11, 13]:
    m = (p-1)//2
    print(f"\n  p={p}, m={m}:")
    for k in range(1, m+1):
        xk = np.cos(np.pi*k/p)
        f = 3/2 - xk**2 - (-1)**k * xk/2
        if k % 2 == 1:
            f_factored = -(2*xk-3)*(xk+1)/2
        else:
            f_factored = -(2*xk+3)*(xk-1)/2
        print(f"    k={k}: f = {f:.10f}, factored = {f_factored:.10f}, match: {abs(f-f_factored)<1e-10}")

# Since x_k ∈ (-1, 1), and specifically x_k = cos(πk/p) with k=1,...,m:
# For odd k: -(2x_k-3)(x_k+1)/2. Since x_k < 3/2 and x_k > -1 (except x_k=cos(π(m)/p)≈0),
#   2x_k-3 < 0 and x_k+1 > 0, so f > 0. ✓
# For even k: -(2x_k+3)(x_k-1)/2. Since x_k < 1 and 2x_k+3 > 0,
#   (x_k-1) < 0, so f > 0. ✓

# Now: prod_{k odd} (-(2x_k-3)(x_k+1)/2) * prod_{k even} (-(2x_k+3)(x_k-1)/2)
# = prod_{k=1}^m 1/2 * prod_{k odd} -(2x_k-3)(x_k+1) * prod_{k even} -(2x_k+3)(x_k-1)
# = (1/2)^m * (-1)^m * prod_{k odd} (2x_k-3)(x_k+1) * prod_{k even} (2x_k+3)(x_k-1)

# This is still complicated. Let me try a completely different approach.
#
# We have the polynomial whose roots are x_k = cos(πk/p):
# The minimal polynomial for cos(2π/p) over Q is the p-th cyclotomic polynomial
# evaluated at (x + 1/x)/2... but this gets complicated.
#
# Simpler: U_{p-1}(x) = prod_{k=1}^{p-1} (x - cos(πk/p)) * leading coefficient
# U_{p-1}(x) = 2^{p-1} * prod_{k=1}^{p-1} (x - cos(πk/p))
# The product of the positive roots gives us part of the story.
# Since the roots come in pairs (cos(πk/p), cos(π(p-k)/p) = -cos(πk/p))... wait no.
# cos(π(p-k)/p) = cos(π - πk/p) = -cos(πk/p). So roots come in ± pairs!
# For k=1,...,m: x_k = cos(πk/p), x_{p-k} = -x_k.
# So prod_{k=1}^{p-1} (x-x_k) = prod_{k=1}^m (x-x_k)(x+x_k) = prod_{k=1}^m (x² - x_k²).
# And U_{p-1}(x) = 2^{p-1} * prod_{k=1}^m (x² - x_k²). But p-1 = 2m, so this works.

# The sum of x_k² for k=1,...,m:
print("\n" + "="*70)
print("SUMS OF cos²(πk/p)")
print("="*70)

for p in [5, 7, 11, 13, 17, 19]:
    m = (p-1)//2
    s1 = sum(np.cos(np.pi*k/p) for k in range(1, m+1))
    s2 = sum(np.cos(np.pi*k/p)**2 for k in range(1, m+1))
    s3 = sum(np.cos(np.pi*k/p)**3 for k in range(1, m+1))
    print(f"  p={p}: sum cos = {s1:.10f}, sum cos² = {s2:.10f} = (p-1)/4 = {(p-1)/4:.1f}?, sum cos³ = {s3:.10f}")

# sum cos²(πk/p) for k=1,...,m:
# cos²θ = (1+cos2θ)/2. sum = m/2 + (1/2)*sum cos(2πk/p) for k=1,...,m.
# sum_{k=1}^m cos(2πk/p) = Re(sum_{k=1}^m ω^k) where ω=e^{2πi/p}.
# Since p is odd, sum_{k=1}^m ω^k = (ω - ω^{m+1})/(1-ω) = ... well for p=2m+1:
# sum_{k=0}^{p-1} ω^k = 0, so sum_{k=1}^{p-1} ω^k = -1.
# sum_{k=1}^m ω^k + sum_{k=m+1}^{p-1} ω^k = -1.
# sum_{k=m+1}^{p-1} = sum_{j=1}^m ω^{p-j} = sum_{j=1}^m conj(ω^j) (since ω^p=1).
# So sum + conj(sum) = -1, i.e., 2*Re(sum) = -1, so Re(sum) = -1/2.
# Therefore sum cos²(πk/p) = m/2 + (-1/2)/2 = m/2 - 1/4 = (2m-1)/4 = (p-2)/4.

print("\n  Expected: sum cos² = (p-2)/4:")
for p in [5, 7, 11, 13, 17, 19]:
    m = (p-1)//2
    s2 = sum(np.cos(np.pi*k/p)**2 for k in range(1, m+1))
    expected = (p-2)/4
    print(f"    p={p}: {s2:.10f} vs {expected:.10f}, match: {abs(s2-expected)<1e-10}")

print("\nDONE.")
