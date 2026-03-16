#!/usr/bin/env python3
"""
rapidity_numbertheory_s116b.py
Creative number theory connections through rapidity = ln(n)/2.

Seven explorations:
1. Rapidity and continued fractions (sqrt(2), phi, e trajectories)
2. Rapidity and quadratic residues (QR vs NQR average rapidity)
3. Rapidity and Bernoulli numbers (spectral zeta at negative integers)
4. Collatz in rapidity (trajectory of n=27)
5. Prime gaps in rapidity (vs Cramer bound)
6. Goldbach in rapidity (decomposition geometry)
7. sqrt rapidity and perfect square detection

kind-pasteur-2026-03-15-S116b
"""

import math
from fractions import Fraction
from collections import defaultdict

def rapidity(n):
    """Rapidity of positive real n = ln(n)/2."""
    if n <= 0:
        return float('-inf')
    return math.log(n) / 2

def fmt(x, d=8):
    return f"{x:>{d+4}.{d}f}"

# =============================================================================
print("=" * 70)
print("RAPIDITY AND NUMBER THEORY: SEVEN EXPLORATIONS")
print("=" * 70)

# =============================================================================
# 1. RAPIDITY AND CONTINUED FRACTIONS
# =============================================================================
print("\n" + "=" * 70)
print("1. RAPIDITY AND CONTINUED FRACTIONS")
print("=" * 70)

def cf_convergents(cf_terms, max_terms=30):
    """Compute convergents p_n/q_n from CF coefficients [a0; a1, a2, ...]."""
    convergents = []
    # p_{-1} = 1, p_{-2} = 0; q_{-1} = 0, q_{-2} = 1
    p_prev2, p_prev1 = 0, 1
    q_prev2, q_prev1 = 1, 0
    for i, a in enumerate(cf_terms[:max_terms]):
        p = a * p_prev1 + p_prev2
        q = a * q_prev1 + q_prev2
        convergents.append((p, q))
        p_prev2, p_prev1 = p_prev1, p
        q_prev2, q_prev1 = q_prev1, q
    return convergents

def cf_sqrt2(n):
    """CF for sqrt(2) = [1; 2, 2, 2, ...]."""
    return [1] + [2] * (n - 1)

def cf_phi(n):
    """CF for golden ratio phi = [1; 1, 1, 1, ...]."""
    return [1] * n

def cf_e(n):
    """CF for e = [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, ...]."""
    terms = [2]
    k = 1
    while len(terms) < n:
        terms.extend([1, 2 * k, 1])
        k += 1
    return terms[:n]

print("\n  Continued fraction convergents and their rapidities:")
print("  rapidity(p_n/q_n) = ln(p_n/q_n)/2")
print()

targets = [
    ("sqrt(2)", cf_sqrt2(20), math.sqrt(2)),
    ("phi    ", cf_phi(20), (1 + math.sqrt(5)) / 2),
    ("e      ", cf_e(20), math.e),
]

for name, cf, true_val in targets:
    convs = cf_convergents(cf)
    true_rap = rapidity(true_val)
    print(f"  {name} = {true_val:.10f}, true rapidity = {true_rap:.10f}")
    print(f"    n  p_n/q_n        rapidity       error from true    ratio_successive")
    prev_err = None
    for i, (p, q) in enumerate(convs[:15]):
        val = p / q
        rap = rapidity(val)
        err = rap - true_rap
        ratio_str = ""
        if prev_err is not None and abs(err) > 1e-15:
            ratio_str = f"  {prev_err / err:+.4f}"
        print(f"    {i:2d}  {p:>8d}/{q:<8d}  {rap:>12.8f}  {err:>+14.2e}{ratio_str}")
        prev_err = err
    print()

# Rapidity differences between successive convergents
print("  RAPIDITY JUMPS between successive convergents:")
print("  delta_n = rapidity(p_{n+1}/q_{n+1}) - rapidity(p_n/q_n)")
print()
for name, cf, true_val in targets:
    convs = cf_convergents(cf)
    print(f"  {name}:")
    for i in range(min(12, len(convs) - 1)):
        p1, q1 = convs[i]
        p2, q2 = convs[i + 1]
        rap1 = rapidity(p1 / q1)
        rap2 = rapidity(p2 / q2)
        delta = rap2 - rap1
        # The jump magnitude should relate to 1/(q_n * q_{n+1})
        expected = 1 / (2 * q1 * q2)  # ~ delta rapidity for |p/q - alpha| ~ 1/(q*q')
        print(f"    delta_{i:2d} = {delta:>+12.8f}   |delta| = {abs(delta):.2e}   "
              f"1/(2*q_n*q_{{n+1}}) = {expected:.2e}   ratio = {abs(delta)/expected:.4f}")
    print()

# Key observation about CF rapidity
print("  KEY OBSERVATION: Rapidity of convergents oscillates around true value.")
print("  For sqrt(2): even convergents approach from below, odd from above.")
print("  The rapidity trajectory is a DAMPED OSCILLATION around ln(alpha)/2.")
print("  Damping rate governed by q_n growth: for phi, q_n = F_n (Fibonacci),")
print("  so delta ~ 1/F_n^2 ~ phi^{-2n}. EXPONENTIAL convergence in rapidity.")
print("  For e, the pattern 1,2k,1 creates a 3-periodic modulation.")

# =============================================================================
# 2. RAPIDITY AND QUADRATIC RESIDUES
# =============================================================================
print("\n" + "=" * 70)
print("2. RAPIDITY AND QUADRATIC RESIDUES")
print("=" * 70)

def quadratic_residues(p):
    """Return set of QRs mod p (excluding 0)."""
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    return sorted(qr)

def non_residues(p):
    """Return set of NQRs mod p."""
    qr = set(quadratic_residues(p))
    return sorted([a for a in range(1, p) if a not in qr])

print("\n  For prime p, QRs = {a : a is a square mod p}.")
print("  rapidity(a) = ln(a)/2 for each residue a in {1, ..., p-1}.")
print()

paley_primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

print(f"  {'p':>4s}  {'QRs':>30s}  {'avg_rap(QR)':>12s}  {'avg_rap(NQR)':>12s}  {'diff':>10s}  {'ratio':>8s}")
print("  " + "-" * 100)

for p in paley_primes:
    qr = quadratic_residues(p)
    nqr = non_residues(p)
    avg_qr = sum(rapidity(a) for a in qr) / len(qr)
    avg_nqr = sum(rapidity(a) for a in nqr) / len(nqr)
    diff = avg_nqr - avg_qr
    ratio = avg_nqr / avg_qr if avg_qr != 0 else float('inf')
    qr_str = str(qr) if len(str(qr)) <= 30 else str(qr)[:27] + "..."
    print(f"  {p:4d}  {qr_str:>30s}  {avg_qr:12.6f}  {avg_nqr:12.6f}  {diff:+10.6f}  {ratio:8.4f}")

print()
print("  DETAILED for Paley primes p = 3, 7, 11:")
for p in [3, 7, 11]:
    qr = quadratic_residues(p)
    nqr = non_residues(p)
    print(f"\n  p = {p}:")
    print(f"    QRs  = {qr}")
    print(f"    NQRs = {nqr}")

    # Sum of rapidities (exact via log product)
    qr_product = 1
    for a in qr:
        qr_product *= a
    nqr_product = 1
    for a in nqr:
        nqr_product *= a

    print(f"    Product of QRs  = {qr_product},  ln(product)/2 = {math.log(qr_product)/2:.8f}")
    print(f"    Product of NQRs = {nqr_product},  ln(product)/2 = {math.log(nqr_product)/2:.8f}")
    print(f"    Sum of QR rapidities  = {sum(rapidity(a) for a in qr):.8f}")
    print(f"    Sum of NQR rapidities = {sum(rapidity(a) for a in nqr):.8f}")

    # Wilson's theorem: product of all residues = (p-1)! = -1 mod p
    # So product(QRs) * product(NQRs) = (p-1)!
    total_product = math.factorial(p - 1)
    print(f"    (p-1)! = {total_product}")
    print(f"    product(QR) * product(NQR) = {qr_product * nqr_product} = (p-1)!? {qr_product * nqr_product == total_product}")

    # Total rapidity = rapidity((p-1)!) = ln((p-1)!)/2
    total_rap = rapidity(total_product)
    print(f"    Total rapidity = ln((p-1)!)/2 = {total_rap:.8f}")
    print(f"    = rapidity(QRs) + rapidity(NQRs) = {sum(rapidity(a) for a in qr) + sum(rapidity(a) for a in nqr):.8f}")

# QR rapidity asymptotic
print("\n  ASYMPTOTIC: As p -> infinity,")
print("    avg rapidity of QRs  ~ (1/|QR|) * sum_{a in QR} ln(a)/2")
print("    avg rapidity of NQRs ~ (1/|NQR|) * sum_{a in NQR} ln(a)/2")
print("  Since |QR| = |NQR| = (p-1)/2 and sum_all ln(a)/2 = ln((p-1)!)/2,")
print("  the TOTAL sum is fixed. The question is: do QRs cluster low or high?")
print()

# Check: are QRs biased toward smaller values?
print("  QR 'smallness' bias (fraction of QRs that are < p/2):")
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73]:
    qr = quadratic_residues(p)
    small = sum(1 for a in qr if a < p / 2)
    frac = small / len(qr)
    print(f"    p={p:3d}: {small}/{len(qr)} = {frac:.4f}   "
          f"{'BELOW 0.5' if frac < 0.5 else 'ABOVE 0.5' if frac > 0.5 else 'EXACT 0.5'}")

# =============================================================================
# 3. RAPIDITY AND BERNOULLI NUMBERS
# =============================================================================
print("\n" + "=" * 70)
print("3. RAPIDITY AND BERNOULLI NUMBERS / SPECTRAL ZETA")
print("=" * 70)

# Bernoulli numbers
def bernoulli_numbers(max_n):
    """Compute Bernoulli numbers B_0, B_1, ..., B_{max_n} using the recursive formula."""
    B = [Fraction(0)] * (max_n + 1)
    B[0] = Fraction(1)
    for m in range(1, max_n + 1):
        s = Fraction(0)
        for k in range(m):
            # binomial(m+1, k)
            binom = 1
            for j in range(k):
                binom = binom * (m + 1 - j) // (j + 1)
            s += binom * B[k]
        B[m] = -s / (m + 1)
    return B

BN = bernoulli_numbers(20)

print("\n  Bernoulli numbers B_n:")
for n in range(17):
    print(f"    B_{n:2d} = {str(BN[n]):>20s}   (= {float(BN[n]):>14.8f})")

print("\n  Riemann zeta at negative integers: zeta(-n) = -B_{n+1}/(n+1)")
print()
for n in range(1, 15):
    zeta_neg = -BN[n + 1] / (n + 1)
    print(f"    zeta(-{n:2d}) = -B_{n+1}/{n+1} = {str(zeta_neg):>15s} = {float(zeta_neg):>14.8f}")

print("\n  OUR SPECTRAL ZETA: zeta_M(s) = sum_{eigenvalues lambda_i} lambda_i^s")
print("  For the n=3 transfer matrix with eigenvalues {-1, -1, 1} (or similar):")
print("  The 'forbidden' values: zeta_M(-3) = 7, zeta_M(-5) = 21")
print()

# For a generic set of eigenvalues, compute spectral zeta
# The n=3 regular tournament transfer matrix has eigenvalues related to
# the cubic lambda^3 - lambda^2 - x*lambda - x = 0
# At x=2 (the tournament point): lambda^3 - lambda^2 - 2*lambda - 2 = 0
# Factor: (lambda^2 - 2)(lambda - 1) = 0 ... no
# Actually lambda^3 - lambda^2 - 2lambda - 2. Try lambda = -1: -1-1+2-2 = -2. No.
# lambda = 2: 8-4-4-2 = -2. No. Let's just compute numerically.

import numpy as np

print("  SPECTRAL ZETA FROM SIMPLE EIGENVALUE SETS:")
print()

# The simplest case: eigenvalues that give zeta(-3)=7, zeta(-5)=21
# Try: eigenvalues {a, b, c} with a^3 + b^3 + c^3 = 7, a^5 + b^5 + c^5 = 21
# Newton's identities: p_k = power sums, e_k = elementary symmetric
# e1 = a+b+c, e2 = ab+ac+bc, e3 = abc
# p1 = e1
# p2 = e1*p1 - 2*e2
# p3 = e1*p2 - e2*p1 + 3*e3 = 7
# p5 = e1*p4 - e2*p3 + e3*p2 = 21

# For eigenvalues {1, lambda, lambda^{-1}} with product 1:
# p_k = 1 + lambda^k + lambda^{-k} = 1 + 2*cosh(k*phi)
# where lambda = e^phi.
# p_3 = 1 + 2*cosh(3*phi) = 7 => cosh(3*phi) = 3 => 3*phi = acosh(3)
# phi = acosh(3)/3 = 1.76275/3 = 0.58758
# Then p_5 = 1 + 2*cosh(5*phi) = 1 + 2*cosh(5*0.58758) = 1 + 2*cosh(2.9379)
phi_spec = math.acosh(3) / 3
p3_check = 1 + 2 * math.cosh(3 * phi_spec)
p5_check = 1 + 2 * math.cosh(5 * phi_spec)
print(f"  If eigenvalues = {{1, e^phi, e^{{-phi}}}} with phi = acosh(3)/3:")
print(f"    phi = {phi_spec:.10f}")
print(f"    e^phi = {math.exp(phi_spec):.10f}")
print(f"    p_3 = 1 + 2*cosh(3*phi) = {p3_check:.10f} (want 7)")
print(f"    p_5 = 1 + 2*cosh(5*phi) = {p5_check:.10f} (want 21)")
print()

# So p_5 != 21 for this set. Let's try a different approach.
# What if the eigenvalues are {1, omega, omega^2} (cube roots of unity)?
# p_k = 1 + omega^k + omega^{2k} = 3 if 3|k, else 0.
# p_3 = 3, p_6 = 3. Not 7.

# Try eigenvalues {a, b} with a + b = s, ab = p:
# p_k(a,b) = a^k + b^k. p_3 = s^3 - 3sp = 7, p_5 = s^5 - 5s^3p + 5sp^2 = 21
# Two equations, two unknowns.

# More interesting: Bernoulli connection to power sums
# sum_{k=1}^{n} k^m = (1/(m+1)) * sum_{j=0}^{m} C(m+1,j) * B_j * n^{m+1-j}
# = (B(n+1)^{m+1} - B^{m+1}) / (m+1) in umbral notation
print("  POWER SUMS via Bernoulli (Faulhaber's formula):")
print("  S_m(n) = 1^m + 2^m + ... + n^m")
print()

def power_sum(m, n):
    """Exact power sum via Bernoulli numbers."""
    result = Fraction(0)
    for j in range(m + 1):
        # binomial(m+1, j)
        binom = Fraction(1)
        for k in range(j):
            binom = binom * (m + 1 - k) / (k + 1)
        result += binom * BN[j] * Fraction(n) ** (m + 1 - j)
    result /= (m + 1)
    return result

print(f"  {'m':>4s}  {'S_m(3)':>10s}  {'S_m(7)':>12s}  {'S_m(11)':>14s}  S_m(3) mod 7  S_m(7) mod 21")
print("  " + "-" * 80)
for m in range(1, 13):
    s3 = power_sum(m, 3)
    s7 = power_sum(m, 7)
    s11 = power_sum(m, 11)
    print(f"  {m:4d}  {str(s3):>10s}  {str(s7):>12s}  {str(s11):>14s}  "
          f"{int(s3) % 7:>11d}  {int(s7) % 21:>12d}")

# Connection to the forbidden numbers
print("\n  FORBIDDEN NUMBERS AND BERNOULLI POWER SUMS:")
print(f"    S_1(3) = 1+2+3 = 6 = 7-1")
print(f"    S_1(6) = 1+...+6 = 21")
print(f"    S_2(3) = 1+4+9 = 14 = 2*7")
print(f"    So 7 = S_1(3)+1, 21 = S_1(6) = 3*7 = T_6 (6th triangular number)")
print(f"    21 = T_6, and T_n = n(n+1)/2 via Bernoulli: -B_2 * C(n+1,2) + ...")
print()

# Rapidity of Bernoulli-related values
print("  RAPIDITY OF POWER SUMS S_m(p) for Paley primes p:")
for p in [3, 7, 11]:
    print(f"    p = {p}:")
    for m in [1, 2, 3, 5]:
        s = int(power_sum(m, p))
        print(f"      S_{m}({p}) = {s:>12d},  rapidity = {rapidity(s):>10.6f}")

# =============================================================================
# 4. COLLATZ IN RAPIDITY
# =============================================================================
print("\n" + "=" * 70)
print("4. COLLATZ TRAJECTORIES IN RAPIDITY SPACE")
print("=" * 70)

def collatz_trajectory(n, max_steps=1000):
    """Return Collatz trajectory starting from n."""
    traj = [n]
    while n != 1 and len(traj) < max_steps:
        if n % 2 == 0:
            n = n // 2
        else:
            n = 3 * n + 1
        traj.append(n)
    return traj

# Key: In rapidity, n/2 step changes rapidity by -ln(2)/2
# 3n+1 step changes rapidity by ln(3)/2 + ln(1 + 1/n)/2
# For large n: rapidity change ~ ln(3)/2 - ln(2)/2 = ln(3/2)/2 per odd step
# (since 3n+1 is always even, the next step is always /2)

ln2_half = math.log(2) / 2
ln3_half = math.log(3) / 2

print(f"\n  Rapidity step sizes:")
print(f"    n -> n/2:    delta_rapidity = -ln(2)/2 = {-ln2_half:.8f}")
print(f"    n -> 3n+1:   delta_rapidity ~ +ln(3)/2 = {ln3_half:.8f} (+ correction O(1/n))")
print(f"    Combined (3n+1)/2: delta ~ ln(3/2)/2 = {math.log(1.5)/2:.8f}")
print()

for start_n in [27, 97, 871, 6171]:
    traj = collatz_trajectory(start_n)
    print(f"  Collatz trajectory for n = {start_n}: {len(traj)} steps to reach 1")

    # Rapidity trajectory
    rap_traj = [rapidity(n) for n in traj]
    max_rap = max(rap_traj)
    max_val = max(traj)
    max_idx = traj.index(max_val)

    print(f"    Start rapidity:   {rap_traj[0]:.6f} (n={start_n})")
    print(f"    Max rapidity:     {max_rap:.6f} (n={max_val} at step {max_idx})")
    print(f"    End rapidity:     {rap_traj[-1]:.6f} (n=1)")
    print(f"    Total descent:    {rap_traj[0] - rap_traj[-1]:.6f}")
    print(f"    Max excursion:    {max_rap - rap_traj[0]:.6f} above start")

    # Count up vs down steps in rapidity
    up_steps = sum(1 for i in range(len(rap_traj) - 1) if rap_traj[i + 1] > rap_traj[i])
    down_steps = sum(1 for i in range(len(rap_traj) - 1) if rap_traj[i + 1] < rap_traj[i])
    avg_up = sum(rap_traj[i + 1] - rap_traj[i] for i in range(len(rap_traj) - 1)
                 if rap_traj[i + 1] > rap_traj[i]) / max(up_steps, 1)
    avg_down = sum(rap_traj[i] - rap_traj[i + 1] for i in range(len(rap_traj) - 1)
                   if rap_traj[i + 1] < rap_traj[i]) / max(down_steps, 1)

    print(f"    Up steps: {up_steps}, avg up = {avg_up:.6f}")
    print(f"    Down steps: {down_steps}, avg down = {avg_down:.6f}")
    print(f"    Ratio down/up steps: {down_steps/max(up_steps,1):.4f}")

    # Show first 40 steps with rapidity
    show_n = min(40, len(traj))
    print(f"    First {show_n} steps:")
    for i in range(show_n):
        step_type = ""
        if i > 0:
            delta = rap_traj[i] - rap_traj[i - 1]
            if traj[i - 1] % 2 == 0:
                step_type = f"  /2   delta={delta:+.6f}"
            else:
                step_type = f"  3n+1 delta={delta:+.6f}"
        print(f"      step {i:3d}: n={traj[i]:>8d}  rapidity={rap_traj[i]:>10.6f}{step_type}")
    if len(traj) > show_n:
        print(f"      ... ({len(traj) - show_n} more steps)")
    print()

# Collatz rapidity drift theorem
print("  COLLATZ RAPIDITY DRIFT THEOREM:")
print("  If the trajectory visits k odd numbers in N total steps,")
print("  then N-k are even steps (dividing by 2), giving:")
print("    net rapidity change = k * ln(3)/2 - (N-k) * ln(2)/2 + corrections")
print("  For the trajectory to reach 1 (rapidity 0 from rapidity ln(n)/2),")
print("  we need: k * ln(3) < (N-k) * ln(2) + ln(n)")
print("  i.e., k/N < ln(2)/ln(3) + ln(n)/(N*ln(3))")
print(f"  ln(2)/ln(3) = {math.log(2)/math.log(3):.10f}")
print("  So asymptotically, the fraction of odd steps must be < 0.6309...")
print()

# Verify for n=27
traj = collatz_trajectory(27)
odd_count = sum(1 for n in traj[:-1] if n % 2 == 1)
total = len(traj) - 1
print(f"  For n=27: {odd_count} odd steps out of {total} = {odd_count/total:.6f}")
print(f"    Compare ln(2)/ln(3) = {math.log(2)/math.log(3):.6f}")

# =============================================================================
# 5. PRIME GAPS IN RAPIDITY
# =============================================================================
print("\n" + "=" * 70)
print("5. PRIME GAPS IN RAPIDITY SPACE")
print("=" * 70)

def sieve_primes(limit):
    """Sieve of Eratosthenes."""
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]

primes = sieve_primes(100000)

print(f"\n  Generated {len(primes)} primes up to {primes[-1]}.")
print()

# Rapidity gaps: gap_n = rapidity(p_{n+1}) - rapidity(p_n) = ln(p_{n+1}/p_n)/2
# Cramer conjecture: max(p_{n+1} - p_n) ~ (ln p_n)^2
# In rapidity: max gap ~ (p_{n+1} - p_n) / (2*p_n) ~ (ln p_n)^2 / (2*p_n)
# But rapidity(p_n) = ln(p_n)/2, so ln(p_n) = 2*R where R = rapidity(p_n)
# Cramer in rapidity: max gap up to rapidity R ~ (2R)^2 / (2 * e^{2R}) = 2R^2 / e^{2R}
# This goes to 0! The rapidity gaps SHRINK exponentially.

# But the actual statement of Cramer is about ARITHMETIC gaps g_n = p_{n+1} - p_n.
# The rapidity gap is g_n / (2 * p_n) approximately.
# So max rapidity gap ~ (ln p_n)^2 / (2 * p_n) by Cramer.

print("  PRIME RAPIDITY GAPS: delta_n = ln(p_{n+1})/2 - ln(p_n)/2 = ln(p_{n+1}/p_n)/2")
print("  Cramer bound on arithmetic gap: g_n <= C * (ln p_n)^2")
print("  Rapidity version: delta_n ~ g_n / (2*p_n) <= C * (ln p_n)^2 / (2*p_n)")
print()

# Collect max rapidity gaps in windows
print(f"  {'Rapidity range':>20s}  {'Max rap gap':>12s}  {'At primes':>20s}  {'Arith gap':>10s}  {'Cramer bnd':>10s}")
print("  " + "-" * 80)

# Track running max rapidity gap
max_gap_records = []
current_max_rgap = 0

for i in range(len(primes) - 1):
    p1 = primes[i]
    p2 = primes[i + 1]
    rgap = rapidity(p2) - rapidity(p1)
    agap = p2 - p1

    if rgap > current_max_rgap:
        current_max_rgap = rgap
        cramer = (math.log(p1) ** 2) / (2 * p1)
        max_gap_records.append((p1, p2, rgap, agap, cramer))

# Show the record-setters
for p1, p2, rgap, agap, cramer in max_gap_records[:25]:
    r1 = rapidity(p1)
    print(f"  rap = {r1:>8.4f} ({p1:>6d})  {rgap:>12.8f}  ({p1:>8d}, {p2:>8d})  {agap:>10d}  {cramer:>10.6f}")

# Asymptotic behavior
print("\n  ASYMPTOTIC: rapidity gaps DECREASE as 1/p on average.")
print("  The rapidity prime 'density' INCREASES: d(rapidity)/dn ~ 1/(2n) for nth number,")
print("  and by PNT the average arithmetic gap at p is ln(p), so avg rapidity gap ~ ln(p)/(2p).")
print()

# Average rapidity gap in ranges
print("  Average rapidity gap by magnitude range:")
ranges = [(2, 100), (100, 1000), (1000, 10000), (10000, 100000)]
for lo, hi in ranges:
    gaps = []
    for i in range(len(primes) - 1):
        if lo <= primes[i] < hi:
            gaps.append(rapidity(primes[i + 1]) - rapidity(primes[i]))
    if gaps:
        avg = sum(gaps) / len(gaps)
        mx = max(gaps)
        # Predicted average: ln(p_mid)/(2*p_mid) where p_mid = sqrt(lo*hi)
        p_mid = math.sqrt(lo * hi)
        pred = math.log(p_mid) / (2 * p_mid)
        print(f"    primes in [{lo:>6d}, {hi:>6d}): avg gap = {avg:.8f}, "
              f"max = {mx:.8f}, predicted avg = {pred:.8f}, ratio = {avg/pred:.4f}")

# =============================================================================
# 6. GOLDBACH IN RAPIDITY
# =============================================================================
print("\n" + "=" * 70)
print("6. GOLDBACH DECOMPOSITIONS IN RAPIDITY SPACE")
print("=" * 70)

prime_set = set(primes)

print("\n  For even n, Goldbach decompositions: n = p + q, both prime.")
print("  In rapidity: (rapidity(p), rapidity(q)) = (ln(p)/2, ln(q)/2).")
print("  Note: rapidity(p) + rapidity(q) = ln(pq)/2, NOT ln(p+q)/2 = ln(n)/2.")
print("  So the points do NOT lie on the line x+y = rapidity(n).")
print("  Instead they satisfy: e^{2x} + e^{2y} = n (the EXPONENTIAL constraint).")
print()

for n in [10, 20, 30, 50, 100, 200, 1000]:
    decomps = []
    for p in primes:
        if p > n // 2:
            break
        q = n - p
        if q in prime_set:
            decomps.append((p, q))

    print(f"  n = {n}: {len(decomps)} Goldbach decompositions")

    if len(decomps) <= 12:
        for p, q in decomps:
            rp = rapidity(p)
            rq = rapidity(q)
            rpq = rapidity(p * q)
            rn = rapidity(n)
            print(f"    ({p:>4d}, {q:>4d}):  rap = ({rp:.6f}, {rq:.6f})  "
                  f"sum = {rp+rq:.6f}  vs rap(n) = {rn:.6f}  "
                  f"rap(pq) = {rpq:.6f}")
    else:
        # Show statistics
        rap_sums = [rapidity(p) + rapidity(q) for p, q in decomps]
        rap_products = [rapidity(p * q) for p, q in decomps]
        rn = rapidity(n)
        print(f"    rap(n) = {rn:.6f}")
        print(f"    rap sum range: [{min(rap_sums):.6f}, {max(rap_sums):.6f}]")
        print(f"    rap product range: [{min(rap_products):.6f}, {max(rap_products):.6f}]")
        print(f"    avg rap sum = {sum(rap_sums)/len(rap_sums):.6f}")
        # Show a few
        for p, q in decomps[:3]:
            rp = rapidity(p)
            rq = rapidity(q)
            print(f"    ({p:>4d}, {q:>4d}):  rap = ({rp:.6f}, {rq:.6f})")
        print(f"    ... and {len(decomps) - 3} more")
    print()

# The rapidity geometry of Goldbach
print("  GEOMETRIC INSIGHT:")
print("  The Goldbach constraint p + q = n in rapidity space is:")
print("    e^{2*phi_p} + e^{2*phi_q} = n")
print("  This is a CONVEX CURVE in (phi_p, phi_q) space.")
print("  At the symmetric point phi_p = phi_q: 2*e^{2*phi} = n, phi = ln(n/2)/2.")
print("  The decomposition (p,q) = (n/2, n/2) would have phi = ln(n/2)/2,")
print("  but n/2 must be prime, which is rare (only n=4).")
print()

# Rapidity spread of Goldbach decompositions
print("  RAPIDITY SPREAD of Goldbach decompositions:")
print("  (difference between most and least balanced decomposition)")
print()
for n in range(4, 102, 2):
    decomps = []
    for p in primes:
        if p > n // 2:
            break
        q = n - p
        if q in prime_set:
            decomps.append((p, q))
    if decomps:
        spreads = [abs(rapidity(p) - rapidity(q)) for p, q in decomps]
        min_spread = min(spreads)  # most balanced
        max_spread = max(spreads)  # least balanced
        most_balanced = decomps[spreads.index(min_spread)]
        if n <= 50 or n % 20 == 0:
            print(f"    n={n:4d}: {len(decomps):3d} decomps, "
                  f"spread [{min_spread:.6f}, {max_spread:.6f}], "
                  f"most balanced = {most_balanced}")

# =============================================================================
# 7. SQRT RAPIDITY AND PERFECT SQUARES
# =============================================================================
print("\n" + "=" * 70)
print("7. SQRT RAPIDITY AND PERFECT SQUARE DETECTION")
print("=" * 70)

print("\n  rapidity(n) = ln(n)/2")
print("  rapidity(n^2) = ln(n^2)/2 = ln(n) = 2 * rapidity(n)")
print("  rapidity(sqrt(n)) = ln(sqrt(n))/2 = ln(n)/4 = rapidity(n)/2")
print()
print("  PERFECT SQUARES have rapidity = 2 * rapidity(sqrt).")
print("  In other words, n is a perfect square iff rapidity(n) is in the set")
print("  {2*rapidity(k) : k >= 1} = {ln(k) : k >= 1}.")
print()
print("  Since rapidity(n) = ln(n)/2, this means n is a perfect square iff")
print("  there exists integer k with ln(n)/2 = ln(k), i.e., n = k^2.")
print("  Tautological? Not quite -- the RAPIDITY LATTICE tells us more.")

print("\n  THE RAPIDITY LATTICE of perfect k-th powers:")
print("  rapidity(n^k) = k * rapidity(n)")
print("  So perfect k-th powers form a SUBLATTICE of rapidity space")
print("  with spacing k times the fundamental.")
print()

print("  n     rapidity    square?  cube?  4th?    rapidity mod ln(2)   mod ln(3)")
print("  " + "-" * 80)
for n in range(1, 41):
    r = rapidity(n)
    is_sq = int(math.isqrt(n)) ** 2 == n
    is_cube = round(n ** (1/3)) ** 3 == n
    is_4th = round(n ** 0.25) ** 4 == n
    r_mod2 = r % math.log(2)
    r_mod3 = r % math.log(3)
    markers = ""
    if is_sq: markers += " SQ"
    if is_cube: markers += " CUBE"
    if is_4th: markers += " 4TH"
    print(f"  {n:3d}   {r:10.6f}   {'Y' if is_sq else '.':>3s}    {'Y' if is_cube else '.':>3s}    "
          f"{'Y' if is_4th else '.':>3s}       {r_mod2:8.6f}        {r_mod3:8.6f}{markers}")

# Rapidity characterization of multiplicative structure
print("\n  MULTIPLICATIVE STRUCTURE IN RAPIDITY:")
print("  The integers form a free abelian monoid under multiplication,")
print("  generated by the primes. In rapidity space, this becomes an")
print("  ADDITIVE monoid generated by {ln(p)/2 : p prime}.")
print()
print("  Perfect square <=> all rapidity components even")
print("  Perfect cube   <=> all rapidity components divisible by 3")
print("  k-th power     <=> all rapidity components divisible by k")
print()

# The rapidity address of special numbers in the project
print("  RAPIDITY ADDRESSES of project-special numbers:")
special = [
    (1, "transitive tournament"),
    (3, "3-cycle, curvature quantum"),
    (5, "n=4 max H"),
    (7, "FORBIDDEN"),
    (9, "3^2 = first square > 1 in our theory"),
    (21, "FORBIDDEN"),
    (42, "2*3*7, product of first forbidden factor set"),
    (45, "n=6 max H"),
    (147, "3*7^2, product of forbidden velocities"),
    (189, "3^3*7, Paley T_7 max H"),
    (661, "n=8 max H (prime!)"),
    (95095, "5*7*11*13*19, Paley T_11"),
]

print(f"  {'n':>8s}  {'rapidity':>10s}  {'2*rapidity':>10s}  {'in lattice?':>12s}  description")
print("  " + "-" * 80)
for n, desc in special:
    r = rapidity(n)
    r2 = 2 * r
    # Check if rapidity is a "nice" multiple of ln(2)/2 or ln(3)/2
    r_in_ln2 = r / (math.log(2) / 2)
    r_in_ln3 = r / (math.log(3) / 2)
    lattice = ""
    if abs(r_in_ln2 - round(r_in_ln2)) < 0.001:
        lattice = f"2^{round(r_in_ln2)}"
    elif abs(r_in_ln3 - round(r_in_ln3)) < 0.001:
        lattice = f"3^{round(r_in_ln3)}"
    print(f"  {n:8d}  {r:10.6f}  {r2:10.6f}  {lattice:>12s}  {desc}")

# DEEPER: Rapidity distance between all pairs of special numbers
print("\n  RAPIDITY DISTANCES between forbidden numbers and key values:")
print("  (rapidity distance = ln(a/b)/2)")
print()
key_vals = [3, 5, 7, 9, 21, 45, 189, 661, 95095]
print(f"  {'':>8s}", end="")
for b in key_vals:
    print(f"  {b:>8d}", end="")
print()
for a in key_vals:
    print(f"  {a:>8d}", end="")
    for b in key_vals:
        d = abs(rapidity(a) - rapidity(b))
        print(f"  {d:>8.4f}", end="")
    print()

# Check which distances are nice multiples of ln(p)/2
print("\n  Notable rapidity distances that are exact ln(p)/2 multiples:")
for i, a in enumerate(key_vals):
    for j, b in enumerate(key_vals):
        if i < j:
            d = abs(rapidity(a) - rapidity(b))
            ratio = a / b if a > b else b / a
            # Check if ratio is a nice fraction
            if ratio == int(ratio):
                print(f"    |rap({a}) - rap({b})| = {d:.8f} = ln({int(ratio)})/2 = rapidity({int(ratio)})")
            else:
                # Check if ratio is p/q for small p,q
                fr = Fraction(a, b) if a > b else Fraction(b, a)
                if fr.denominator <= 20:
                    print(f"    |rap({a}) - rap({b})| = {d:.8f} = ln({fr})/2")

# =============================================================================
# SYNTHESIS
# =============================================================================
print("\n" + "=" * 70)
print("SYNTHESIS: RAPIDITY AS NUMBER-THEORETIC COORDINATE")
print("=" * 70)

print("""
  1. CONTINUED FRACTIONS: Convergent rapidities form a damped oscillation
     around the true rapidity. The damping rate is governed by the CF
     partial quotients: phi (all 1s) has geometric damping phi^{-2n},
     sqrt(2) has damping (1+sqrt(2))^{-2n}, e has 3-periodic modulation.

  2. QUADRATIC RESIDUES: QRs and NQRs have nearly equal average rapidity
     (splitting of ln((p-1)!)/2 into two halves). The small bias comes from
     the correlation between being a QR and being small. For Paley primes,
     this bias encodes the tournament structure.

  3. BERNOULLI NUMBERS: Power sums S_m(n) = sum k^m connect to Bernoulli
     numbers via Faulhaber's formula. The forbidden 7 = S_1(3)+1 and
     21 = S_1(6) = T_6 are triangular-number adjacent. The spectral zeta
     zeta_M(-s) = sum lambda_i^{-s} is expressible via Newton-Bernoulli
     when the eigenvalues come from a polynomial char equation.

  4. COLLATZ: The Collatz map in rapidity space is a biased random walk
     with drift -ln(2)/2 (even step) or +ln(3)/2 (odd step). The critical
     ratio ln(2)/ln(3) = 0.6309... governs the fraction of odd steps.
     For n=27 (111 steps), the rapidity trajectory shows dramatic excursions
     before the final descent to 0.

  5. PRIME GAPS: In rapidity space, prime gaps DECREASE as O(ln(p)/(2p)).
     Cramer's conjecture translates to: max rapidity gap at rapidity R is
     O(R^2 * e^{-2R}), which vanishes doubly-exponentially. The primes
     become DENSE in rapidity space. Record rapidity gaps decrease monotonically
     in the long run.

  6. GOLDBACH: The Goldbach constraint p+q=n becomes e^{2x}+e^{2y}=n in
     rapidity space -- a convex curve, not a line. The most balanced
     decomposition (closest rapidities) corresponds to p,q near n/2.
     The rapidity spread measures how "asymmetric" the decomposition is.

  7. SQRT RAPIDITY: Perfect k-th powers form sublattices of rapidity space
     with spacing k times the prime rapidity generators. The multiplicative
     structure of Z becomes purely additive in rapidity. The forbidden numbers
     7 and 21 are NOT in any power sublattice (both have rapidity that is
     not a simple multiple of any prime rapidity).

  OVERARCHING THEME: Rapidity = ln(n)/2 converts multiplicative number theory
  into additive geometry. The fundamental theorem of arithmetic becomes:
  every rapidity is a unique sum of prime rapidities. This is the LOG of
  the integers, but normalized by 1/2 to match the Cayley/Lorentz structure.
""")

print("=" * 70)
print("END OF RAPIDITY NUMBER THEORY EXPLORATIONS")
print("=" * 70)
