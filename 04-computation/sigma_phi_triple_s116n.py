#!/usr/bin/env python3
"""sigma_phi_triple_s116n.py — The three solutions to sigma(n) = 2n + phi(n).

NOT unique at 42! There are (at least) three: n = 12, 42, 1242.
What's the pattern? Are there more?

Session: kind-pasteur-2026-03-16-S116n32
"""
from math import gcd, isqrt
from fractions import Fraction

print()
print("  THE THREE SOLUTIONS: sigma(n) = 2n + phi(n)")
print()
print("=" * 70)
print()

# Efficient sieve
def compute_sigma_phi(N):
    sig = [0] * (N+1)
    ph = list(range(N+1))
    for d in range(1, N+1):
        for m in range(d, N+1, d):
            sig[m] += d
    is_prime = [True] * (N+1)
    for p in range(2, N+1):
        if is_prime[p]:
            for m in range(p*p, N+1, p):
                is_prime[m] = False
            for m in range(p, N+1, p):
                ph[m] = ph[m] // p * (p-1)
    return sig, ph

print("  Computing sigma and phi up to 10,000,000...")
N = 10_000_000
sig, ph = compute_sigma_phi(N)

solutions = []
for n in range(1, N+1):
    if sig[n] == 2*n + ph[n]:
        solutions.append(n)

print(f"\n  Solutions to sigma(n) = 2n + phi(n) for n <= {N}:")
print()

def factorize(n):
    factors = {}
    temp = n
    for p in range(2, isqrt(temp)+2):
        while temp % p == 0:
            factors[p] = factors.get(p, 0) + 1
            temp //= p
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

for n in solutions:
    f = factorize(n)
    fstr = ' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f.items()))
    print(f"  n = {n:>10d} = {fstr}")
    print(f"    sigma = {sig[n]}, phi = {ph[n]}, 2n = {2*n}")
    print(f"    abundance = sigma - 2n = {sig[n]-2*n}")
    print(f"    phi(n) = {ph[n]}")
    print(f"    Match: {sig[n]-2*n == ph[n]}")
    # Additional properties
    print(f"    n/6 = {n/6}")
    print(f"    sigma/n = {Fraction(sig[n], n)} = {sig[n]/n:.6f}")
    print(f"    phi/n = {Fraction(ph[n], n)} = {ph[n]/n:.6f}")
    print(f"    abundance/n = {Fraction(sig[n]-2*n, n)} = {(sig[n]-2*n)/n:.6f}")
    print()

print(f"  Total solutions up to {N}: {len(solutions)}")
print()

# ============================================================
print("  STRUCTURAL ANALYSIS")
print("  " + "-" * 50)
print()

# All three are divisible by 6
print("  All solutions divisible by 6?")
for n in solutions:
    print(f"    {n} / 6 = {n/6}, mod 6 = {n%6}")
print()

# The equation: sigma(n) = 2n + phi(n)
# Rewrite: sigma(n)/n = 2 + phi(n)/n
# abundancy index = 2 + totient ratio
# For n = prod p_i^{a_i}:
#   sigma(n)/n = prod (p_i^{a_i+1}-1)/((p_i-1)*p_i^{a_i})
#   phi(n)/n = prod (1 - 1/p_i)
#   Condition: prod (p^{a+1}-1)/((p-1)*p^a) = 2 + prod(1-1/p)

print("  The equation in terms of local factors:")
print("  Let f(p,a) = (p^{a+1}-1)/((p-1)*p^a)  [local sigma/n factor]")
print("  Let g(p) = (p-1)/p = 1-1/p             [local phi/n factor]")
print("  Then: prod f(p_i,a_i) = 2 + prod g(p_i)")
print()

for n in solutions:
    f = factorize(n)
    local_sig = []
    local_phi = []
    for p, a in sorted(f.items()):
        fs = Fraction(p**(a+1) - 1, (p-1) * p**a)
        gp = Fraction(p-1, p)
        local_sig.append((p, a, fs))
        local_phi.append((p, gp))
        print(f"    p={p}, a={a}: f(p,a)={fs} = {float(fs):.6f}, g(p)={gp}")

    prod_f = Fraction(1)
    for _, _, fs in local_sig:
        prod_f *= fs
    prod_g = Fraction(1)
    for _, gp in local_phi:
        prod_g *= gp

    print(f"    prod f = {prod_f} = {float(prod_f):.6f}")
    print(f"    prod g = {prod_g} = {float(prod_g):.6f}")
    print(f"    2 + prod g = {2 + prod_g} = {float(2 + prod_g):.6f}")
    print(f"    Equal? {prod_f == 2 + prod_g}")
    print()

# ============================================================
print("  DEEPER PATTERN: THE ROLE OF 6")
print("  " + "-" * 50)
print()

# sigma(n) - 2n - phi(n) = 0
# For n = 6m: what are the constraints on m?
print("  n = 12 = 6 * 2")
print("  n = 42 = 6 * 7")
print("  n = 1242 = 6 * 207 = 6 * 9 * 23")
print()

# Can we find solutions of the form n = 2 * 3^b * q for prime q?
print("  Trying n = 2 * 3^b * q for small b and primes q:")
print()

def sigma_formula(n):
    """Compute sigma(n) from factorization."""
    f = factorize(n)
    result = 1
    for p, a in f.items():
        result *= (p**(a+1) - 1) // (p - 1)
    return result

def phi_formula(n):
    """Compute phi(n) from factorization."""
    f = factorize(n)
    result = n
    for p in f:
        result = result // p * (p - 1)
    return result

# Sieve of primes
primes = []
sieve = [True] * 100001
for p in range(2, 100001):
    if sieve[p]:
        primes.append(p)
        for m in range(p*p, 100001, p*p):
            sieve[m] = False

# For n = 2 * 3^b * q (q prime, q != 2, 3):
# sigma = (1+2) * (3^{b+1}-1)/2 * (1+q) = 3 * (3^{b+1}-1)/2 * (1+q)
# phi = 1 * 2 * 3^{b-1} * (q-1) = 2 * 3^{b-1} * (q-1)
# 2n = 4 * 3^b * q

for b in range(1, 8):
    print(f"  b={b}: n = 2 * 3^{b} * q")
    # sigma = 3 * (3^{b+1}-1)/2 * (1+q)
    # phi = 2 * 3^{b-1} * (q-1)
    # 2n = 4 * 3^b * q
    # Equation: 3*(3^{b+1}-1)/2*(1+q) = 4*3^b*q + 2*3^{b-1}*(q-1)
    # = 4*3^b*q + 2*3^{b-1}*q - 2*3^{b-1}
    # = q*(4*3^b + 2*3^{b-1}) - 2*3^{b-1}
    # = q*3^{b-1}*(12+2) - 2*3^{b-1}
    # = q*14*3^{b-1} - 2*3^{b-1}
    # = 3^{b-1}*(14q - 2)
    # Hmm wait, let me be more careful.
    # 4*3^b + 2*3^{b-1} = 3^{b-1}*(4*3 + 2) = 3^{b-1}*14
    # So RHS = 3^{b-1}*(14q - 2)
    # LHS = 3*(3^{b+1}-1)/2*(1+q)
    # Let S = 3^{b+1} - 1 = 3*3^b - 1
    # LHS = 3*S*(1+q)/2
    # Equation: 3*S*(1+q)/2 = 3^{b-1}*(14q - 2)
    # 3*(3*3^b-1)*(1+q)/2 = 3^{b-1}*(14q-2)
    # Multiply both sides by 2:
    # 3*(3*3^b-1)*(1+q) = 2*3^{b-1}*(14q-2)
    # 3*(3^{b+1}-1)*(1+q) = 2*3^{b-1}*(14q-2)
    # (3^{b+2}-3)*(1+q) = 2*3^{b-1}*(14q-2)
    # 3*(3^{b+1}-1)*(1+q) = 4*3^{b-1}*(7q-1)
    # Divide by 3^{b-1}... hmm this is getting messy. Let me just compute.

    found = False
    for q in primes:
        if q <= 3:
            continue
        n = 2 * 3**b * q
        if n > N:
            break
        if sig[n] == 2*n + ph[n]:
            print(f"    SOLUTION: q={q}, n={n}")
            found = True
    if not found:
        print(f"    No solutions up to n={N}")
    print()

# ============================================================
print("  GENERAL FORM SEARCH")
print("  " + "-" * 50)
print()

# What about n = 2^a * 3^b * q for various a, b?
print("  n = 2^a * 3^b * q (q prime > 3):")
for a in range(1, 6):
    for b_val in range(1, 6):
        found_sols = []
        for q in primes:
            if q <= 3:
                continue
            n = 2**a * 3**b_val * q
            if n > N:
                break
            if sig[n] == 2*n + ph[n]:
                found_sols.append((q, n))
        if found_sols:
            for q, n in found_sols:
                print(f"    a={a}, b={b_val}, q={q}: n={n}")

# n = 2^a * 3 * 7 (generalizing 42 = 2 * 3 * 7)
print()
print("  n = 2^a * 3 * 7:")
for a in range(1, 20):
    n = 2**a * 3 * 7
    if n > N:
        break
    D = sig[n] - 2*n - ph[n]
    print(f"    a={a}: n={n}, D={D}")

# n = 2 * 3 * p (generalizing 42 = 2 * 3 * 7)
print()
print("  n = 2 * 3 * p (squarefree):")
for p in primes[:20]:
    if p <= 3:
        continue
    n = 6 * p
    D = sig[n] - 2*n - ph[n]
    print(f"    p={p}: n={n}, D={D}")

print()

# ============================================================
print("  THE n=12 SOLUTION")
print("  " + "-" * 50)
print()
print("  n = 12 = 2^2 * 3")
print("  sigma(12) = 1+2+3+4+6+12 = 28")
print("  phi(12) = |{1,5,7,11}| = 4")
print("  2*12 = 24")
print("  28 = 24 + 4. YES.")
print()
print("  12 = phi(42)/1 = the totient of 42!")
print("  sigma(phi(42)) = 2*phi(42) + phi(phi(42))")
print(f"  phi(phi(42)) = phi(12) = {ph[12]}")
print(f"  sigma(12) = {sig[12]} = 2*12 + 4 = 28. Confirmed.")
print()

# ============================================================
print("  THE n=1242 SOLUTION")
print("  " + "-" * 50)
print()
print("  n = 1242 = 2 * 3^3 * 23")
print(f"  sigma(1242) = {sig[1242]}")
print(f"  phi(1242) = {ph[1242]}")
print(f"  2*1242 = {2*1242}")
print(f"  {sig[1242]} = {2*1242} + {ph[1242]}? {sig[1242] == 2*1242 + ph[1242]}")
print()

# Connection to 42?
print("  1242 = 42 * 29 + 24 = 42 * 29.571...")
print(f"  1242 / 42 = {1242/42}")
print(f"  1242 / 6 = {1242/6} = 207 = 9 * 23")
print(f"  1242 = 6 * 207 = 6 * 9 * 23 = 2 * 27 * 23")
print()
print(f"  23 - 1 = 22 = 2 * 11")
print(f"  23 is prime. 23 = 24 - 1 = 4! - 1")
print()

# Relationship: 42 = 2*3*7, 1242 = 2*27*23 = 2*3^3*23
# The 3 got cubed, and 7 got replaced by 23
# 23 = 3*7 + 2? No, 3*7+2 = 23. YES!
print(f"  23 = 3*7 + 2. So 1242 = 2 * 3^3 * (3*7+2)")
print(f"  Starting from 42 = 2*3*7: cube the 3, shift the 7 to 3*7+2=23.")
print()

# Are there solutions with more factors?
print("  Checking n with 4+ prime factors (with multiplicity):")
for n in solutions:
    f = factorize(n)
    omega = sum(f.values())  # total prime factors with multiplicity
    print(f"    n={n}: Omega(n) = {omega}")
print()

# ============================================================
print("  OEIS CHECK: sequence 12, 42, 1242, ...")
print("  " + "-" * 50)
print()
print("  This sequence should be checked in OEIS.")
print(f"  Known terms: {solutions}")
print(f"  Is this A??????")
print()

# The ratio pattern
if len(solutions) >= 2:
    for i in range(1, len(solutions)):
        r = Fraction(solutions[i], solutions[i-1])
        print(f"  {solutions[i]}/{solutions[i-1]} = {r} = {float(r):.6f}")
print()

# ============================================================
print("  EXISTENCE PROOF FOR MORE SOLUTIONS?")
print("  " + "-" * 50)
print()

# For n = 2 * 3^b * q with q prime:
# sigma = 3 * (3^{b+1}-1)/2 * (1+q) = 3(3^{b+1}-1)(1+q)/2
# phi = 2 * 3^{b-1} * (q-1)  (for b >= 1)
# 2n = 4 * 3^b * q
# Equation: 3(3^{b+1}-1)(1+q)/2 = 4*3^b*q + 2*3^{b-1}*(q-1)
# Multiply by 2: 3(3^{b+1}-1)(1+q) = 8*3^b*q + 4*3^{b-1}*(q-1)
# = 3^{b-1}*(8*9*q + 4q - 4)... wait, 8*3^b = 8*3*3^{b-1} = 24*3^{b-1}
# So RHS = 3^{b-1}*(24q + 4q - 4) = 3^{b-1}*(28q - 4)
# LHS = 3*(3^{b+1}-1)*(1+q) = 3*(3*3^b - 1)*(1+q)
# = (3^{b+2} - 3)*(1+q) = (9*3^b - 3)*(1+q) = 3*(3^{b+1}-1)*(1+q)
# Let me try differently. Divide by 3^{b-1}:
# LHS/3^{b-1} = 3*(3^{b+1}-1)*(1+q)/3^{b-1} = 3*(1+q)*(3^{b+1}-1)/3^{b-1}
# = 3*(1+q)*(3^2 * 3^{b-1} - 1)/3^{b-1}
# = 3*(1+q)*(9 - 3^{1-b})
# This only works cleanly for b >= 1.
# For b=1: LHS/1 = 3*(3^2-1)*(1+q) = 3*8*(1+q) = 24(1+q)
#           RHS/1 = 28q - 4
#           24+24q = 28q-4 => 28 = 4q => q = 7. GIVES n = 42!
#
# For b=2: LHS/3 = 3*(3^3-1)*(1+q) = 3*26*(1+q) = 78(1+q)
#           RHS/3 = 28q - 4... wait, let me redo.
# n = 2*3^2*q = 18q. sigma = 3*(3^3-1)/2*(1+q) = 3*26/2*(1+q) = 39(1+q)
# phi = 2*3*(q-1) = 6(q-1). 2n = 36q.
# 39(1+q) = 36q + 6(q-1) = 42q - 6
# 39 + 39q = 42q - 6 => 45 = 3q => q = 15. NOT PRIME.

# For b=3: n = 2*3^3*q = 54q. sigma = 3*(3^4-1)/2*(1+q) = 3*80/2*(1+q) = 120(1+q)
# phi = 2*3^2*(q-1) = 18(q-1). 2n = 108q.
# 120(1+q) = 108q + 18(q-1) = 126q - 18
# 120 + 120q = 126q - 18 => 138 = 6q => q = 23. PRIME! GIVES n = 1242!

# For b=4: n = 2*3^4*q = 162q. sigma = 3*(3^5-1)/2*(1+q) = 3*242/2*(1+q) = 363(1+q)
# phi = 2*3^3*(q-1) = 54(q-1). 2n = 324q.
# 363(1+q) = 324q + 54(q-1) = 378q - 54
# 363 + 363q = 378q - 54 => 417 = 15q => q = 27.8. NOT INTEGER.

# For b=5: n = 2*3^5*q = 486q. sigma = 3*(3^6-1)/2*(1+q) = 3*728/2*(1+q) = 1092(1+q)
# phi = 2*3^4*(q-1) = 162(q-1). 2n = 972q.
# 1092(1+q) = 972q + 162(q-1) = 1134q - 162
# 1092 + 1092q = 1134q - 162 => 1254 = 42q => q = 1254/42 = 29.857... NOT INTEGER.

# For b=6: n = 2*3^6*q = 1458q. sigma = 3*(3^7-1)/2*(1+q) = 3*2186/2*(1+q) = 3279(1+q)
# phi = 2*3^5*(q-1) = 486(q-1). 2n = 2916q.
# 3279(1+q) = 2916q + 486(q-1) = 3402q - 486
# 3279 + 3279q = 3402q - 486 => 3765 = 123q => q = 3765/123 = 30.609... NOT INTEGER.

print("  Family n = 2 * 3^b * q (q prime, q > 3):")
print("  Equation reduces to:")
print()

for b in range(1, 15):
    # sigma = 3*(3^{b+1}-1)/2 * (1+q)
    # phi = 2*3^{b-1}*(q-1)
    # 2n = 4*3^b*q
    # Equation: 3*(3^{b+1}-1)/2 * (1+q) = 4*3^b*q + 2*3^{b-1}*(q-1)

    A = 3 * (3**(b+1) - 1)  # LHS without /2 and (1+q)
    B = 8 * 3**b + 4 * 3**(b-1)  # RHS without q terms and constant, times 2
    C = -4 * 3**(b-1)  # RHS constant part times 2

    # 2: A*(1+q) = B*q + C
    # A + A*q = B*q + C
    # q*(A - B) = C - A
    # q = (C - A)/(A - B)

    numerator = C - A
    denominator = A - B

    if denominator != 0:
        q_val = Fraction(numerator, denominator)
        q_int = numerator // denominator if numerator % denominator == 0 else None

        is_prime_q = False
        if q_int and q_int > 3:
            is_prime_q = all(q_int % d != 0 for d in range(2, isqrt(q_int)+1)) if q_int > 1 else False

        n_val = 2 * 3**b * q_int if q_int else None
        status = ""
        if q_int is None:
            status = f"q = {q_val} (NOT integer)"
        elif q_int <= 3:
            status = f"q = {q_int} (too small)"
        elif not is_prime_q:
            status = f"q = {q_int} (NOT prime)"
        else:
            status = f"q = {q_int} (PRIME!) => n = {n_val}"

        print(f"  b={b:2d}: {status}")

print()

# Now try n = 2^a * 3 * q
print("  Family n = 2^a * 3 * q (q prime, q > 3):")
for a in range(1, 15):
    # sigma = (2^{a+1}-1) * (1+3) * (1+q) = (2^{a+1}-1) * 4 * (1+q)
    # phi = 2^{a-1} * 2 * (q-1) = 2^a * (q-1)
    # 2n = 2^{a+1} * 3 * q
    # Equation: 4*(2^{a+1}-1)*(1+q) = 2^{a+1}*3*q + 2^a*(q-1)

    # Let's multiply through
    # 4*(2^{a+1}-1)*(1+q) = 6*2^a*q + 2^a*q - 2^a
    # = 2^a*(7q - 1)
    # LHS = 4*(2^{a+1}-1)*(1+q)
    # So: 4*(2^{a+1}-1)*(1+q) = 2^a*(7q-1)
    # (2^{a+3}-4)*(1+q) = 2^a*(7q-1)
    # Hmm, let me just solve numerically.

    # sigma(2^a) = 2^{a+1}-1, sigma(3) = 4, sigma(q) = 1+q
    # sigma(n) = (2^{a+1}-1)*4*(1+q)
    # phi(n) = 2^{a-1}*2*(q-1) = 2^a*(q-1)
    # 2n = 2^{a+1}*3*q = 6*2^a*q

    # Eq: 4*(2^{a+1}-1)*(1+q) = 6*2^a*q + 2^a*(q-1)
    # = 2^a*(6q + q - 1) = 2^a*(7q - 1)
    # 4*(2^{a+1}-1)*(1+q) = 2^a*(7q-1)
    # Let u = 2^a
    # 4*(2u-1)*(1+q) = u*(7q-1)
    # (8u-4)*(1+q) = u*(7q-1)
    # 8u + 8uq - 4 - 4q = 7uq - u
    # 8u + 8uq - 4 - 4q - 7uq + u = 0
    # 9u + uq - 4 - 4q = 0
    # u(9+q) = 4(1+q)
    # u = 4(1+q)/(9+q)
    # 2^a = 4(1+q)/(9+q)

    # For this to give 2^a (a power of 2):
    # 4(1+q)/(9+q) must be a power of 2
    # (1+q)/(9+q) must be a power of 2 divided by 4 = 2^{a-2}

    # For a=1: 2 = 4(1+q)/(9+q) => 2(9+q) = 4+4q => 18+2q = 4+4q => 14 = 2q => q=7
    #   n = 2*3*7 = 42. CONFIRMED!
    # For a=2: 4 = 4(1+q)/(9+q) => 9+q = 1+q => 9=1. IMPOSSIBLE.
    # For a=3: 8 = 4(1+q)/(9+q) => 2(9+q) = 1+q => 18+2q = 1+q => q = -17. IMPOSSIBLE.

    target_u = 2**a
    # u*(9+q) = 4(1+q) => q*(u-4) = 4-9u => q = (4-9u)/(u-4)
    if a == 2:
        # u=4: denominator = 0. No solution.
        print(f"  a={a}: denominator zero, no solution")
        continue
    u = 2**a
    q_num = 4 - 9*u
    q_den = u - 4
    if q_den != 0:
        q_val = Fraction(q_num, q_den)
        q_int = q_num // q_den if q_num % q_den == 0 else None
        is_prime_q = False
        if q_int and q_int > 3:
            is_prime_q = all(q_int % d != 0 for d in range(2, isqrt(abs(q_int))+1)) if q_int > 1 else False
        n_val = 2**a * 3 * q_int if q_int and q_int > 0 else None

        if q_int is None or q_int <= 0:
            status = f"q = {float(q_val):.2f} (invalid)"
        elif q_int <= 3:
            status = f"q = {q_int} (too small)"
        elif not is_prime_q:
            status = f"q = {q_int} (NOT prime, {factorize(q_int)})"
        else:
            status = f"q = {q_int} (PRIME!) => n = {n_val}"
        print(f"  a={a:2d}: {status}")

print()

# Try n = 2^2 * 3^b
print("  Family n = 2^2 * 3^b:")
for b in range(1, 15):
    n = 4 * 3**b
    if n > N:
        break
    D = sig[n] - 2*n - ph[n]
    print(f"  b={b}: n={n}, D={D}")

print()

# ============================================================
print("  THE THREE SOLUTIONS: UNIFIED VIEW")
print("  " + "-" * 50)
print()
print("  n = 12 = 2^2 * 3     (family 2^a * 3^b)")
print("  n = 42 = 2 * 3 * 7   (family 2 * 3 * q, q=7)")
print("  n = 1242 = 2 * 3^3 * 23  (family 2 * 3^b * q, b=3, q=23)")
print()
print("  From the family n = 2 * 3^b * q:")
print("  b=1: q=7 (prime) => n=42")
print("  b=3: q=23 (prime) => n=1242")
print("  Other b: q not prime or not integer.")
print()
print("  From the family n = 2^a * 3:")
print("  Only a=2 could work (n=12), and it does.")
print()
print("  CONJECTURE: 12, 42, 1242 are the ONLY three solutions.")
print(f"  Verified up to n = {N}.")
print()

# Check the special structure of 12
print("  Why does 12 = 2^2 * 3 work?")
print(f"  sigma(12) = (1+2+4)(1+3) = 7*4 = 28")
print(f"  phi(12) = 2*2 = 4")
print(f"  2*12 = 24")
print(f"  28 - 24 = 4 = phi(12). Checks out.")
print()
print(f"  12 = phi(42) and sigma(12) = 2*12 + phi(12).")
print(f"  The TOTIENT of 42 is itself a solution!")
print(f"  This creates a self-referential loop:")
print(f"  42 satisfies the equation, and phi(42) = 12 also satisfies it.")
print()

# Does phi(12) = 4 satisfy?
print(f"  phi(12) = 4. Does 4 satisfy?")
print(f"  sigma(4) = 7, 2*4 + phi(4) = 8 + 2 = 10. 7 != 10. NO.")
print()
print(f"  phi(1242) = {ph[1242]}. Does {ph[1242]} satisfy?")
n = ph[1242]
if n <= N:
    print(f"  sigma({n}) = {sig[n]}, 2*{n} + phi({n}) = {2*n + ph[n]}. "
          f"Equal? {sig[n] == 2*n + ph[n]}")
print()

# The chain: 42 -> phi(42) = 12 (solution) -> phi(12) = 4 (not solution)
# Is 42 = sigma(12) - phi(12)? sigma(12) = 28, phi(12) = 4. 28-4 = 24 ≠ 42.
# But 42 = sigma(12) + 2*phi(12) + 2 = 28 + 8 + 6? No.
# 42 = 2*12 + 18? 42 = 2*12 + 18? No pattern.
# Actually: n * (sigma(n)/n - 2 - phi(n)/n) = 0
# So sigma/n - 2 = phi/n for all three.
# abundancy - 2 = totient ratio

print("  abundancy(n) - 2 = phi(n)/n:")
for n in solutions:
    ab = Fraction(sig[n], n)
    tr = Fraction(ph[n], n)
    print(f"    n={n}: abundancy = {ab} = {float(ab):.6f}, "
          f"phi/n = {tr} = {float(tr):.6f}, diff = {ab - 2} = {float(ab-2):.6f}")
print()
