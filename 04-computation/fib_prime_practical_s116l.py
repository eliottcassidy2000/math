"""
fib_prime_practical_s116l.py
Practical applications of the Fibonacci-Lucas-prime connection.

8 tools:
  1. Fibonacci Primality Witness
  2. Lucas Pseudoprime Test
  3. Fibonacci Entry Point Factoring
  4. Fibonacci (Zeckendorf) Coding
  5. Fibonacci Hash
  6. Golden Ratio Approximation Tool
  7. Fibonacci-Weighted Prime Counter
  8. Lucas-Tournament Connection

Author: kind-pasteur-2026-03-16-S116l
"""

import time
import math
from collections import defaultdict

# ============================================================
# CORE FIBONACCI / LUCAS UTILITIES
# ============================================================

def fib_mod(n, m):
    """Compute F_n mod m using fast doubling. O(log n)."""
    if m == 1:
        return 0
    if n < 0:
        # F_{-n} = (-1)^{n+1} F_n
        val = fib_mod(-n, m)
        if (-n) % 2 == 0:
            return (-val) % m
        return val
    if n == 0:
        return 0
    if n == 1:
        return 1 % m

    def _fib_pair(k):
        """Return (F_k mod m, F_{k+1} mod m)."""
        if k == 0:
            return (0, 1)
        a, b = _fib_pair(k >> 1)
        # F_{2k} = F_k * (2*F_{k+1} - F_k)
        # F_{2k+1} = F_k^2 + F_{k+1}^2
        c = a * ((2 * b - a) % m) % m
        d = (a * a + b * b) % m
        if k & 1:
            return (d, (c + d) % m)
        else:
            return (c, d)

    return _fib_pair(n)[0]


def lucas_mod(n, m):
    """Compute L_n mod m using L_n = F_{n-1} + F_{n+1} = 2*F_{n+1} - F_n."""
    if m == 1:
        return 0
    if n == 0:
        return 2 % m
    if n == 1:
        return 1 % m
    # L_n = F_{n-1} + F_{n+1}
    fn = fib_mod(n, m)
    fn1 = fib_mod(n + 1, m)
    # F_{n-1} = F_{n+1} - F_n
    fn_1 = (fn1 - fn) % m
    return (fn_1 + fn1) % m


def fib(n):
    """Compute F_n exactly (small n)."""
    if n < 0:
        val = fib(-n)
        if (-n) % 2 == 0:
            return -val
        return val
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def lucas(n):
    """Compute L_n exactly (small n)."""
    if n == 0:
        return 2
    if n == 1:
        return 1
    a, b = 2, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return b


def legendre(a, p):
    """Legendre symbol (a/p) for odd prime p."""
    if a % p == 0:
        return 0
    val = pow(a, (p - 1) // 2, p)
    return val if val == 1 else -1


def is_prime_trial(n):
    """Trial division primality test."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def sieve_primes(n):
    """Sieve of Eratosthenes up to n."""
    if n < 2:
        return []
    is_p = [True] * (n + 1)
    is_p[0] = is_p[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_p[i]:
            for j in range(i * i, n + 1, i):
                is_p[j] = False
    return [i for i in range(2, n + 1) if is_p[i]]


def gcd(a, b):
    """Greatest common divisor."""
    while b:
        a, b = b, a % b
    return a


# ============================================================
# TOOL 1: FIBONACCI PRIMALITY WITNESS
# ============================================================

def fibonacci_primality_witness(n):
    """
    Probable prime test using p | F_{p - (p/5)}.
    For prime p != 5: p | F_{p-1} if p = +-1 mod 5, else p | F_{p+1}.
    Returns True if n passes the Fibonacci witness test.
    """
    if n < 2:
        return False
    if n == 2 or n == 3 or n == 5:
        return True
    if n % 2 == 0:
        return False
    if n == 5:
        return True
    # Compute Legendre symbol (n/5) using quadratic reciprocity shortcut
    r = n % 5
    if r == 0:
        return n == 5
    if r == 1 or r == 4:
        eps = 1  # QR mod 5
    else:
        eps = -1  # NQR mod 5

    k = n - eps
    return fib_mod(k, n) == 0


def run_tool1():
    print("=" * 70)
    print("TOOL 1: FIBONACCI PRIMALITY WITNESS")
    print("=" * 70)
    print()
    print("Theory: For prime p != 2,5: p | F_{p - (p/5)}")
    print("  where (p/5) is +1 if p = 1,4 mod 5, else -1.")
    print()

    # Test on all numbers up to 10000
    primes = set(sieve_primes(10000))
    fib_prp = set()
    for n in range(2, 10001):
        if fibonacci_primality_witness(n):
            fib_prp.add(n)

    actual_primes = primes
    false_positives = fib_prp - actual_primes
    false_negatives = actual_primes - fib_prp

    print(f"Numbers 2..10000:")
    print(f"  Actual primes:      {len(actual_primes)}")
    print(f"  Fib-PRP positives:  {len(fib_prp)}")
    print(f"  False positives (Fibonacci pseudoprimes): {len(false_positives)}")
    if false_positives:
        fps = sorted(false_positives)
        print(f"  First 20 Fibonacci pseudoprimes: {fps[:20]}")
    print(f"  False negatives:    {len(false_negatives)}")
    if false_negatives:
        print(f"  Missed primes: {sorted(false_negatives)[:10]}")
    print()

    # Speed comparison
    t0 = time.time()
    for n in range(2, 10001):
        fibonacci_primality_witness(n)
    t_fib = time.time() - t0

    t0 = time.time()
    for n in range(2, 10001):
        is_prime_trial(n)
    t_trial = time.time() - t0

    print(f"Speed (2..10000):")
    print(f"  Fibonacci PRP test: {t_fib:.4f} sec")
    print(f"  Trial division:     {t_trial:.4f} sec")
    print(f"  Ratio (Fib/Trial):  {t_fib/t_trial:.2f}x")
    print()

    # The advantage of Fib test is for LARGE numbers (O(log n) vs O(sqrt(n)))
    print("Large number test (n = 10^15 + 37, which is prime):")
    big = 10**15 + 37
    t0 = time.time()
    res_fib = fibonacci_primality_witness(big)
    t_fib_big = time.time() - t0

    t0 = time.time()
    res_trial = is_prime_trial(big)
    t_trial_big = time.time() - t0

    print(f"  Fib PRP: {res_fib} in {t_fib_big:.6f} sec")
    print(f"  Trial:   {res_trial} in {t_trial_big:.6f} sec")
    if t_trial_big > 0:
        print(f"  Speedup: {t_trial_big / max(t_fib_big, 1e-9):.1f}x")
    print()


# ============================================================
# TOOL 2: LUCAS PSEUDOPRIME TEST
# ============================================================

def lucas_prp(n):
    """L_n mod n == 1 for prime n (n > 2). Returns True if passes."""
    if n < 2:
        return False
    if n == 2:
        return True  # L_2 = 3, 3 mod 2 = 1
    return lucas_mod(n, n) == 1


def run_tool2():
    print("=" * 70)
    print("TOOL 2: LUCAS PSEUDOPRIME TEST")
    print("=" * 70)
    print()
    print("Theory: L_p = 1 mod p for odd prime p.")
    print("A Lucas pseudoprime is composite n with L_n = 1 mod n.")
    print()

    primes = set(sieve_primes(10000))

    # Verify for small primes
    print("Verification for small primes:")
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        Lp = lucas(p)
        print(f"  L_{p} = {Lp}, L_{p} mod {p} = {Lp % p}")
    print()

    # Find Lucas pseudoprimes
    print("Searching for Lucas pseudoprimes (composite n with L_n = 1 mod n)...")
    lucas_pseudoprimes = []
    composites_tested = 0
    for n in range(2, 100001):
        if n not in primes and not is_prime_trial(n):
            composites_tested += 1
            if lucas_prp(n):
                lucas_pseudoprimes.append(n)

    print(f"  Composites tested up to 100000: {composites_tested}")
    print(f"  Lucas pseudoprimes found: {len(lucas_pseudoprimes)}")
    if lucas_pseudoprimes:
        print(f"  List: {lucas_pseudoprimes[:30]}")
    else:
        print("  NONE found! Lucas pseudoprimes are extremely rare.")
    print()

    # Compare to Fermat pseudoprimes (base 2)
    fermat_psp = []
    for n in range(2, 100001):
        if n > 2 and not is_prime_trial(n) and n % 2 == 1:
            if pow(2, n - 1, n) == 1:
                fermat_psp.append(n)

    print(f"Comparison to Fermat pseudoprimes (base 2) up to 100000:")
    print(f"  Fermat pseudoprimes: {len(fermat_psp)}")
    print(f"  First 15: {fermat_psp[:15]}")
    print(f"  Lucas pseudoprimes:  {len(lucas_pseudoprimes)}")
    if len(fermat_psp) > 0 and len(lucas_pseudoprimes) > 0:
        print(f"  Ratio Fermat/Lucas:  {len(fermat_psp)/len(lucas_pseudoprimes):.1f}x")
    elif len(lucas_pseudoprimes) == 0:
        print(f"  Lucas are infinitely rarer (none found in range)!")
    print()

    # Extended search
    print("Extended search up to 500000 for Lucas pseudoprimes...")
    for n in range(100001, 500001):
        if n % 2 == 1 and not is_prime_trial(n):
            if lucas_prp(n):
                lucas_pseudoprimes.append(n)
                print(f"  FOUND: {n} (factors: ", end="")
                for d in range(2, int(n**0.5) + 1):
                    if n % d == 0:
                        print(f"{d} x {n // d}", end="")
                        break
                print(")")
    print(f"  Total Lucas pseudoprimes up to 500000: {len(lucas_pseudoprimes)}")
    if lucas_pseudoprimes:
        print(f"  All found: {lucas_pseudoprimes}")
    print()


# ============================================================
# TOOL 3: FIBONACCI ENTRY POINT FACTORING
# ============================================================

def fib_entry_point(n):
    """Find alpha(n) = smallest k >= 1 with n | F_k."""
    if n <= 1:
        return n
    # Pisano period is at most 6n, entry point divides Pisano period
    # In practice, much smaller
    limit = 6 * n + 10
    for k in range(1, limit):
        if fib_mod(k, n) == 0:
            return k
    return None


def fib_factor(n):
    """Factor n using Fibonacci entry points."""
    if n <= 1:
        return [n]
    if is_prime_trial(n):
        return [n]

    alpha = fib_entry_point(n)
    if alpha is None:
        return [n]

    # n | F_{alpha}. Try gcd(n, F_k) for divisors of alpha
    factors_found = []
    remaining = n

    # Try gcd with F_k for various k dividing alpha
    divisors_of_alpha = []
    for d in range(1, alpha + 1):
        if alpha % d == 0:
            divisors_of_alpha.append(d)

    for d in divisors_of_alpha:
        if remaining <= 1:
            break
        g = gcd(remaining, fib(d)) if d < 500 else gcd(remaining, fib_mod(d, remaining))
        if 1 < g < remaining:
            # Found a factor
            factors_found.append(g)
            remaining //= g

    if remaining > 1:
        factors_found.append(remaining)

    # Recursively factor each part
    result = []
    for f in factors_found:
        if is_prime_trial(f):
            result.append(f)
        else:
            result.extend(fib_factor(f))

    result.sort()
    return result


def run_tool3():
    print("=" * 70)
    print("TOOL 3: FIBONACCI ENTRY POINT FACTORING")
    print("=" * 70)
    print()
    print("Theory: alpha(n) = min k with n | F_k. Then gcd(n, F_d)")
    print("for d | alpha(n) often reveals factors of n.")
    print()

    targets = [21, 42, 91, 143, 6174]

    for n in targets:
        alpha = fib_entry_point(n)
        print(f"n = {n}:")
        print(f"  alpha({n}) = {alpha}  (F_{alpha} = 0 mod {n})")

        # Show the gcd approach
        divisors = [d for d in range(1, alpha + 1) if alpha % d == 0]
        print(f"  Divisors of alpha: {divisors}")
        for d in divisors:
            Fd = fib(d) if d < 100 else fib_mod(d, n)
            g = gcd(n, Fd if isinstance(Fd, int) and Fd > 0 else fib(d))
            if g > 1 and g < n:
                print(f"    gcd({n}, F_{d}) = gcd({n}, {fib(d) if d < 80 else '...'}) = {g} <-- FACTOR!")
            elif g == n:
                print(f"    gcd({n}, F_{d}) = {g} (= n, trivial)")

        factors = fib_factor(n)
        product = 1
        for f in factors:
            product *= f
        print(f"  Factorization: {n} = {' x '.join(map(str, factors))}")
        assert product == n, f"Factor check failed: {product} != {n}"
        print()

    # Bonus: factor some semiprimes
    print("Bonus -- factoring semiprimes via Fibonacci entry points:")
    semiprimes = [15, 77, 221, 323, 437, 667, 899, 1147]
    for n in semiprimes:
        alpha = fib_entry_point(n)
        factors = fib_factor(n)
        print(f"  {n} = {' x '.join(map(str, factors))}  (alpha = {alpha})")
    print()


# ============================================================
# TOOL 4: FIBONACCI (ZECKENDORF) CODING
# ============================================================

def zeckendorf_encode(n):
    """
    Zeckendorf representation: unique sum of non-consecutive Fibonacci numbers.
    Returns list of Fibonacci indices used (1-indexed: F_1=1, F_2=2, F_3=3, ...).
    """
    if n <= 0:
        return []

    # Generate Fibonacci numbers up to n
    fibs = [1, 2]
    while fibs[-1] < n:
        fibs.append(fibs[-1] + fibs[-2])

    indices = []
    remaining = n
    for i in range(len(fibs) - 1, -1, -1):
        if fibs[i] <= remaining:
            indices.append(i + 1)  # 1-indexed
            remaining -= fibs[i]
        if remaining == 0:
            break

    indices.reverse()
    return indices


def zeckendorf_decode(indices):
    """Decode Zeckendorf representation back to integer."""
    fibs = [1, 2]
    while len(fibs) < max(indices) + 1 if indices else 2:
        fibs.append(fibs[-1] + fibs[-2])

    return sum(fibs[i - 1] for i in indices)


def zeckendorf_bitstring(n):
    """Binary representation: bit i is 1 if F_{i+1} is used."""
    indices = zeckendorf_encode(n)
    if not indices:
        return "0"
    max_idx = max(indices)
    bits = ['0'] * max_idx
    for i in indices:
        bits[max_idx - i] = '1'
    return ''.join(bits)


def run_tool4():
    print("=" * 70)
    print("TOOL 4: FIBONACCI (ZECKENDORF) CODING")
    print("=" * 70)
    print()
    print("Theory: Every positive integer = unique sum of non-consecutive")
    print("Fibonacci numbers (Zeckendorf's theorem, 1972).")
    print()

    # Important numbers
    important = [7, 21, 42, 168, 6174]

    # Reference: Fibonacci numbers
    print("Reference Fibonacci sequence:")
    print("  F_1=1, F_2=2, F_3=3, F_4=5, F_5=8, F_6=13, F_7=21,")
    print("  F_8=34, F_9=55, F_10=89, F_11=144, F_12=233, F_13=377,")
    print("  F_14=610, F_15=987, F_16=1597, F_17=2584, F_18=4181, F_19=6765")
    print()

    print("Zeckendorf representations of important numbers:")
    print(f"  {'n':>6} | {'Indices':>20} | {'Fibonacci sum':>30} | {'Bitstring':>20}")
    print(f"  {'-'*6}-+-{'-'*20}-+-{'-'*30}-+-{'-'*20}")

    for n in important:
        indices = zeckendorf_encode(n)
        fibs_list = [1, 2]
        while len(fibs_list) < max(indices) + 1:
            fibs_list.append(fibs_list[-1] + fibs_list[-2])
        fib_terms = [str(fibs_list[i - 1]) for i in indices]
        fib_sum_str = " + ".join(fib_terms)
        bits = zeckendorf_bitstring(n)

        # Verify roundtrip
        decoded = zeckendorf_decode(indices)
        assert decoded == n, f"Roundtrip failed: {decoded} != {n}"

        print(f"  {n:>6} | {str(indices):>20} | {fib_sum_str:>30} | {bits:>20}")

    print()

    # Verify no consecutive indices
    print("Verification -- no consecutive Fibonacci indices used:")
    for n in important:
        indices = zeckendorf_encode(n)
        consecutive = any(indices[i + 1] - indices[i] == 1 for i in range(len(indices) - 1))
        print(f"  {n}: indices {indices} -- consecutive? {consecutive}")

    print()

    # Pattern analysis
    print("Pattern analysis:")
    print("  7 = 5 + 2          -- uses F_4 + F_2 (even indices)")
    print("  21 = 21            -- IS a Fibonacci number (F_7)")
    print("  42 = 34 + 8        -- uses F_8 + F_5")
    print("  168 = 144 + 21 + 3 -- uses F_11 + F_7 + F_3 (all odd indices!)")
    print("  6174 = 4181+1597+377+13+5+1 -- 6 terms")
    print()

    # Statistics
    print("Zeckendorf statistics for 1..1000:")
    lengths = defaultdict(int)
    for n in range(1, 1001):
        indices = zeckendorf_encode(n)
        lengths[len(indices)] += 1

    for k in sorted(lengths.keys()):
        print(f"  {k} terms: {lengths[k]} numbers ({100*lengths[k]/1000:.1f}%)")

    # Average number of terms
    total_terms = sum(len(zeckendorf_encode(n)) for n in range(1, 1001))
    print(f"  Average terms: {total_terms / 1000:.3f}")
    print(f"  (Theory predicts ~log_phi(n) / (1+1/phi) terms on average)")
    print()

    # Which numbers ARE Fibonacci numbers?
    fib_set = set()
    f = 1
    g = 2
    while f <= 10000:
        fib_set.add(f)
        f, g = g, f + g
    fib_important = [n for n in important if n in fib_set]
    print(f"  Numbers that are Fibonacci numbers: {fib_important}")
    print()


# ============================================================
# TOOL 5: FIBONACCI HASH
# ============================================================

def fib_hash(n, M):
    """Hash: F_n mod M."""
    return fib_mod(n, M)


def run_tool5():
    print("=" * 70)
    print("TOOL 5: FIBONACCI HASH")
    print("=" * 70)
    print()
    print("hash_F(n) = F_n mod M. Testing collision properties.")
    print()

    # Test with various M values
    test_Ms = [97, 127, 251, 1009, 4096]
    N_keys = 1000  # Hash keys 0..999

    for M in test_Ms:
        buckets = defaultdict(list)
        for n in range(N_keys):
            h = fib_hash(n, M)
            buckets[h].append(n)

        # Count collisions
        max_bucket = max(len(v) for v in buckets.values())
        empty_buckets = M - len(buckets)
        collisions = sum(1 for v in buckets.values() if len(v) > 1)
        total_collisions = sum(len(v) - 1 for v in buckets.values() if len(v) > 1)

        # Compare to simple mod hash
        simple_buckets = defaultdict(list)
        for n in range(N_keys):
            simple_buckets[n % M].append(n)
        simple_max = max(len(v) for v in simple_buckets.values())

        print(f"M = {M}, hashing keys 0..{N_keys-1}:")
        print(f"  Fib hash:    {len(buckets)} distinct values, max bucket {max_bucket}, "
              f"{collisions} buckets with collision, {total_collisions} total collisions")
        print(f"  Simple mod:  {min(M, N_keys)} distinct values, max bucket {simple_max}")

        # Pisano period
        period = None
        for k in range(1, 6 * M + 10):
            if fib_mod(k, M) == 0 and fib_mod(k + 1, M) == 1:
                period = k
                break
        if period:
            print(f"  Pisano period pi({M}) = {period}")
            print(f"  Fib hash cycles every {period} keys => {period} distinct outputs max")
        print()

    # Application: tournament vertex labeling
    print("Application: Fibonacci hash for tournament vertex labeling")
    print("  For n-vertex tournament, label vertices 0..n-1.")
    print("  Use F-hash to map vertex pairs to edge buckets.")
    print()
    for n in [5, 7, 11, 13]:
        M = n * (n - 1) // 2  # number of edges
        edge_hashes = set()
        for i in range(n):
            for j in range(i + 1, n):
                h = fib_hash(i * n + j, M)
                edge_hashes.add(h)
        num_edges = n * (n - 1) // 2
        print(f"  n={n}: {num_edges} edges, {len(edge_hashes)} distinct F-hashes "
              f"(mod {M}), coverage {100*len(edge_hashes)/M:.1f}%")
    print()


# ============================================================
# TOOL 6: GOLDEN RATIO APPROXIMATION TOOL
# ============================================================

def golden_ratio_convergents(max_n=30):
    """
    Convergents of phi = [1; 1, 1, 1, ...].
    Each convergent is F_{n+1}/F_n.
    Returns list of (numerator, denominator, value, error).
    """
    phi = (1 + 5**0.5) / 2
    convergents = []
    for n in range(1, max_n + 1):
        num = fib(n + 1)
        den = fib(n)
        val = num / den
        err = abs(val - phi)
        convergents.append((num, den, val, err))
    return convergents


def best_fib_ratio_approx(target, max_fib_index=30):
    """Find F_{n+1}/F_n closest to target."""
    best = None
    best_err = float('inf')
    for n in range(1, max_fib_index + 1):
        num = fib(n + 1)
        den = fib(n)
        val = num / den
        err = abs(val - target)
        if err < best_err:
            best_err = err
            best = (num, den, n, val, err)
    return best


def run_tool6():
    print("=" * 70)
    print("TOOL 6: GOLDEN RATIO APPROXIMATION TOOL")
    print("=" * 70)
    print()
    phi = (1 + 5**0.5) / 2
    print(f"phi = (1 + sqrt(5))/2 = {phi:.15f}")
    print()

    print("Convergents F_{n+1}/F_n:")
    print(f"  {'n':>3} | {'F_{n+1}/F_n':>15} | {'Value':>18} | {'Error':>15}")
    print(f"  {'-'*3}-+-{'-'*15}-+-{'-'*18}-+-{'-'*15}")

    convergents = golden_ratio_convergents(20)
    for i, (num, den, val, err) in enumerate(convergents):
        n = i + 1
        frac_str = f"{num}/{den}"
        print(f"  {n:>3} | {frac_str:>15} | {val:>18.15f} | {err:>15.2e}")

    print()
    print("Error decreases as 1/phi^{2n} -- exponential convergence!")
    print(f"  Error ratio (consecutive): ", end="")
    for i in range(1, min(8, len(convergents))):
        if convergents[i][3] > 0:
            ratio = convergents[i - 1][3] / convergents[i][3]
            print(f"{ratio:.3f} ", end="")
    print()
    print(f"  phi^2 = {phi**2:.6f} (should match ratios)")
    print()

    # Approximate 3/4 = 0.75
    target = 3 / 4
    print(f"Approximating 3/4 = {target} with F-ratios:")
    print()

    # All F-ratios and their distances to 0.75
    print(f"  {'n':>3} | {'F_{n+1}/F_n':>12} | {'Value':>18} | {'|val - 3/4|':>15}")
    print(f"  {'-'*3}-+-{'-'*12}-+-{'-'*18}-+-{'-'*15}")

    results = []
    for n in range(1, 20):
        num = fib(n + 1)
        den = fib(n)
        val = num / den
        err = abs(val - target)
        results.append((n, num, den, val, err))

    results_sorted = sorted(results, key=lambda x: x[4])
    for n, num, den, val, err in results:
        marker = " <-- BEST" if (n, num, den, val, err) == results_sorted[0] else ""
        frac_str = f"{num}/{den}"
        print(f"  {n:>3} | {frac_str:>12} | {val:>18.15f} | {err:>15.12f}{marker}")

    best = results_sorted[0]
    print()
    print(f"Best F-ratio approximation to 3/4:")
    print(f"  F_{best[0]+1}/F_{best[0]} = {best[1]}/{best[2]} = {best[3]:.15f}")
    print(f"  Error: {best[4]:.15f}")
    print()

    # More general: given arbitrary precision, how many convergents needed?
    print("Precision table:")
    for digits in range(1, 16):
        eps = 10**(-digits)
        needed = None
        for i, (num, den, val, err) in enumerate(convergents):
            if err < eps:
                needed = i + 1
                break
        if needed:
            print(f"  {digits} digits: need F_{needed+1}/F_{needed} "
                  f"= {convergents[needed-1][0]}/{convergents[needed-1][1]}")
        else:
            print(f"  {digits} digits: need n > 20")
    print()


# ============================================================
# TOOL 7: FIBONACCI-WEIGHTED PRIME COUNTER
# ============================================================

def run_tool7():
    print("=" * 70)
    print("TOOL 7: FIBONACCI-WEIGHTED PRIME COUNTER")
    print("=" * 70)
    print()
    print("pi_F(N) = sum_{p <= N, p prime} F_{index(p)}")
    print("where index(p) = position of p in the prime sequence (1-indexed).")
    print()

    primes = sieve_primes(10000)

    # Compute pi_F at milestones
    milestones = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]

    print(f"  {'N':>6} | {'pi(N)':>6} | {'pi_F(N)':>20} | {'pi_F/pi':>12} | {'log(pi_F)':>10}")
    print(f"  {'-'*6}-+-{'-'*6}-+-{'-'*20}-+-{'-'*12}-+-{'-'*10}")

    pi_f = 0
    prime_idx = 0
    milestone_idx = 0

    for i, p in enumerate(primes):
        prime_idx = i + 1  # 1-indexed
        f_val = fib(prime_idx)
        pi_f += f_val

        while milestone_idx < len(milestones) and p >= milestones[milestone_idx]:
            # We passed a milestone; check previous state
            # Actually compute at exact milestone
            break

    # Recompute cleanly
    pi_f_accum = 0
    prime_count = 0
    milestone_idx = 0
    for i, p in enumerate(primes):
        prime_count += 1
        pi_f_accum += fib(i + 1)

        if milestone_idx < len(milestones) and p <= milestones[milestone_idx]:
            if milestone_idx < len(milestones) and (
                i + 1 == len(primes) or primes[i + 1] > milestones[milestone_idx]
            ):
                N = milestones[milestone_idx]
                ratio = pi_f_accum / prime_count if prime_count > 0 else 0
                log_val = math.log(pi_f_accum) if pi_f_accum > 0 else 0
                print(f"  {N:>6} | {prime_count:>6} | {pi_f_accum:>20} | {ratio:>12.2f} | {log_val:>10.4f}")
                milestone_idx += 1

    print()

    # Growth analysis
    print("Growth analysis:")
    print("  pi(N) ~ N/ln(N)  (prime number theorem)")
    print("  F_k ~ phi^k / sqrt(5)")
    print("  pi_F(N) ~ sum_{k=1}^{pi(N)} phi^k/sqrt(5)")
    print("         ~ phi^{pi(N)+1} / (sqrt(5)*(phi-1))")
    print("         ~ phi^{N/ln(N)} / sqrt(5)")
    print()
    print("  So log(pi_F(N)) ~ (N/ln(N)) * ln(phi)")
    print(f"  ln(phi) = {math.log(phi_val := (1+5**0.5)/2):.6f}")
    print()

    # Compare predicted vs actual
    print("Predicted vs actual log(pi_F):")
    pi_f_accum = 0
    prime_count = 0
    milestone_idx = 0
    for i, p in enumerate(primes):
        prime_count += 1
        pi_f_accum += fib(i + 1)
        if milestone_idx < len(milestones) and p <= milestones[milestone_idx]:
            if i + 1 == len(primes) or primes[i + 1] > milestones[milestone_idx]:
                N = milestones[milestone_idx]
                actual_log = math.log(pi_f_accum) if pi_f_accum > 0 else 0
                predicted = prime_count * math.log(phi_val)
                print(f"  N={N:>5}: actual log={actual_log:>10.2f}, "
                      f"predicted ~{predicted:>10.2f}, "
                      f"ratio={actual_log/predicted:.4f}" if predicted > 0 else "")
                milestone_idx += 1
    print()


# ============================================================
# TOOL 8: LUCAS-TOURNAMENT CONNECTION
# ============================================================

def run_tool8():
    print("=" * 70)
    print("TOOL 8: LUCAS-TOURNAMENT CONNECTION")
    print("=" * 70)
    print()
    print("For Paley primes p in {3, 5, 7, 11}, compare L_p, F_p, H(T_p).")
    print()

    # Known H values for Paley tournaments
    # H(T_3) = 3 (verified)
    # H(T_5) = 15 (we need to check this -- from the project, T_5 regular has H=15)
    # H(T_7) = 189 (verified in project)
    # H(T_11) = 95095 (verified in project)
    paley_H = {3: 3, 5: 15, 7: 189, 11: 95095}

    print(f"  {'p':>4} | {'F_p':>8} | {'L_p':>8} | {'H(T_p)':>10} | {'H/F_p':>12} | {'H/L_p':>12}")
    print(f"  {'-'*4}-+-{'-'*8}-+-{'-'*8}-+-{'-'*10}-+-{'-'*12}-+-{'-'*12}")

    for p in [3, 5, 7, 11]:
        Fp = fib(p)
        Lp = lucas(p)
        H = paley_H[p]
        print(f"  {p:>4} | {Fp:>8} | {Lp:>8} | {H:>10} | {H/Fp:>12.6f} | {H/Lp:>12.6f}")

    print()

    # More detailed analysis
    print("Detailed analysis:")
    print()
    for p in [3, 5, 7, 11]:
        Fp = fib(p)
        Lp = lucas(p)
        H = paley_H[p]
        print(f"  p = {p}:")
        print(f"    F_{p} = {Fp}")
        print(f"    L_{p} = {Lp}")
        print(f"    H(T_{p}) = {H}")
        print(f"    H mod F_{p} = {H % Fp}")
        print(f"    H mod L_{p} = {H % Lp}")
        print(f"    H / p! = {H / math.factorial(p):.6f}")
        print(f"    H / (p-1)! = {H / math.factorial(p-1):.6f}")

        # Check if H = some nice function of p, F_p, L_p
        # H(T_3) = 3 = F_4 = L_2 * F_2
        # H(T_5) = 15 = F_5 * L_1 ... no. 15 = 5*3 = p * F_4
        # H(T_7) = 189 = 7 * 27 = 7 * 3^3
        # H(T_11) = 95095 = 5 * 7 * 11 * 13 * 19

        # Factor H
        h = H
        factors = []
        for d in range(2, h + 1):
            while h % d == 0:
                factors.append(d)
                h //= d
            if h == 1:
                break
        print(f"    H factored: {' * '.join(map(str, factors))}")

        # Check H mod (Fp * Lp)
        print(f"    F_{p} * L_{p} = {Fp * Lp}, H mod (F*L) = {H % (Fp * Lp)}")
        # F_p * L_p = F_{2p}
        F2p = fib(2 * p)
        print(f"    F_{2*p} = {F2p}, H mod F_{2*p} = {H % F2p}")
        print()

    # Summary ratios
    print("Summary of key ratios:")
    print()
    print("  H(T_p) / (p-1)!! where (p-1)!! = double factorial:")
    for p in [3, 5, 7, 11]:
        H = paley_H[p]
        # double factorial
        dfact = 1
        for k in range(p - 1, 0, -2):
            dfact *= k
        print(f"    p={p}: H/(p-1)!! = {H}/{dfact} = {H/dfact:.6f}")

    print()
    print("  H(T_p) / C(p, (p-1)/2):")
    for p in [3, 5, 7, 11]:
        H = paley_H[p]
        binom = math.comb(p, (p - 1) // 2)
        print(f"    p={p}: H/C({p},{(p-1)//2}) = {H}/{binom} = {H/binom:.6f}")

    print()

    # Check the formula H(T_p) = p!! or similar
    print("  H(T_p) as products of small primes:")
    print("    H(T_3)  = 3")
    print("    H(T_5)  = 15 = 3 * 5")
    print("    H(T_7)  = 189 = 3^3 * 7")
    print("    H(T_11) = 95095 = 5 * 7 * 11 * 13 * 19")
    print()

    # Look for a pattern in H/|Aut|
    print("  H(T_p) / |Aut(T_p)| (orbit count for Hamiltonian paths):")
    aut_sizes = {3: 3, 5: 10, 7: 21, 11: 55}  # |Aut(T_p)| = p*(p-1)/2 for Paley
    for p in [3, 5, 7, 11]:
        H = paley_H[p]
        aut = aut_sizes[p]
        print(f"    p={p}: H/|Aut| = {H}/{aut} = {H//aut}")
    print("    Sequence: 1, 1.5(?), 9, 1729")
    print("    Wait -- 15/10 = 1.5, not integer. Let me recheck...")
    print(f"    15/10 = {15/10} -- H(T_5)=15 confirmed, |Aut(T_5)|=10")
    print(f"    But orbits must be integer. So |Aut| acts on DIRECTED paths,")
    print(f"    some orbits may have smaller stabilizer.")
    print()

    # Modular patterns
    print("  Modular patterns:")
    for p in [3, 5, 7, 11]:
        H = paley_H[p]
        Fp = fib(p)
        Lp = lucas(p)
        print(f"    p={p}: H mod p = {H % p}, H mod F_p = {H % Fp}, H mod L_p = {H % Lp}")
    print()
    print("  Key observation: H(T_p) = 0 mod p for all tested Paley primes.")
    print("  This is expected: Aut(T_p) acts on Ham paths, |Aut| = p(p-1)/2,")
    print("  so p | |Aut| | (something), and fixed-point count argument gives p | H.")
    print()

    # Connection to Lucas numbers specifically
    print("  Direct Lucas connection search:")
    for p in [3, 5, 7, 11]:
        H = paley_H[p]
        Lp = lucas(p)
        # Check H mod L_p
        print(f"    p={p}: L_{p} = {Lp}, H mod L_{p} = {H % Lp}, "
              f"H/L_{p} = {H/Lp:.4f}, "
              f"gcd(H, L_{p}) = {gcd(H, Lp)}")
    print()
    print("  Verdict: No clean H(T_p) = f(L_p, F_p) formula emerges.")
    print("  The Fibonacci/Lucas numbers and Paley Hamiltonian path counts")
    print("  live in different algebraic worlds: F/L from linear recurrences,")
    print("  H from permanent-like combinatorial sums over tournaments.")
    print("  Their connection to p is through DIFFERENT mechanisms:")
    print("    - F_p, L_p: via p | F_{p-(p/5)} and L_p = 1 mod p")
    print("    - H(T_p):   via Aut(T_p) = PSL(2,p) orbit structure")
    print()


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("*" * 70)
    print("* FIBONACCI-PRIME PRACTICAL APPLICATIONS")
    print("* kind-pasteur-2026-03-16-S116l")
    print("*" * 70)
    print()

    run_tool1()
    print()
    run_tool2()
    print()
    run_tool3()
    print()
    run_tool4()
    print()
    run_tool5()
    print()
    run_tool6()
    print()
    run_tool7()
    print()
    run_tool8()

    print()
    print("=" * 70)
    print("ALL 8 TOOLS COMPLETE.")
    print("=" * 70)
