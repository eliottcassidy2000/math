# Independent audit of THM-4887, part 5: the comparison claim in the status block.
#  pi(2^32) - pi(2^31) by a segmented odd-only sieve (numpy); C(N, 31) vs 10^213; product bound 2^992 < 2^1000;
#  X^(2/3) at X = 2^1000; the corollary's actual value there; and how far the elementary product construction
#  (products of k distinct primes from [2^31, 2^32], all certified by THM-4875 (D), joined by THM-1862) stays ahead
#  of the corollary X^(2/3)/(13 ln X + 121).
import math, time
import numpy as np

def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]

def count_primes_range(lo, hi, seg=1 << 24):
    """number of primes in [lo, hi), lo odd, hi even, lo > sqrt(hi)."""
    assert lo % 2 == 1 and hi % 2 == 0
    ps = small_primes(int(math.isqrt(hi)) + 1)[1:]          # odd primes up to sqrt(hi)
    total = 0
    start = lo
    while start < hi:
        end = min(start + 2 * seg, hi)                       # odd numbers start, start+2, ..., < end
        L = (end - start) // 2
        mark = np.ones(L, dtype=bool)
        for p in ps:
            p = int(p)
            # first odd multiple of p that is >= start (and >= p*p)
            f = max(p * p, ((start + p - 1) // p) * p)
            if f % 2 == 0: f += p
            if f >= end: continue
            mark[(f - start) // 2::p] = False
        total += int(mark.sum())
        start = end
    return total

t0 = time.time()
N = count_primes_range((1 << 31) + 1, 1 << 32)
print(f"pi(2^32) - pi(2^31) = {N} by segmented sieve ({time.time()-t0:.0f}s); claimed 98,182,656: {N == 98182656}")
try:
    import sympy
    a, b = sympy.primepi(2 ** 32), sympy.primepi(2 ** 31)
    print(f"sympy.primepi: pi(2^32) = {a}, pi(2^31) = {b}, difference {a - b}")
except Exception as e:
    print("sympy unavailable:", e)

c31 = math.comb(N, 31)
print(f"C(N, 31) = {c31:.6e} = 10^{math.log10(c31):.3f}; > 10^213: {c31 > 10 ** 213}")
print(f"(2^32)^31 = 2^{32*31} < 2^1000: {32*31 < 1000};  X^(2/3) at X = 2^1000 = 10^{2000/3*math.log10(2):.3f}")
lnX = 1000 * math.log(2)
cor = 2000 / 3 * math.log10(2) - math.log10(13 * lnX + 121)
print(f"corollary X^(2/3)/(13 ln X + 121) at X = 2^1000 = 10^{cor:.3f}")

# distinct-prime construction, made continuous in X: for k >= 1, a set of k distinct primes from [2^31, 2^32) times an
# odd certified cofactor y < 2^31 (y = 1 allowed; y != 7, 21) is a Hamiltonian-path count (THM-4875 (D) + THM-1862),
# and different (set, y) give different integers (the primes >= 2^31 of the product are exactly the set).  Products
# lie below 2^(32k) y.  So for X in [2^(32k), 2^(32k+32)):  A(X) >= C(N, k) * #{odd y <= min(X/2^(32k), 2^31 - 1),
# y not in {7, 21}}.  Find the first X (on a grid of log2 X with step 1/64, plus each interval's endpoints) where this
# falls below the corollary X^(2/3)/(13 ln X + 121).
def log10_comb(n, k):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(10)
def log10_cor(log2X):
    lx = log2X * math.log(2)
    return (2 / 3) * log2X * math.log10(2) - math.log10(13 * lx + 121)
def log10_odd_count(Y):
    c = (int(Y) + 1) // 2 - (1 if Y >= 7 else 0) - (1 if Y >= 21 else 0)
    return math.log10(c) if c > 0 else -1e9
first_fail = None
for k in range(1, 2000):
    lc = log10_comb(N, k)
    for j in range(0, 32 * 64):
        l2 = 32 * k + j / 64
        Y = min(2.0 ** (j / 64), 2 ** 31 - 1)
        if lc + log10_odd_count(Y) < log10_cor(l2):
            first_fail = l2; break
    if first_fail is not None: break
print(f"products (k distinct primes in [2^31,2^32)) x (odd cofactor < 2^31) stay above the corollary for every X < 2^{first_fail:.2f}"
      f" (k = {k} at the first failure)")
print(f"at the first failure: construction 10^{lc + log10_odd_count(min(2.0 ** (first_fail - 32 * k), 2**31 - 1)):.1f}, corollary 10^{log10_cor(first_fail):.1f}")
