#!/usr/bin/env python3
"""
S628 — computational SPEEDUPS for LRC, from our structure (THM-411/412/415/418, HYP-2245).
Each fast routine is paired with a brute reference and a correctness check + benchmark.

Convention (canon): n runners, speeds = n-1 distinct positive ints, gap 1/n,
M(S) = max_t min_i ||v_i t||.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, time, random

# ---------------------------------------------------------------- primes
def isprime(k):
    if k < 2: return False
    i = 2
    while i*i <= k:
        if k % i == 0: return False
        i += 1
    return True
def least_prime_ge(k):
    while not isprime(k): k += 1
    return k

# ============================================================ (A) good multiplier
def good_mult_brute(S, p):
    """O(p*n): scan every multiplier a, check all runners off band {0,+-1}."""
    for a in range(1, p):
        if all(min((a*v) % p, p-((a*v) % p)) >= 2 for v in S):
            return a
    return None

def good_mult_fast(S, p):
    """O(n): forbidden a = union_i {+-v_i^{-1} mod p}; pick any a not forbidden.
       Correct for PRIME p with p∤v_i (THM-418: band {0,+-1} forbids exactly a=+-v_i^{-1})."""
    bad = bytearray(p)                      # bad[a]=1 if forbidden
    bad[0] = 1
    for v in S:
        r = v % p
        if r == 0:                          # runner stuck at 0 -> no good a on this shell
            return None
        inv = pow(r, p-2, p)                # modular inverse, p prime
        bad[inv] = 1; bad[p-inv] = 1
    for a in range(1, p):
        if not bad[a]:
            return a
    return None

# ============================================================ (B) gap via shell witnesses
def norm(x):
    f = x - (x.numerator // x.denominator)
    return f if f <= Fr(1, 2) else 1-f

def gap_brute(speeds):
    """O(n^2 * vmax): all corner + crossing critical times (the old exact gap)."""
    V = [abs(v) for v in speeds]; cands = set()
    for i in range(len(V)):
        vi = V[i]
        for k in range(0, 2*vi+1):
            t = Fr(2*k+1, 2*vi)
            if 0 < t <= Fr(1, 2): cands.add(t)
        for j in range(i):
            vj = V[j]
            for d in (vi+vj, abs(vi-vj)):
                if d == 0: continue
                kk = 1
                while Fr(kk, d) <= Fr(1, 2):
                    cands.add(Fr(kk, d)); kk += 1
    best = Fr(0)
    for t in cands:
        m = min(norm(v*t) for v in V)
        if m > best: best = m
    return best

def gap_shells(speeds, n, Mmax=None):
    """O(sum_{m<=2n-1} phi(m) * n): max of min_i||v_i a/m|| over shell witnesses a/m, m<=2n-1.
       Exact for tight/near-tight configs (witness on the (Z/m)* orbit, THM-411/415); a lower
       bound in general. Integer-only (dist via residues)."""
    if Mmax is None: Mmax = 2*n-1
    best = Fr(0)
    for m in range(2, Mmax+1):
        for a in range(1, m//2+1):
            if gcd(a, m) != 1: continue
            md = min(min((a*v) % m, m-((a*v) % m)) for v in speeds)
            val = Fr(md, m)
            if val > best: best = val
    return best

# ============================================================ (C) fast looseness certificate
def is_loose_fast(S, n):
    """fast lower bound on M via prime-shell dodge: M >= 2/p, p=least prime>=2n with p∤v_i.
       Returns (lower_bound_Fraction, witness_a, p)."""
    p = least_prime_ge(2*n)
    while any(v % p == 0 for v in S):
        p = least_prime_ge(p+1)
    a = good_mult_fast(S, p)
    return (Fr(2, p) if a is not None else Fr(0), a, p)

# ============================================================ benchmarks
if __name__ == "__main__":
    rng = random.Random(0)
    print("(A) GOOD-MULTIPLIER  fast O(n)  vs  brute O(p*n)")
    mism = 0; tb = tf = 0.0
    for _ in range(2000):
        n = rng.randint(6, 30); p = least_prime_ge(2*n)
        S = rng.sample(range(1, p), n-1)
        t0 = time.perf_counter(); b = good_mult_brute(S, p); tb += time.perf_counter()-t0
        t0 = time.perf_counter(); f = good_mult_fast(S, p); tf += time.perf_counter()-t0
        # both should agree on existence; fast may return a different valid a
        if (b is None) != (f is None): mism += 1
        if f is not None and not all(min((f*v) % p, p-((f*v) % p)) >= 2 for v in S): mism += 1
    print(f"   existence/validity mismatches: {mism}/2000")
    print(f"   brute {tb*1000:.1f} ms   fast {tf*1000:.1f} ms   speedup x{tb/tf:.1f}")

    print("\n(B) GAP  shell-restricted (<=2n-1)  vs  brute crossings  (correctness + speed)")
    eq = 0; tot = 0; tb = ts = 0.0; ge = 0
    for _ in range(400):
        n = rng.randint(5, 9); R = rng.randint(2*n, 3*n)
        S = sorted(rng.sample(range(1, R), n-1))
        t0 = time.perf_counter(); gb = gap_brute(S); tb += time.perf_counter()-t0
        t0 = time.perf_counter(); gs = gap_shells(S, n); ts += time.perf_counter()-t0
        tot += 1
        if gs == gb: eq += 1
        if gs <= gb: ge += 1
    print(f"   shell-gap == brute-gap: {eq}/{tot}  (shell-gap <= brute always: {ge}/{tot})")
    print(f"   brute {tb*1000:.1f} ms   shell {ts*1000:.1f} ms   speedup x{tb/ts:.1f}")
    print("   (shell-gap is EXACT when the optimum is on the (Z/m)* orbit; else a fast lower bound)")

    print("\n(C) LOOSENESS certificate  is_loose_fast (O(n))  — sample")
    for S, n in [([14,2,3,4,5,6,7,8,9,10,11,12,13], 14), ([1,2,3,4,5,6], 7)]:
        lb, a, p = is_loose_fast(S, n)
        print(f"   n={n} S={S[:4]}...: M >= {lb}={float(lb):.4f} (>1/n? {lb>Fr(1,n)}) via a/{p}={a}/{p}")
