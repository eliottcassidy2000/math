# The prime 2 for strong tournaments.  s(n)/n! = sum over strong classes of 1/|Aut| is a 2-adic integer (Moon: |Aut| odd),
# so v2(s(n)) >= v2(n!), and s(n)/2^v2(n!) = #strong classes (mod 2).  The EGF T(x) = sum 2^C(n,2) x^n/n! is congruent to
# Phi_3(x) = 1 + x + x^2 modulo 2, hence the parity of the class count is [n != 2 mod 3].  This script computes
# e(n) = v2(s(n)) - v2(n!) exactly from Wright's recurrence truncated mod 2^N (the terms 2^C(n-k,2) die quickly),
# and tests: e(n) = 0 iff n != 2 mod 3; e(n) = 2 for n = 5 mod 6; e(n) = 3 + v2((n-2)/6) for n = 2 mod 6.
# mac-mini-2026-10-21-fanopsp.
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
N = NMAX + 64                                  # working precision 2^N  (v2(n!) < n)
MOD = 1 << N
def v2(x):
    return (x & -x).bit_length() - 1 if x else 10**9
s = [0] * (NMAX + 1); fact_v2 = [0] * (NMAX + 1)
binom = [1]                                    # row n of Pascal's triangle mod 2^N
for n in range(1, NMAX + 1):
    binom = [1] + [(binom[i - 1] + binom[i]) % MOD for i in range(1, n)] + [1]
    acc = (1 << (n * (n - 1) // 2)) % MOD if n * (n - 1) // 2 < N else 0
    for k in range(1, n):
        e = (n - k) * (n - k - 1) // 2
        if e >= N: continue
        acc -= binom[k] * s[k] << e
    s[n] = acc % MOD
    fact_v2[n] = fact_v2[n - 1] + v2(n)
bad = []
for n in range(1, NMAX + 1):
    e = v2(s[n]) - fact_v2[n]
    if n % 3 != 2: pred = 0
    elif n % 2 == 1: pred = 2
    else: pred = 3 + v2((n - 2) // 6) if n > 2 else None
    if n == 2: continue                        # s(2) = 0
    if e != pred: bad.append((n, e, pred))
print(f"n <= {NMAX}: e(n) = v2(s(n)) - v2(n!) checked; failures: {len(bad)}")
for b in bad[:20]: print("  n=%d e=%d predicted=%s" % b)
print("first values e(1..40):", [v2(s[n]) - fact_v2[n] if s[n] else None for n in range(1, 41)])
