#!/usr/bin/env python3
"""
untouchable_census_goldbach_dependence_kps_20261025.py  (kind-pasteur-2026-10-25-S4)

Exact census of untouchable numbers (positive integers m that are not s(n) = sigma(n) - n for
any n >= 2) up to N, with:
  * a complete case split proving that every n with s(n) <= N is examined (no sieve factor
    heuristics): n = p*k with p = spf(n), k = n/p; since 1 and k are proper divisors of a composite
    n, s(n) >= 1 + k, so k <= N - 1.  Cases: p <= C (n <= C*N, direct sigma sieve); p > C and k = 1
    (n prime, s = 1); p > C and k prime (n = p*q or p^2: Goldbach convolution / p + 1);
    p > C and k composite with spf(k) >= p (explicit enumeration, few pairs).
  * controls against A057978 (counts <= 10^4, 10^5, 10^6 are 1212, 13863, 150232) and A005114
    (first terms), plus an independent brute-force path for N = 3000 (all n <= (N-1)^2).
  * the odd untouchables (expected: only 5), the Goldbach-dependence statistic: odd m whose only
    sources are semiprimes p*q (so that their touchability rests on m-1 = p+q with p != q),
  * residue tables (mod 2,3,4,5,6,7,12,14,15) of the untouchables <= N,
  * the Chowla shift: m is s'-touchable (s' = sigma - n - 1) iff m+1 is s-touchable.

Usage: python3 this.py [N]   (default 10**6).  Needs numpy.
"""
import sys
import time
import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
C = 10                      # sieve covers every n with spf(n) <= C, i.e. n <= C*(N-1) <= L
L = C * N
t0 = time.time()
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("" if detail == "" else "  [" + str(detail) + "]"))


def sigma_sieve(limit):
    """sigma(n) for 0 <= n <= limit (sigma(0) = 0)."""
    s = np.zeros(limit + 1, dtype=np.int64)
    for d in range(1, limit + 1):
        s[d::d] += d
    return s


def spf_sieve(limit):
    spf = np.zeros(limit + 1, dtype=np.int64)
    spf[1] = 1
    for i in range(2, limit + 1):
        if spf[i] == 0:
            spf[i::i][spf[i::i] == 0] = i
    return spf


def touchable_set(N, C):
    """Boolean array touched[0..N]; touched[m] iff m = s(n) for some n >= 2.  Also returns the
    array odd_nonsemiprime_touch[0..N] (odd m having a source n that is NOT p*q with p<q odd primes
    nor p^2), and the number of (p,k) pairs enumerated in the composite case."""
    L = C * N
    sig = sigma_sieve(L)
    n_arr = np.arange(L + 1, dtype=np.int64)
    s_arr = sig - n_arr
    s_arr[0] = -1
    s_arr[1] = -1                      # n = 1 has s = 0; excluded (n >= 2)
    touched = np.zeros(N + 1, dtype=bool)
    nonsemi = np.zeros(N + 1, dtype=bool)   # odd m touched by a non-semiprime, non-prime-square n
    # --- case (i): all n <= L (covers every n with spf(n) <= C), plus any n <= L at all
    mask = (s_arr >= 1) & (s_arr <= N)
    touched[s_arr[mask]] = True
    # which n <= L are semiprimes p*q (p<q) or prime squares?  Use spf and the cofactor.
    spf = spf_sieve(L)
    cof = np.zeros(L + 1, dtype=np.int64)
    cof[2:] = n_arr[2:] // spf[2:]
    is_prime = np.zeros(L + 1, dtype=bool)
    is_prime[2:] = (spf[2:] == n_arr[2:])
    semi_or_sq = np.zeros(L + 1, dtype=bool)
    semi_or_sq[2:] = is_prime[cof[2:]]          # n = p * q with q prime (includes p = q and p = 2)
    semi_or_sq[is_prime] = False
    # odd m sources: n with s(n) odd.  non-semiprime sources among n <= L:
    mask2 = mask & ~semi_or_sq & ~is_prime
    vals = s_arr[mask2]
    nonsemi[vals[(vals % 2 == 1)]] = True
    # --- case (ii): n prime -> s = 1
    touched[1] = True
    # --- case (iii): n = p*q, p <= q primes, p > C (p <= C already covered since n <= C*q <= C*N
    #     when q <= N; but we simply add ALL p <= q via the convolution; duplicates are harmless)
    P = is_prime[: N + 1].astype(np.float64)
    size = 1
    while size < 2 * (N + 1):
        size *= 2
    F = np.fft.rfft(P, size)
    conv = np.fft.irfft(F * F, size)[: N + 1]
    r = np.rint(conv).astype(np.int64)          # ordered pairs (p,q), p+q = e, both prime
    e = np.arange(N + 1)
    half_prime = np.zeros(N + 1, dtype=bool)
    half_prime[e % 2 == 0] = is_prime[(e[e % 2 == 0] // 2)]
    distinct = r - half_prime.astype(np.int64)  # ordered pairs with p != q
    gold = np.zeros(N + 1, dtype=bool)
    gold[1:] = distinct[:-1] > 0                # m = e + 1 with e = p + q, p != q
    touched |= gold
    # n = p^2: s = p + 1 (p <= N - 1)
    sq = np.zeros(N + 1, dtype=bool)
    pr = np.nonzero(is_prime[: N])[0]
    sq[pr + 1] = True
    touched |= sq
    # --- case (iv): p > C, k composite, spf(k) >= p.  Enumerate k <= N - 1 composite with spf(k) > C.
    pairs = 0
    primes_N = np.nonzero(is_prime[: N + 1])[0]
    sigN = sig[: N + 1]
    spfN = spf[: N + 1]
    ks = np.nonzero((spfN > C) & (~is_prime[: N + 1]) & (np.arange(N + 1) >= 4))[0]
    import bisect
    plist = primes_N.tolist()
    sig_k_list = sigN.tolist()
    spf_k_list = spfN.tolist()
    for k in ks.tolist():
        if k > N - 1:
            continue
        sp = spf_k_list[k]
        hi = bisect.bisect_right(plist, sp)
        lo = bisect.bisect_right(plist, C)
        sk = sig_k_list[k]
        for p in plist[lo:hi]:
            n = p * k
            if n <= L:
                continue                         # already in case (i)
            pairs += 1
            if k % p == 0:
                a = 0
                kk = k
                while kk % p == 0:
                    kk //= p
                    a += 1
                sn = sk * (p ** (a + 2) - 1) // (p ** (a + 1) - 1) - n
            else:
                sn = sk * (p + 1) - n
            if 1 <= sn <= N:
                touched[sn] = True
                if sn % 2 == 1:
                    nonsemi[sn] = True
    return touched, nonsemi, pairs, is_prime[: N + 1], distinct


touched, nonsemi, pairs, is_prime_N, distinct = touchable_set(N, C)
unt = np.nonzero(~touched[2:])[0] + 2           # untouchable m in [2, N]
print("N = %d, sieve limit L = %d, composite-case pairs = %d, time %.1fs" % (N, L, pairs, time.time() - t0))
print("untouchables <= N:", len(unt))
first = unt[:40].tolist()
print("first 40:", first)
A005114_first = [2, 5, 52, 88, 96, 120, 124, 146, 162, 188, 206, 210, 216, 238, 246, 248, 262, 268, 276, 288,
                 290, 292, 304, 306, 322, 324, 326, 336, 342, 372, 406, 408, 426, 430, 448, 472, 474, 498, 516, 518]
check("U1 first 40 untouchables equal A005114", first == A005114_first, first[:12])
A057978 = {10 ** 1: 2, 10 ** 2: 5, 10 ** 3: 89, 10 ** 4: 1212, 10 ** 5: 13863, 10 ** 6: 150232, 10 ** 7: 1574973}
for X in sorted(A057978):
    if X <= N:
        cnt = int(np.sum(unt <= X))
        check("U2 count of untouchables <= %d equals A057978 (%d)" % (X, A057978[X]), cnt == A057978[X], cnt)
odd_unt = unt[unt % 2 == 1].tolist()
check("U3 the only odd untouchable <= N is 5", odd_unt == [5], odd_unt[:10])

# independent brute force for a small range: all n <= (M-1)^2 (complete by the k <= M-1 argument)
M = 3000
sig_small = sigma_sieve((M - 1) ** 2)
n_small = np.arange((M - 1) ** 2 + 1, dtype=np.int64)
s_small = sig_small - n_small
touched_small = np.zeros(M + 1, dtype=bool)
msk = (n_small >= 2) & (s_small >= 1) & (s_small <= M)
touched_small[s_small[msk]] = True
unt_small = (np.nonzero(~touched_small[2:])[0] + 2).tolist()
check("U4 independent brute force (all n <= (M-1)^2, M=3000) agrees with the case-split census",
      unt_small == unt[unt <= M].tolist(), (len(unt_small), int(np.sum(unt <= M))))

# Goldbach dependence of odd touchables
odd_m = np.arange(3, N + 1, 2)
gold_dep = odd_m[(~nonsemi[odd_m]) & touched[odd_m]]          # odd, touchable, no non-semiprime source
has_gold = np.zeros(N + 1, dtype=bool)
has_gold[1:] = distinct[:-1] > 0
odd_touch = odd_m[touched[odd_m]]
print("odd m in [3,N]: %d; touchable: %d; with a non-semiprime source: %d; Goldbach-dependent (only p*q sources): %d"
      % (len(odd_m), len(odd_touch), int(np.sum(nonsemi[odd_m])), len(gold_dep)))
print("share of odd touchables that are Goldbach-dependent: %.4f" % (len(gold_dep) / len(odd_touch)))
print("first 30 Goldbach-dependent odd m:", gold_dep[:30].tolist())
print("largest 10 Goldbach-dependent odd m <= N:", gold_dep[-10:].tolist())
# every odd touchable m >= 9 is covered by Goldbach (distinct) in this range?
odd_ge9 = odd_m[odd_m >= 9]
check("U5 every odd m in [9, N] has m-1 = p+q with p != q primes (distinct Goldbach, this range only)",
      bool(np.all(has_gold[odd_ge9])))
check("U6 odd m in {1,3,7} touchable (s(p)=1, s(4)=3, s(8)=7) and 5 not", touched[1] and touched[3] and touched[7] and not touched[5])

# residue tables of untouchables
for q in (2, 3, 4, 5, 6, 7, 12, 14, 15):
    counts = [int(np.sum(unt % q == r)) for r in range(q)]
    tot = len(unt)
    print("untouchables <= N by residue mod %2d: %s" % (q, counts))
# even untouchables mod 4 and mod 6 shares
ev = unt[unt % 2 == 0]
print("even untouchables: %d; == 0 mod 4: %d; == 2 mod 4: %d; == 0 mod 6: %d; == 2 mod 6: %d; == 4 mod 6: %d"
      % (len(ev), int(np.sum(ev % 4 == 0)), int(np.sum(ev % 4 == 2)), int(np.sum(ev % 6 == 0)), int(np.sum(ev % 6 == 2)), int(np.sum(ev % 6 == 4))))
print("density of untouchables <= N: %.5f" % (len(unt) / N))

# Chowla shift: s'-untouchables = s-untouchables - 1
chowla_unt = (unt - 1).tolist()
print("first 20 Chowla-untouchables (m not of the form sigma(n)-n-1):", chowla_unt[:20])
check("U7 Chowla shift: m is s'-touchable iff m+1 is s-touchable (definitional; 4 = 5-1 is the only even one <= N-1)",
      [m for m in chowla_unt if m % 2 == 0] == [4])

bad = [n for n, ok in CHECKS if not ok]
print("%d checks, %d failed; total time %.1fs" % (len(CHECKS), len(bad), time.time() - t0))
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL CHECKS PASSED")
