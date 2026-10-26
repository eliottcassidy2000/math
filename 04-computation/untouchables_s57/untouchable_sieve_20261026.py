#!/usr/bin/env python3
"""Exact untouchable numbers up to N and the multiplicity |s^{-1}(k)| of every even k <= N (opus-2026-10-26-S57).

s(m) = sigma(m) - m.  Parity: s(m) is odd iff (m odd and not a square) or (m even and m a square or twice a square).
So an even k >= 2 has preimages only among
  (i)  even m with s(m) = k; since m/2 is a proper divisor, s(m) >= m/2 + 1, so m <= 2k - 2;
  (ii) odd squares q^2: s(q^2) even.  q prime gives s = q + 1 (k - 1 prime); q composite gives s(q^2) > q^2/p_min >= q^1.5,
       so only q <= N^(2/3) matter.
Odd k: k = 1 + p + q = s(pq) for distinct primes p, q whenever k - 1 is a sum of two distinct primes (verified strong
Goldbach for even numbers up to 4e18, Oliveira e Silva et al.), so the odd untouchables <= N are exactly {5} (checked
independently here for k <= N by a small-prime search).
Usage: python untouchable_sieve_20261026.py N [segment]
Writes mult_even_N.npy (uint8 multiplicities, capped at 255, index j <-> k = 2j) and a summary.
"""
import sys, os, time, math, json
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))


@nb.njit(cache=False)
def primes_upto(L):
    s = np.ones(L + 1, np.uint8)
    s[:2] = 0
    for i in range(2, int(L ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = 0
    return np.nonzero(s)[0].astype(np.int64)


@nb.njit(cache=False)
def even_segment(lo, cnt, primes, N, mult):
    """even m = lo + 2i, i < cnt (lo even).  Increment mult[s/2] for s = s(m) even, s <= N."""
    rem = np.empty(cnt, np.int64)
    sig = np.ones(cnt, np.int64)
    for i in range(cnt):
        rem[i] = lo + 2 * i
    # p = 2
    for i in range(cnt):
        pk = 1
        s = 1
        while rem[i] % 2 == 0:
            rem[i] //= 2
            pk *= 2
            s += pk
        sig[i] = s
    hi = lo + 2 * cnt
    for t in range(1, primes.shape[0]):
        p = primes[t]
        if p * p >= hi:
            break
        # first m >= lo with m = 0 mod 2p
        step = 2 * p
        first = ((lo + step - 1) // step) * step
        for m in range(first, hi, step):
            i = (m - lo) // 2
            pk = 1
            s = 1
            while rem[i] % p == 0:
                rem[i] //= p
                pk *= p
                s += pk
            sig[i] *= s
    for i in range(cnt):
        m = lo + 2 * i
        if m < 2:
            continue
        r = rem[i]
        sg = sig[i]
        if r > 1:
            sg *= (r + 1)
        sv = sg - m
        if sv >= 2 and sv <= N and (sv & 1) == 0:
            j = sv // 2
            if mult[j] < 255:
                mult[j] += 1


@nb.njit(cache=False)
def odd_prime_squares(N, mult):
    """k = q + 1 for odd primes q <= N - 1 (preimage q^2): segmented sieve of primes up to N."""
    L = int(N ** 0.5) + 2
    base = primes_upto(L)
    seg = 1 << 24
    lo = 3
    while lo < N:
        hi = min(lo + seg, N)
        flags = np.ones(hi - lo, np.uint8)
        for t in range(base.shape[0]):
            p = base[t]
            if p * p >= hi:
                break
            st = max(p * p, ((lo + p - 1) // p) * p)
            for m in range(st, hi, p):
                flags[m - lo] = 0
        for i in range(hi - lo):
            q = lo + i
            if flags[i] and (q & 1) == 1:
                k = q + 1
                if k <= N:
                    j = k // 2
                    if mult[j] < 255:
                        mult[j] += 1
        lo = hi


@nb.njit(cache=False)
def sigma_small(n):
    s = 1
    x = n
    p = 2
    while p * p <= x:
        if x % p == 0:
            pk = 1
            t = 1
            while x % p == 0:
                x //= p
                pk *= p
                t += pk
            s *= t
        p += 1
    if x > 1:
        s *= (1 + x)
    return s


@nb.njit(cache=False)
def odd_composite_squares(N, mult):
    Q = int(N ** (2.0 / 3.0)) + 2
    for q in range(9, Q + 1, 2):
        # composite odd q only
        isprime = True
        d = 3
        while d * d <= q:
            if q % d == 0:
                isprime = False
                break
            d += 2
        if isprime:
            continue
        # sigma(q^2) via factorisation of q
        s = 1
        x = q
        p = 3
        while p * p <= x:
            if x % p == 0:
                e = 0
                while x % p == 0:
                    x //= p
                    e += 1
                t = 1
                pk = 1
                for _ in range(2 * e):
                    pk *= p
                    t += pk
                s *= t
            p += 2
        if x > 1:
            s *= (1 + x + x * x)
        sv = s - q * q
        if sv >= 2 and sv <= N and (sv & 1) == 0:
            j = sv // 2
            if mult[j] < 255:
                mult[j] += 1


@nb.njit(cache=False)
def hist_u8(a):
    h = np.zeros(256, np.int64)
    for i in range(a.shape[0]):
        h[a[i]] += 1
    return h


def run(N, seg=1 << 22):
    t0 = time.time()
    primes = primes_upto(int(math.isqrt(2 * N)) + 2)
    mult = np.zeros(N // 2 + 1, np.uint8)
    M = 2 * N - 2
    lo = 2
    while lo <= M:
        cnt = min(seg, (M - lo) // 2 + 1)
        even_segment(lo, cnt, primes, N, mult)
        lo += 2 * cnt
    odd_prime_squares(N, mult)
    odd_composite_squares(N, mult)
    return mult, time.time() - t0


def odd_untouchables(N):
    """Odd k <= N not touched by s(pq) = 1+p+q (distinct primes), s(4)=3, s(8)=7, s(p)=1: search small p."""
    L = min(N, 10 ** 7)
    pr = primes_upto(L)
    isp = np.zeros(L + 1, np.uint8)
    isp[pr] = 1
    bad = []
    for k in range(1, L + 1, 2):
        if k in (1, 3, 7):
            continue
        n = k - 1
        ok = False
        for p in pr:
            if 2 * p >= n:
                break
            if isp[n - p] and n - p != p:
                ok = True
                break
        if not ok:
            bad.append(k)
    return bad, L


if __name__ == '__main__':
    N = int(float(sys.argv[1]))
    mult, dt = run(N)
    bc = hist_u8(mult[1:])                      # j >= 1 <-> k = 2j >= 2 (no int64 temporaries at N = 1e10)
    even_untouch = int(bc[0])
    oddbad, L = odd_untouchables(N)
    total = even_untouch + len(oddbad)
    print(f'N = {N}: even untouchables {even_untouch}, odd untouchables (k <= {L} checked) {oddbad}, '
          f'total untouchables {total} ({total / N:.6f} of [1,N]); sieve {dt:.0f}s', flush=True)
    big = os.environ.get('S57_BIGDIR', HERE)   # large arrays go to the scratchpad, not the repo
    np.save(os.path.join(big, f'mult_even_{N}.npy'), mult)
    # counts at powers of ten (and 2^e) for comparison with OEIS A070015
    marks = {}
    x = 10
    while x <= N:
        marks[str(x)] = int(np.count_nonzero(mult[1:x // 2 + 1] == 0)) + (1 if x >= 5 else 0)
        x *= 10
    print('untouchables <= 10^e:', marks, flush=True)
    dist = {int(v): int(c) for v, c in enumerate(bc) if c}
    json.dump({'N': N, 'counts_at_powers_of_ten': marks, 'even_untouchables': even_untouch, 'odd_untouchables_checked_to': L, 'odd_untouchables': oddbad,
               'total': total, 'share': total / N, 'multiplicity_distribution_even_k': dist},
              open(os.path.join(HERE, f'untouchables_{N}.json'), 'w'), indent=1)
