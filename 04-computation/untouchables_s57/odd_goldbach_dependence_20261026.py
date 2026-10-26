#!/usr/bin/env python3
"""How much of odd touchability rests on Goldbach alone, to 10^8 (opus-2026-10-26-S57; extends kind-pasteur's 10^6
statistic, untouchable_census_goldbach_dependence_kps_20261025.py).

Odd k >= 3 has preimages only among odd non-squares and even squares / twice squares (parity of sigma).  Split:
  (S) semiprimes n = pq, p < q odd primes: s(pq) = 1 + p + q, i.e. k - 1 is a sum of two distinct odd primes;
  (O) every other source: odd non-square n = p*m with p = spf(n) <= spf(m), m odd composite (then s(n) >= 1 + m, so
      m <= k - 1); even n = u^2 or 2u^2 (s(n) >= n/2 + 1, so n <= 2k - 2); n = 4 (s = 3).
k is Goldbach-dependent if it is touched by (S) and not by (O).  We count, for N = 10^6, 10^7, 10^8: odd touchables,
those touched by (O), the Goldbach-dependent ones, and those whose only source is a single semiprime.
Usage: python odd_goldbach_dependence_20261026.py N
"""
import sys, os, json, math, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))


@nb.njit(cache=False)
def spf_sieve(L):
    spf = np.zeros(L + 1, np.uint32)
    for i in range(2, L + 1):
        if spf[i] == 0:
            spf[i] = i
            if i <= L // i:
                for j in range(i * i, L + 1, i):
                    if spf[j] == 0:
                        spf[j] = i
    return spf


@nb.njit(cache=False)
def sigma_from_spf(m, spf):
    s = 1
    while m > 1:
        p = np.int64(spf[m])
        pk = 1
        t = 1
        while m % p == 0:
            m //= p
            pk *= p
            t += pk
        s *= t
    return s


@nb.njit(cache=False)
def other_sources(N, spf, cov):
    """cov[k] += 1 for odd k <= N hit by a non-semiprime source (O)."""
    # odd non-square n = p*m, m odd composite, p prime <= spf(m)
    for m in range(9, N, 2):
        q = np.int64(spf[m])
        if q == m:
            continue                      # m prime
        sm = sigma_from_spf(m, spf)
        # exponent of q in m and sigma(q^e)
        e = 0
        x = m
        while x % q == 0:
            x //= q
            e += 1
        sq = 0
        pk = 1
        for _ in range(e + 1):
            sq += pk
            pk *= q
        sq_next = sq + pk                   # sigma(q^(e+1))
        # p < q: s(pm) = sigma(m) + p s(m) increases with p, so the loop may stop at the first k > N
        p = 3
        while p < q:
            if spf[p] == p:
                k = (1 + p) * sm - p * m
                if k > N:
                    break
                if (k & 1) == 1 and k >= 3:
                    cov[k] += 1
            p += 2
        # p = q is separate: s(qm) can be smaller than s(pm) for p < q (the hostile n = 7^2 * 229, k = 1889)
        k = (sm // sq) * sq_next - q * m
        if 3 <= k <= N and (k & 1) == 1:      # odd k <=> n = q*m is not a square
            cov[k] += 1
    # even squares and twice squares, n <= 2N
    u = 1
    while 2 * u * u <= 2 * N:
        for n in (2 * u * u, 4 * u * u):
            if n <= 2 * N and n >= 4:
                k = sigma_from_spf(n, spf) - n if n <= spf.shape[0] - 1 else -1
                if 3 <= k <= N and (k & 1) == 1:
                    cov[k] += 1
        u += 1


@nb.njit(cache=False)
def semiprime_reps(N, spf, rep):
    """rep[k] = number of representations k - 1 = p + q, p < q odd primes (capped at 2)."""
    primes = []
    for i in range(3, N + 1, 2):
        if spf[i] == i:
            primes.append(i)
    P = np.array(primes, np.int64)
    for k in range(9, N + 1, 2):
        e = k - 1
        c = 0
        for t in range(P.shape[0]):
            p = P[t]
            if 2 * p >= e:
                break
            if spf[e - p] == e - p:
                c += 1
                if c >= 2:
                    break
        rep[k] = c


if __name__ == '__main__':
    N = int(float(sys.argv[1]))
    t0 = time.time()
    spf = spf_sieve(2 * N + 2)
    cov = np.zeros(N + 1, np.int32)
    other_sources(N, spf, cov)
    rep = np.zeros(N + 1, np.int8)
    semiprime_reps(N, spf, rep)
    out = {'N': N}
    for M in sorted(set([10 ** 6, 10 ** 7, 10 ** 8, N])):
        if M > N:
            continue
        ks = np.arange(3, M + 1, 2)
        o = cov[ks] > 0
        s = rep[ks] > 0
        touch = o | s | (ks == 3) | (ks == 7)
        dep = s & ~o
        single = (rep[ks] == 1) & ~o
        unt = ks[~touch]
        row = {'odd_touchable': int(touch.sum()), 'touched_by_other_sources': int(o.sum()),
               'goldbach_dependent': int(dep.sum()), 'share_of_odd_touchables': float(dep.sum() / touch.sum()),
               'single_semiprime_only': int(single.sum()), 'single_semiprime_only_list': [int(x) for x in ks[single][:20]],
               'odd_untouchables': [int(x) for x in unt[:10]],
               'largest_goldbach_dependent': int(ks[dep].max()) if dep.any() else None}
        out[str(M)] = row
        print(f'odd k <= {M}: touchable {row["odd_touchable"]}, by non-semiprime sources {row["touched_by_other_sources"]}, '
              f'Goldbach-dependent {row["goldbach_dependent"]} ({100 * row["share_of_odd_touchables"]:.3f}%), '
              f'single-semiprime-only {row["single_semiprime_only"]}, odd untouchables {row["odd_untouchables"]} '
              f'({time.time() - t0:.0f}s)', flush=True)
    json.dump(out, open(os.path.join(HERE, f'odd_goldbach_dependence_{N}.json'), 'w'), indent=1)
