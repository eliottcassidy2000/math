#!/usr/bin/env python3
"""How often Pollack-Pomerance's key assumption fails at computable sizes (opus-2026-10-26-S57).

PP sec. 3.1 assume that s maps T_a = {n : gcd(n, A_y) = a} into T_a, "asymptotically true if n > e^(e^y)".  For
n = a m (m coprime to A_y), s(n) = sigma(a) sigma(m) - a m lies in T_a when every small prime power divides sigma(m)
to high order; the assumption fails at the prime l whenever l does not divide sigma(n).  Landau-Wirsing: the n <= X
with l not dividing sigma(n) have relative density ~ C_l (log X)^(-1/(l-1)) (no prime q = -1 mod l to an odd power,
...).  We tabulate f_l(X) = share of n in (X/2, X] (all n, and even n) with l not dividing sigma(n), l = 3, 5, 7, 11, 13,
for X = 2^10 .. 2^E, and the fitted exponent of log X.
Usage: python pp_assumption_failure_20261026.py [E]
"""
import sys, os, json, math, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from untouchable_sieve_20261026 import primes_upto

LS = np.array([3, 5, 7, 11, 13], np.int64)


@nb.njit(cache=False)
def sigma_segment(lo, cnt, primes, out):
    rem = np.empty(cnt, np.int64)
    for i in range(cnt):
        rem[i] = lo + i
        out[i] = 1
    hi = lo + cnt
    for t in range(primes.shape[0]):
        p = primes[t]
        if p * p >= hi:
            break
        first = ((lo + p - 1) // p) * p
        for m in range(first, hi, p):
            i = m - lo
            pk = 1
            s = 1
            while rem[i] % p == 0:
                rem[i] //= p
                pk *= p
                s += pk
            out[i] *= s
    for i in range(cnt):
        if rem[i] > 1:
            out[i] *= (rem[i] + 1)


@nb.njit(cache=False)
def tally(lo, cnt, sg, LS, tot, bad, tot_e, bad_e):
    for i in range(cnt):
        n = lo + i
        ev = (n & 1) == 0
        tot[0] += 1
        if ev:
            tot_e[0] += 1
        for t in range(LS.shape[0]):
            if sg[i] % LS[t] != 0:
                bad[t] += 1
                if ev:
                    bad_e[t] += 1


if __name__ == '__main__':
    E = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    primes = primes_upto(int(math.isqrt(2 ** E)) + 2)
    rows = []
    t0 = time.time()
    for e in range(10, E + 1):
        lo, hi = 2 ** (e - 1) + 1, 2 ** e
        tot = np.zeros(1, np.int64); bad = np.zeros(LS.shape[0], np.int64)
        tot_e = np.zeros(1, np.int64); bad_e = np.zeros(LS.shape[0], np.int64)
        L = lo
        while L <= hi:
            cnt = min(1 << 22, hi - L + 1)
            sg = np.empty(cnt, np.int64)
            sigma_segment(L, cnt, primes, sg)
            tally(L, cnt, sg, LS, tot, bad, tot_e, bad_e)
            L += cnt
        f = [bad[t] / tot[0] for t in range(LS.shape[0])]
        fe = [bad_e[t] / tot_e[0] for t in range(LS.shape[0])]
        rows.append({'e': e, 'f_all': f, 'f_even': fe})
        print(f'2^{e:>2}: share with l not dividing sigma(n), l = 3,5,7,11,13: all ' + ' '.join(f'{x:.4f}' for x in f) +
              '   even ' + ' '.join(f'{x:.4f}' for x in fe) + f'  ({time.time() - t0:.0f}s)', flush=True)
    # fitted exponent of log X between the first and last windows
    fit = {}
    for t, l in enumerate(LS):
        a, b = rows[0], rows[-1]
        la, lb = math.log(2 ** a['e']), math.log(2 ** b['e'])
        fit[int(l)] = {'fitted_exponent_all': math.log(b['f_all'][t] / a['f_all'][t]) / math.log(lb / la),
                       'landau_exponent': -1 / (l - 1)}
    print('fitted exponent of log X (first to last window) vs -1/(l-1):', fit)
    json.dump({'rows': rows, 'fit': fit}, open(os.path.join(HERE, 'pp_assumption_failure.json'), 'w'), indent=1)
