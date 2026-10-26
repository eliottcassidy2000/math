#!/usr/bin/env python3
"""Two different departures of the real map s from Pollack-Pomerance's random-map model (opus-2026-10-26-S57,
written after the independent audit).

(a) PP's literal assumption: s maps T_a = {n : gcd(n, A_y) = a} into T_a, i.e. gcd(s(n), A_y) = gcd(n, A_y).
    At a prime l with l^e || A_y it fails iff min(v_l(s(n)), e) != min(v_l(n), e).  For l not dividing n this is
    l | s(n), i.e. sigma(n) = n (mod l).
(b) Residue neutrality: the random-map model treats the nonzero residues of s(n) mod l alike.  The real map does so
    on inputs with l | sigma(n) (there s(n) = -n mod l); the complementary share is f_l = share(l not dividing
    sigma(n)) ~ C_l (log X)^(-1/(l-1)) (Ramanujan-Watson, Rankin, Scourfield).
For n in a dyadic window (all n, and even n): per-prime failure of (a) (capped exponents for y = 13 and uncapped),
joint failure of (a) at y = 13 and y = 53, and f_l.
Usage: python pp_class_preservation_20261026.py e
"""
import sys, os, json, math, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from untouchable_sieve_20261026 import primes_upto
from pp_assumption_failure_20261026 import sigma_segment

PS = np.array([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53], np.int64)


@nb.njit(cache=False)
def vl(x, l):
    c = 0
    while x % l == 0 and x != 0:
        x //= l
        c += 1
    return c


@nb.njit(cache=False)
def tally(lo, cnt, sg, PS, e13, e53, out_all, out_even):
    # columns: [count, per-prime uncapped fail x16, per-prime capped(y=13) fail x6, joint13, joint53, f_l (l not | sigma) x16]
    for i in range(cnt):
        n = lo + i
        s = sg[i] - n
        ev = (n & 1) == 0
        row = np.zeros(1 + 16 + 6 + 2 + 16, np.int64)
        row[0] = 1
        j13 = False
        j53 = False
        for t in range(PS.shape[0]):
            l = PS[t]
            a = vl(n, l)
            b = vl(s, l) if s != 0 else 99
            if a != b:
                row[1 + t] = 1
            if t < 6:
                ec = e13[t]
                if min(a, ec) != min(b, ec):
                    row[17 + t] = 1
                    j13 = True
            ec5 = e53[t]
            if min(a, ec5) != min(b, ec5):
                j53 = True
            if sg[i] % l != 0:
                row[25 + t] = 1
        row[23] = 1 if j13 else 0
        row[24] = 1 if j53 else 0
        for c in range(row.shape[0]):
            out_all[c] += row[c]
            if ev:
                out_even[c] += row[c]


if __name__ == '__main__':
    e = int(sys.argv[1]) if len(sys.argv) > 1 else 27
    e13 = np.array([3, 2, 1, 1, 1, 1], np.int64)
    e53 = np.array([int(math.floor(math.log(53) / math.log(p) + 1e-12)) for p in PS], np.int64)
    primes = primes_upto(int(math.isqrt(2 ** e)) + 2)
    lo, hi = 2 ** (e - 1) + 1, 2 ** e
    out_all = np.zeros(41, np.int64)
    out_even = np.zeros(41, np.int64)
    t0 = time.time()
    L = lo
    while L <= hi:
        cnt = min(1 << 22, hi - L + 1)
        sg = np.empty(cnt, np.int64)
        sigma_segment(L, cnt, primes, sg)
        tally(L, cnt, sg, PS, e13, e53, out_all, out_even)
        L += cnt
    res = {'window': [lo, hi]}
    for name, o in (('all', out_all), ('even', out_even)):
        N = o[0]
        res[name] = {
            'per_prime_uncapped_fail': {int(PS[t]): o[1 + t] / N for t in range(16)},
            'per_prime_capped_y13_fail': {int(PS[t]): o[17 + t] / N for t in range(6)},
            'joint_fail_y13': o[23] / N, 'joint_fail_y53': o[24] / N,
            'f_l_not_dividing_sigma': {int(PS[t]): o[25 + t] / N for t in range(16)}}
        r = res[name]
        print(f'{name} n in (2^{e - 1}, 2^{e}]: joint failure of gcd(s(n),A_y) = gcd(n,A_y): y=13 {r["joint_fail_y13"]:.3f}, '
              f'y=53 {r["joint_fail_y53"]:.3f}', flush=True)
        print('   per prime l = 2,3,5,7,11,13: literal failure (capped, y=13) ' +
              ' '.join(f'{r["per_prime_capped_y13_fail"][int(l)]:.3f}' for l in PS[:6]) +
              ';  uncapped ' + ' '.join(f'{r["per_prime_uncapped_fail"][int(l)]:.3f}' for l in PS[:6]) +
              ';  f_l (l not | sigma) ' + ' '.join(f'{r["f_l_not_dividing_sigma"][int(l)]:.3f}' for l in PS[:6]), flush=True)
    print(f'({time.time() - t0:.0f}s)')
    json.dump(res, open(os.path.join(HERE, f'pp_class_preservation_2e{e}.json'), 'w'), indent=1)
