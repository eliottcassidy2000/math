#!/usr/bin/env python3
"""The Pollack-Pomerance constant Delta, evaluated precisely (opus-2026-10-26-S57).

PP Conjecture 1.4 (Trans. AMS Ser. B 3 (2016)): the nonaliquot numbers have density
    Delta = lim_{y->oo} (1/log y) sum_{a<=y, 2|a} a^{-1} e^{-a/s(a)}.
A log-average equals the limit of the natural-density averages when the latter exists, so
    Delta = (1/2) E[ e^{-1/h} ],   h = s(a)/a = sigma(a)/a - 1,  a a natural-density-random even integer.
(1) exact local averages of e^{-a/s(a)} over even a in windows [10^k, 10^(k+1)) up to N, and the cumulative PP
    expression (1/log y) * sum_{a<=y even} e^{-a/s(a)}/a at y = 10^k (all even a >= 2 included);
(2) the limit E[e^{-1/h}] from the independent-valuation model: v_2(a) - 1 ~ Geom(1/2), v_p ~ P(v_p = j) = (1-1/p)p^-j
    for odd p (Davenport / Erdos-Wintner: sigma(a)/a has this limiting law); Monte Carlo with exact small primes.
(3) the empirical local untouchable density per window from mult_even_N.npy, for comparison.
Usage: python pp_delta_20261026.py N bigdir [mc_samples]
"""
import sys, os, json, math, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from untouchable_sieve_20261026 import primes_upto
from untouchable_profile_20261026 import s_even_segment


@nb.njit(cache=False)
def window_sums(lo, cnt, sv, mult, edges, acc_f, acc_n, acc_u, acc_harm):
    for i in range(cnt):
        k = lo + 2 * i
        s = sv[i]
        f = math.exp(-k / s) if s > 0 else 0.0
        acc_harm[0] += f / k
        # window index: largest e with edges[e] <= k
        e = 0
        while e + 1 < edges.shape[0] and edges[e + 1] <= k:
            e += 1
        acc_f[e] += f
        acc_n[e] += 1
        if mult[k // 2] == 0:
            acc_u[e] += 1


@nb.njit(cache=False)
def mc_model(nsamp, primes, seed):
    np.random.seed(seed)
    tot = 0.0
    tot2 = 0.0
    for t in range(nsamp):
        # p = 2: v >= 1, P(v = j) = 2^-j
        v = 1
        while np.random.random() < 0.5:
            v += 1
        r = 2.0 - 0.5 ** v          # sigma(2^v)/2^v = 2 - 2^-v
        for q in range(1, primes.shape[0]):
            p = primes[q]
            # P(v_p >= 1) = 1/p ; then geometric: P(v_p >= j+1 | v_p >= j) = 1/p
            if np.random.random() < 1.0 / p:
                vp = 1
                while np.random.random() < 1.0 / p:
                    vp += 1
                r *= (1.0 - (1.0 / p) ** (vp + 1)) / (1.0 - 1.0 / p)
        h = r - 1.0
        f = math.exp(-1.0 / h)
        tot += f
        tot2 += f * f
    m = tot / nsamp
    return m, math.sqrt(max(tot2 / nsamp - m * m, 0.0) / nsamp)


if __name__ == '__main__':
    N = int(float(sys.argv[1]))
    bigdir = sys.argv[2]
    mcs = int(float(sys.argv[3])) if len(sys.argv) > 3 else 20000000
    out = {'N': N}
    mult = np.load(os.path.join(bigdir, f'mult_even_{N}.npy'))
    primes = primes_upto(int(math.isqrt(2 * N)) + 2)
    edges = [2]
    x = 10
    while x <= N:
        edges.append(x)
        x *= 10
    edges = np.array(edges, np.int64)
    acc_f = np.zeros(len(edges))
    acc_n = np.zeros(len(edges), np.int64)
    acc_u = np.zeros(len(edges), np.int64)
    acc_harm = np.zeros(1)
    cum = {}
    seg = 1 << 22
    lo = 2
    t0 = time.time()
    while lo <= N:
        cnt = min(seg, (N - lo) // 2 + 1)
        # stop segments exactly at powers of ten for the cumulative expression
        nxt = None
        for e in edges[1:]:
            if lo < e <= lo + 2 * (cnt - 1):
                nxt = e
                break
        if nxt is not None:
            cnt = (nxt - lo) // 2 + 1
        sv = np.empty(cnt, np.int64)
        s_even_segment(lo, cnt, primes, sv)
        window_sums(lo, cnt, sv, mult, edges, acc_f, acc_n, acc_u, acc_harm)
        hi = lo + 2 * (cnt - 1)
        if hi in set(int(e) for e in edges[1:]):
            cum[int(hi)] = acc_harm[0] / math.log(hi)
        lo += 2 * cnt
    out['cumulative_pp_expression'] = cum
    print('cumulative PP expression (1/log y) sum_{even a <= y} e^{-a/s(a)}/a:',
          {k: round(v, 6) for k, v in cum.items()}, f'({time.time() - t0:.0f}s)', flush=True)
    rows = []
    for e in range(len(edges)):
        if acc_n[e] == 0:
            continue
        loc_pp = acc_f[e] / acc_n[e] / 2
        loc_un = acc_u[e] / acc_n[e] / 2
        rows.append((int(edges[e]), int(acc_n[e]), loc_pp, loc_un))
        print(f'   window from {edges[e]:>11}: n_even={acc_n[e]:>11}  local PP (1/2)E[e^-a/s(a)] = {loc_pp:.6f}   '
              f'local untouchable density (even untouchables / all integers) = {loc_un:.6f}', flush=True)
    out['windows'] = rows
    t0 = time.time()
    small = primes_upto(2000)
    m, se = mc_model(mcs, small, 12345)
    out['model_E_exp_minus_1_over_h'] = m
    out['model_se'] = se
    out['Delta_model'] = m / 2
    print(f'independent-valuation model (primes <= 2000, {mcs} samples): E[e^-1/h] = {m:.6f} +- {se:.6f}  =>  '
          f'Delta = {m / 2:.6f} +- {se / 2:.6f}  ({time.time() - t0:.0f}s)', flush=True)
    json.dump(out, open(os.path.join(HERE, f'pp_delta_{N}.json'), 'w'), indent=1)
