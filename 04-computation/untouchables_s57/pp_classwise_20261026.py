#!/usr/bin/env python3
"""Class-by-class test of the Pollack-Pomerance model (opus-2026-10-26-S57).

PP Section 3.1: for n in T_a = {n : gcd(n, A_y) = a} (a even), the number of preimages is modelled as Poisson with mean
lambda_a = a/s(a), so P(untouchable) = exp(-a/s(a)).  We test, for even k in [N/10, N] grouped by a = gcd(k, A) with
A = lcm(1..10) = 2520:  the empirical mean multiplicity versus a/s(a); the empirical zero share versus exp(-a/s(a)) and
versus exp(-empirical mean) (separating a wrong mean from a non-Poisson law); also the variance/mean.
Usage: python pp_classwise_20261026.py N bigdir
"""
import sys, os, json, math
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
A = 2520


def sigma(n):
    s, d = 0, 1
    while d * d <= n:
        if n % d == 0:
            s += d
            if d * d != n:
                s += n // d
        d += 1
    return s


@nb.njit(cache=False)
def classwise(mult, lo_j, hi_j, A, idx_of, cnt, sm, sm2, zr):
    for j in range(lo_j, hi_j):
        k = 2 * j
        a = k
        b = A
        while b:
            a, b = b, a % b
        c = idx_of[a]
        m = mult[j]
        cnt[c] += 1
        sm[c] += m
        sm2[c] += m * m
        if m == 0:
            zr[c] += 1


if __name__ == '__main__':
    N = int(float(sys.argv[1]))
    mult = np.load(os.path.join(sys.argv[2], f'mult_even_{N}.npy'))
    divs = [d for d in range(2, A + 1, 2) if A % d == 0]
    idx_of = np.full(A + 1, -1, np.int64)
    for i, d in enumerate(divs):
        idx_of[d] = i
    nc = len(divs)
    cnt = np.zeros(nc, np.int64)
    sm = np.zeros(nc, np.float64)
    sm2 = np.zeros(nc, np.float64)
    zr = np.zeros(nc, np.int64)
    classwise(mult, N // 20, N // 2 + 1, A, idx_of, cnt, sm, sm2, zr)
    rows = []
    tot = cnt.sum()
    pp_total = 0.0
    emp_total = 0.0
    for i, d in enumerate(divs):
        if cnt[i] == 0:
            continue
        lam_pp = d / (sigma(d) - d)
        mean = sm[i] / cnt[i]
        var = sm2[i] / cnt[i] - mean ** 2
        z = zr[i] / cnt[i]
        rows.append({'a': d, 'n': int(cnt[i]), 'weight': cnt[i] / tot, 'lambda_PP': lam_pp, 'mean': mean,
                     'var_over_mean': var / mean if mean > 0 else None, 'zero': z, 'zero_PP': math.exp(-lam_pp),
                     'zero_poisson_emp_mean': math.exp(-mean)})
        pp_total += cnt[i] * math.exp(-lam_pp)
        emp_total += zr[i]
    rows.sort(key=lambda r: -r['weight'])
    print(f'even k in [{N // 10}, {N}], classes a = gcd(k, {A}); top classes by weight:')
    print('     a   weight  lambda_PP  mean_mult  var/mean  zero_emp  zero_PP  e^-mean')
    for r in rows[:24]:
        print(f'{r["a"]:>6} {r["weight"]:8.4f} {r["lambda_PP"]:9.4f} {r["mean"]:10.4f} {r["var_over_mean"]:9.3f} '
              f'{r["zero"]:9.4f} {r["zero_PP"]:8.4f} {r["zero_poisson_emp_mean"]:8.4f}', flush=True)
    wm_pp = sum(r['weight'] * r['lambda_PP'] for r in rows)
    wm_emp = sum(r['weight'] * r['mean'] for r in rows)
    print(f'weighted mean multiplicity: PP {wm_pp:.4f}, empirical {wm_emp:.4f}')
    print(f'even zero share: PP class model {pp_total / tot:.4f}, empirical {emp_total / tot:.4f}')
    corr_mean = np.corrcoef([r['lambda_PP'] for r in rows], [r['mean'] for r in rows])[0, 1]
    print(f'correlation across classes of lambda_PP with empirical mean: {corr_mean:.4f}')
    json.dump({'N': N, 'A': A, 'rows': rows, 'weighted_mean_PP': wm_pp, 'weighted_mean_emp': wm_emp,
               'zero_share_PP_class_model': pp_total / tot, 'zero_share_emp': emp_total / tot},
              open(os.path.join(HERE, f'pp_classwise_{N}.json'), 'w'), indent=1)
