#!/usr/bin/env python3
"""Exact incoming-Q identity audit and the conditional covariance/Pell model.

No sampling hypothesis for Q, or Collatz rate, is inferred from these checks.
The finite Q universe is 1..2^16; the covariance model has rational controls
in [-1,1], endpoint controls and an inadmissible-input hostile control.
"""
import argparse
import json
from fractions import Fraction as F


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# Exact Q(sqrt(2)) arithmetic, independent of any decimal root calculation.
def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def neg(x):
    return (-x[0], -x[1])


def mul(x, y):
    return (x[0] * y[0] + 2 * x[1] * y[1],
            x[0] * y[1] + x[1] * y[0])


def div(x, y):
    norm = y[0] ** 2 - 2 * y[1] ** 2
    require(norm != 0, 'field denominator')
    z = mul(x, (y[0], -y[1]))
    return (z[0] / norm, z[1] / norm)


def scalar(x):
    return (F(x), F(0))


def covariance_step(r):
    require(-1 <= r <= 1, 'correlation outside admissible interval')
    return (r - 1) / (3 - r)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    N = 2 ** 16
    q = [0] * (N + 1)
    q[1] = q[2] = 1
    checks = 0
    d = [F(0)] * (N + 1)
    d[1], d[2] = F(1, 2), F(0)
    for n in range(3, N + 1):
        p, v = n - q[n - 1], n - q[n - 2]
        require(1 <= p < n and 1 <= v < n, 'finite Q source rank')
        q[n] = q[p] + q[v]
        d[n] = F(q[n], n) - F(1, 2)
        require(F(p, n) == F(1, 2) + F(1, 2 * n)
                - (1 - F(1, n)) * d[n - 1], 'p/n correction')
        require(F(v, n) == F(1, 2) + F(1, n)
                - (1 - F(2, n)) * d[n - 2], 'q/n correction')
        source = F(p, n) * d[p] + F(v, n) * d[v]
        require(d[n] + (1 - F(1, n)) * d[n - 1] / 2
                + (1 - F(2, n)) * d[n - 2] / 2
                == source + F(3, 4 * n), 'normalized Q identity')
        require(n + 1 - q[n - 1] == p + 1, 'adjacent source identity')
        checks += 5
    require(F(2, 3) != F(1, 2) - (1 - F(1, 3)) * d[2],
            'old p/n formula hostile at n=3')
    require(q[1536] == 808 and q[1536] != 1024, 'dyadic spike hostile')
    checks += 2

    rm = (F(1), F(-1))
    rp = (F(1), F(1))
    lam = (F(3), F(-2))
    require(mul(rm, rm) == lam, 'fixed correlation squared is Pell inverse')
    require(mul(lam, (F(3), F(2))) == scalar(1), 'Pell norm')
    require(div(add(rm, scalar(-1)), add(scalar(3), neg(rm))) == rm,
            'field fixed point')
    checks += 3
    for numerator in range(-16, 17):
        r = F(numerator, 16)
        nxt = covariance_step(r)
        require(F(-1, 2) <= nxt <= 0, 'invariant correlation interval')
        c = r / 2
        g0 = (3 - r) / 2  # source variance normalized to 1
        g1 = nxt * g0
        g2 = -(g1 + g0) / 2
        require(g0 == -g1 / 2 - g2 / 2 + 1 - c / 2,
                'first covariance equation')
        require(g1 == -g0 / 2 - g1 / 2 + c,
                'second covariance equation')
        require(g2 == -g1 / 2 - g0 / 2, 'third covariance equation')
        # Exact two-step conjugacy to multiplication by 3+2 sqrt(2).
        s = 1 - r
        require(1 - covariance_step(nxt) == (3 * s + 4) / (2 * s + 3),
                'two-step Pell conjugacy')
        cr0 = div(add(scalar(r), neg(rm)), add(scalar(r), neg(rp)))
        cr1 = div(add(scalar(nxt), neg(rm)), add(scalar(nxt), neg(rp)))
        require(cr1 == mul(lam, cr0), 'exact cross-ratio contraction')
        checks += 6
    try:
        covariance_step(F(2))
    except RuntimeError:
        checks += 1
    else:
        raise RuntimeError('inadmissible covariance accepted')

    p, denominator = 0, 1
    r, variance = F(0), F(1)
    reports = []
    for k in range(1, 1025):
        variance *= (3 - r) / 4
        r = covariance_step(r)
        p, denominator = p - denominator, -p + 3 * denominator
        require(r == F(p, denominator), 'matrix versus scalar iteration')
        require(variance == F(denominator, 4 ** k),
                'variance product versus matrix entry')
        checks += 2
        if k in (1, 2, 3, 4, 20, 1024):
            reports.append({'generation': k,
                            'correlation_exact': str(r),
                            'variance_ratio_exact': str(variance)})
    result = {'status': 'FINITE-EXACT Q controls; exact conditional covariance model controls',
              'Q_universe': [1, N], 'checks': checks,
              'special_Q': {str(n): q[n] for n in (20, 1024, 1536)},
              'old_weight_formula_first_hostile': 3,
              'map': '(r-1)/(3-r)',
              'two_step_matrix_in_s_equals_1_minus_r': [[3, 4], [2, 3]],
              'cross_ratio_multiplier': '3-2sqrt(2)',
              'Q_sampling_hypotheses': 'OPEN; not certified by this script',
              'model_iterates': reports}
    with open(args.output, 'w') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'model_iterates'}, indent=2))


if __name__ == '__main__':
    main()
