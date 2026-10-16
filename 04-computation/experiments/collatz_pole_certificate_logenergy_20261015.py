#!/usr/bin/env python3
"""Arb interval certificate for the CONDITIONAL Collatz Mercator pole under LOG-GROWING dyadic energy budgets.

Extends 04-computation/experiments/collatz_pole_certificate_20261013c.py (renewal_shapes_oct11; Theorem B3 of
collatz_amplitude_boundary_20261013c.md) without changing its finite numerator: the same exact block moments through
N = 2^34, the same root disk, the same quadrature and R-tail.  Only the future E-tail is replaced.

Hypothesis family (UNPROVED), for every dyadic X = 2^(34+j), j >= 0:
    (9b)   integral_X^(2X) H(floor(x))^2 dx <= X^2 (a + b j).
Under (9b), Cauchy-Schwarz on each dyadic block gives
    integral_N^inf |H(x)| x^(-sigma-1) dx <= sqrt((1 - 2^(-2 sigma - 1))/(2 sigma + 1)) N^(-eta) sum_j sqrt(a + b j) q^j,
with eta = sigma - 1/2, q = 2^(-eta).  The j-sum is evaluated exactly in Arb for j < J0 and bounded above for j >= J0 by
sqrt(a + b) sum_(j >= J0) j q^j = sqrt(a + b) q^J0 (J0 - (J0 - 1) q)/(1 - q)^2   (a + b j <= (a + b) j for j >= 1).
(a, b) = (2, 0) reproduces the original constant-budget interval.  Session opus-2026-10-15-S37.
Usage: python collatz_pole_certificate_logenergy_20261015.py --moments M.npz --census C.json --output out.json
"""
import argparse, json, hashlib, time
import numpy as np
import mpmath as mp
from flint import arb, acb, ctx


def need(x, msg):
    if not x:
        raise RuntimeError(msg)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--moments', required=True)
    p.add_argument('--census', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    N = 2 ** 34
    B = 2 ** 17
    d = np.load(args.moments)
    v = d['moments']
    need(int(d['N']) == N and int(d['B']) == B, 'stored universe')
    census = json.load(open(args.census))
    need(census['N'] == N, 'census population')
    need(int(v[:, 0].sum()) + 1 == census['rows'][-1]['S'], 'moment total')
    hN = abs(census['rows'][-1]['residues'][2][0])
    ctx.prec = 160
    mp.mp.dps = 85
    theta = mp.log(3, 2)
    defect = mp.mpf('26.5') * theta - mp.nint(mp.mpf('26.5') * theta)
    w = 0j
    for j in range(1, 121):
        w = mp.findroot(lambda z: mp.exp(-z) + mp.exp((theta - 1) * z + 2j * mp.pi * defect * j / 120) - 2, w)
    zz = 1 + (w + 2j * mp.pi * mp.mpf('26.5')) / mp.log(2)
    c = acb(str(zz.real), str(zz.imag))
    radius = arb('1e-30')
    l2 = arb(2).log()
    l15 = (arb(3) / 2).log()

    def D(s):
        return 1 + acb(2) ** (-s) + acb(arb(3) / 2) ** s / 3

    def Dp(s):
        return -l2 * acb(2) ** (-s) + l15 * acb(arb(3) / 2) ** s / 3
    M = l2 * l2 * arb(2) ** (-c.real + radius) + l15 * l15 * (arb(3) / 2) ** (c.real + radius) / 3
    need(D(c).abs_upper() + M * radius * radius / 2 < Dp(c).abs_lower() * radius, 'Rouche root disk')
    need(Dp(c).abs_lower() > M * radius, 'simple root')
    s = acb(arb(str(zz.real), '1e-30'), arb(str(zz.imag), '1e-30'))
    sigma = s.real
    L = (2 * N - 1) // 3
    J = L // 2
    start_long = (J + 1 + B - 1) // B
    start_short = (L + 1 + B - 1) // B
    longsum = acb(0)
    shortsum = acb(0)
    started = time.monotonic()
    for j in range(start_long, len(v)):
        lo = j * B
        center = arb(2 * lo + B - 1) / 2
        sm = int(v[j, 0])
        mu1 = arb(2 * int(v[j, 1]) - (B - 1) * sm) / 2
        power = acb(center) ** (-s)
        val = power * (sm - s * mu1 / center)
        longsum += val
        if j >= start_short:
            shortsum += val
        if j % 32768 == 0:
            print('certified blocks', j, 'seconds', round(time.monotonic() - started, 1), flush=True)
    last = acb(N) ** (-s)
    longsum += last
    shortsum += last

    def quadrature_error(cut, start):
        lo = start * B
        partial = (lo - cut - 1) * arb(cut + 1) ** (-sigma)
        sum2 = arb(B) * (B * B - 1) / 12
        derivative_sum = arb(lo) ** (-sigma - 2) + arb(lo) ** (-sigma - 1) / (B * (sigma + 1))
        return partial + (s * (s + 1)).abs_upper() * sum2 * derivative_sum / 2
    errlong = quadrature_error(J, start_long)
    errshort = quadrature_error(L, start_short)
    numerator = acb(arb(2) / 3) ** s * (shortsum + acb(2) ** (-s) * longsum)
    error = (arb(2) / 3) ** sigma * (errshort + arb(2) ** (-sigma) * errlong)
    numerator += acb(arb(0, error.upper()), arb(0, error.upper()))
    residue = -acb(arb(3) / 2) ** s * numerator / Dp(s)
    eta = sigma - arb(1) / 2
    rtail = s.abs_upper() / (2 * sigma) * (1 - arb(1) / (2 * N)) ** (-sigma - 1) * arb(N) ** (-sigma)
    scale = 2 * arb(2 ** 29) ** (sigma - 1)
    q = arb(2) ** (-eta)
    block = ((1 - arb(2) ** (-2 * sigma - 1)) / (2 * sigma + 1)).sqrt()
    J0 = 400
    rows = []
    for (a, b) in [(2, 0), (2, 0.25), (3, 0.5), (4, 1), (8, 2), (16, 4), (32, 8), (64, 16)]:
        A = arb(str(a))
        Bc = arb(str(b))
        S = arb(0)
        qj = arb(1)
        for jj in range(J0):
            S += (A + Bc * jj).sqrt() * qj
            qj *= q
        tail = (A + Bc).sqrt() * q ** J0 * (J0 - (J0 - 1) * q) / (1 - q) ** 2
        S += tail
        etail = arb(hN) / 3 * arb(N) ** (-sigma) + s.abs_upper() / 3 * arb(N) ** (-eta) * block * S
        rho_e = (arb(3) / 2) ** sigma / Dp(s).abs_lower() * (etail + rtail)
        lo_ = (scale * (residue.abs_lower() - rho_e)).lower()
        hi_ = (scale * (residue.abs_upper() + rho_e)).upper()
        rows.append(dict(a=a, b=b, hypothesis=f'for every j >= 0: integral_X^(2X) H^2 dx <= X^2 ({a} + {b} j), X = 2^(34+j)',
                         normalized_residue_lower=str(lo_), normalized_residue_upper=str(hi_),
                         nonzero_certified=bool(lo_ > 0)))
        print(f'(a, b) = ({a}, {b}): 2|r| 2^(29(sigma-1)) in [{float(lo_):.4f}, {float(hi_):.4f}]  nonzero: {bool(lo_ > 0)}', flush=True)
    report = dict(status='CONDITIONAL: rigorously interval-checked implications; the energy budgets are UNPROVED',
                  N=N, B=B, root_disk_center=str(c), residue_before_future_tail=str(residue),
                  hypotheses=rows, moments_sha256=hashlib.sha256(open(args.moments, 'rb').read()).hexdigest(),
                  session='opus-2026-10-15-S37')
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)


if __name__ == '__main__':
    main()
