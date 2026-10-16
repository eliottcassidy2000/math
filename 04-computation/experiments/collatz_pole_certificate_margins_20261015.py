#!/usr/bin/env python3
"""Budget margins of the CONDITIONAL Collatz Mercator-pole certificate at a general census size N = 2^LOGN.

Bisects (in Arb, each test rigorous) for the largest constant budget a, the largest slope b with a = 2, and the largest
c in (4c, c) that still certify a nonzero residue.  Formula as in collatz_pole_certificate_general_N_20261015.py:

Same formula as collatz_pole_certificate_logenergy_20261015.py (itself the log-budget extension of renewal_shapes_oct11's
collatz_pole_certificate_20261013c.py, Theorem B3 of collatz_amplitude_boundary_20261013c.md): exact block moments
through N, the same root disk, quadrature and R-tail; the boundary term t(N) N^(-s) carries the sign t(N) = (-1)^LOGN
(the N = 2^34 script hard-codes t(N) = +1); the energy hypothesis is needed only for dyadic X >= N.  With exact data
through 2^37 (collatz_boundary_moments_2p37_20261015.py) the hypothesis is weakened from X >= 2^34 to X >= 2^37.
Certificates at two census sizes give nested intervals, but this is implied by the exact data in between (Lemma B2:
the finite numerators differ by that data, which the census energies bound far below the difference of the radii),
so nesting checks the code, not the budget (second S37 referee, E1).

Hypothesis family (UNPROVED), for every dyadic X = N 2^j, j >= 0:
    (9b)   integral_X^(2X) H(floor(x))^2 dx <= X^2 (a + b j).
Under (9b), Cauchy-Schwarz on each dyadic block gives
    integral_N^inf |H(x)| x^(-sigma-1) dx <= sqrt((1 - 2^(-2 sigma - 1))/(2 sigma + 1)) N^(-eta) sum_j sqrt(a + b j) q^j,
with eta = sigma - 1/2, q = 2^(-eta).  The j-sum is evaluated exactly in Arb for j < J0 and bounded above for j >= J0 by
sqrt(a + b) sum_(j >= J0) j q^j = sqrt(a + b) q^J0 (J0 - (J0 - 1) q)/(1 - q)^2   (a + b j <= (a + b) j for j >= 1).
(a, b) = (2, 0) reproduces the original constant-budget interval.  Session opus-2026-10-15-S37.
Usage: python collatz_pole_certificate_margins_20261015.py --moments M.npz --census C.json --output out.json
(census JSON from collatz_boundary_moments_2p37_20261015.py: keys N, S, H, t_N)
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
    d = np.load(args.moments)
    v = d['moments']
    N = int(d['N'])
    B = int(d['B'])
    need(B == 2 ** 17 and N % B == 0 and len(v) == N // B, 'stored universe')
    census = json.load(open(args.census))
    need(census['N'] == N, 'census population')
    tN = int(census['t_N'])
    logN = N.bit_length() - 1
    need(N == 1 << logN and tN == (-1) ** logN, 't(2^k) = (-1)^k')
    need(int(v[:, 0].sum()) + tN == census['S'], 'moment total')
    hN = abs(int(census['H']))
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
    last = tN * acb(N) ** (-s)
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

    def lower(a, b):
        A = arb(str(a))
        Bc = arb(str(b))
        S = arb(0)
        qj = arb(1)
        for jj in range(J0):
            S += (A + Bc * jj).sqrt() * qj
            qj *= q
        S += (A + Bc).sqrt() * q ** J0 * (J0 - (J0 - 1) * q) / (1 - q) ** 2
        etail = arb(hN) / 3 * arb(N) ** (-sigma) + s.abs_upper() / 3 * arb(N) ** (-eta) * block * S
        rho_e = (arb(3) / 2) ** sigma / Dp(s).abs_lower() * (etail + rtail)
        return float((scale * (residue.abs_lower() - rho_e)).lower())

    def bisect(f, lo, hi, it=60):
        for _ in range(it):
            mid = (lo + hi) / 2
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        return lo
    amax = bisect(lambda a: lower(a, 0), 2.0, 1e7)
    bmax = bisect(lambda b: lower(2, b), 0.0, 1e7)
    cmax = bisect(lambda c: lower(4 * c, c), 0.5, 1e7)
    rows = dict(largest_constant_budget=amax, largest_slope_with_a2=bmax, largest_c_for_4c_c=cmax)
    print(f'N = 2^{logN}: largest certified constant budget E_X <= a X^2 (X >= N): a = {amax:.3f}; '
          f'largest slope with a = 2: b = {bmax:.3f}; largest (4c, c): c = {cmax:.4f}', flush=True)
    report = dict(status='CONDITIONAL: rigorously interval-checked implications; the energy budgets are UNPROVED',
                  N=N, B=B, root_disk_center=str(c), residue_before_future_tail=str(residue),
                  margins=rows, moments_sha256=hashlib.sha256(open(args.moments, 'rb').read()).hexdigest(),
                  session='opus-2026-10-15-S37', t_N=tN, H_N=int(census['H']))
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)


if __name__ == '__main__':
    main()
