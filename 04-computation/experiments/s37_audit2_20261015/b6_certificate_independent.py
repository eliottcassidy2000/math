#!/usr/bin/env python3
"""S37 audit 2, script b6: independent re-implementation of the conditional residue certificate (Theorem B3 of
collatz_amplitude_boundary_20261013c.md) at N = 2^34 and N = 2^37, at the Mercator zero and at other zeros of D.

Independent choices (not the session's code):
 * root enclosure: own Rouche check on a disk of radius 1e-25 (Arb, 200 bits);
 * block sums expanded about the block START lo (not the centre): sum t(n) n^-s = lo^-s [m0 - (s/lo) m1] + err,
   |err| <= |s(s+1)|/2 * sum_(k<B) k^2 * lo^(-sigma-2), with m0, m1 the stored moments (m1 = sum (n - lo) t(n));
 * the two partial blocks after J = floor(L/2) and L = floor((2N-1)/3) use exact t(n) (full trajectories, numba with
   exact Python fallback) and the same expansion about their first point;
 * own tail bounds: E-tail |H(N)| N^-sigma/3 + |s|/3 N^-eta block(sigma) sum_j sqrt(a+bj) q^j (sum by mpmath nsum, no
   J0 cut), R-tail |s|/(2 sigma) (N - 1/2)^-sigma (from |(n-1/2)^-s - n^-s| <= |s|/2 (n-1/2)^(-sigma-1));
 * margins in closed form (constant and (4c, c)) or by bisection on nsum (slope).
Then: (i) intervals and margins at the Mercator zero (N = 2^34, 2^37) versus the session's JSON numbers;
(ii) the nested-interval "consistency test": the exact census energies on [2^34, 2^37) bound |C_(2^37) - C_(2^34)| by
     Cauchy-Schwarz, which makes nesting automatic; the margin is computed, and also for hypothetical e_k = 2;
(iii) the four other "certified nonzero" zeros and two uncertified ones at N = 2^37;
(iv) sanity check at zeros with Re > 1, where F is holomorphic and the residue must vanish: the finite numerator must lie
     inside the unconditional tail bound, and the certificate must not certify a nonzero residue.
Usage: python b6_certificate_independent.py moments37.npz census_ext.json
"""
import sys, json, math, time
import numpy as np
import numba as nb
import mpmath as mp
from flint import arb, acb, ctx

ctx.prec = 200
mp.mp.dps = 50
MOM = np.load(sys.argv[1])['moments']
CJ = json.load(open(sys.argv[2]))
B = 1 << 17
SUMK2 = (B - 1) * B * (2 * B - 1) // 6


@nb.njit
def tpar_full(n):
    # parity of tau(n) by the full trajectory; -1 if a value would exceed 2^63
    x = np.uint64(n)
    s = 0
    G = np.uint64(1) << np.uint64(63)
    while x != np.uint64(1):
        if x & np.uint64(1):
            if x >= G:
                return -1
            x = x + ((x + np.uint64(1)) >> np.uint64(1))
        else:
            x = x >> np.uint64(1)
        s += 1
    return s & 1


def t_exact(n):
    p = tpar_full(n)
    if p < 0:
        s = 0
        while n != 1:
            n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
            s += 1
        p = s & 1
    return 1 - 2 * p


def Dm(s):
    return 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3


def enclose(guess):
    z = mp.findroot(Dm, mp.mpc(*guess))
    c = acb(mp.nstr(z.real, 45), mp.nstr(z.imag, 45))
    r = arb('1e-25')
    l2 = arb(2).log(); l15 = (arb(3) / 2).log()
    D = lambda s: 1 + acb(2) ** (-s) + acb(arb(3) / 2) ** s / 3
    Dp = lambda s: -l2 * acb(2) ** (-s) + l15 * acb(arb(3) / 2) ** s / 3
    M2 = l2 * l2 * arb(2) ** (-c.real + r) + l15 * l15 * (arb(3) / 2) ** (c.real + r) / 3
    ok = (D(c).abs_upper() + M2 * r * r / 2) < (Dp(c).abs_lower() * r)
    assert ok, 'Rouche'
    s = acb(arb(mp.nstr(z.real, 45), '1e-25'), arb(mp.nstr(z.imag, 45), '1e-25'))
    return z, s, Dp(s)


def partial_sum(a, b, s):
    # sum_(a <= n <= b) t(n) n^-s with exact t, expanded about a; returns (acb value, error bound arb)
    if b < a:
        return acb(0), arb(0)
    m0 = 0; m1 = 0; k2 = 0
    for n in range(a, b + 1):
        t = t_exact(n)
        m0 += t; m1 += (n - a) * t; k2 += (n - a) ** 2
    val = acb(a) ** (-s) * (m0 - s * arb(m1) / a)
    err = (s * (s + 1)).abs_upper() / 2 * k2 * arb(a) ** (-s.real - 2)
    return val, err


def numerator(s, N, tN, moments):
    L = (2 * N - 1) // 3
    J = L // 2
    nb_ = N // B
    jl = -(-(J + 1) // B); js = -(-(L + 1) // B)
    sig = s.real
    long_, short_ = acb(0), acb(0)
    errl, errs = arb(0), arb(0)
    fac = (s * (s + 1)).abs_upper() / 2 * SUMK2
    for j in range(jl, nb_):
        lo = j * B
        m0 = int(moments[j, 0]); m1 = int(moments[j, 1])
        v = acb(lo) ** (-s) * (m0 - s * arb(m1) / lo)
        e = fac * arb(lo) ** (-sig - 2)
        long_ += v; errl += e
        if j >= js:
            short_ += v; errs += e
    pl, pel = partial_sum(J + 1, jl * B - 1, s)
    ps, pes = partial_sum(L + 1, js * B - 1, s)
    last = tN * acb(N) ** (-s)
    long_ += pl + last; short_ += ps + last
    errl += pel; errs += pes
    num = acb(arb(2) / 3) ** s * (short_ + acb(2) ** (-s) * long_)
    err = (arb(2) / 3) ** sig * (errs + arb(2) ** (-sig) * errl)
    return num + acb(arb(0, err.upper()), arb(0, err.upper())), err


def tails(z, N, hN):
    sig = z.real
    eta = sig - mp.mpf(1) / 2
    q = mp.power(2, -eta)
    block = mp.sqrt((1 - mp.power(2, -2 * sig - 1)) / (2 * sig + 1))
    hterm = mp.mpf(hN) * mp.power(N, -sig) / 3
    rtail = abs(z) / (2 * sig) * mp.power(N - mp.mpf(1) / 2, -sig)
    pref = abs(z) / 3 * mp.power(N, -eta) * block
    return dict(sig=sig, q=q, block=block, hterm=hterm, rtail=rtail, pref=pref)


def Ssum(a, b, q):
    if b == 0:
        return mp.sqrt(a) / (1 - q)
    return mp.nsum(lambda j: mp.sqrt(a + b * j) * q ** j, [0, mp.inf])


def to_mp(x):
    return mp.mpf(x.str(40, radius=False))


def cabs_mid(z):
    return abs(mp.mpc(to_mp(z.real.mid()), to_mp(z.imag.mid())))


def certificate(z, s, Dps, N, tN, hN, moments, budgets):
    num, err = numerator(s, N, tN, moments)
    r = -acb(arb(3) / 2) ** s * num / Dps
    kap = (arb(3) / 2) ** s.real / Dps.abs_lower()
    T = tails(z, N, hN)
    rlo = to_mp(r.abs_lower().lower())
    rhi = to_mp(r.abs_upper().upper())
    k = to_mp(kap.upper())
    out = []
    for (a, b) in budgets:
        rad = k * (T['hterm'] + T['pref'] * Ssum(a, b, T['q']) + T['rtail'])
        out.append(((a, b), rlo - rad, rhi + rad))
    return dict(r=r, num=num, rlo=rlo, rhi=rhi, kappa=k, T=T, rows=out)


hN37, hN34 = abs(-67763), 26869
N37, N34 = mp.mpf(2) ** 37, mp.mpf(2) ** 34
mom34 = MOM[: (1 << 34) // B]
print('(i) Mercator zero')
t0 = time.time()
z0, s0, Dp0 = enclose(('0.99820922076525681364', '240.24783014069124432530'))
scale = 2 * mp.power(2, 29 * (z0.real - 1))
budgets = [(2, 0), (4, 1), (8, 2), (16, 4), (32, 8), (64, 16)]
res = {}
for lab, N, tN, hN, mom in (('2^34', 1 << 34, 1, hN34, mom34), ('2^37', 1 << 37, -1, hN37, MOM)):
    c = certificate(z0, s0, Dp0, N, tN, hN, mom, budgets)
    res[lab] = c
    print(f'  N = {lab}: residue before tail {c["r"]}; |r| in [{mp.nstr(c["rlo"], 10)}, {mp.nstr(c["rhi"], 10)}]')
    print('    normalized intervals: ' + '; '.join(f'{ab}: [{mp.nstr(scale * lo, 5)}, {mp.nstr(scale * hi, 5)}]'
                                               for ab, lo, hi in c['rows']))
    T = c['T']
    target = (c['rlo'] / c['kappa'] - T['hterm'] - T['rtail']) / T['pref']
    amax = (target * (1 - T['q'])) ** 2
    cmax = (target / mp.nsum(lambda j: mp.sqrt(4 + j) * T['q'] ** j, [0, mp.inf])) ** 2
    lo_b, hi_b = mp.mpf(0), mp.mpf(1e6)
    for _ in range(80):
        mid = (lo_b + hi_b) / 2
        if Ssum(2, mid, T['q']) < target:
            lo_b = mid
        else:
            hi_b = mid
    print(f'    margins: largest constant a = {mp.nstr(amax, 7)}, largest slope b (a = 2) = {mp.nstr(lo_b, 7)}, '
          f'largest (4c, c): c = {mp.nstr(cmax, 7)}')
print(f'  ({time.time() - t0:.0f} s)')

print('(ii) nested-interval test: is nesting automatic given the exact energies on [2^34, 2^37)?')
c34, c37 = res['2^34'], res['2^37']
num_diff = abs(mp.mpc(to_mp(c37['num'].real.mid()), to_mp(c37['num'].imag.mid()))
               - mp.mpc(to_mp(c34['num'].real.mid()), to_mp(c34['num'].imag.mid())))
T34, T37 = c34['T'], c37['T']
sig = z0.real
E = {k: int(CJ['octave_sums_H_H2'][str(k)][1]) for k in (34, 35, 36)}
cs = mp.mpf(0)
cs2 = mp.mpf(0)
for k in (34, 35, 36):
    X = mp.mpf(2) ** k
    l2n = mp.sqrt(X ** (-2 * sig - 1) * (1 - mp.power(2, -2 * sig - 1)) / (2 * sig + 1))
    cs += mp.sqrt(E[k]) * l2n
    cs2 += mp.sqrt(2) * X * l2n
pref = abs(z0) / 3
bound_actual = (hN34 * N34 ** (-sig) + hN37 * N37 ** (-sig)) / 3 + pref * cs + abs(z0) / (2 * sig) * (N34 - 0.5) ** (-sig)
bound_e2 = (hN34 * N34 ** (-sig) + hN37 * N37 ** (-sig)) / 3 + pref * cs2 + abs(z0) / (2 * sig) * (N34 - 0.5) ** (-sig)
tail34 = T34['hterm'] + T34['pref'] * Ssum(2, 0, T34['q']) + T34['rtail']
tail37 = T37['hterm'] + T37['pref'] * Ssum(2, 0, T37['q']) + T37['rtail']
print(f'  |C(2^37) - C(2^34)| (numerators at s0) = {mp.nstr(num_diff, 6)}; Cauchy-Schwarz bound with the exact energies '
      f'e_34..36 = {mp.nstr(bound_actual, 6)}; same bound if every e_k were 2: {mp.nstr(bound_e2, 6)}')
print(f'  tail widths (numerator units, budget (2,0)): T(2^34) = {mp.nstr(tail34, 6)}, T(2^37) = {mp.nstr(tail37, 6)}, '
      f'difference {mp.nstr(tail34 - tail37, 6)}')
print(f'  nesting margin (T34 - T37) - bound: exact energies {mp.nstr(tail34 - tail37 - bound_actual, 6)}; with e_k = 2: '
      f'{mp.nstr(tail34 - tail37 - bound_e2, 6)} (the 2^34 tail counts the H(2^37) boundary term once, the bound twice)')

print('(iii) other zeros at N = 2^37, budget (2, 0): |r| intervals (fit values from the residue table for comparison)')
others = [('0.98494', '720.7401'), ('0.94973', '1572.5561'), ('0.86345', '131.1562'), ('0.69392', '22.2532'),
          ('0.92158', '1092.0834'), ('0.99317', '2533.5254')]
for g in others:
    t1 = time.time()
    z, s, Dps = enclose(g)
    c = certificate(z, s, Dps, 1 << 37, -1, hN37, MOM, [(2, 0)])
    (_, lo, hi), = c['rows']
    print(f'  rho = {mp.nstr(z.real, 7)} + {mp.nstr(z.imag, 9)} i: |r_N| = {mp.nstr((c["rlo"] + c["rhi"]) / 2, 6)}, '
          f'tail radius {mp.nstr(hi - c["rhi"], 6)}; |r| in [{mp.nstr(lo, 5)}, {mp.nstr(hi, 5)}], '
          f'nonzero certified: {lo > 0} ({time.time() - t1:.0f} s)')

print('(iv) sanity check at zeros with Re > 1 (F holomorphic there, residue 0):')
for g in (('2.3029', '23.604'), ('1.2', '30.0'), ('1.05', '250.0')):
    try:
        z, s, Dps = enclose(g)
    except Exception as ex:
        print(f'  start {g}: {ex}')
        continue
    if z.real <= 1:
        print(f'  start {g} converged to {mp.nstr(z, 8)} (Re <= 1), skipped')
        continue
    c = certificate(z, s, Dps, 1 << 37, -1, hN37, MOM, [(2, 0)])
    (_, lo, hi), = c['rows']
    sigz = z.real
    # unconditional tail of E + R for sigma > 1: |E - E_N| <= (2/3) sum_(n>N) n^-sigma <= (2/3) N^(1-sigma)/(sigma-1)
    unc = mp.mpf(2) / 3 * mp.power(N37, 1 - sigz) / (sigz - 1) + abs(z) / (2 * sigz) * (N37 - 0.5) ** (-sigz)
    numabs = cabs_mid(c['num'])
    print(f'  rho = {mp.nstr(z, 8)}: |C_N(rho)| = {mp.nstr(numabs, 4)} <= unconditional tail {mp.nstr(unc, 4)}: {numabs <= unc}; '
          f'conditional |r| interval [{mp.nstr(lo, 4)}, {mp.nstr(hi, 4)}] certifies nonzero: {lo > 0}')
