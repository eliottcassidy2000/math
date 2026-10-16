#!/usr/bin/env python3
"""S37 audit, part 6: independent (mpmath, 60 digits) check of the HYP-9283(B) log-budget tail.

Under E_X <= X^2 (a + b j), X = 2^(34+j), Cauchy-Schwarz per dyadic block gives
  int_N^inf |H| x^(-sigma-1) dx <= block * N^(-eta) * S(a,b),   S(a,b) = sum_j sqrt(a + b j) q^j,  q = 2^(-eta),
  block = sqrt((1 - 2^(-2 sigma - 1))/(2 sigma + 1)),  eta = sigma - 1/2.
The future residue radius is rho(a,b) = (3/2)^sigma/|D'(s0)| * (|H(N)| N^-sigma/3 + |s0|/3 N^-eta block S + Rtail).
Anchoring on the session's (2,0) row, every other row must satisfy
  lo(a,b) = lo(2,0) - scale*(rho(a,b) - rho(2,0)),  hi(a,b) = hi(2,0) + scale*(rho(a,b) - rho(2,0)).
S(a,b) is summed here to convergence (no tail bound), so this also checks the script's J0 = 400 tail estimate.
Also: (i) closed form of the original constant budget, (ii) the largest certifiable budgets, (iii) the clock-mesh
translation C(P_b) <= X^2 (a' + b' j)  =>  E_X <= X^2 (4(a'+b'j) + 4/27 + O(X^-1/2)), checked for (15, 4).
"""
import json, os
import mpmath as mp

mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results'))
J = json.load(open(os.path.join(RES, 'collatz_pole_certificate_logenergy_20261015.json')))
s0 = mp.mpc('0.9982092207652568136445639331675041376460', '240.2478301406912443253007903875287137674')
D = lambda s: 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3
s0 = mp.findroot(D, s0)
Dp = mp.diff(D, s0)
sig = s0.real
eta = sig - mp.mpf(1) / 2
q = mp.power(2, -eta)
N = mp.mpf(2) ** 34
hN = 26869
block = mp.sqrt((1 - mp.power(2, -2 * sig - 1)) / (2 * sig + 1))
scale = 2 * mp.power(2, 29 * (sig - 1))
rtail = abs(s0) / (2 * sig) * mp.power(1 - 1 / (2 * N), -sig - 1) * mp.power(N, -sig)
print(f's0 = {mp.nstr(s0, 25)}, |D(s0)| = {mp.nstr(abs(D(s0)), 3)}, |D\'(s0)| = {mp.nstr(abs(Dp), 12)}, q = {mp.nstr(q, 12)}')


def Ssum(a, b):
    return mp.nsum(lambda j: mp.sqrt(a + b * j) * q ** j, [0, mp.inf])


def Sscript(a, b, J0=400):
    s = mp.fsum(mp.sqrt(a + b * j) * q ** j for j in range(J0))
    return s + mp.sqrt(a + b) * q ** J0 * (J0 - (J0 - 1) * q) / (1 - q) ** 2


def rho(a, b, S=None):
    S = Ssum(a, b) if S is None else S
    et = mp.mpf(hN) / 3 * mp.power(N, -sig) + abs(s0) / 3 * mp.power(N, -eta) * block * S
    return mp.power(mp.mpf(3) / 2, sig) / abs(Dp) * (et + rtail)


rows = {(r['a'], r['b']): (mp.mpf(r['normalized_residue_lower'].split('+/-')[0].strip('[ ')),
                           mp.mpf(r['normalized_residue_upper'].split('+/-')[0].strip('[ '))) for r in J['hypotheses']}
lo20, hi20 = rows[(2, 0)]
r20 = rho(2, 0)
print('row (a,b): session [lo, hi] | predicted from (2,0) anchor with exact S | S exact vs script tail-bounded S')
worst = 0
for (a, b), (lo, hi) in rows.items():
    r = rho(mp.mpf(a), mp.mpf(b))
    plo = lo20 - scale * (r - r20)
    phi = hi20 + scale * (r - r20)
    worst = max(worst, abs(plo - lo), abs(phi - hi))
    print(f'  ({a}, {b}): [{mp.nstr(lo, 6)}, {mp.nstr(hi, 6)}] | [{mp.nstr(plo, 6)}, {mp.nstr(phi, 6)}] | '
          f'{mp.nstr(Ssum(mp.mpf(a), mp.mpf(b)), 12)} vs {mp.nstr(Sscript(mp.mpf(a), mp.mpf(b)), 12)}')
print(f'max |predicted - session| over rows: {mp.nstr(worst, 3)}')
# original constant-budget closed form (amplitude note (11a)): sqrt(2) block / (1 - q)
print(f'(2,0): S = sqrt(2)/(1-q) closed form {mp.nstr(mp.sqrt(2) / (1 - q), 15)} vs nsum {mp.nstr(Ssum(2, 0), 15)}')
# |r| bounds implied by the (2,0) anchor
rlo = lo20 / scale + r20
rhi = hi20 / scale - r20
print(f'|r| bounds implied by the (2,0) row: [{mp.nstr(rlo, 8)}, {mp.nstr(rhi, 8)}]; future radius rho(2,0) = {mp.nstr(r20, 8)}')


def lower(a, b):
    return scale * (rlo - rho(a, b))


# largest certifiable constant budget and largest (4c, c) budget
f_const = lambda a: lower(a, 0)
a_max = mp.findroot(f_const, 100)
f_lin = lambda c: lower(4 * c, c)
c_max = mp.findroot(f_lin, 20)
print(f'largest constant budget E_X <= a X^2 still certified: a = {mp.nstr(a_max, 6)}; '
      f'largest (4c, c): c = {mp.nstr(c_max, 6)} (i.e. ({mp.nstr(4 * c_max, 5)}, {mp.nstr(c_max, 5)}))')
f_slope = lambda b: lower(2, b)
print(f'largest slope with a = 2: b = {mp.nstr(mp.findroot(f_slope, 30), 6)}')
# clock-mesh translation for (a', b') = (15, 4)
ok = True
worst_margin = None
for j in range(0, 400):
    X = mp.mpf(2) ** (34 + j)
    Sx = 15 + 4 * j + mp.mpf(1) / 27 + 1 / (3 * X)                 # S_X/X^2 upper bound (even k worst case)
    Ex = (2 * mp.sqrt(Sx) + 1 / mp.sqrt(X)) ** 2                    # E_X/X^2 <= (2 sqrt(S_X/X^2) + X^-1/2)^2
    margin = 64 + 16 * j - Ex
    worst_margin = margin if worst_margin is None else min(worst_margin, margin)
    ok &= margin > 0
print(f'clock-mesh (15, 4): E_X/X^2 <= (2 sqrt(15 + 4j + 1/27 + 1/(3X)) + X^-1/2)^2 <= 64 + 16 j for j = 0..399: {ok}; '
      f'minimum margin {mp.nstr(worst_margin, 6)} (at j = 0)')
