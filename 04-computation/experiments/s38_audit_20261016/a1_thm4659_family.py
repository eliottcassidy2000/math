#!/usr/bin/env python3
"""S38 audit, item 1: THM-4659 (zeros of D(s) = 1 + 2^(-s) + (3/2)^s/3 near Re s = 1).

(1) the identity D(1 + i pi q/L + w) = Phi(w, eps), eps = q theta - p (p even, q odd), at random complex w;
(2) the Taylor coefficients of the branch w(eps), w(0) = 0, by Cauchy integrals on circles in the eps-plane
    (a_1 = i pi/kappa, a_2 = -C, parity of the coefficients), the nearest branch point (double root of Phi) and the
    radius of convergence;
(3) the quoted zeros (84/53, 1054/665, 50508/31867, 8/5 continuation, the odd multiples (3162,1995), (5270,3325)) and
    the ratios (1 - Re rho)/eps^2;
(4) every zero with Re > 0.98 in collatz_D_zeros_and_boundary_residues_20261016.json identified with an explicit pair
    (p even, q odd) and compared with Newton on Phi; frequency offsets from the half-integers;
(5) statement 5 / title: alias pairs on the grids with P = 1024 and P = 1000 points per octave, detunings;
(6) convergents of log2 3 with p even, q odd.
mpmath, 60 digits.  Usage: python a1_thm4659_family.py zeros.json
"""
import sys, json, math, os, random
import mpmath as mp

mp.mp.dps = 60
ZJ = sys.argv[1]
L = mp.log(2)
M = mp.log(mp.mpf(3) / 2)
K = L - M
TH = mp.log(3) / L
C = mp.pi ** 2 * L ** 2 / K ** 3
HERE = os.path.dirname(os.path.abspath(__file__))
out = {'C': mp.nstr(C, 15), 'pi_over_kappa': mp.nstr(mp.pi / K, 15)}
print('C =', mp.nstr(C, 15), ' pi/kappa =', mp.nstr(mp.pi / K, 15))


def D(s):
    return 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3


def Dp(s):
    return -L * mp.power(2, -s) + M * mp.power(mp.mpf(3) / 2, s) / 3


def Phi(w, e):
    return 1 - mp.power(2, -w) / 2 - mp.exp(1j * mp.pi * e) * mp.power(mp.mpf(3) / 2, w) / 2


def dPhi(w, e):
    return L * mp.power(2, -w) / 2 - M * mp.exp(1j * mp.pi * e) * mp.power(mp.mpf(3) / 2, w) / 2


def wroot(e, w0=None):
    if w0 is None:
        w0 = 1j * mp.pi * e / K - C * e ** 2
    return mp.findroot(lambda w: Phi(w, e), w0, tol=mp.mpf(10) ** (-50))


def eps(p, q):
    return q * TH - p


# (1) identity
rng = random.Random(1)
idmax = 0
for (p, q) in [(84, 53), (1054, 665), (50508, 31867), (8, 5), (3162, 1995), (2192, 1383)]:
    e = eps(p, q)
    for _ in range(5):
        w = mp.mpc(rng.uniform(-0.6, 0.6), rng.uniform(-3, 3))
        dlt = abs(D(1 + 1j * mp.pi * q / L + w) - Phi(w, e))
        idmax = max(idmax, dlt)
out['identity_max_abs_difference'] = mp.nstr(idmax, 5)
print('(1) max |D(1 + i pi q/L + w) - Phi(w, eps)| over 30 random (w, pair):', mp.nstr(idmax, 5))

# (2) Taylor coefficients of w(eps) by Cauchy integrals on |eps| = h
coefs = {}
for h in (mp.mpf('0.004'), mp.mpf('0.008')):
    Kp = 64
    vals = []
    prev = mp.mpc(0)
    # continue along the radius first, then around the circle
    for j in range(1, 21):
        prev = wroot(h * j / 20, prev)
    for k in range(Kp):
        e = h * mp.exp(2j * mp.pi * k / Kp)
        prev = wroot(e, prev)
        vals.append(prev)
    a = []
    for n in range(0, 9):
        s = sum(vals[k] * mp.exp(-2j * mp.pi * k * n / Kp) for k in range(Kp)) / Kp / h ** n
        a.append(s)
    coefs[str(h)] = a
a = coefs['0.004']
print('(2) Taylor coefficients a_n of w(eps) (h = 0.004 and 0.008):')
for n in range(0, 9):
    print(f'   a_{n} = {mp.nstr(a[n], 12)}   (h = 0.008: {mp.nstr(coefs["0.008"][n], 12)})')
out['taylor_a1'] = mp.nstr(a[1], 15)
out['taylor_a2'] = mp.nstr(a[2], 15)
out['taylor_a3'] = mp.nstr(a[3], 12)
out['taylor_a4'] = mp.nstr(a[4], 12)
print('   a_1 - i pi/kappa =', mp.nstr(a[1] - 1j * mp.pi / K, 5), '; a_2 + C =', mp.nstr(a[2] + C, 5))
# branch point: double root of Phi
u = M / (L + M)
v = L / (L + M)
wb = -mp.log(2 * u) / L
ce = 2 * v * mp.exp(-M * wb)        # = e^(i pi eps)
eb = -1j * mp.log(ce) / mp.pi
out['branch_point_eps'] = mp.nstr(eb, 12)
out['branch_point_w'] = mp.nstr(wb, 12)
print('   double root of Phi at w =', mp.nstr(wb, 10), '(Re s =', mp.nstr(1 + wb, 8), ') for eps =', mp.nstr(eb, 10),
      '; |a_n|^(1/n) for n = 4..8:', [mp.nstr(abs(a[n]) ** (mp.mpf(1) / n), 6) for n in range(4, 9)],
      ' 1/|eps_b| =', mp.nstr(1 / abs(eb), 6))
# confirm the branch from w(0) = 0 actually reaches the branch point along eps = -i y
prev = mp.mpc(0)
for j in range(1, 200):
    y = abs(eb) * j / 200
    prev = wroot(-1j * y, prev)
out['w_at_0.995_branch_point'] = mp.nstr(prev, 10)
print('   branch w(-i y) at y = 0.995 |eps_b|:', mp.nstr(prev, 10), '(-> double root', mp.nstr(wb, 8), ')')

# (3) quoted zeros
res3 = []
for (p, q, label) in [(84, 53, 's0'), (1054, 665, 'rho*'), (50508, 31867, '50508/31867'), (3162, 1995, '3x'),
                      (5270, 3325, '5x'), (2192, 1383, '1383')]:
    e = eps(p, q)
    w = wroot(e)
    rho = 1 + 1j * mp.pi * q / L + w
    rho = mp.findroot(D, rho)
    row = dict(p=p, q=q, label=label, eps=mp.nstr(e, 12), rho=mp.nstr(rho, 15),
               one_minus_re=mp.nstr(1 - rho.real, 10), ratio=mp.nstr((1 - rho.real) / e ** 2, 10),
               C_eps2=mp.nstr(C * e ** 2, 10), im_minus_pertur=mp.nstr(rho.imag - (mp.pi * q / L + mp.pi * e / K), 6),
               absDp=mp.nstr(abs(Dp(rho)), 10))
    res3.append(row)
    print('(3)', row)
# 8/5 continuation along real eps
e85 = eps(8, 5)
prev = mp.mpc(0)
path = []
for j in range(1, 401):
    prev = wroot(e85 * j / 400, prev)
rho85 = mp.findroot(D, 1 + 1j * mp.pi * 5 / L + prev)
guess = 1 + 1j * mp.pi * 5 / L + 1j * mp.pi * e85 / K - C * e85 ** 2
direct = mp.findroot(lambda w: Phi(w, e85), 1j * mp.pi * e85 / K - C * e85 ** 2)
rho85d = 1 + 1j * mp.pi * 5 / L + direct
row = dict(eps=mp.nstr(e85, 10), continuation=mp.nstr(rho85, 12), newton_from_perturbative_guess=mp.nstr(rho85d, 12),
           C_eps2=mp.nstr(C * e85 ** 2, 6), one_minus_re=mp.nstr(1 - rho85.real, 6), freq=mp.nstr(rho85.imag * L / (2 * mp.pi), 8))
res3.append(row)
print('(3) 8/5:', row)
out['quoted'] = res3

# (4) the zeros with Re > 0.98 in the JSON
Z = json.load(open(ZJ))['zeros_and_residues']
fam = [r for r in Z if r['re'] > 0.98]
rows4 = []
maxoff = 0
for r in fam:
    rho_j = mp.mpc(r['re'], r['im'])
    q = int(mp.nint(rho_j.imag * L / mp.pi))
    p = int(mp.nint(q * TH / 2) * 2)
    e = eps(p, q)
    try:
        w = wroot(e)
        rho = mp.findroot(D, 1 + 1j * mp.pi * q / L + w)
        dist = abs(rho - rho_j)
    except Exception:
        rho, dist = None, None
    off = float(rho_j.imag * L / (2 * mp.pi) - mp.mpf(q) / 2)
    maxoff = max(maxoff, abs(off))
    rows4.append(dict(im=r['im'], re=r['re'], q=q, p=p, q_odd=(q % 2 == 1), p_even=(p % 2 == 0), eps=float(e),
                      dist_to_family_zero=float(dist) if dist is not None else None, freq_offset=off,
                      ratio=float((1 - rho_j.real) / e ** 2), abs_r=r['abs_r']))
out['family'] = rows4
print(f'(4) {len(fam)} zeros with Re > 0.98 below height {json.load(open(ZJ))["T"]}:')
for x in rows4:
    print(f"   {x['re']:.8f} + {x['im']:10.4f}i  (p, q) = ({x['p']}, {x['q']})  odd/even {x['q_odd']}/{x['p_even']}  "
          f"eps = {x['eps']:+.6e}  |rho - rho_(p,q)| = {x['dist_to_family_zero']:.2e}  freq - q/2 = {x['freq_offset']:+.5f}  "
          f"(1-Re)/eps^2 = {x['ratio']:.2f}  |r_hat| = {x['abs_r']:.4f}")
out['max_freq_offset'] = maxoff
print('   max |freq - q/2| =', maxoff)

# (5) aliasing: pairs q + q' = 2 P k
for P in (1024, 1000):
    pairs = []
    qs = {x['q']: x for x in rows4}
    for x in rows4:
        for kk in (1, 2, 3, 4):
            q2 = 2 * P * kk - x['q']
            if q2 in qs and q2 > x['q']:
                y = qs[q2]
                det = (mp.mpf(x['im']) + mp.mpf(y['im'])) * L / (2 * mp.pi) - P * kk   # f + f' - P k (cycles/octave)
                pairs.append((x['im'], y['im'], kk, float(det), x['eps'] + y['eps']))
    out[f'alias_pairs_P{P}'] = pairs
    print(f'(5) P = {P}: alias pairs among the family members (Im, Im\', k, detuning cycles/octave, eps + eps\'):')
    for pr in pairs:
        print('   ', pr)
    print(f'   2 P theta = {mp.nstr(2 * P * TH, 12)}; distance to the nearest even integer '
          f'{mp.nstr(2 * P * TH - 2 * mp.nint(P * TH), 8)}')
out['phase_drift_P1024_over_2p22_2p40_cycles'] = None

# (6) convergents
cf = []
x = TH
for _ in range(25):
    a0 = int(mp.floor(x))
    cf.append(a0)
    x = 1 / (x - a0)
h0, h1, k0, k1 = 1, cf[0], 0, 1
conv = [(cf[0], 1)]
for a0 in cf[1:]:
    h0, h1 = h1, a0 * h1 + h0
    k0, k1 = k1, a0 * k1 + k0
    conv.append((h1, k1))
matched = [(p_, q_, float(eps(p_, q_))) for p_, q_ in conv if p_ % 2 == 0 and q_ % 2 == 1]
out['continued_fraction'] = cf
out['convergents'] = conv[:16]
out['parity_matched_convergents'] = matched[:8]
print('(6) CF of log2 3:', cf[:16])
print('    convergents:', conv[:14])
print('    parity matched (p even, q odd):', matched[:7])
json.dump(out, open(os.path.join(HERE, 'a1_thm4659_family.json'), 'w'), indent=1, default=str)
