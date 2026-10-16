#!/usr/bin/env python3
"""S38 audit, item 5 (and THM-4660 Lemmas A, C): an independent enumeration of the zeros of
D(s) = 1 + 2^(-s) + (3/2)^s/3 by the phase torus, compared with collatz_D_zeros_and_boundary_residues_20261016.json.

D(sigma + i t) = 0  iff  1 + A + B = 0 with |A| = a(sigma) = 2^(-sigma), |B| = b(sigma) = (3/2)^sigma/3,
arg A = -t ln 2, arg B = t ln(3/2) (mod 2 pi).  For each sigma with a triangle (a + b >= 1, |a - b| <= 1) the phases are
(phi_A, phi_B) = +-(-arccos cA, arccos cB), cA = (b^2 - 1 - a^2)/(2a), cB = (a^2 - 1 - b^2)/(2b).  Writing
t = (-phi_A + 2 pi k1)/L = (phi_B + 2 pi k2)/M (L = ln 2, M = ln 3/2) gives the scalar equation
psi(sigma) := -M phi_A(sigma) - L phi_B(sigma) = 2 pi (L k2 - M k1).
psi is strictly monotone on each sigma-interval used (checked), so every zero is one lattice pair (k1, k2, sign):
the count of zeros in a box is a count of lattice points (float64 with a margin report; no argument principle).
Outputs: the zero list for 1/2 < Re < 1, 0 < Im < 20000 against the JSON; |D'(rho)| as a function of Re rho (Lemma A);
zeros per unit height in [1/4, 3] (Lemma C); the near-1/2 family, |D(1/2 + it)| minima and the multiplier; the
almost-period 2 pi 665/ln 2; gaps between consecutive zeros.
Usage: python a2_zeros_torus_count.py zeros.json
"""
import sys, json, math, os
import numpy as np
from scipy.optimize import brentq

ZJ = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
L, M = math.log(2.0), math.log(1.5)
KAP = L - M
T = 20000.0


def ab(s):
    return 2.0 ** (-s), 1.5 ** s / 3.0


def phases(s):
    a, b = ab(s)
    cA = (b * b - 1 - a * a) / (2 * a)
    cB = (a * a - 1 - b * b) / (2 * b)
    cA = min(1.0, max(-1.0, cA))
    cB = min(1.0, max(-1.0, cB))
    return -math.acos(cA), math.acos(cB)          # configuration '+'


def psi(s):
    pA, pB = phases(s)
    return -M * pA - L * pB


def Dval(s):
    return 1 + np.exp(-s * L) + np.exp(s * M) / 3


def Dp(s):
    return -L * np.exp(-s * L) + M * np.exp(s * M) / 3


def enumerate_zeros(s_lo, s_hi, Tmax):
    # monotonicity of psi on [s_lo, s_hi]
    ss = np.linspace(s_lo, s_hi, 20001)
    ps = np.array([psi(x) for x in ss])
    dps = np.diff(ps)
    mono = bool((dps < 0).all() or (dps > 0).all())
    p_lo, p_hi = psi(s_lo), psi(s_hi)
    lo, hi = min(p_lo, p_hi), max(p_lo, p_hi)
    zs = []
    margin = 1.0
    kmax = int(Tmax * L / (2 * math.pi)) + 3
    for sign in (+1, -1):
        for k1 in range(-2, kmax + 1):
            # value v = 2 pi (L k2 - M k1) must lie in sign*[lo, hi]
            e1 = (sign * lo / (2 * math.pi) + M * k1) / L
            e2 = (sign * hi / (2 * math.pi) + M * k1) / L
            for k2 in range(int(math.floor(min(e1, e2))) - 1, int(math.ceil(max(e1, e2))) + 2):
                v = 2 * math.pi * (L * k2 - M * k1)
                w = sign * v
                if not (lo < w < hi):
                    if lo - 1e-9 < w < hi + 1e-9:
                        margin = 0.0
                    continue
                margin = min(margin, w - lo, hi - w)
                sg = brentq(lambda x: psi(x) - w, s_lo, s_hi, xtol=1e-15, rtol=1e-15)
                pA, pB = phases(sg)
                pA, pB = sign * pA, sign * pB
                t = (-pA + 2 * math.pi * k1) / L
                if 0 < t < Tmax:
                    zs.append(complex(sg, t))
    zs.sort(key=lambda z: z.imag)
    return np.array(zs), mono, margin


out = {}
Z, mono, margin = enumerate_zeros(0.5, 1.0, T)
# polish with Newton in complex128 and check D
Zp = Z.copy()
for _ in range(5):
    Zp = Zp - Dval(Zp) / Dp(Zp)
out['n_zeros_half_to_one_below_20000'] = int(len(Z))
out['psi_monotone_on_half_one'] = mono
out['min_margin_of_lattice_values'] = margin
out['max_newton_correction'] = float(np.abs(Zp - Z).max())
out['max_abs_D_at_torus_zeros'] = float(np.abs(Dval(Z)).max())
print(f'torus enumeration: {len(Z)} zeros with 1/2 < Re < 1, 0 < Im < {T:.0f}; psi monotone {mono}; '
      f'min margin {margin:.3e}; max Newton correction {np.abs(Zp - Z).max():.2e}', flush=True)
J = json.load(open(ZJ))
ZJl = np.array([complex(r['re'], r['im']) for r in J['zeros_and_residues']])
dist = np.array([np.abs(Zp - z).min() for z in ZJl])
dist2 = np.array([np.abs(ZJl - z).min() for z in Zp])
out['json_n'] = int(len(ZJl))
out['json_max_dist_to_torus_zero'] = float(dist.max())
out['torus_max_dist_to_json_zero'] = float(dist2.max())
print(f'JSON list: {len(ZJl)} zeros; max distance JSON -> torus {dist.max():.2e}; torus -> JSON {dist2.max():.2e}', flush=True)
# per-500 rectangle counts against the JSON argument-principle table
ap = J['argument_principle']
mis = 0
for (a_, b_, w_, c_) in ap:
    c2 = int(((Zp.imag > a_) & (Zp.imag < b_)).sum())
    if c2 != c_ or abs(w_ - c2) > 0.01:
        mis += 1
out['rectangles_disagreeing_with_json_table'] = mis
# Lemma A: |D'| as a function of sigma
def absDp_sigma(s):
    pA, pB = phases(s)
    a, b = ab(s)
    A = a * complex(math.cos(pA), math.sin(pA))
    B = b * complex(math.cos(pB), math.sin(pB))
    return abs(-L * A + M * B)
grid1 = np.linspace(0.5, 1.0, 50001)
v1 = np.array([absDp_sigma(x) for x in grid1])
out['min_absDp_over_sigma_in_half_one'] = [float(v1.min()), float(grid1[v1.argmin()])]
out['absDp_at_zeros_min'] = float(np.abs(Dp(Zp)).min())
out['kappa_over_2'] = KAP / 2
# feasible sigma range and global inf over 0 <= Re <= 3
feas = [x for x in np.linspace(-1, 3.2, 42001) if (sum(ab(x)) >= 1 and abs(ab(x)[0] - ab(x)[1]) <= 1)]
feas = np.array(feas)
gl = np.array([absDp_sigma(x) for x in feas if 0 <= x <= 3])
fs = np.array([x for x in feas if 0 <= x <= 3])
out['feasible_sigma_range'] = [float(feas.min()), float(feas.max())]
out['inf_absDp_0_le_sigma_le_3'] = [float(gl.min()), float(fs[gl.argmin()])]
print('Lemma A: min |D\'| over sigma in [1/2, 1] =', out['min_absDp_over_sigma_in_half_one'], '(kappa/2 =', KAP / 2, ');',
      'over zeros listed:', out['absDp_at_zeros_min'], '; inf over 0 <= Re <= 3:', out['inf_absDp_0_le_sigma_le_3'],
      '; feasible Re range', out['feasible_sigma_range'], flush=True)
# Lemma C: zeros in [1/4, 3] per unit height
Za, _, _ = enumerate_zeros(0.25, 1.0, T)
Zb, monob, _ = enumerate_zeros(2.0 + 1e-12, 3.0 - 1e-12, T)
allz = np.sort(np.concatenate([Za.imag, Zb.imag]))
mx = 0
j = 0
for i in range(len(allz)):
    while allz[i] - allz[j] > 1.0:
        j += 1
    mx = max(mx, i - j + 1)
out['zeros_in_quarter_to_one'] = int(len(Za))
out['zeros_in_two_to_three'] = int(len(Zb))
out['psi_monotone_on_2_3'] = monob
out['max_zeros_per_unit_height_in_[1/4,3]'] = int(mx)
gaps = np.diff(np.sort(Zp.imag))
out['min_gap_Im_half_one'] = float(gaps.min())
out['min_gap_Im_quarter_three'] = float(np.diff(allz).min())
print(f'Lemma C: zeros with Im < {T:.0f}: {len(Za)} in [1/4, 1), {len(Zb)} in (2, 3); max per unit height {mx}; '
      f'min gap of Im in (1/2,1): {gaps.min():.3f}, in [1/4,3]: {np.diff(allz).min():.4f}', flush=True)
# near-1/2 family
near = Zp[Zp.real < 0.502]
out['near_half'] = [[float(z.real), float(z.imag)] for z in near]
print('zeros with Re < 0.502:', [(round(z.real, 6), round(z.imag, 2)) for z in near])
sp = np.diff(near.imag)
out['near_half_spacings'] = [float(x) for x in sp]
mins = []
for z in near:
    ts = np.linspace(z.imag - 0.05, z.imag + 0.05, 200001)
    v = np.abs(Dval(0.5 + 1j * ts))
    i = v.argmin()
    mins.append([float(ts[i]), float(v[i]), float(math.sqrt(1.5) / (3 * v[i]))])
out['min_abs_D_on_half_line_near_them'] = mins
print('min |D(1/2 + it)| and multiplier sqrt(3/2)/(3|D|) near them:', [(round(a, 3), f'{b:.3e}', round(c, 0)) for a, b, c in mins])
# almost period 2 pi 665/ln 2
P665 = 2 * math.pi * 665 / L
dd = []
for z in Zp:
    w = z + 1j * P665
    if w.imag < T:
        dd.append(np.abs(Zp - w).min())
out['shift_6028_median_distance'] = float(np.median(dd))
out['shift_6028_max_distance'] = float(np.max(dd))
out['shift_6028_n'] = len(dd)
print(f'shift by 2 pi 665/ln 2 = {P665:.4f}: median distance to the nearest zero {np.median(dd):.5f}, max {np.max(dd):.4f}, '
      f'{len(dd)} zeros', flush=True)
out['density_per_1000'] = len(Z) / (T / 1000)
json.dump(out, open(os.path.join(HERE, 'a2_zeros_torus_count.json'), 'w'), indent=1)
np.save(os.path.join(HERE, 'a2_zeros_half_one.npy'), Zp)
