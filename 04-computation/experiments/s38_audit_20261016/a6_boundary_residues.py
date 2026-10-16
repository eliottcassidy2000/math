#!/usr/bin/env python3
"""S38 audit, item 5: the boundary residue estimates r_hat(rho) = -(3/2)^rho C_N(rho)/D'(rho) at the zeros of D.

C_N(s) = (2/3)^s [(F_N - F_L) + 2^(-s)(F_N - F_J)], L = floor((2N-1)/3), J = floor(L/2) (formula (4) of
collatz_amplitude_boundary_20261013c.md), with F_N - F_X evaluated from the exact signed block moments (blocks of 2^17,
first-order Taylor about the block centre; the partial block after the cut is omitted, as in the S38 float script).
Own implementation from the formula, not from the S38 code.  Checks:
(a) agreement with the residues of collatz_D_zeros_and_boundary_residues_20261016.json at all 295 zeros;
(b) stability of r_hat in N (N = 2^26 .. 2^40) at rho* (1054/665), s_0 (84/53) and other zeros;
(c) the conditional radius under the budget (2, 0) (same formula as the Arb certificate, in floats) at every zero for
    N = 2^40, and its scaling in N at rho* (predicted 2^42, 2^43, 2^44);
(d) the vertical-line medians of |C_N(sigma + it)| and the crossing with the tail bound at sigma = 0.85;
(e) max |r_hat| over the 295 zeros and over the Mercator family (Re > 0.98).
Usage: python a6_boundary_residues.py moments.npz zeros.json
"""
import sys, json, math, os
import numpy as np

MOM, ZJ = sys.argv[1:3]
HERE = os.path.dirname(os.path.abspath(__file__))
L2, L15 = math.log(2.0), math.log(1.5)
B = 1 << 17
v = np.load(MOM)['moments']
Nfull = len(v) * B
sm_all = v[:, 0].astype(np.float64)
mu1_all = (2 * v[:, 1].astype(np.float64) - (B - 1) * sm_all) / 2      # sum (n - centre) t
lo_all = np.arange(len(v), dtype=np.float64) * B
cen_all = lo_all + (B - 1) / 2
lc_all = np.log(cen_all)
J = json.load(open(ZJ))
rows = J['zeros_and_residues']
Z = np.array([complex(r['re'], r['im']) for r in rows])


def Dp(s):
    return -L2 * np.exp(-s * L2) + L15 * np.exp(s * L15) / 3


def shell_sums(s, N):
    """returns (F_N - F_L, F_N - F_J) from blocks starting at or after the cut (block-Taylor), plus t(N) N^(-s)"""
    logN = N.bit_length() - 1
    tN = (-1) ** logN
    Lc = (2 * N - 1) // 3
    Jc = Lc // 2
    jl = (Jc + 1 + B - 1) // B
    js = (Lc + 1 + B - 1) // B
    jN = N // B
    sl = slice(jl, jN)
    val = np.exp(-s * lc_all[sl]) * (sm_all[sl] - s * mu1_all[sl] / cen_all[sl])
    cs = np.cumsum(val[::-1])[::-1]
    last = tN * np.exp(-s * math.log(N))
    long_ = cs[0] + last
    short_ = (cs[js - jl] if js - jl < len(cs) else 0) + last
    return short_, long_


def CN(s, N):
    a, b = shell_sums(s, N)
    return np.exp(s * math.log(2 / 3)) * (a + np.exp(-s * L2) * b)


def rhat(s, N):
    return -np.exp(s * L15) * CN(s, N) / Dp(s)


def radius20(s, N, HN):
    sig = s.real
    eta = sig - 0.5
    q = 2 ** (-eta)
    Ssum = math.sqrt(2) / (1 - q)
    block = math.sqrt((1 - 2 ** (-2 * sig - 1)) / (2 * sig + 1))
    etail = abs(HN) / 3 * N ** (-sig) + abs(s) / 3 * N ** (-eta) * block * Ssum
    rtail = abs(s) / (2 * sig) * (1 - 1 / (2 * N)) ** (-sig - 1) * N ** (-sig)
    return 1.5 ** sig / abs(Dp(s)) * (etail + rtail)


out = {}
N40 = 1 << 40
HN40 = -55691
R = np.array([rhat(z, N40) for z in Z])
Rj = np.array([complex(r['r_re'], r['r_im']) for r in rows])
rel = np.abs(R - Rj) / np.maximum(np.abs(Rj), 1e-300)
out['a_max_rel_diff_vs_json'] = float(rel.max())
out['a_max_abs_diff_vs_json'] = float(np.abs(R - Rj).max())
print(f'(a) 295 residues: max relative difference to the JSON {rel.max():.2e}, max absolute {np.abs(R - Rj).max():.2e}', flush=True)
# (c) radii at N = 2^40
rad = np.array([radius20(z, N40, HN40) for z in Z])
meaningful = np.abs(R) > rad
out['c_n_zeros_with_abs_rhat_gt_radius20'] = int(meaningful.sum())
out['c_zeros_with_abs_rhat_gt_radius20'] = [[float(z.real), float(z.imag), float(abs(r)), float(rr)] for z, r, rr, m in zip(Z, R, rad, meaningful) if m]
out['c_median_radius_over_abs_rhat'] = float(np.median(rad / np.abs(R)))
out['c_radius_over_abs_rhat_by_re'] = {}
for lo_, hi_ in ((0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 0.98), (0.98, 1.0)):
    m = (Z.real > lo_) & (Z.real <= hi_)
    if m.any():
        out['c_radius_over_abs_rhat_by_re'][f'{lo_}-{hi_}'] = [int(m.sum()), float(np.median(rad[m] / np.abs(R[m])))]
print(f'(c) N = 2^40, budget (2,0): {meaningful.sum()} of 295 zeros have |r_hat| > radius; median radius/|r_hat| = '
      f'{np.median(rad / np.abs(R)):.3g}; by Re band (count, median ratio): {out["c_radius_over_abs_rhat_by_re"]}', flush=True)
# (e) max |r_hat|
out['e_max_abs_rhat'] = [float(np.abs(R).max()), [float(Z[np.abs(R).argmax()].real), float(Z[np.abs(R).argmax()].imag)]]
fam = Z.real > 0.98
out['e_family_abs_rhat_range'] = [float(np.abs(R[fam]).min()), float(np.abs(R[fam]).max())]
out['e_family_outside_0.005_0.02'] = int(((np.abs(R[fam]) < 0.005) | (np.abs(R[fam]) > 0.02)).sum())
print('(e) max |r_hat| =', out['e_max_abs_rhat'], '; Mercator family |r_hat| range', out['e_family_abs_rhat_range'],
      '; members outside [0.005, 0.02]:', out['e_family_outside_0.005_0.02'], 'of', int(fam.sum()), flush=True)
# (b) stability in N
sel = {'rho*': 0.99999921 + 3014.0202j, 's0': 0.99820922 + 240.2478j, '720.74': 0.98494 + 720.7401j,
       '1572.56': 0.94973 + 1572.556j, '22.25': 0.69392 + 22.2532j, '9042.06': 0.99999289 + 9042.0605j,
       '6268.29': 0.99805786 + 6268.2882j}
# H at 2^k from the census json values is needed only for the tiny |H_N| term; use |H_N| <= 2^(k/2+1) bound
stab = {}
for name, z0 in sel.items():
    z = Z[np.argmin(np.abs(Z - z0))]
    vals = []
    for k in range(26, 41):
        N = 1 << k
        r = rhat(z, N)
        vals.append([k, float(abs(r)), float(np.angle(r)), float(radius20(z, N, 2 ** (k / 2 + 1)))])
    stab[name] = dict(zero=[float(z.real), float(z.imag)], rows=vals)
    print(f'(b) {name} {z:.8f}: k, |r_hat|, arg, radius(2,0):', flush=True)
    for row in vals:
        print(f'      2^{row[0]}: {row[1]:.6f}  {row[2]:+.4f}  {row[3]:.5f}')
out['b_stability'] = stab
zs = Z[np.argmin(np.abs(Z - (0.99999921 + 3014.0202j)))]
ab = [abs(rhat(zs, 1 << k)) for k in range(32, 41)]
out['b_rhostar_abs_rhat_2p32_2p40_min_max'] = [float(min(ab)), float(max(ab)), float((max(ab) - min(ab)) / np.mean(ab))]
pred = {k: radius20(zs, float(2 ** k), 2 ** (k / 2 + 1)) for k in (40, 41, 42, 43, 44, 45)}
out['c_rhostar_radius_by_N'] = pred
print('(b) rho*: |r_hat| over N = 2^32..2^40 min/max/relative spread:', out['b_rhostar_abs_rhat_2p32_2p40_min_max'])
print('(c) rho*: radius under (2,0) at N = 2^40..2^45:', {k: round(x, 5) for k, x in pred.items()}, flush=True)
# (d) vertical lines (same t-grid as the S38 script), recursion in t
ts = np.arange(10.0, J['T'], 37.0)
N = N40
Lc = (2 * N - 1) // 3
Jc = Lc // 2
jl = (Jc + 1 + B - 1) // B
js = (Lc + 1 + B - 1) // B
sl = slice(jl, N // B)
lc = lc_all[sl]
cen = cen_all[sl]
smv = sm_all[sl]
muv = mu1_all[sl]
ms = np.zeros(len(lc), bool)
ms[js - jl:] = True
base = np.exp(-1j * ts[0] * lc)
step = np.exp(-1j * 37.0 * lc)
sigs = (0.55, 0.7, 0.85, 1.0)
amp = {sg: np.exp(-sg * lc) for sg in sigs}
vals = {sg: [] for sg in sigs}
z = base.copy()
for k, t in enumerate(ts):
    if k % 64 == 0:
        z = np.exp(-1j * t * lc)          # re-anchor to limit rounding drift
    for sg in sigs:
        s = sg + 1j * t
        term = amp[sg] * z * (smv - s * muv / cen)
        lsum = term.sum() + np.exp(-s * math.log(N))
        ssum = term[ms].sum() + np.exp(-s * math.log(N))
        vals[sg].append(abs(np.exp(s * math.log(2 / 3)) * (ssum + np.exp(-s * L2) * lsum)))
    z = z * step
bands = ((10, 1000), (1000, 4000), (4000, 10000), (10000, J['T']))
dres = {}
for sg in sigs:
    vv = np.array(vals[sg])
    q = [float(np.median(vv[(ts >= a) & (ts < b)])) for a, b in bands]
    dres[str(sg)] = q
    print(f'(d) median |C_N({sg} + it)| by band: {[round(x, 4) for x in q]}  (S38: {J["vertical_lines"][str(sg)]["median_abs_numer_by_height_band"]})', flush=True)
out['d_vertical_medians'] = dres
# tail bound of the numerator at sigma = 0.85 under (2,0): |s|/3 N^(-eta) block S
sg = 0.85
eta = sg - 0.5
tb = lambda t: abs(complex(sg, t)) / 3 * N ** (-eta) * math.sqrt((1 - 2 ** (-2 * sg - 1)) / (2 * sg + 1)) * math.sqrt(2) / (1 - 2 ** (-eta))
vv = np.array(vals[0.85])
cross = [float(t) for t, x in zip(ts, vv) if tb(t) > x]
out['d_tail_bound_sigma085_at_t200'] = tb(200.0)
out['d_first_t_where_tail_bound_exceeds_numer_sigma085'] = cross[0] if cross else None
out['d_fraction_of_t_above_200_where_bound_exceeds'] = float(np.mean([tb(t) > x for t, x in zip(ts, vv) if t >= 200]))
print('(d) sigma = 0.85: tail bound at t = 200:', round(tb(200.0), 5), '; first grid t with bound > |C_N|:',
      out['d_first_t_where_tail_bound_exceeds_numer_sigma085'], '; fraction of t >= 200 with bound > |C_N|:',
      round(out['d_fraction_of_t_above_200_where_bound_exceeds'], 3), flush=True)
json.dump(out, open(os.path.join(HERE, 'a6_boundary_residues.json'), 'w'), indent=1)
