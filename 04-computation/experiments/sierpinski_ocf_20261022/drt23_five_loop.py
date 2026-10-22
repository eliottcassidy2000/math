# Five loops for the doubly regular 23-tournaments (HYP-9371 follow-up).  Sign-network engine (loop_engine5.py), DRT
# two-label mode: per core, W2[z, t] over {P, Shat} edge labels at 64 Cauchy points on |z| = 1.5 and at z = 1; per class,
# Gamma2[t] ({P, Shat} homomorphism numbers, loop_engine4.Network); X(z) = sum_cores <W2(z), Gamma2>/|G|.
# Classes: the 13 of drt23_classes.json and the doubled-Paley class of drt23_doubled.json (14), relative to P_23.
#   excess 4: all 1076 cores (cores with <= 4 vertices must be universal: checked);
#   excess 5: the 2046 cores with V <= 5 (cores_e5_v5.json); only its z^10 coefficient is complete (larger excess-5
#             cores have >= 11 edges, so they start at z^11).
# The decisive identity: at z^10 only excess 3, 4 and 5 can differ between doubly regular tournaments of one order
# (excess >= 6 has <= 4 blocks at L = 10), so for every class
#     Delta c_10 (exact, forest numbers) = delta_10 (four loops) + Delta X_4^(10) + Delta X_5^(10).
# Output: drt23_five_loop.out, drt23_five_loop.json.  mac-mini-2026-10-22-loop4.
import sys, json, time
from fractions import Fraction
from math import factorial
import numpy as np
sys.path.insert(0, '.')
import loop_engine4 as LE4
import loop_engine5 as LE5
from drt_switch import code_to_S

n = 23; K = 64; R = 1.5
zs = np.concatenate([R * np.exp(2j * np.pi * np.arange(K) / K), [1.0]])
ops = LE5.ops_drt(n, zs)
classes = json.load(open('drt23_classes.json'))
dbl = json.load(open('drt23_doubled.json'))
labels = [c['label'] for c in classes] + ['D23_dbl']
codes = [c['code'] for c in classes] + [dbl['code']]
TT5 = [c['TT5'] for c in classes] + [dbl['TT5']]
Hs = [int(c['H']) for c in classes] + [int(dbl['H'])]
iP = labels.index('P_23')
Ss = [code_to_S(c) for c in codes]
assert all(LE4.check_drt(S) for S in Ss)
P = np.eye(n) - np.ones((n, n)) / n
Bs = [np.array([P, S / n]) for S in Ss]
NC = len(Ss)

def run(excess, vmax=None, vmin=1):
    cores = LE5.load_cores(excess, vmax)
    vals = {}                                              # V -> (NC, K+1)
    t0 = time.time()
    for i, (V, edges, aut) in enumerate(cores):
        if V < vmin: continue
        W2 = LE5.SignCore(V, edges, aut).contract('drt', ops).reshape(len(zs), -1) / aut
        net = LE4.Network(V, edges)
        G = np.array([net.contract(B).reshape(-1) for B in Bs])      # (NC, 2^E)
        vals.setdefault(V, np.zeros((NC, len(zs)), dtype=complex))
        vals[V] += G @ W2.T
        if (i + 1) % 200 == 0: print(f"   excess {excess}: {i + 1}/{len(cores)} cores, {time.time() - t0:.0f}s", flush=True)
    return vals

def coeffs(v):
    """z-series coefficients L = 0..40 from the K Cauchy values; value at z = 1"""
    c = np.array([(np.mean(v[:K] * np.exp(-2j * np.pi * np.arange(K) * L / K)) / R ** L).real for L in range(41)])
    return c, v[K].real

out = {'labels': labels, 'TT5': TT5}
for excess, vmax, vmin in ((3, None, 1), (4, None, 1), (5, 5, 5)):
    t0 = time.time()
    vals = run(excess, vmax, vmin)
    tot = sum(vals.values())
    ncores = sum(1 for c in LE5.load_cores(excess, vmax) if c[0] >= vmin)
    if vmin == 1:
        small = sum(v for V, v in vals.items() if V <= 4)
        cs_small = [coeffs(small[k] - small[iP])[0] for k in range(NC)]
        print(f"excess {excess}: {ncores} cores in {time.time() - t0:.0f}s; cores with V <= 4 universal: max class difference "
              f"{max(np.max(np.abs(c)) for c in cs_small):.1e}", flush=True)
    else:
        print(f"excess {excess}: the {ncores} cores with V = {vmin} in {time.time() - t0:.0f}s (only differences are meaningful)", flush=True)
    res = []
    for k in range(NC):
        c, v1 = coeffs(tot[k]); cP, v1P = coeffs(tot[iP])
        res.append({'coeffs': c.tolist(), 'value': v1, 'dcoeffs': (c - cP).tolist(), 'dvalue': v1 - v1P,
                    'byV': {str(V): (coeffs(vv[k] - vv[iP])[1]) for V, vv in vals.items()}})
    out[f'X{excess}'] = res
    cP, v1P = coeffs(tot[iP])
    if vmin == 1: print(f"   P_23: X_{excess}(1) = {v1P:+.10e}" + (f", n^{excess} X_{excess} = {n ** excess * v1P:+.4f}" if vmax is None else '') +
          "; series " + ' '.join(f"z^{L}:{cP[L]:+.3e}" for L in range(6, 17, 2)), flush=True)

# exact log Phi coefficients (forest numbers) for all 14 classes
ex = {r['label']: [int(x) for x in r['N']] for r in json.load(open('drt23_forest_exact.json'))}
ex['D23_dbl'] = [int(x) for x in dbl['N']]
def logc(N, LM=24):
    phi = [Fraction(N[L], n ** L) if L < n else Fraction(0) for L in range(LM + 1)]
    c = [Fraction(0)] * (LM + 1)
    for L in range(1, LM + 1):
        s = L * phi[L]
        for j in range(1, L): s -= j * c[j] * phi[L - j]
        c[L] = s / L
    return c
cex = {lab: logc(ex[lab]) for lab in labels}
print("\nThe z^10 identity: Delta c_10 (exact) = delta_10 (4 loops) + Delta X_4^(10) + Delta X_5^(10) (V = 5 cores)")
print(f"{'class':8s} {'TT5':>5s} {'exact Delta c_10':>17s} {'delta_10':>13s} {'Delta X_4^(10)':>15s} {'Delta X_5^(10)':>15s} {'sum':>17s} {'exact - sum':>12s}")
rows = []
for k, lab in enumerate(labels):
    if k == iP: continue
    dc10 = float(cex[lab][10] - cex['P_23'][10]); dc8 = float(cex[lab][8] - cex['P_23'][8])
    d3 = out['X3'][k]['dcoeffs']; d4 = out['X4'][k]['dcoeffs']; d5 = out['X5'][k]['dcoeffs']
    s = d3[10] + d4[10] + d5[10]
    rows.append((lab, TT5[k], dc10, d3[10], d4[10], d5[10], s, dc8, d3[8], d4[8], d5[8]))
    print(f"{lab:8s} {TT5[k]:5d} {dc10:+.10e} {d3[10]:+.6e} {d4[10]:+.8e} {d5[10]:+.8e} {s:+.10e} {dc10 - s:+.1e}")
print("   (z^8 check: exact Delta c_8 - delta_8, and the excess-4/5 z^8 coefficients, which must vanish: max "
      f"{max(abs(r[7] - r[8]) for r in rows):.1e}, {max(abs(r[9]) for r in rows):.1e}, {max(abs(r[10]) for r in rows):.1e})")
json.dump(out, open('drt23_five_loop.json', 'w'))
