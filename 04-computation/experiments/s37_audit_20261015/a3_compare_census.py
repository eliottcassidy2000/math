#!/usr/bin/env python3
"""S37 audit, part 3: compare the independent 2^34 census (a2_census34.json/.npz) with
  (i)  the session's census 05-knowledge/results/collatz_octave_memory_census_20261015.json/.npz,
  (ii) the earlier repository census 05-knowledge/results/collatz_discrepancy_20261013c.json (C++, 2026-10-13),
  (iii) the session's analysis table (collatz_octave_memory_analysis_20261015.json) and the note's section-2 table.
Also recomputes from the audit census: the transport e_4X vs (e + 2c + delta)/4, its exact residual R_X/(16X^2),
kappa = c/sqrt(e delta), the descent bound, and the 'agreement to 4 decimals' claim.
"""
import json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results'))
A = json.load(open(os.path.join(HERE, 'a2_census34.json')))
AZ = np.load(os.path.join(HERE, 'a2_census34.npz'))
SJ = json.load(open(os.path.join(RES, 'collatz_octave_memory_census_20261015.json')))
SZ = np.load(os.path.join(RES, 'collatz_octave_memory_census_20261015.npz'))
OLD = json.load(open(os.path.join(RES, 'collatz_discrepancy_20261013c.json')))
AN = json.load(open(os.path.join(RES, 'collatz_octave_memory_analysis_20261015.json')))

print(f"audit census: S(2^34) = {A['S_end']}, H(2^34) = {A['H_end']}; session: {SJ['S_end']}, {SJ['H_end']}")
# (ii) earlier census at powers of two
oldrows = {r['k']: r for r in OLD['rows']}
mis = []
for k in range(0, 35):
    s, h = A['pow2'][str(k)]
    r = oldrows[k]
    if r['S'] != s or r['residues'][2][0] != h:
        mis.append(k)
print(f"earlier repository census (2026-10-13 C++): S(2^k), H(2^k) for k = 0..34: mismatches at k = {mis}")
print(f"  earlier census row 34: S = {oldrows[34]['S']}, H_2 = {oldrows[34]['residues'][2][0]} "
      f"(the amplitude note states 24953018 and 26869)")
# earlier census prefix sums of squares -> dyadic energies over [X, 2X)
octA = {o['k']: o for o in A['octaves']}
maxrel_old = 0.0
for k in range(0, 34):
    X = 1 << k
    sq_blk = int(oldrows[k + 1]['residues'][2][3])     # old field = sum of H^2 over (2^k, 2^(k+1)]
    H_lo = oldrows[k]['residues'][2][0]
    H_hi = oldrows[k + 1]['residues'][2][0]
    E_old = sq_blk + H_lo ** 2 - H_hi ** 2               # -> [X, 2X)
    E_new = int(octA[k]['sumH2'])
    if E_old != E_new:
        maxrel_old = max(maxrel_old, abs(E_old - E_new) / max(1, E_new))
        print(f'  k = {k}: energy mismatch old {E_old} vs audit {E_new}')
print(f"  exact dyadic energies E_X, k = 0..33 (old census block sums over (X,2X] converted to [X,2X)): "
      f"max relative difference {maxrel_old}")

# (i) session census octave table
octS = {o['k']: o for o in SJ['octaves']}
worst = {}
for k in range(0, 34):
    a, s = octA[k], octS[k]
    X = 1 << k
    pairs = [('e', a['e'], s['energy_over_X2'][2]),
             ('mean_part', a['mean_part'], s['top_mean_part']),
             ('refinement', a['refinement'], s['top_innovation_total']),
             ('mean_over_sqrtX', a['sumH'] / X / math.sqrt(X), s['mean_over_sqrtX'][2])]
    if 'delta' in a and 'delta_energy_over_X2' in s:
        pairs += [('delta', a['delta'], s['delta_energy_over_X2']), ('cross', a['cross'], s['H_delta_cross_over_X2'])]
    for name, x, y in pairs:
        d = abs(x - y) / max(abs(x), 1e-300)
        if d > worst.get(name, (0, None))[0]:
            worst[name] = (d, k)
print('session census octave table vs audit census, k = 0..33 (delta, cross: k = 0..30): worst relative differences')
for name, (d, k) in worst.items():
    print(f'  {name:16s} {d:.3e} (at k = {k})')
print(f"  delta/cross rows present in session JSON: k <= {max(k for k in octS if 'delta_energy_over_X2' in octS[k])}")

# (i) session grid samples
gA, gS = AZ['grid'], SZ['grid']
print(f'grid identical: {np.array_equal(gA, gS)} ({len(gA)} points)')
C243 = AZ['C243']
S_A = C243.sum(axis=1)
C3 = np.stack([C243[:, r::3].sum(axis=1) for r in range(3)], axis=1)
C9 = np.stack([C243[:, r::9].sum(axis=1) for r in range(9)], axis=1)
C27 = np.stack([C243[:, r::27].sum(axis=1) for r in range(27)], axis=1)
print(f"  S: {np.array_equal(S_A, SZ['S'])}; C3: {np.array_equal(C3, SZ['C3'])}; C9: {np.array_equal(C9, SZ['C9'])}; "
      f"C27: {np.array_equal(C27, SZ['C27'])}")
print(f"  prefix sums of H at grid: {np.array_equal(AZ['prefix_sum'], SZ['prefix_sum'][:, 2])}; "
      f"prefix sums of H^2 (float64) max rel diff: {np.max(np.abs(AZ['prefix_sq'] - SZ['prefix_sq'][:, 2]) / np.maximum(1, AZ['prefix_sq'])):.2e}")

# (iii) transport, kappa, descent from the audit census
print('audit-census transport table: k, e_X, mean part, refinement, delta, kappa, e_4X, (e+2c+delta)/4, diff, '
      'descent bound, 4-decimal match')
bad4 = []
kap = []
for k in range(10, 31):
    a = octA[k]
    e, dl, cx = a['e'], a['delta'], a['cross']
    e4 = octA[k + 2]['e']
    pred = (e + 2 * cx + dl) / 4
    kp = cx / math.sqrt(e * dl)
    kap.append((kp, k))
    desc = ((math.sqrt(e) + math.sqrt(dl)) / 2 + math.sqrt(3 / 2 ** k)) ** 2
    m4 = round(e4, 4) == round(pred, 4)
    if not m4:
        bad4.append(k)
    print(f'  {k:2d} {e:.4f} {a["mean_part"]:.4f} {a["refinement"]:.4f} {dl:.4f} {kp:+.3f} {e4:.4f} {pred:.4f} '
          f'{e4 - pred:+.6f} {desc:.4f} {m4}')
print(f'  octaves 2^10..2^30 where e_4X and the transport prediction differ at 4 decimals: {bad4}')
print(f'  kappa range: min {min(kap)[0]:+.3f} (k = {min(kap)[1]}), max {max(kap)[0]:+.3f} (k = {max(kap)[1]})')
# compare with session analysis JSON table
anrows = {r['k']: r for r in AN['octaves']}
mx = 0
for k in range(10, 31):
    r = anrows[k]
    mx = max(mx, abs(r['kappa'] - kap[k - 10][0]), abs(r['transport_prediction'] - (octA[k]['e'] + 2 * octA[k]['cross'] + octA[k]['delta']) / 4))
print(f'  session analysis table vs audit (kappa, prediction): max abs diff {mx:.2e}')
# note's section-2 table
note = {18: (0.274, 0.222, 0.053, 1.604, 0.80, 0.733, 0.807), 20: (0.733, 0.701, 0.033, 0.213, 0.65, 0.365, 0.436),
        22: (0.365, 0.335, 0.031, 4.537, 0.84, 1.769, 1.872), 24: (1.769, 1.659, 0.110, 0.982, 0.60, 1.084, 1.348),
        26: (1.084, 0.982, 0.102, 2.835, -0.92, 0.173, 1.857), 28: (0.173, 0.028, 0.145, 1.181, -0.36, 0.256, 0.564),
        30: (0.256, 0.163, 0.094, 0.954, 0.20, 0.353, 0.550)}
nbad = 0
for k, vals in note.items():
    a = octA[k]
    e, dl, cx = a['e'], a['delta'], a['cross']
    mine = (e, a['mean_part'], a['refinement'], dl, cx / math.sqrt(e * dl), octA[k + 2]['e'],
            ((math.sqrt(e) + math.sqrt(dl)) / 2 + math.sqrt(3 / 2 ** k)) ** 2)
    for i, (x, y) in enumerate(zip(vals, mine)):
        if abs(x - y) > (0.0051 if i == 4 else 0.0016):     # kappa is printed to 2 decimals
            nbad += 1
            print(f'  note table k = {k}: {x} vs audit {y:.4f}')
print(f'  note section-2 table entries off by more than rounding: {nbad}')
# energy statistics claims
es = [(octA[k]['e'], k) for k in range(10, 34)]
print(f'  dyadic e_X, k = 10..33: max {max(es)[0]:.6f} at k = {max(es)[1]}; octaves with e_X < 0.26: '
      f'{[k for e, k in es if e < 0.26]}')
low_mean = [(k, round(octA[k]['mean_part'], 4)) for k in range(14, 34) if octA[k]['mean_part'] < 0.05]
print(f'  octaves 2^14..2^33 with mean part < 0.05 (no deterministic 0.26 floor): {low_mean}')
