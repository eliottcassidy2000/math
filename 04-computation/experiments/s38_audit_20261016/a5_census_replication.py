#!/usr/bin/env python3
"""S38 audit, item 3: replication checks of the dense 2^40 census against the S37 2^40 census (same kernels, separate run).

(a) all grid points shared with the S37 64-per-octave census: class sums mod 243, H, prefix sums of H and H^2 (strings);
(b) exact octave sums of H and H^2 (k = 20..39) of both runs;
(c) the full signed sub-block moments (2^23 rows x 2: sum t and sum (n - lo) t over every block of 2^17) of both runs;
(d) internal consistency of the resume state dense40_state.json (class counts below 2^40, H(2^40 - 1), sum H, sum H^2)
    against dense40.npz / dense40.json and S(2^40), H(2^40);
(e) S, H at the dense grid points re-derived from the stored class sums of dense40_full.npz.
Usage: python a5_census_replication.py SCRATCH_S38_DIR S37_MOMENTS_NPZ REPO_ROOT
"""
import sys, os, json, hashlib
import numpy as np

SD, M37, ROOT = sys.argv[1:4]
out = {}
d = np.load(os.path.join(SD, 'dense40.npz'))
dj = json.load(open(os.path.join(SD, 'dense40.json')))
st = json.load(open(os.path.join(SD, 'dense40_state.json')))
c37 = np.load(os.path.join(ROOT, '05-knowledge/results/collatz_census_extension_2p40_20261015.npz'))
j37 = json.load(open(os.path.join(ROOT, '05-knowledge/results/collatz_census_extension_2p40_20261015.json')))
full = np.load(os.path.join(SD, 'dense40_full.npz'))

g = np.array(d['grid'])
S = np.array(d['S'])
H = np.array(d['H'])
C3 = np.array(d['C3'])
p1 = [int(x) for x in d['preH1_str']]
p2 = [int(x) for x in d['preH2_str']]
gf = np.array(full['grid'])
cf = np.array(full['classsums'])
assert np.array_equal(g, gf)
# (e) S, C3, H from the class sums
S_e = cf.sum(axis=1)
C3_e = np.stack([cf[:, a::3].sum(axis=1) for a in range(3)], axis=1)
H_e = 3 * C3_e[:, 2] - S_e
out['e_S_from_classsums'] = bool(np.array_equal(S_e, S))
out['e_C3_from_classsums'] = bool(np.array_equal(C3_e, C3))
out['e_H_from_classsums_equals_chained_H'] = bool(np.array_equal(H_e, H))
out['n_dense_points'] = int(len(g))
# (a) shared points
g0 = np.array(c37['grid'])
idx = np.searchsorted(g, g0)
ok = bool(np.array_equal(g[idx], g0))
out['a_all_S37_points_in_dense_grid'] = ok
out['a_classsums_equal'] = bool(np.array_equal(cf[idx], np.array(c37['classsums'])))
out['a_H_equal'] = bool(np.array_equal(H[idx], np.array(c37['H'])))
q1 = [int(x) for x in c37['preH1_str']]
q2 = [int(x) for x in c37['preH2_str']]
out['a_preH1_equal'] = all(p1[int(i)] == q for i, q in zip(idx, q1))
out['a_preH2_equal'] = all(p2[int(i)] == q for i, q in zip(idx, q2))
out['a_n_shared'] = int(len(g0))
# (b) octave sums
od, o37 = dj['octave_sums_H_H2'], j37['octave_sums_H_H2']
out['b_octave_sums_equal'] = all(od[k] == o37[k] for k in o37) and set(od) == set(o37)
out['b_S_H_2p40'] = [dj['S'], dj['H'], j37['S'], j37['H']]
# (c) moments
m40 = np.load(os.path.join(SD, 'dense40_moments.npz'))
mm = np.load(M37)
A, B = np.array(m40['moments']), np.array(mm['moments'])
out['c_shapes'] = [list(A.shape), list(B.shape)]
out['c_block0_dense'] = int(m40['block0'])
eq = (A == B)
out['c_rows_differing'] = int((~eq.all(axis=1)).sum())
out['c_moments_identical'] = bool(eq.all())
out['c_sha256_dense40_moments_array'] = hashlib.sha256(A.tobytes()).hexdigest()
out['c_sha256_s37_moments_array'] = hashlib.sha256(B.tobytes()).hexdigest()
out['c_sum_t_plus_tN'] = int(A[:, 0].sum()) + 1
# (d) state
cnt = [int(x) for x in st['class_counts_below']]
N = 1 << 40
tN = 1
S_below = sum(cnt)
H_last = 3 * sum(cnt[2::3]) - S_below
cnt2 = list(cnt)
cnt2[N % 243] += tN
S_N = sum(cnt2)
H_N = 3 * sum(cnt2[2::3]) - S_N
out['d_state_S_at_2k_from_counts'] = [S_N, int(st['S_at_2k'])]
out['d_state_H_at_2k_from_counts'] = [H_N, int(st['H_at_2k'])]
out['d_state_H_last_from_counts'] = [H_last, int(st['H_last'])]
out['d_state_S_equals_moment_total'] = (S_N == out['c_sum_t_plus_tN'])
# sum H below 2^40 = prefix through 2^20 - H(2^20) + octave sums 20..39
i20 = int(np.nonzero(g == (1 << 20))[0][0])
sH = p1[i20] - int(H[i20]) + sum(int(od[str(k)][0]) for k in range(20, 40))
sH2 = p2[i20] - int(H[i20]) ** 2 + sum(int(od[str(k)][1]) for k in range(20, 40))
out['d_sumH_below'] = [sH, int(st['sumH_below']), sH == int(st['sumH_below'])]
out['d_sumH2_below'] = [sH2, int(st['sumH2_below']), sH2 == int(st['sumH2_below'])]
# last dense grid point below 2^40 against the state: walk-free check S(g_last) <= consistency only
out['d_last_grid_point'] = int(g[-1])
for k, v in out.items():
    print(k, v)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a5_census_replication.json'), 'w'), indent=1)
