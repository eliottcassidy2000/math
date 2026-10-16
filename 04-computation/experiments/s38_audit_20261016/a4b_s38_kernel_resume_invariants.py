#!/usr/bin/env python3
"""S38 audit, item 3: the S38 census code itself, at small scale.

(1) The descent kernel tval of collatz_census_extension_2p40_20261015.py (used by the dense script), with a parity
    table to 2^LOW = 2^30, on the path-record starting values in [2^37, 2^42) and on every point of their trajectories
    above 2^40 (the only inputs that exercise the two-limb 128-bit branch), against Python big-integer trajectories.
(2) Fault injection into the per-block invariant check block_invariants of the dense script: single bit flips in each
    accumulator (class counts, register total, sub-block sums of t, first moments, H increment, sums of hp and hp^2);
    which flips are detected.
(3) Resume logic: the dense script run fresh to 2^26 and split (2^24, then resumed to 2^26) as subprocesses
    (NUMBA_NUM_THREADS = 2, LOW = 22, 1024 points per octave), compared with each other and with the independent census
    of a4_independent_census34.py (a4_census34_grid.npz) at every grid point in [2^24, 2^26).
Usage: python a4b_s38_kernel_resume_invariants.py REPO_ROOT TMPDIR
"""
import sys, os, json, subprocess, importlib.util, random
import numpy as np

ROOT, TMP = sys.argv[1:3]
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(ROOT, '04-computation', 'experiments')
os.makedirs(TMP, exist_ok=True)
out = {}
# (1) kernel
spec = importlib.util.spec_from_file_location('c40', os.path.join(EXP, 'collatz_census_extension_2p40_20261015.py'))
c40 = importlib.util.module_from_spec(spec)
argv = sys.argv
sys.argv = ['x', '30', '26']
spec.loader.exec_module(c40)
sys.argv = argv
LOW = 30
P, patched = c40.ext.parity_bits(LOW)
A, Cc = c40.jump_tables()
low = 1 << LOW


def tpar(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return 1 - 2 * (s & 1)


def maxval(n):
    m = n
    while n >= low:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        m = max(m, n)
    return m


records = [204430613247, 231913730799, 272025660543, 446559217279, 567839862631, 871673828443, 2674309547647,
           3716509988199]
tests = []
for r in records:
    x = r
    while x >= (1 << 37) and len(tests) < 6000:
        if x < (1 << 63):
            tests.append(x)
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
bad = 0
twolimb = 0
for n in tests:
    if maxval(n) >= 2 ** 64:
        twolimb += 1
    got = c40.tval(P, n, low, A, Cc, c40.JGUARD, c40.GUARD)
    if got != tpar(n):
        bad += 1
out['kernel_tests'] = len(tests)
out['kernel_tests_with_values_beyond_2p64'] = twolimb
out['kernel_mismatches'] = bad
print(f'(1) tval (LOW = 30) on {len(tests)} record-trajectory points >= 2^37 below 2^63 '
      f'({twolimb} of them pass values >= 2^64): {bad} mismatches', flush=True)
# (2) fault injection into block_invariants
spec2 = importlib.util.spec_from_file_location('dense', os.path.join(EXP, 'collatz_census_dense_resumable_20261016.py'))
dense = importlib.util.module_from_spec(spec2)
sys.argv = ['x', '26', '22', '1024', os.path.join(TMP, 'unused')]
spec2.loader.exec_module(dense)
sys.argv = argv
P22, _ = c40.ext.parity_bits(22)
cnt, mom, hs1, hs2, hend, mx, badb, tott = dense.heavy_range(P22, 3, 7, 1 << 20, 1 << 17, 243, 1 << 22, A, Cc, c40.JGUARD, c40.GUARD)
nsub = 8
base = dense.block_invariants(cnt, mom, hend, tott, nsub)
det = {}
rng = random.Random(7)
for name in ('cnt', 'tott', 'mom0', 'mom1', 'hend', 'hs1', 'hs2'):
    hits = 0
    for trial in range(20):
        c_, m_, he_, tt_ = cnt.copy(), mom.copy(), hend.copy(), tott.copy()
        bit = np.int64(1) << np.int64(rng.randrange(0, 62))
        i = rng.randrange(0, cnt.shape[0])
        if name == 'cnt':
            c_[i, rng.randrange(243)] ^= bit
        elif name == 'tott':
            tt_[i] ^= bit
        elif name == 'mom0':
            m_[i * nsub + rng.randrange(nsub), 0] ^= bit
        elif name == 'mom1':
            m_[i * nsub + rng.randrange(nsub), 1] ^= bit
        elif name == 'hend':
            he_[i] ^= bit
        # hs1, hs2 are not inputs of block_invariants: a flip there can never be detected
        inv = dense.block_invariants(c_, m_, he_, tt_, nsub)
        hits += int(inv.any())
    det[name] = f'{hits}/20 detected'
out['baseline_inconsistent_blocks'] = int(base.sum())
out['fault_injection'] = det
print('(2) block_invariants: baseline', int(base.sum()), 'inconsistent; single-bit flips detected:', det, flush=True)
# (3) resume test
env = dict(os.environ, NUMBA_NUM_THREADS='2')
script = os.path.join(EXP, 'collatz_census_dense_resumable_20261016.py')
fa = os.path.join(TMP, 'fresh26')
fb = os.path.join(TMP, 'split24')
fc = os.path.join(TMP, 'resumed26')
for args in (['26', '22', '1024', fa], ['24', '22', '1024', fb], ['26', '22', '1024', fc, '--resume', fb + '_state.json']):
    r = subprocess.run([sys.executable, script] + args, env=env, cwd=ROOT, capture_output=True, text=True)
    print('   run', ' '.join(args[:3]), args[3][-10:], 'rc', r.returncode, (r.stdout.strip().splitlines() or [''])[-2:], flush=True)
    if r.returncode:
        print(r.stderr[-2000:])
A_ = np.load(fa + '.npz')
C_ = np.load(fc + '.npz')
ga, gc = np.array(A_['grid']), np.array(C_['grid'])
sel = (ga >= (1 << 24)) & (ga < (1 << 26))
same = bool(np.array_equal(ga[sel], gc) and np.array_equal(np.array(A_['S'])[sel], np.array(C_['S'])) and
            np.array_equal(np.array(A_['H'])[sel], np.array(C_['H'])) and
            np.array_equal(np.array(A_['C3'])[sel], np.array(C_['C3'])) and
            list(np.array(A_['preH2_str'])[sel]) == list(np.array(C_['preH2_str'])) and
            list(np.array(A_['preH1_str'])[sel]) == list(np.array(C_['preH1_str'])))
sa = json.load(open(fa + '_state.json'))
sc = json.load(open(fc + '_state.json'))
out['resume_grid_values_equal_fresh'] = same
out['resume_final_state_equal_fresh'] = (sa == sc)
mine = np.load(os.path.join(HERE, 'a4_census34_grid.npz'))
gm = np.array(mine['grid'])
idx = np.searchsorted(gm, gc)
eq_mine = bool(np.array_equal(gm[idx], gc) and np.array_equal(np.array(mine['S'])[idx], np.array(C_['S'])) and
               np.array_equal(np.array(mine['H'])[idx], np.array(C_['H'])) and
               list(np.array(mine['preH2_str'])[idx]) == list(np.array(C_['preH2_str'])))
out['resumed_equal_independent_census'] = eq_mine
ma = np.load(fa + '_moments.npz')
mc = np.load(fc + '_moments.npz')
out['resumed_moments_equal_fresh_rows'] = bool(np.array_equal(np.array(ma['moments'])[(1 << 24) >> 17:], np.array(mc['moments'])))
print('(3) resume: grid values of the resumed run equal the fresh run:', same, '; final states equal:', sa == sc,
      '; equal to the independent census:', eq_mine, '; moments equal:', out['resumed_moments_equal_fresh_rows'], flush=True)
json.dump(out, open(os.path.join(HERE, 'a4b_s38_kernel_resume_invariants.json'), 'w'), indent=1)
