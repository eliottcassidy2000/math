# Audit S83, step 2: exact |S_n(theta)| / 2^((n-1)/2) at the five irrationals and the energy E_n / |X_n|, n = 12..NMAX
# (every n, not only even n), with the independent engine. Also a cross-check of the session's float64 phase method.
import sys, time, math
import numpy as np
from fractions import Fraction
from oa_engine import g_block, S_exact, energy, IRR

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 26
names = list(IRR)
print('n | ' + ' '.join(f'{k:>8}' for k in names) + ' | E_n/|X_n|   distinct/|X_n|   max mult   (float64-method max abs diff of ratios)')
for n in range(12, NMAX + 1):
    t0 = time.time()
    vals = g_block(n)
    # sanity: H odd and <= 2^n
    assert vals.min() >= 1 and vals.max() <= (1 << n) and np.all(vals & 1)
    ratios = []
    for k in names:
        s = S_exact(vals, IRR[k])
        ratios.append(abs(s) / 2 ** ((n - 1) / 2))
    # session-style float64 phases for comparison
    diffs = []
    for k, r in zip(names, ratios):
        th = float(IRR[k]); acc = 0j
        for i in range(0, vals.size, 1 << 22):
            acc += np.exp(2j * np.pi * th * vals[i:i + (1 << 22)].astype(np.float64)).sum()
        diffs.append(abs(abs(acc) / 2 ** ((n - 1) / 2) - r))
    E, D, mx = energy(vals, n)
    X = 1 << (n - 1)
    print(f'{n:2d} | ' + ' '.join(f'{r:8.3f}' for r in ratios) + f' | {E / X:.4f}   {D / X:.4f}   {mx}   ({max(diffs):.1e})   [{time.time() - t0:.1f}s]', flush=True)
    del vals
