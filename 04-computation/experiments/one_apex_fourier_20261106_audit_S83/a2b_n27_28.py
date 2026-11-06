# Audit S83, step 2b: n = 27, 28 (fixed-point phases, numba energy; ~2 GB peak at n = 28).
import time
from oa_engine import g_block, S_exact, energy, IRR
names = list(IRR)
for n in (27, 28):
    t0 = time.time()
    vals = g_block(n)
    ratios = [abs(S_exact(vals, IRR[k])) / 2 ** ((n - 1) / 2) for k in names]
    E, D, mx = energy(vals, n)
    X = 1 << (n - 1)
    print(f'{n:2d} | ' + ' '.join(f'{k}={r:.3f}' for k, r in zip(names, ratios)) + f' | E/|X| = {E / X:.4f}  distinct/|X| = {D / X:.4f}  max mult {mx}  [{time.time() - t0:.1f}s]', flush=True)
    del vals
