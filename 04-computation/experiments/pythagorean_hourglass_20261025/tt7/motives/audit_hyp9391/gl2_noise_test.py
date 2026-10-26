# Audit of HYP-9391 (c): could a GL(2) motive (|x| <= 2) plus the projection error produce the observed exceedances?
# Model: x = 2 cos(theta) (Sato-Tate SU(2)) + e, e ~ N(0, s^2) with s = 0.19 (max white-noise SD of the fitted known
# part at q >= 2000, part_c_analysis.out) and the pessimistic s = 0.3.  Compare the expected counts of |x| > 2, 2.5, 3
# among the 464 primes >= 2000 with the observed counts for m1..m4 (fit on [307, 9999], mE.npy).
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(sys.argv[1]); q = d[:, 0]; m = d[:, 1:]
sel = q >= 2000; n = int(sel.sum())
rng = np.random.default_rng(1)
N = 4_000_000
th = rng.uniform(0, np.pi, N * 3)
acc = rng.uniform(0, 1, N * 3) < np.sin(th) ** 2
x0 = 2 * np.cos(th[acc])[:N]
out = [f"{n} primes >= 2000"]
for s in (0.19, 0.3):
    x = x0 + rng.normal(0, s, len(x0))
    out.append(f"GL(2) + N(0, {s}^2): expected counts among {n}: |x|>2: {n*np.mean(np.abs(x)>2):.1f}, |x|>2.5: {n*np.mean(np.abs(x)>2.5):.3f}, |x|>3: {n*np.mean(np.abs(x)>3):.2e}")
for j in range(4):
    a = np.abs(m[sel, j])
    out.append(f"observed m{j+1}: |x|>2: {(a>2).sum()}, |x|>2.5: {(a>2.5).sum()}, |x|>3: {(a>3).sum()}, max {a.max():.3f}")
print('\n'.join(out))
open(os.path.join(HERE, 'gl2_noise_test.out'), 'w').write('\n'.join(out) + '\n')
