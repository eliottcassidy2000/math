# AUDIT task 2 extension: the 6-vertex excess-3 cores (9 edges) first contribute at z^10, and the validation in
# excess3_validate.out reaches z^12 only at n = 7, i.e. weight-tensor entries with <= 3 multi-segment edges.  Here the
# engine series of X_3 is compared with the exact blow-up scaling up to z^16 on the non-circulant regular 7-tournament
# (gentourng class 1, as in audit2/c_dumbbell_n7_L16.py), with the share of the 6-vertex cores printed per order.
import sys, time, resource
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE); sys.path.insert(0, HERE + '/audit2')
import loop_engine4 as LE4
from b_mscaling import mscaling
from common import skew, code_to_A, regular_classes
S = skew(code_to_A(regular_classes(7)[1])); n = 7; Sh = S / n
print("tournament:", regular_classes(7)[1], " max |x_k| =", np.max(np.abs(np.linalg.eigvals(Sh).imag)))
t = time.time(); X = mscaling(S, 16, verbose=False); print(f"m-scaling to z^16: {time.time() - t:.0f}s, maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.0f} MB", flush=True)
t = time.time(); co = LE4.series(Sh, 3, 16, K=64, r=1.0, per_core=True); print(f"engine series (K=64, r=1.0): {time.time() - t:.0f}s, maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.0f} MB", flush=True)
six = [i for i, (V, e, a) in enumerate(LE4.load_cores(3)) if V == 6]
worst = 0
for L in range(2, 17, 2):
    ex = X[L].get(3, 0); tot = co[L].sum(); s6 = co[L][six].sum()
    worst = max(worst, abs(tot - float(ex)) / max(1e-300, abs(float(ex)) if ex else 1))
    print(f"   L={L:2d}: exact X_3^(L) {float(ex):+.15e}  engine {tot:+.15e}  diff {tot - float(ex):+.1e}   (6-vertex cores: {s6:+.4e})")
print(f"max relative deviation (absolute where exact = 0): {worst:.1e}")
