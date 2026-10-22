# Audit2 (C): the dumbbell term of THM-4757 (4) on the non-circulant (non-vertex-transitive) regular 7-tournament.
# For every regular tournament d_2 = 2 diag(Shat^2) is constant and G 1 = 0, so the dumbbell has no z^L term below
# L = 14.  Method M (blow-up scaling, b_mscaling.py) gives the exact excess-1 cluster sums through z^16; compare with the
# closed formula (exact rationals).  Columns: L, X0 exact?, X1 exact?, X1, [figure-eight, theta, dumbbell, total].
import sys, time, resource
sys.path.insert(0, '.')
from b_mscaling import mscaling
from common import skew, X1_series, X0_series, code_to_A, regular_classes
codes = regular_classes(7)
S = skew(code_to_A(codes[1]))
print('tournament (gentourng code):', codes[1])
t = time.time()
X = mscaling(S, 16, verbose=False)
print(f'time {time.time() - t:.1f}s, maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.0f} MB')
x1 = X1_series(S, 16, parts=True); x0 = X0_series(S, 16)
for L in range(2, 17, 2):
    print(f"L={L:2d} X0 exact: {X[L].get(0, 0) == x0[L]}  X1 exact: {X[L].get(1, 0) == x1[L][3]}  X1 = {float(X[L].get(1, 0)):+.12e}  "
          f"[F8 {float(x1[L][0]):+.4e}, theta {float(x1[L][1]):+.4e}, dumbbell {float(x1[L][2]):+.4e}]  exact dumbbell = {x1[L][2]}")
