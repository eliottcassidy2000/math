# Audit2 (C): the two-loop formula on a generic Eulerian skew matrix (entries 0, +-1, zero row sums, unequal degrees;
# c_two_loop_check.eulerian_skew(7, 1, 3)) through z^K (argv[1]; z^20 in c_eulerian_L20.out, peak RSS 2.4 GB) by the blow-up scaling method M: here the dumbbell (with loops of
# lengths 2, 6, 10, ...) enters from z^6, so long dumbbell, theta and figure-eight kernels are all exercised.
import sys, time, resource
sys.path.insert(0, '.')
from b_mscaling import mscaling
from common import X0_series, X1_series
from c_two_loop_check import eulerian_skew
M = eulerian_skew(7, 1, 3)
print('degrees', (M != 0).sum(axis=1).tolist())
K = int(sys.argv[1]) if len(sys.argv) > 1 else 16
t = time.time()
X = mscaling(M, K, verbose=False)
print(f'time {time.time() - t:.1f}s, maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.0f} MB')
x0 = X0_series(M, K); x1 = X1_series(M, K, parts=True)
ok = True
for L in range(2, K + 1, 2):
    e0 = X[L].get(0, 0); e1 = X[L].get(1, 0); good = (e0 == x0[L]) and (e1 == x1[L][3]); ok &= good
    print(f"L={L:2d} X0 exact: {e0 == x0[L]}  X1 exact: {e1 == x1[L][3]}  X1 = {float(e1):+.10e}  [F8 {float(x1[L][0]):+.3e}, theta {float(x1[L][1]):+.3e}, dumbbell {float(x1[L][2]):+.3e}]")
print('ALL EXACT' if ok else 'MISMATCH')
