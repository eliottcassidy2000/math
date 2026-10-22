# Audit2 (C): set-partition brute force of the excess-1 cluster sums at z^10 and z^12 for the regular 9-tournament
# gentourng class 3 (non-vertex-transitive), compared exactly with the THM-4757 (4) formula.  Output: c_partition_L12.out
import sys, time
sys.path.insert(0, '.')
from b_partitions import cluster_sums
from common import rot_A, skew, X1_series, X0_series, code_to_A, regular_classes
import numpy as np
codes = regular_classes(9)
A = code_to_A(codes[3]); S = skew(A)
t = time.time()
cs = cluster_sums(S, 12, regular=True, emax=1, odd=False, Lmin=10)
print('time', time.time() - t)
x1 = X1_series(S, 12, parts=True)
for L in (10, 12):
    print(L, cs[L].get(1, 0) == x1[L][3], float(cs[L].get(1, 0)), [float(v) for v in x1[L]])
