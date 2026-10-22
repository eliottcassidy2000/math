# Audit2 (C): the session's own test cases (R_15, P_19; z^4, z^6), now exactly: brute-force excess-0/1 cluster sums
# (set partitions, b_partitions.py) against the closed formulas, as exact rationals, printed to 12 digits.
import sys
sys.path.insert(0, '.')
from b_partitions import cluster_sums
from common import rot_A, paley_A, skew, X0_series, X1_series
for label, A in (('R_15', rot_A(15)), ('P_19', paley_A(19))):
    S = skew(A)
    cs = cluster_sums(S, 8, regular=True, emax=1, odd=False, Lmin=2, verbose=False)
    x0 = X0_series(S, 8); x1 = X1_series(S, 8, parts=True)
    for L in (2, 4, 6, 8):
        e0 = cs[L].get(0, 0); e1 = cs[L].get(1, 0)
        print(f"{label} z^{L}: excess-0 brute {float(e0):+.12e} formula {float(x0[L]):+.12e} exact={e0 == x0[L]};  "
              f"excess-1 brute {float(e1):+.12e} formula {float(x1[L][3]):+.12e} exact={e1 == x1[L][3]}", flush=True)
