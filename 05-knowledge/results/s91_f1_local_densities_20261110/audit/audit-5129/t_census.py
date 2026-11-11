# audit-5129: test of the fast order-10 counter of step 7 against the routines of acore (validated literally in
# step 1), and a timing estimate.
import time
import numpy as np
import acore as A
from s7_order10 import census_one, pass1

rng = np.random.default_rng(7)
n = 10
LOG = np.zeros(1024, np.int64)
for i in range(10):
    LOG[1 << i] = i
adj = np.zeros(n, np.int64)
hp = np.zeros((1 << n, n), np.int64); sp = np.zeros((1 << n, n), np.int64); Hs = np.zeros(1 << n, np.int64)
sp2 = np.zeros((1024, 10), np.int64); Hs2 = np.zeros(1024, np.int64); Uout = np.zeros(10, np.int64)
st = np.zeros(4, np.int64)
bad = 0
for it in range(300):
    code = int(rng.integers(0, 1 << 45))
    A.adj_from_code(code, n, adj)
    A.tables(adj, n, hp, sp, Hs)
    D = census_one(adj, n, sp2, Hs2, LOG, Uout)
    if D != Hs[(1 << n) - 1]:
        bad += 1
    for y in range(n):
        A.root_stats(n, hp, sp, Hs, y, st)
        if st[2] + st[3] != Uout[y]:
            bad += 1
print("census_one against tables/root_stats on 300 random tournaments of order 10: %d mismatches" % bad)
# literal check on one tournament of order 10 (10! permutations)
perms = A.perm_array(10)
E = np.zeros(n, np.int64); Ue = np.zeros(n, np.int64); Uo = np.zeros(n, np.int64)
Hl = A.literal_stats(adj, n, perms, E, Ue, Uo)
D = census_one(adj, n, sp2, Hs2, LOG, Uout)
print("literal (10! permutations) on the last one: H %d vs %d ; U equal at all roots: %s" %
      (Hl, D, bool(np.all(Ue + Uo == Uout))))
t0 = time.time()
for it in range(20000):
    census_one(adj, n, sp2, Hs2, LOG, Uout)
dt = (time.time() - t0) / 20000
print("time per tournament of order 10: %.1f microseconds" % (dt * 1e6))
