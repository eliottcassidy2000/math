import time, sys
import knight_paths2 as K
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF
eng = K.Engine(6)
eng.seed(20)
path = [0, 6*1 + 2]
EP = [eng.eid[(min(path), max(path))]]
cand = [k for k in range(len(eng.E)) if k not in EP]
pool, seen = eng.pool_for(EP); print("pool", len(pool))
def mhs(pool, cand):
    w = WCNF()
    for C in pool:
        w.append([k + 1 for k in C])
    for k in cand:
        w.append([-(k + 1)], weight=1)
    with RC2(w, solver='g4') as rc2:
        m = rc2.compute()
        return rc2.cost, sorted(v - 1 for v in m if v > 0)
t = time.time(); cost, B = mhs(pool, cand); print("RC2 cost", cost, time.time() - t, B)
