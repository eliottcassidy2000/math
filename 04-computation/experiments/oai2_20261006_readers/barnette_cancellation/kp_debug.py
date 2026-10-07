import time, sys
sys.argv = ['x', '6', '1']
import knight_paths2 as K
eng = K.Engine(6)
eng.seed(20)
path = [0, 6*1 + 2]
EP = [eng.eid[(min(path), max(path))]]
cand = [k for k in range(len(eng.E)) if k not in EP]
t = time.time(); pool, seen = eng.pool_for(EP); print("pool", len(pool), time.time() - t)
for lb in range(1, 8):
    t = time.time(); st, B = K.master(cand, pool, size=lb); print("master lb", lb, st, time.time() - t, B)
    if st == "OK":
        t = time.time(); c = K.find_hc(eng.E, eng.nV, set(B), EP); print(" oracle", type(c), time.time() - t)
        break
