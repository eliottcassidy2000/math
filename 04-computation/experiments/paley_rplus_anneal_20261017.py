#!/usr/bin/env python3
"""On rotation survivors (R+, THM-4666 1c), does any non-rotation order beat the natural Paley order? (S40)

For each prime: the Farey sweep gives the best rotation order (the natural order when p is in R+); it is polished
(insertion + window DP) and a batch of simulated-annealing runs from random starts is added, using the S40 audit's
search routines (s40_audit_20261017/audit_search.py: polish, sa_batch). Any order beating the natural one is a
certificate (excess rechecked directly).
Usage: python paley_rplus_anneal_20261017.py RESTARTS STEPS_PER_VERTEX p1 p2 ...
Session opus-2026-10-17-S40.
"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 's40_audit_20261017'))
from paley_class_number_ranking_20261017 import kronecker_table, class_number
from paley_rotation_farey_sweep_20261017 import sweep, rotation_order, excess_of
import audit_search as AS


def main():
    R = int(sys.argv[1]); steps = int(sys.argv[2]); primes = [int(a) for a in sys.argv[3:]]
    rows = []
    for p in primes:
        t0 = time.time()
        chi8 = kronecker_table(-p)
        chi = np.array(chi8, dtype=np.int64)
        h = class_number(-p)
        nat = p * h
        best, ba, bb, bc, bd, Efin, total = sweep(chi8, p, h)
        rot = rotation_order(p, ba + bc, bb + bd) if best > nat else np.arange(p, dtype=np.int64)
        achi = AS.leg_table(p)
        pol = AS.polish(rot.copy(), achi, p, 14)
        res, orders = AS.sa_batch(achi, p, R, steps * p, 3.0, 0.15, 14, 1000 * p)
        b = int(np.argmax(res))
        cands = [(int(best), rot), (int(pol), None), (int(res[b]), orders[b].copy())]
        top = max(c[0] for c in cands)
        if res[b] == top:
            e_check = int(excess_of(np.array(orders[b], dtype=np.int64), chi8, p))
            assert e_check == top
        row = dict(p=p, h=h, natural=nat, best_rotation=int(best), polished_rotation=int(pol), sa_best=int(res[b]),
                   sa_median=int(np.median(res)), restarts=R, steps=steps * p, best=top, beaten=bool(top > nat),
                   seconds=round(time.time() - t0, 1))
        rows.append(row)
        print(row, flush=True)
    json.dump(rows, open(os.path.join(HERE, '..', '..', '05-knowledge', 'results', 'paley_rplus_anneal_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
