#!/usr/bin/env python3
"""
paley127_baseline_kps2.py — kind-pasteur-2026-06-09-S2, Branch I coda

Head-to-head at order 127: trans(T127) = 15 (tower core, computed in
tower_step_structure_kps2). Here: trans(Paley_107), trans(Paley_127), and
5 random 127-vertex tournaments for baseline.

Output: 05-knowledge/results/paley127_baseline_kps2.out
"""
import sys, time
import numpy as np

sys.setrecursionlimit(100000)
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\04-computation')
from erdos_moser_trans_tower_kps2 import TransSolver, masks_of, trans_of, paley, random_T

OUT = open(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\paley127_baseline_kps2.out',
           'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def main():
    t_start = time.time()
    w('=== paley127_baseline_kps2 — order-127 head-to-head ===')
    w('reference: trans(T127 tower core) = 15 (tower_step_structure_kps2)')
    w('')
    for p in (107, 127):
        t0 = time.time()
        tr, sv = trans_of(paley(p))
        w(f'trans(Paley_{p}) = {tr}   (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
    w('')
    rng = np.random.default_rng(20260609)
    vals = []
    for trial in range(5):
        if time.time() - t_start > 2000:
            w(f'random trial {trial}: SKIPPED (budget)')
            continue
        t0 = time.time()
        tr, sv = trans_of(random_T(127, rng))
        vals.append(tr)
        w(f'trans(random_127 #{trial}) = {tr}   (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
    if vals:
        w(f'random n=127 summary: {sorted(vals)}  mean={sum(vals)/len(vals):.2f}  '
          f'[2*log2(127)+1 = {2*np.log2(127)+1:.2f}]')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
