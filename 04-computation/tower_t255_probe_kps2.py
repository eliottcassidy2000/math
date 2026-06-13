#!/usr/bin/env python3
"""
tower_t255_probe_kps2.py — kind-pasteur-2026-06-09-S2, Branch I stretch test

Two-step law observed on the Mersenne tower: trans(level k) = 2*trans(level k-2) + 1
  (2, 3, 5, 7, 11, 15 at n = 3, 7, 15, 31, 63, 127 — matches A027383 pattern:
   odd k: 2^{(k+1)/2} - 1, even k: 3*2^{k/2-1} - 1).
Prediction: trans(T255) = 2*trans(T63) + 1 = 23.

Bracket via has_TT decisions, budget-limited.

Output: 05-knowledge/results/tower_t255_probe_kps2.out
"""
import sys, time
import numpy as np

sys.setrecursionlimit(100000)
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\04-computation')
from skew_doubling_core_kps1 import normalize_first_row, core_tournament, is_skew_hadamard
from erdos_moser_trans_tower_kps2 import TransSolver, masks_of, double_S

OUT = open(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\tower_t255_probe_kps2.out',
           'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def main():
    t_start = time.time()
    w('=== tower_t255_probe_kps2 — trans(T255) bracket (prediction: 23) ===')
    S = np.array([[1]], dtype=np.int64)
    while S.shape[0] < 256:
        S = double_S(S)
    assert is_skew_hadamard(S)
    T255 = core_tournament(normalize_first_row(S))
    w(f'T255 built: {T255.shape[0]} vertices, regular score check: '
      f'{sorted(set(int(x) for x in T255.sum(axis=1)))}')
    mk, n = masks_of(T255)
    sv = TransSolver(mk)
    full = (1 << n) - 1
    last = None
    for k in (23, 24):
        t0 = time.time()
        try:
            res = sv.has_TT(full, k)
        except MemoryError:
            w(f'has_TT(T255, {k}): MemoryError after {time.time()-t0:.1f}s '
              f'(nodes={sv.nodes})')
            break
        w(f'has_TT(T255, {k}) = {res}   (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
        if res:
            last = k
        if time.time() - t_start > 2400:
            w('budget reached')
            break
    if last == 23:
        # check 24 result if available to declare exact
        pass
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
