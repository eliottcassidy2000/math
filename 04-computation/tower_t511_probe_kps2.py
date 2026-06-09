#!/usr/bin/env python3
"""
tower_t511_probe_kps2.py — kind-pasteur-2026-06-09-S2, Branch I stretch test 2

Two-step law trans(level k) = 2*trans(level k-2) + 1 now confirmed at
n = 3..255 (trans = 2, 3, 5, 7, 11, 15, 23 = A027383). Prediction:
trans(T511) = 2*15 + 1 = 31. Bracket via has_TT(511, 31) and has_TT(511, 32).

Output: 05-knowledge/results/tower_t511_probe_kps2.out
"""
import sys, time
import numpy as np

sys.setrecursionlimit(100000)
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\04-computation')
from skew_doubling_core_kps1 import normalize_first_row, core_tournament, is_skew_hadamard
from erdos_moser_trans_tower_kps2 import TransSolver, masks_of, double_S

OUT = open(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\tower_t511_probe_kps2.out',
           'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def main():
    t_start = time.time()
    w('=== tower_t511_probe_kps2 — trans(T511) bracket (prediction: 31) ===')
    S = np.array([[1]], dtype=np.int64)
    while S.shape[0] < 512:
        S = double_S(S)
    assert is_skew_hadamard(S)
    T511 = core_tournament(normalize_first_row(S))
    w(f'T511 built: {T511.shape[0]} vertices, score set: '
      f'{sorted(set(int(x) for x in T511.sum(axis=1)))}')
    mk, n = masks_of(T511)
    sv = TransSolver(mk)
    full = (1 << n) - 1
    for k in (31, 32):
        t0 = time.time()
        try:
            res = sv.has_TT(full, k)
            w(f'has_TT(T511, {k}) = {res}   (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
        except MemoryError:
            w(f'has_TT(T511, {k}): MemoryError (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
            break
        if time.time() - t_start > 3000:
            w('budget reached')
            break
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
