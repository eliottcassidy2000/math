#!/usr/bin/env python3
"""
t511_direct_capped_kps2.py — kind-pasteur-2026-06-09-S2, Branch I final probe

Direct has_TT(T511, 31) with memo capped at 12M entries and node budget,
after the uncapped attempt died. Certified already: trans(T511) in [30, 46]
(witness lift gives 30; upper from trans(D(T255))+1 <= 2*23+1 = 47 minus 1).
Prediction: 31.

Output: 05-knowledge/results/t511_direct_capped_kps2.out
"""
import sys, time
import numpy as np

sys.setrecursionlimit(100000)
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\04-computation')
from skew_doubling_core_kps1 import normalize_first_row, core_tournament
from erdos_moser_trans_tower_kps2 import masks_of, double_S

OUT = open(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\t511_direct_capped_kps2.out',
           'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

CAP = 12_000_000
BUDGET_NODES = 220_000_000

class Capped:
    def __init__(self, masks):
        self.masks = masks
        self.dmemo = {}
        self.nodes = 0

    def has_TT(self, S, k):
        if k <= 0:
            return True
        if S.bit_count() < k:
            return False
        key = (S, k)
        m = self.dmemo.get(key)
        if m is not None:
            return m
        self.nodes += 1
        if self.nodes > BUDGET_NODES:
            raise TimeoutError
        items = []
        T = S
        mk = self.masks
        while T:
            b = T & -T
            u = b.bit_length() - 1
            T ^= b
            d = (mk[u] & S).bit_count()
            if d >= k - 1:
                items.append((d, u))
        items.sort(reverse=True)
        res = False
        for d, u in items:
            if self.has_TT(mk[u] & S, k - 1):
                res = True
                break
        if len(self.dmemo) < CAP:
            self.dmemo[key] = res
        return res

def main():
    t_start = time.time()
    w('=== t511_direct_capped_kps2 — has_TT(T511, 31) with capped memo ===')
    S = np.array([[1]], dtype=np.int64)
    while S.shape[0] < 512:
        S = double_S(S)
    T511 = core_tournament(normalize_first_row(S))
    mk, n = masks_of(T511)
    sv = Capped(mk)
    full = (1 << n) - 1
    for k in (31, 32):
        t0 = time.time()
        try:
            res = sv.has_TT(full, k)
            w(f'has_TT(T511, {k}) = {res}   (nodes={sv.nodes}, memo={len(sv.dmemo)}, '
              f'{time.time()-t0:.1f}s)')
            if not res:
                break
        except TimeoutError:
            w(f'has_TT(T511, {k}): NODE BUDGET EXCEEDED (nodes={sv.nodes}, '
              f'{time.time()-t0:.1f}s) — undecided')
            break
        except MemoryError:
            w(f'has_TT(T511, {k}): MemoryError (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
            break
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
