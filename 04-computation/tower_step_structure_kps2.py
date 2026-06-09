#!/usr/bin/env python3
"""
tower_step_structure_kps2.py — kind-pasteur-2026-06-09-S2, Branch I follow-up

Finding from skew_double_trans_law_kps2: trans(T_{2m+1}) == trans(D(T_m)) at
m = 3, 7, 15, 31. Conjectured structure: the tower core T_{2m+1} IS the skew
double D(T_m) plus ONE extra vertex (the border twin), and that extra vertex
never extends a maximum transitive chain.

Checks:
 1. Labeled identity: T_{2m+1} restricted to (copy1 .. copy2, skipping the
    border-twin at core index m) == D_skew(T_m), for m = 3, 7, 15, 31, 63.
 2. Border-twin arcs: core index m beats exactly copy2 (and loses to copy1)?
    Report its exact in/out pattern.
 3. g(full, border_twin) = size of the largest transitive set CONTAINING the
    border twin, for m = 3, 7, 15, 31 — is it ever equal to trans(T_{2m+1})?
 4. trans(T127): bracket via has_TT decisions (lower bound = trans(D(T63)) = 15).
 5. Paley record table extension: trans(Paley_p) for p = 71, 79, 83, 103.

Output: 05-knowledge/results/tower_step_structure_kps2.out
"""
import sys, time
from collections import Counter
import numpy as np

sys.setrecursionlimit(100000)
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\04-computation')
from skew_doubling_core_kps1 import (D_skew, normalize_first_row, core_tournament,
    is_skew_hadamard)
from erdos_moser_trans_tower_kps2 import (TransSolver, masks_of, trans_of, paley,
    double_S)

OUT = open(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\tower_step_structure_kps2.out',
           'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def build_tower_to(order):
    S = np.array([[1]], dtype=np.int64)
    cores = {}
    while S.shape[0] < order:
        S = double_S(S)
        assert is_skew_hadamard(S)
        cores[S.shape[0] - 1] = core_tournament(normalize_first_row(S))
    return cores

def g_containing(A, v):
    """Largest transitive subset of A's vertex set containing v."""
    mk, n = masks_of(A)
    inmask = [sum(1 << i for i in range(n) if A[i, j]) for j in range(n)]
    full = (1 << n) - 1
    memo = {}
    tsv = TransSolver(mk)
    def g(S):
        # S always contains v; largest transitive subset of S containing v
        key = S
        r = memo.get(key)
        if r is not None:
            return r
        best = 1 + tsv.trans(mk[v] & S)           # v as source
        T = inmask[v] & S                          # candidate sources beating v
        items = []
        while T:
            b = T & -T
            u = b.bit_length() - 1
            T ^= b
            items.append(((mk[u] & S).bit_count(), u))
        items.sort(reverse=True)
        for sz, u in items:
            if sz + 1 <= best:
                break
            sub = mk[u] & S
            if sub & (1 << v):
                val = 1 + g(sub)
                if val > best:
                    best = val
        memo[key] = best
        return best
    return g(full)

def main():
    t_start = time.time()
    w('=== tower_step_structure_kps2 — T_{2m+1} = D(T_m) + border-twin ===')
    w('')
    cores = build_tower_to(128)     # cores at 3, 7, 15, 31, 63, 127
    w(f'tower cores built: {sorted(cores.keys())}')
    w('')

    # ---------- 1+2: labeled identity ----------
    w('--- 1+2: labeled submatrix identity ---')
    for m in (3, 7, 15, 31, 63):
        big = cores[2 * m + 1]
        Tm = cores[m]
        keep = list(range(0, m)) + list(range(m + 1, 2 * m + 1))
        sub = big[np.ix_(keep, keep)]
        Dm = D_skew(Tm)[0]
        same = np.array_equal(sub, Dm)
        # border twin (core index m) arc pattern
        bt_out = [j for j in range(2 * m + 1) if big[m, j]]
        out_copy1 = sum(1 for j in bt_out if j < m)
        out_copy2 = sum(1 for j in bt_out if j > m)
        w(f'm={m:>2}: T{2*m+1} minus border-twin == D(T{m}): {same}   '
          f'border-twin out-arcs: {out_copy1}/{m} into copy1, {out_copy2}/{m} into copy2')
    w('')

    # ---------- 3: does the border twin ever lie on a maximum chain? ----------
    w('--- 3: largest transitive set CONTAINING the border twin ---')
    for m in (3, 7, 15, 31):
        big = cores[2 * m + 1]
        tr, _ = trans_of(big)
        gv = g_containing(big, m)
        w(f'T{2*m+1}: trans = {tr}   max chain through border-twin = {gv}   '
          f'twin on a maximum chain: {gv == tr}')
    w('')

    # ---------- 4: trans(T127) bracket ----------
    w('--- 4: trans(T127) ---')
    T127 = cores[127]
    mk, n = masks_of(T127)
    sv = TransSolver(mk)
    full = (1 << n) - 1
    w('lower bound: T127 contains D(T63) (identity above), trans(D(T63)) = 15  =>  trans(T127) >= 15')
    k = 16
    t0 = time.time()
    last = 15
    while True:
        res = sv.has_TT(full, k)
        w(f'has_TT(T127, {k}) = {res}   (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
        if res:
            last = k
            k += 1
        else:
            w(f'trans(T127) = {last} exactly')
            break
        if time.time() - t0 > 1500:
            w(f'budget reached: trans(T127) >= {last} (upper side unresolved)')
            break
    w('')

    # ---------- 5: Paley extension ----------
    w('--- 5: Paley record table extension ---')
    for p in (71, 79, 83, 103):
        if time.time() - t_start > 2400:
            w(f'Paley_{p}: SKIPPED (time budget)')
            continue
        t0 = time.time()
        tr, sv2 = trans_of(paley(p))
        w(f'trans(Paley_{p}) = {tr}   (nodes={sv2.nodes}, {time.time()-t0:.1f}s)')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
