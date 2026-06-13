#!/usr/bin/env python3
"""
trans_tower_erdos_moser_kps2.py — kind-pasteur-2026-06-09-S2 (branch I, run by session owner
after two agent API failures)

Largest transitive subtournament (trans) of:
  - the skew-Sylvester Mersenne tower cores T7, T15, T31, T63 (THM-448)
  - Paley_31 control, random controls at n=15, 31, 63
  - D(T) vs T for all iso classes n=3..5  (HYP-2357: trans(D(T)) = trans(T)+1?)
  - the per-level tower growth (THM-455 part 3)

Algorithm: trans(C) = 1 + max_{v in C} trans(C ∩ N+(v)), memoized on frozensets, with
branch-and-bound pruning (1 + |C| <= best). For DRTs candidate sets halve each level
(|N+ ∩ N+| = (n-3)/4), so this is fast.

Anchors: trans(transitive_n) = n; trans(C3) = 2; trans(Paley T7) = 3 (unique largest
TT4-free tournament).

Output: 05-knowledge/results/trans_tower_erdos_moser_kps2.out
"""
import sys, itertools, random
import numpy as np
sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import (all_tournaments, iso_classes, M_of, A_of, scores,
                                     D_skew, border, normalize_first_row, core_tournament,
                                     is_DRT)

def tower_core(k):
    """Core tournament of the skew tower of order 2^k (on 2^k - 1 vertices)."""
    S = np.array([[1]], dtype=np.int64)
    for _ in range(k):
        n = S.shape[0]
        I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
    return core_tournament(normalize_first_row(S))

def paley(p, residues=None):
    if residues is None:
        residues = {(x * x) % p for x in range(1, p)}
    A = np.zeros((p, p), dtype=np.int64)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in residues:
                A[i, j] = 1
    return A

def trans_exact(A):
    """Largest transitive subtournament via memoized recursion with B&B pruning."""
    n = A.shape[0]
    out = [frozenset(j for j in range(n) if A[i, j]) for i in range(n)]
    memo = {}
    best_holder = [0]

    def rec(C, depth):
        if not C:
            return 0
        key = C
        if key in memo:
            return memo[key]
        best = 0
        # order candidates by descending out-degree within C (better pruning)
        cand = sorted(C, key=lambda v: -len(out[v] & C))
        for v in cand:
            if 1 + len(out[v] & C) <= best:
                continue  # cannot beat current best within this call
            sub = rec(out[v] & C, depth + 1)
            if 1 + sub > best:
                best = 1 + sub
        memo[key] = best
        return best

    return rec(frozenset(range(n)), 0)

def main():
    out = open('05-knowledge/results/trans_tower_erdos_moser_kps2.out', 'w', encoding='utf-8')
    def w(s=''):
        out.write(s + '\n'); out.flush(); print(s)

    w('=== trans_tower_erdos_moser_kps2 — Erdős–Moser on the Mersenne tower ===\n')

    # anchors
    tr_n = np.triu(np.ones((6, 6), dtype=np.int64), 1)
    C3 = np.array([[0,1,0],[0,0,1],[1,0,0]], dtype=np.int64)
    w(f'anchor trans(transitive_6) = {trans_exact(tr_n)}   (expect 6)')
    w(f'anchor trans(C3)           = {trans_exact(C3)}   (expect 2)')
    P7 = paley(7)
    w(f'anchor trans(Paley_7)      = {trans_exact(P7)}   (expect 3)')
    w('')

    # tower cores
    w('--- tower cores ---')
    vals = {}
    for k in (3, 4, 5, 6):
        T = tower_core(k)
        n = T.shape[0]
        t = trans_exact(T)
        vals[n] = t
        w(f'T{n}: trans = {t}   (DRT={is_DRT(T)})')
    w('')

    # Paley 31 control + random controls
    w('--- controls ---')
    P31 = paley(31)
    w(f'trans(Paley_31) = {trans_exact(P31)}')
    rng = random.Random(20260609)
    for n in (15, 31, 63):
        rvals = []
        for _ in range(10):
            A = np.zeros((n, n), dtype=np.int64)
            for i in range(n):
                for j in range(i + 1, n):
                    if rng.random() < 0.5:
                        A[i, j] = 1
                    else:
                        A[j, i] = 1
            rvals.append(trans_exact(A))
        w(f'random n={n}: trans values {sorted(rvals)}')
    w('')

    # HYP-2357: trans(D(T)) vs trans(T), exhaustive iso classes n=3..5
    w('--- trans(D(T)) vs trans(T), iso classes n=3..5 (HYP-2357) ---')
    law_plus1 = 0; total = 0; excess = []
    for n in (3, 4, 5):
        for idx, A in enumerate(iso_classes(n)):
            t = trans_exact(A)
            td = trans_exact(D_skew(A)[0])
            total += 1
            if td == t + 1:
                law_plus1 += 1
            else:
                excess.append((n, idx, t, td))
            w(f'n={n} idx={idx} scores={scores(A)}: trans(T)={t}  trans(D(T))={td}  '
              f'{"= t+1" if td == t+1 else "EXCEPTION"}')
    w(f'law trans(D(T)) = trans(T)+1 holds in {law_plus1}/{total} classes')
    if excess:
        w(f'EXCEPTIONS: {excess}')
    w('')

    # bordered doubling growth (the tower step n -> 2n+1)
    w('--- tower-step growth: trans(T_{2m-1}) vs trans(T_{m-1}) ---')
    seq = [(7, vals[7]), (15, vals[15]), (31, vals[31]), (63, vals[63])]
    for (n1, t1), (n2, t2) in zip(seq, seq[1:]):
        w(f'T{n1} -> T{n2}: trans {t1} -> {t2}   (step +{t2 - t1})')
    w('')
    w('=== done ===')
    out.close()

if __name__ == '__main__':
    main()
