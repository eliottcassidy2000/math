#!/usr/bin/env python3
"""
verify_trans_tower_kps2.py — independent verification of trans_tower_erdos_moser_kps2 results.

(1) SINK recursion trans(C) = 1 + max_v trans(C ∩ N-(v)) — different search space than the
    source recursion (different memo sets), must agree.
(2) Explicit witness extraction + validation for T63 (TT11) and T31 (TT7).
(3) BRUTE FORCE subset enumeration for all D(T) at n=3..5 (2^10 subsets each) — fully
    independent check of the 17/18 law and the idx10 exception.
Output: 05-knowledge/results/verify_trans_tower_kps2.out
"""
import sys, itertools
import numpy as np
sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import iso_classes, D_skew, scores, is_DRT, normalize_first_row, core_tournament
from trans_tower_erdos_moser_kps2 import tower_core, paley, trans_exact

def trans_sink(A):
    n = A.shape[0]
    inn = [frozenset(j for j in range(n) if A[j, i]) for i in range(n)]
    memo = {}
    def rec(C):
        if not C: return 0
        if C in memo: return memo[C]
        best = 0
        for v in sorted(C, key=lambda v: -len(inn[v] & C)):
            if 1 + len(inn[v] & C) <= best: continue
            s = rec(inn[v] & C)
            if 1 + s > best: best = 1 + s
        memo[C] = best
        return best
    return rec(frozenset(range(n)))

def witness(A):
    """Extract one maximum transitive set (source recursion with reconstruction)."""
    n = A.shape[0]
    out = [frozenset(j for j in range(n) if A[i, j]) for i in range(n)]
    memo = {}
    def rec(C):
        if not C: return (0, [])
        if C in memo: return memo[C]
        best = (0, [])
        for v in sorted(C, key=lambda v: -len(out[v] & C)):
            if 1 + len(out[v] & C) <= best[0]: continue
            s, chain = rec(out[v] & C)
            if 1 + s > best[0]: best = (1 + s, [v] + chain)
        memo[C] = best
        return best
    return rec(frozenset(range(n)))

def is_transitive_chain(A, chain):
    return all(A[chain[i], chain[j]] for i in range(len(chain)) for j in range(i + 1, len(chain)))

def brute_trans(A):
    n = A.shape[0]
    for k in range(n, 0, -1):
        for sub in itertools.combinations(range(n), k):
            # transitive iff acyclic iff some ordering works: tournaments — sort by score within sub
            S = A[np.ix_(sub, sub)]
            order = sorted(range(k), key=lambda v: -S[v].sum())
            if all(S[order[i], order[j]] for i in range(k) for j in range(i + 1, k)):
                return k
    return 0

def main():
    out = open('05-knowledge/results/verify_trans_tower_kps2.out', 'w', encoding='utf-8')
    def w(s=''):
        out.write(s + '\n'); out.flush(); print(s)

    w('=== verify_trans_tower_kps2 ===\n')
    # (1) sink vs source on tower + Paley
    for k in (3, 4, 5, 6):
        T = tower_core(k)
        a, b = trans_exact(T), trans_sink(T)
        w(f'T{2**k - 1}: source={a} sink={b} agree={a == b}')
    P31 = paley(31)
    w(f'Paley_31: source={trans_exact(P31)} sink={trans_sink(P31)}')
    w('')
    # (2) witnesses
    for k in (5, 6):
        T = tower_core(k)
        val, chain = witness(T)
        ok = is_transitive_chain(T, chain)
        w(f'T{2**k-1}: witness TT{val} chain={chain} valid={ok}')
    w('')
    # (3) brute force doubles n=3..5
    w('--- brute force trans(D(T)) n=3..5 ---')
    mism = 0
    for n in (3, 4, 5):
        for idx, A in enumerate(iso_classes(n)):
            t_b = brute_trans(A)
            td_b = brute_trans(D_skew(A)[0])
            t_f = trans_exact(A)
            td_f = trans_exact(D_skew(A)[0])
            tag = 'OK' if (t_b, td_b) == (t_f, td_f) else 'MISMATCH'
            if tag == 'MISMATCH': mism += 1
            w(f'n={n} idx={idx}: brute ({t_b},{td_b}) vs fast ({t_f},{td_f}) {tag}'
              + ('   <-- EXCEPTION CLASS' if (n, idx) == (5, 10) else ''))
    w(f'mismatches: {mism}')
    # exception class detail: exhibit the TT5 in D(idx10)
    A = iso_classes(5)[10]
    Dd = D_skew(A)[0]
    val, chain = witness(Dd)
    w(f'exception class n=5 idx10: D(T) TT{val} witness {chain} valid={is_transitive_chain(Dd, chain)}')
    w(f'  base T arcs (i beats j): {[ (i,j) for i in range(5) for j in range(5) if A[i,j] ]}')
    w('')
    w('=== done ===')
    out.close()

if __name__ == '__main__':
    main()
