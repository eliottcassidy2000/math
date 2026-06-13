#!/usr/bin/env python3
"""
drt_tower31_paley_kps1.py — kind-pasteur-2026-06-09-S1, Branch B step 4.

T31 = core of the order-32 skew tower (seed S_2, doubling S -> [[S,S],[S-2I,2I-S]]).
Question: is T31 isomorphic to Paley T_31 (QR circulant)?

Probes (cheapest decisive first):
  1. DRT verification for both; c3; M^2 = J - 31I.
  2. Per-vertex invariants: out-neighborhood subtournament B_v = T[N+(v)] —
     (sorted scores of B_v, c3(B_v), sorted per-vertex 3-cycle counts in B_v).
     Constancy across v (vertex-transitivity necessary condition) + multiset equality.
  3. Triple distribution |N+(u) cap N+(v) cap N+(w)| over all C(31,3)=4495 triples.
  4. Quadruple distribution over all C(31,4)=31465 quadruples.
  5. rank(A) mod 2, 3, 5.
  6. Circulant DRT(31) enumeration: all 2^15 = 32768 antisymmetric S-sets of Z_31;
     which are DRT difference sets? (expected: exactly QR and NQR = Paley + Paley^op).
  7. H of the out-neighborhood subtournament B_0 (15-vertex Ham-path DP) for both.
  8. Explicit iso backtracking T31 -> Paley31 with bitmask pruning (budgeted).
  9. |Aut(T31)| (budgeted backtracking); Paley31 has |Aut| = 31*15 = 465 (verify).

Output tee'd to 05-knowledge/results/drt_tower31_paley_kps1.out
"""
import sys, time
from math import comb
from collections import Counter
import numpy as np

sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import (H_count, M_of, scores, is_skew_hadamard,
    normalize_first_row, core_tournament, is_DRT)

OUT = open('05-knowledge/results/drt_tower31_paley_kps1.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def double_S(S):
    n = S.shape[0]
    I = np.eye(n, dtype=np.int64)
    return np.block([[S, S], [S - 2 * I, 2 * I - S]])

def c3_of(A):
    n = A.shape[0]
    s = A.sum(axis=1)
    return comb(n, 3) - sum(comb(int(x), 2) for x in s)

def rank_mod(Ain, p):
    Amat = (Ain.astype(np.int64) % p).copy()
    n = Amat.shape[0]; r = 0
    for c in range(n):
        piv = None
        for i in range(r, n):
            if Amat[i, c] % p:
                piv = i; break
        if piv is None:
            continue
        Amat[[r, piv]] = Amat[[piv, r]]
        inv = pow(int(Amat[r, c]), p - 2, p)
        Amat[r] = (Amat[r] * inv) % p
        for i in range(n):
            if i != r and Amat[i, c]:
                Amat[i] = (Amat[i] - Amat[i, c] * Amat[r]) % p
        r += 1
    return r

def triple_dist(A):
    n = A.shape[0]; cnt = Counter()
    for u in range(n):
        for v in range(u + 1, n):
            vec = A[u] * A[v]
            vals = A @ vec
            for t in range(v + 1, n):
                cnt[int(vals[t])] += 1
    return cnt

def quad_dist(A):
    n = A.shape[0]; cnt = Counter()
    for u in range(n):
        for v in range(u + 1, n):
            pre2 = A[u] * A[v]
            for x in range(v + 1, n):
                pre3 = pre2 * A[x]
                vals = A @ pre3
                for t in range(x + 1, n):
                    cnt[int(vals[t])] += 1
    return cnt

def out_sub(A, v):
    idx = np.flatnonzero(A[v])
    return A[np.ix_(idx, idx)]

def per_vertex_inv(A):
    """For each v: (sorted scores of B_v, c3(B_v), sorted per-vertex cyc-counts in B_v)."""
    invs = []
    for v in range(A.shape[0]):
        B = out_sub(A, v)
        m = B.shape[0]
        cycper = sorted(int(B[t] @ B @ B[:, t]) for t in range(m))
        invs.append((tuple(sorted(scores(B))), c3_of(B), tuple(cycper)))
    return invs

def iso_bt(A, B, count_all=True, budget=None, record_perms=False, time_limit=None):
    n = A.shape[0]
    Arows = [[int(A[i, j]) for j in range(n)] for i in range(n)]
    Brow = [sum(int(B[i, j]) << j for j in range(n)) for i in range(n)]
    res = {'count': 0, 'nodes': 0, 'aborted': False, 'perms': [], 'found': False}
    perm = [0] * n
    t0 = time.time()
    def bt(k, usedmask):
        if k == n:
            res['count'] += 1; res['found'] = True
            if record_perms:
                res['perms'].append(perm.copy())
            return not count_all
        Ak = Arows[k]
        need = 0
        for j in range(k):
            if Ak[j]:
                need |= 1 << perm[j]
        for img in range(n):
            b = 1 << img
            if usedmask & b:
                continue
            res['nodes'] += 1
            if budget is not None and res['nodes'] > budget:
                res['aborted'] = True; return True
            if time_limit is not None and res['nodes'] % 1048576 == 0 \
                    and time.time() - t0 > time_limit:
                res['aborted'] = True; return True
            if (Brow[img] & usedmask) == need:
                perm[k] = img
                if bt(k + 1, usedmask | b):
                    return True
        return False
    bt(0, 0)
    return res

def main():
    t_start = time.time()
    w('=== drt_tower31_paley_kps1 — is the tower T31 the Paley tournament? ===')
    w('')
    # build tower to order 32
    S = np.array([[1]], dtype=np.int64)
    for k in range(1, 6):
        S = double_S(S)
    assert S.shape[0] == 32
    w(f'order 32 skew-Hadamard: {is_skew_hadamard(S)}')
    T31 = core_tournament(normalize_first_row(S))
    w(f'T31: scores set = {sorted(set(scores(T31)))}  DRT = {is_DRT(T31)}  c3 = {c3_of(T31)}')
    M = M_of(T31)
    J = np.ones((31, 31), dtype=np.int64)
    w(f'T31: M^2 == J - 31I exactly: {np.array_equal(M @ M, J - 31 * np.eye(31, dtype=np.int64))}')

    # Paley 31
    QR = sorted({(x * x) % 31 for x in range(1, 31)})
    w(f'QR mod 31 ({len(QR)}): {QR}')
    P = np.zeros((31, 31), dtype=np.int64)
    for i in range(31):
        for j in range(31):
            if i != j and (j - i) % 31 in QR:
                P[i, j] = 1
    w(f'Paley31: scores set = {sorted(set(scores(P)))}  DRT = {is_DRT(P)}  c3 = {c3_of(P)}')
    w('')

    # per-vertex invariants
    w('--- per-vertex out-neighborhood invariants ---')
    invT = per_vertex_inv(T31)
    invP = per_vertex_inv(P)
    cT, cP = Counter(invT), Counter(invP)
    w(f'T31:   {len(cT)} distinct per-vertex invariant values; counts = {sorted(cT.values(), reverse=True)}')
    w(f'Paley: {len(cP)} distinct per-vertex invariant values; counts = {sorted(cP.values(), reverse=True)}')
    w(f'T31 per-vertex invariants CONSTANT across v: {len(cT) == 1}  (necessary for vertex-transitivity)')
    w(f'multisets equal (T31 vs Paley): {cT == cP}')
    if len(cT) <= 4:
        for val, mult in cT.most_common():
            w(f'  T31 inv x{mult}: scores={val[0][:6]}... c3(B_v)={val[1]} cycper={val[2][:6]}...')
    for val, mult in cP.most_common():
        w(f'  Paley inv x{mult}: c3(B_v)={val[1]}  scores(B_v)={val[0]}')
    if len(cT) <= 4:
        for val, mult in cT.most_common():
            w(f'  T31  inv x{mult}: c3(B_v)={val[1]}  scores(B_v)={val[0]}')
    w('')

    # triple + quad distributions
    w('--- triple / quadruple common-out-neighborhood distributions ---')
    t0 = time.time()
    tdT, tdP = triple_dist(T31), triple_dist(P)
    w(f'T31   triple dist: {dict(sorted(tdT.items()))}')
    w(f'Paley triple dist: {dict(sorted(tdP.items()))}')
    w(f'triple dists EQUAL: {tdT == tdP}   ({time.time()-t0:.1f}s)')
    t0 = time.time()
    qdT, qdP = quad_dist(T31), quad_dist(P)
    w(f'T31   quad dist: {dict(sorted(qdT.items()))}')
    w(f'Paley quad dist: {dict(sorted(qdP.items()))}')
    w(f'quad dists EQUAL: {qdT == qdP}   ({time.time()-t0:.1f}s)')
    w('')

    # ranks
    w('--- small-prime ranks of adjacency ---')
    for p in (2, 3, 5):
        w(f'rank mod {p}: T31 = {rank_mod(T31, p)}   Paley = {rank_mod(P, p)}')
    w('')

    # circulant DRT(31) enumeration
    w('--- circulant DRT(31) enumeration: 2^15 antisymmetric S-sets of Z_31 ---')
    t0 = time.time()
    pairs = [(d, 31 - d) for d in range(1, 16)]
    IDX = (np.arange(31)[None, :] - np.arange(31)[:, None]) % 31  # IDX[d,x] = (x-d)%31
    found_sets = []
    for bits in range(1 << 15):
        s = np.zeros(31, dtype=np.int64)
        for kk, (a, b) in enumerate(pairs):
            if (bits >> kk) & 1:
                s[a] = 1
            else:
                s[b] = 1
        SR = s[IDX]            # row d = s shifted by d
        corrs = SR @ s         # corr[d] = |N+(0) cap N+(d)| for circulant
        inner = corrs[1:]
        if inner.min() == inner.max():
            found_sets.append(frozenset(int(x) for x in np.flatnonzero(s)))
    w(f'DRT difference sets found: {len(found_sets)}   ({time.time()-t0:.1f}s)')
    QRset = frozenset(QR)
    NQRset = frozenset(range(1, 31)) - QRset
    for fs in found_sets:
        tag = 'QR (Paley)' if fs == QRset else ('NQR (Paley^op)' if fs == NQRset else 'OTHER (!)')
        w(f'  S = {sorted(fs)}  -> {tag}')
    w(f'=> every circulant DRT(31) is Paley or Paley^op: {all(fs in (QRset, NQRset) for fs in found_sets)}')
    w('')

    # H of out-neighborhood subtournament
    w('--- H of the out-neighborhood subtournament B_0 (15 vertices) ---')
    t0 = time.time()
    hT = H_count(out_sub(T31, 0))
    w(f'H(B_0(T31))   = {hT}   ({time.time()-t0:.1f}s)')
    t0 = time.time()
    hP = H_count(out_sub(P, 0))
    w(f'H(B_0(Paley)) = {hP}   ({time.time()-t0:.1f}s)')
    # a second vertex of T31 in a different invariant class if any, else vertex 1
    t0 = time.time()
    hT1 = H_count(out_sub(T31, 1))
    w(f'H(B_1(T31))   = {hT1}   ({time.time()-t0:.1f}s)')
    t0 = time.time()
    hT8 = H_count(out_sub(T31, 8))
    w(f'H(B_8(T31))   = {hT8}   ({time.time()-t0:.1f}s)')
    w('')

    # explicit iso attempt
    w('--- explicit iso backtracking T31 -> Paley31 (budget 80M nodes / 420s) ---')
    t0 = time.time()
    r = iso_bt(T31, P, count_all=False, budget=80_000_000, time_limit=420)
    verdict = ('ISOMORPHIC (explicit map found)' if r['found']
               else ('INCONCLUSIVE (budget exhausted)' if r['aborted']
                     else 'NON-ISOMORPHIC (search space exhausted)'))
    w(f'result: {verdict}  nodes={r["nodes"]}  time={time.time()-t0:.1f}s')
    w('')

    # Aut counts
    w('--- automorphism counts (budgeted) ---')
    t0 = time.time()
    rP = iso_bt(P, P, count_all=True, budget=80_000_000, time_limit=300)
    w(f'|Aut(Paley31)| = {rP["count"]}{" (LOWER BOUND, aborted)" if rP["aborted"] else ""} '
      f' nodes={rP["nodes"]}  time={time.time()-t0:.1f}s  [expected 465]')
    t0 = time.time()
    rT = iso_bt(T31, T31, count_all=True, budget=80_000_000, time_limit=300, record_perms=True)
    w(f'|Aut(T31)| = {rT["count"]}{" (LOWER BOUND, aborted)" if rT["aborted"] else ""} '
      f' nodes={rT["nodes"]}  time={time.time()-t0:.1f}s')
    if not rT['aborted']:
        # orbits
        parent = list(range(31))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for pp in rT['perms']:
            for v in range(31):
                a, b = find(v), find(pp[v])
                if a != b:
                    parent[a] = b
        orb = sorted(Counter(find(v) for v in range(31)).values(), reverse=True)
        w(f'T31 orbit sizes: {orb}  vertex-transitive = {orb == [31]}')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
