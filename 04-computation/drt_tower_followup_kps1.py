#!/usr/bin/env python3
"""
drt_tower_followup_kps1.py — kind-pasteur-2026-06-09-S1, Branch B follow-up.

Leads from drt_mersenne_tower_kps1 + drt_tower31_paley_kps1:
  L1. H(B_0(T31)) = 198335025 = H(T15) exactly  ->  is link(v)=T[N+(v)] of T31 iso to T15?
      (tower self-similarity: vertex links reproduce the previous tower level)
  L2. |Aut(T7)| = |Aut(T15)| = |Aut(T31)| = 21  ->  is the group F_21 (Frobenius) at every
      level?  Element-order profiles: F_21 = {1:1, 3:14, 7:6}; Z_21 = {1:1, 3:2, 7:6, 21:12}.
      Orbit partitions explicit (binary-index pattern of fixed points 0,1,3,7,...?).
  L3. Link-H spectrum: H(B_v) for ALL v of T15 and T31 — splits vertices finer than
      (scores, c3, cycle-counts) invariants; align with Aut orbits.
  L4. Is B_8(T31) (link-H 197147697 != H(T15)) a SECOND DRT(15)?  Is Paley31's link a DRT?
  L5. T63: per-vertex link triple-dist coloring, |Aut(T63)| (color-pruned backtracking,
      budgeted), orbit partition, element orders.  Predict |Aut| = 21, fixed points 7.
  L6. Self-converse probes for T31, T63 (T15 was NOT self-converse).
  L7. B_0(T63) iso T31?  (next level of self-similarity)

Output tee'd to 05-knowledge/results/drt_tower_followup_kps1.out
"""
import sys, time
from math import comb, gcd
from collections import Counter
import numpy as np

sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import (H_count, M_of, scores, is_skew_hadamard,
    normalize_first_row, core_tournament, is_DRT)

OUT = open('05-knowledge/results/drt_tower_followup_kps1.out', 'w', encoding='utf-8')
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

def triple_dist(A):
    n = A.shape[0]; cnt = Counter()
    for u in range(n):
        for v in range(u + 1, n):
            vec = A[u] * A[v]
            vals = A @ vec
            for t in range(v + 1, n):
                cnt[int(vals[t])] += 1
    return cnt

def out_sub(A, v):
    idx = np.flatnonzero(A[v])
    return A[np.ix_(idx, idx)]

def in_sub(A, v):
    idx = np.flatnonzero(A[:, v])
    return A[np.ix_(idx, idx)]

def iso_bt(A, B, colA=None, colB=None, count_all=True, budget=None,
           record_perms=False, time_limit=None):
    """Backtracking iso A->B; optional vertex colors (must match); rare colors first."""
    n = A.shape[0]
    if colA is None:
        colA = [0] * n
    if colB is None:
        colB = [0] * n
    freq = Counter(colA)
    order = sorted(range(n), key=lambda v: (freq[colA[v]], colA[v], v))
    Ar = A[np.ix_(order, order)]
    colAr = [colA[v] for v in order]
    Arows = [[int(Ar[i, j]) for j in range(n)] for i in range(n)]
    Brow = [sum(int(B[i, j]) << j for j in range(n)) for i in range(n)]
    res = {'count': 0, 'nodes': 0, 'aborted': False, 'perms': [], 'found': False}
    perm = [0] * n
    t0 = time.time()
    def bt(k, usedmask):
        if k == n:
            res['count'] += 1; res['found'] = True
            if record_perms:
                # translate back: order[k] -> perm[k]
                full = [0] * n
                for kk in range(n):
                    full[order[kk]] = perm[kk]
                res['perms'].append(full)
            return not count_all
        Ak = Arows[k]
        ck = colAr[k]
        need = 0
        for j in range(k):
            if Ak[j]:
                need |= 1 << perm[j]
        for img in range(n):
            b = 1 << img
            if usedmask & b:
                continue
            if colB[img] != ck:
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

def orbit_partition(n, perms):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for p in perms:
        for v in range(n):
            a, b = find(v), find(p[v])
            if a != b:
                parent[a] = b
    groups = {}
    for v in range(n):
        groups.setdefault(find(v), []).append(v)
    return sorted(groups.values(), key=lambda g: (-len(g), g))

def perm_order(p):
    n = len(p); seen = [False] * n; o = 1
    for v in range(n):
        if seen[v]:
            continue
        l = 0; x = v
        while not seen[x]:
            seen[x] = True; x = p[x]; l += 1
        o = o * l // gcd(o, l)
    return o

def main():
    t_start = time.time()
    w('=== drt_tower_followup_kps1 — self-similarity, F_21, T63 ===')
    w('')
    S = np.array([[1]], dtype=np.int64)
    cores = {}
    for k in range(1, 7):
        S = double_S(S)
        cores[S.shape[0]] = core_tournament(normalize_first_row(S))
    T7, T15, T31, T63 = cores[8], cores[16], cores[32], cores[64]
    P31 = np.zeros((31, 31), dtype=np.int64)
    QR = sorted({(x * x) % 31 for x in range(1, 31)})
    for i in range(31):
        for j in range(31):
            if i != j and (j - i) % 31 in QR:
                P31[i, j] = 1

    # ---------- L2: group structure ----------
    w('--- L2: Aut group structure at each level ---')
    for name, T in (('T7', T7), ('T15', T15), ('T31', T31)):
        r = iso_bt(T, T, count_all=True, record_perms=True)
        orders = Counter(perm_order(p) for p in r['perms'])
        orb = orbit_partition(T.shape[0], r['perms'])
        w(f'{name}: |Aut| = {r["count"]}  element orders = {dict(sorted(orders.items()))}')
        w(f'   orbits: {[(len(g), g) for g in orb]}')
    w('   [F_21 profile = {1:1, 3:14, 7:6};  Z_21 profile = {1:1, 3:2, 7:6, 21:12}]')
    w('')

    # ---------- L3: link-H spectra ----------
    w('--- L3: link-H spectrum H(T[N+(v)]) per vertex ---')
    for name, T in (('T15', T15), ('T31', T31)):
        vals = {}
        for v in range(T.shape[0]):
            h = H_count(out_sub(T, v))
            vals.setdefault(h, []).append(v)
        w(f'{name}: link-H values:')
        for h, vs in sorted(vals.items(), reverse=True):
            w(f'   H = {h}: vertices {vs}  (count {len(vs)})')
    hP = {}
    for v in range(31):
        h = H_count(out_sub(P31, v))
        hP.setdefault(h, []).append(v)
    w(f'Paley31: link-H values: {[(h, len(vs)) for h, vs in sorted(hP.items())]} (VT => constant)')
    w('')

    # ---------- L1 + L4: links as tournaments ----------
    w('--- L1/L4: are T31 links copies of T15 / second DRT(15)? ---')
    for v in (0, 1, 8):
        B = out_sub(T31, v)
        drt = is_DRT(B)
        td = triple_dist(B)
        r = iso_bt(B, T15, count_all=False, budget=30_000_000)
        autB = iso_bt(B, B, count_all=True)
        w(f'B_{v}(T31): DRT={drt}  H={H_count(B)}  |Aut|={autB["count"]}  '
          f'triple dist={dict(sorted(td.items()))}')
        w(f'   iso to T15: {r["found"]} (nodes={r["nodes"]}, aborted={r["aborted"]})')
    # in-neighborhood of 0
    Bin = in_sub(T31, 0)
    r = iso_bt(Bin, T15, count_all=False, budget=30_000_000)
    w(f'B-_0(T31) (IN-link): DRT={is_DRT(Bin)}  H={H_count(Bin)}  iso to T15: {r["found"]}')
    # Paley link
    BP = out_sub(P31, 0)
    autBP = iso_bt(BP, BP, count_all=True)
    w(f'B_0(Paley31): DRT={is_DRT(BP)}  H={H_count(BP)}  |Aut|={autBP["count"]}  '
      f'triple dist={dict(sorted(triple_dist(BP).items()))}')
    rPL = iso_bt(BP, T15, count_all=False, budget=30_000_000)
    w(f'   iso to T15: {rPL["found"]}')
    # T15 links vs T7
    w('T15 links vs T7:')
    for v in (0, 8):
        B = out_sub(T15, v)
        r = iso_bt(B, T7, count_all=False)
        w(f'   B_{v}(T15): DRT={is_DRT(B)}  H={H_count(B)}  iso to T7(=Paley7): {r["found"]}')
    w('')

    # ---------- L6: self-converse ----------
    w('--- L6: self-converse probes ---')
    r = iso_bt(T31, np.ascontiguousarray(T31.T), count_all=False, budget=50_000_000, time_limit=240)
    w(f'T31 self-converse: found={r["found"]}  aborted={r["aborted"]}  nodes={r["nodes"]}')
    w('')

    # ---------- L5: T63 ----------
    w('--- L5: T63 — coloring, Aut, orbits ---')
    w(f'T63: c3 = {c3_of(T63)}  M^2==J-63I: '
      f'{np.array_equal(M_of(T63) @ M_of(T63), np.ones((63,63),dtype=np.int64) - 63*np.eye(63,dtype=np.int64))}')
    t0 = time.time()
    link_sig = []
    for v in range(63):
        B = out_sub(T63, v)
        td = tuple(sorted(triple_dist(B).items()))
        link_sig.append(td)
    sigmap = {}
    colT63 = []
    for s_ in link_sig:
        if s_ not in sigmap:
            sigmap[s_] = len(sigmap)
        colT63.append(sigmap[s_])
    cc = Counter(colT63)
    w(f'link triple-dist coloring: {len(sigmap)} classes, sizes {sorted(cc.values(), reverse=True)}'
      f'   ({time.time()-t0:.1f}s)')
    t0 = time.time()
    r63 = iso_bt(T63, T63, colA=colT63, colB=colT63, count_all=True,
                 record_perms=True, budget=200_000_000, time_limit=420)
    w(f'|Aut(T63)| = {r63["count"]}{" (LOWER BOUND — aborted)" if r63["aborted"] else ""}'
      f'  nodes={r63["nodes"]}  time={time.time()-t0:.1f}s')
    if not r63['aborted']:
        orders = Counter(perm_order(p) for p in r63['perms'])
        orb = orbit_partition(63, r63['perms'])
        w(f'element orders = {dict(sorted(orders.items()))}')
        w(f'orbit sizes = {[len(g) for g in orb]}')
        w(f'fixed points = {sorted(g[0] for g in orb if len(g) == 1)}')
        w(f'orbits: {[(len(g), g) for g in orb]}')
    # T63 self-converse with colors
    t0 = time.time()
    link_sig_op = []
    T63op = np.ascontiguousarray(T63.T)
    for v in range(63):
        B = out_sub(T63op, v)
        td = tuple(sorted(triple_dist(B).items()))
        link_sig_op.append(td)
    if Counter(link_sig) != Counter(link_sig_op):
        w('T63 self-converse: NO (link triple-dist multisets of T63 and T63^op differ)')
    else:
        colOp = [sigmap.get(s_, -1) for s_ in link_sig_op]
        rsc = iso_bt(T63, T63op, colA=colT63, colB=colOp, count_all=False,
                     budget=100_000_000, time_limit=240)
        w(f'T63 self-converse: found={rsc["found"]}  aborted={rsc["aborted"]}  nodes={rsc["nodes"]}')
    w('')

    # ---------- L7: B_0(T63) iso T31? ----------
    w('--- L7: next-level self-similarity ---')
    for v in (0, 1):
        B = out_sub(T63, v)
        w(f'B_{v}(T63): DRT={is_DRT(B)}  triple dist={dict(sorted(triple_dist(B).items()))}')
        r = iso_bt(B, T31, count_all=False, budget=50_000_000, time_limit=240)
        w(f'   iso to T31: {r["found"]}  (nodes={r["nodes"]}, aborted={r["aborted"]})')
    tdT31 = triple_dist(T31)
    w(f'(reference T31 triple dist = {dict(sorted(tdT31.items()))})')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
