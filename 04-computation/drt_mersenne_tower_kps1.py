#!/usr/bin/env python3
"""
drt_mersenne_tower_kps1.py — kind-pasteur-2026-06-09-S1, Branch B (HYP-2333, THM-448)

The DRT/Mersenne tower: seed S_2 = [[1,1],[-1,1]] (= border of the 1-vertex tournament),
iterate S -> [[S, S], [S - 2I, 2I - S]] to orders 4, 8, 16, 32, 64.

Parts:
 1. Tower verification: skew-Hadamard at each order; core tournament on 2^k - 1
    vertices is a DRT.  T7 (order 8 core) vs Paley T_7 iso re-check.
 2. T15 deep dive: scores, c3, M^2 = J - 15I, eigenvalues, |Aut| (backtracking with
    bitmask pruning), orbits / vertex-transitivity, self-converse probe,
    rank of A mod 2/3/5, triple-intersection distribution.
 3. Circulant check: enumerate all 2^7 = 128 antisymmetric S-sets of Z_15, find every
    circulant DRT(15), compare to T15 by invariants then explicit iso search.
 4. H comparisons: H(T15) exact, transitive T15 (=1), 5 random 15-tournaments,
    plus any circulant DRT(15) found.
 5. Walsh structure: S16_tower vs Sylvester H_16 row correlations; the exact
    correction matrix C = S_tower - H_sylvester (same recursion, seed 0);
    nnz(C_N) = N log2 N claim; per-row Hamming/dot profile; N=32, 64 summaries.

Output tee'd to 05-knowledge/results/drt_mersenne_tower_kps1.out
"""
import sys, time, itertools
from math import comb, factorial
from collections import Counter
import numpy as np

sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import (H_count, M_of, A_of, scores, is_skew_hadamard,
    normalize_first_row, core_tournament, is_DRT, is_iso)

OUT = open('05-knowledge/results/drt_mersenne_tower_kps1.out', 'w', encoding='utf-8')
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
    """Multiset of |N+(u) cap N+(v) cap N+(t)| over all unordered triples."""
    n = A.shape[0]; cnt = Counter()
    for u in range(n):
        for v in range(u + 1, n):
            vec = A[u] * A[v]
            vals = A @ vec
            for t in range(v + 1, n):
                cnt[int(vals[t])] += 1
    return cnt

def iso_bt(A, B, count_all=True, budget=None, record_perms=False, time_limit=None):
    """Backtracking iso search A -> B with bitmask adjacency pruning.
    Returns dict: count, nodes, aborted, perms (if recorded), found."""
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

def orbits_from_perms(n, perms):
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
    groups = Counter(find(v) for v in range(n))
    return sorted(groups.values(), reverse=True)

def main():
    t_start = time.time()
    w('=== drt_mersenne_tower_kps1 — Branch B: the DRT/Mersenne tower (HYP-2333, THM-448) ===')
    w('')

    # ---------- PART 1: the tower ----------
    w('--- PART 1: bordered tower from S_2 = [[1,1],[-1,1]], doubling S -> [[S,S],[S-2I,2I-S]] ---')
    S = np.array([[1]], dtype=np.int64)
    cores = {}
    for k in range(1, 7):
        S = double_S(S)
        order = S.shape[0]
        sh = is_skew_hadamard(S)
        row0_allones = bool((S[0] == 1).all())
        Sn = normalize_first_row(S)
        normalize_noop = np.array_equal(S, Sn)
        T = core_tournament(Sn)
        m = order - 1
        drt = is_DRT(T)
        sc = sorted(set(scores(T)))
        cores[order] = T
        w(f'order {order:>3}: skew-Hadamard={sh}  row0-all-ones={row0_allones} '
          f'(normalize no-op={normalize_noop})  core T{m}: DRT={drt}  score set={sc}')
    w('')

    T3, T7, T15, T31, T63 = cores[4], cores[8], cores[16], cores[32], cores[64]

    # T7 vs Paley re-check
    P7 = np.zeros((7, 7), dtype=np.int64)
    for i in range(7):
        for j in range(7):
            if (j - i) % 7 in (1, 2, 4):
                P7[i, j] = 1
    w(f'T7 (order-8 core) iso to Paley T_7: {is_iso(T7, P7)}')
    w(f'T3 (order-4 core) = C_3: scores={scores(T3)}  c3={c3_of(T3)}')
    w('')

    # ---------- PART 2: T15 deep dive ----------
    w('--- PART 2: T15 deep dive ---')
    sc15 = scores(T15)
    w(f'scores(T15) = {sc15}  (7-regular: {set(sc15) == {7}})')
    c3 = c3_of(T15)
    w(f'c3(T15) = {c3}  [note: ANY 7-regular 15-tournament has c3 = C(15,3) - 15*C(7,2) = {comb(15,3) - 15*comb(7,2)};')
    w(f'        max possible for n=15 (n=3 mod 4, DRT-attained) — NOT 105]')
    M15 = M_of(T15)
    J = np.ones((15, 15), dtype=np.int64)
    w(f'M^2 == J - 15I exactly: {np.array_equal(M15 @ M15, J - 15 * np.eye(15, dtype=np.int64))}')
    ev = np.linalg.eigvals(M15.astype(float))
    ev_sorted = sorted(ev, key=lambda z: (round(z.imag, 6), round(z.real, 6)))
    evc = Counter((round(z.real, 6), round(z.imag, 6)) for z in ev_sorted)
    w(f'eigenvalues of M (multiplicity): {dict(evc)}')
    w(f'sqrt(15) = {np.sqrt(15):.6f}')

    t0 = time.time()
    aut = iso_bt(T15, T15, count_all=True, record_perms=True)
    orb = orbits_from_perms(15, aut['perms'])
    w(f'|Aut(T15)| = {aut["count"]}  (search nodes={aut["nodes"]}, {time.time()-t0:.1f}s)')
    w(f'orbit sizes = {orb}  vertex-transitive = {orb == [15]}')

    t0 = time.time()
    sconv = iso_bt(T15, T15.T.copy(), count_all=False, budget=20_000_000)
    w(f'T15 self-converse (iso to T15^op): {sconv["found"]}  '
      f'(nodes={sconv["nodes"]}, aborted={sconv["aborted"]}, {time.time()-t0:.1f}s)')

    r2, r3, r5 = rank_mod(T15, 2), rank_mod(T15, 3), rank_mod(T15, 5)
    w(f'rank(A) mod 2,3,5 = {r2}, {r3}, {r5}')
    td15 = triple_dist(T15)
    w(f'triple |N+ cap N+ cap N+| distribution (455 triples): {dict(sorted(td15.items()))}')
    w('')

    # ---------- PART 3: circulant check on Z_15 ----------
    w('--- PART 3: circulant DRT(15) enumeration (all 128 antisymmetric S-sets of Z_15) ---')
    pairs = [(d, 15 - d) for d in range(1, 8)]
    drt_sets = []
    for bits in range(128):
        Sset = set()
        for kk, (a, b) in enumerate(pairs):
            Sset.add(a if (bits >> kk) & 1 else b)
        Ac = np.zeros((15, 15), dtype=np.int64)
        for i in range(15):
            for j in range(15):
                if i != j and (j - i) % 15 in Sset:
                    Ac[i, j] = 1
        if is_DRT(Ac):
            drt_sets.append((sorted(Sset), Ac))
    w(f'circulant DRT(15) count (labeled S-sets): {len(drt_sets)}')
    circ_reps = []   # iso-class representatives among circulant DRTs
    for Sset, Ac in drt_sets:
        placed = False
        for rep in circ_reps:
            r = iso_bt(Ac, rep[1], count_all=False, budget=10_000_000)
            if r['found']:
                rep[0].append(Sset); placed = True; break
        if not placed:
            circ_reps.append([[Sset], Ac])
    w(f'circulant DRT(15) iso classes: {len(circ_reps)}')
    for irep, rep in enumerate(circ_reps):
        Ac = rep[1]
        td = triple_dist(Ac)
        w(f'  class {irep}: {len(rep[0])} S-sets, e.g. S={rep[0][0]}')
        w(f'    triple dist: {dict(sorted(td.items()))}  '
          f'rank mod 2,3,5 = {rank_mod(Ac,2)},{rank_mod(Ac,3)},{rank_mod(Ac,5)}')
        autc = iso_bt(Ac, Ac, count_all=True)
        w(f'    |Aut| = {autc["count"]}')
        r = iso_bt(T15, Ac, count_all=False, budget=20_000_000)
        w(f'    T15 iso to this circulant: {r["found"]}  (nodes={r["nodes"]}, aborted={r["aborted"]})')
        rep.append(r['found'])
    if not drt_sets:
        w('  NO circulant DRT exists on Z_15 — T15 is automatically non-circulant.')
    w('')

    # ---------- PART 4: H comparisons ----------
    w('--- PART 4: H comparisons at n=15 ---')
    t0 = time.time()
    H15 = H_count(T15)
    w(f'H(T15) = {H15}   ({time.time()-t0:.1f}s)')
    # transitive
    Tr = np.triu(np.ones((15, 15), dtype=np.int64), 1)
    t0 = time.time()
    Htr = H_count(Tr)
    w(f'H(transitive T15) = {Htr}   ({time.time()-t0:.1f}s)')
    rng = np.random.default_rng(42)
    rand_hs = []
    for trial in range(5):
        A = np.zeros((15, 15), dtype=np.int64)
        for i in range(15):
            for j in range(i + 1, 15):
                if rng.integers(2):
                    A[i, j] = 1
                else:
                    A[j, i] = 1
        t0 = time.time()
        h = H_count(A)
        rand_hs.append(h)
        w(f'H(random 15-tournament #{trial}) = {h}   ({time.time()-t0:.1f}s)')
    w(f'random sample: min={min(rand_hs)} max={max(rand_hs)} mean={sum(rand_hs)/5:.1f}')
    # circulant DRT H values (per iso class)
    for irep, rep in enumerate(circ_reps):
        t0 = time.time()
        hc = H_count(rep[1])
        w(f'H(circulant DRT(15) class {irep}) = {hc}   ({time.time()-t0:.1f}s)'
          + ('   [iso to T15]' if rep[2] else '   [NOT iso to T15]'))
    w(f'scale: 15!/2^14 = {factorial(15) // 2**14}  '
      f'(Alon max-H order of magnitude Theta(n!/2^(n-1)))')
    w(f'H(T15)/ (15!/2^14) = {H15 / (factorial(15) / 2**14):.4f}')
    w('')

    # ---------- PART 5: Walsh structure ----------
    w('--- PART 5: Walsh indexing — S_tower vs Sylvester H ---')
    # rebuild towers of both kinds, track correction C = S - W (seed C_1 = 0)
    Stow = np.array([[1]], dtype=np.int64)
    Wsyl = np.array([[1]], dtype=np.int64)
    for k in range(1, 7):
        n = Stow.shape[0]
        I = np.eye(n, dtype=np.int64)
        Stow = np.block([[Stow, Stow], [Stow - 2 * I, 2 * I - Stow]])
        Wsyl = np.block([[Wsyl, Wsyl], [Wsyl, -Wsyl]])
        N = Stow.shape[0]
        C = Stow - Wsyl
        nnz = int((C != 0).sum())
        vals = sorted(set(C.flatten().tolist()))
        w(f'order {N:>3}: nnz(C = S_tower - H_sylvester) = {nnz}  '
          f'(claim N*log2(N) = {N * k}: {nnz == N * k}; claim N^2/2 = {N * N // 2}: {nnz == N * N // 2})'
          f'  C values in {vals}')
        if N == 16:
            S16, H16 = Stow.copy(), Wsyl.copy()
        if N == 32:
            S32, H32 = Stow.copy(), Wsyl.copy()
        if N == 64:
            S64, H64 = Stow.copy(), Wsyl.copy()
    # verify Wsyl == kron-Sylvester
    H2 = np.array([[1, 1], [1, -1]], dtype=np.int64)
    Hk = np.array([[1]], dtype=np.int64)
    for _ in range(4):
        Hk = np.kron(H2, Hk)
    w(f'W tower (order 16) == Sylvester H2^kron4: {np.array_equal(H16, Hk)}')
    w('')
    w('16x16 row-correlation profile: D = S16 @ H16 (H16 symmetric; D[i,j] = <S16 row i, Walsh char j>)')
    D = S16 @ H16
    w(f'{"row":>3} {"diag dot":>9} {"max|dot|":>9} {"argmax":>7} {"argmax==row":>12} {"Hamming to best":>16}')
    diag_eq_max = 0
    for i in range(16):
        row = D[i]
        mx = int(np.max(np.abs(row)))
        am = int(np.argmax(np.abs(row)))
        dd = int(row[i])
        ham = (16 - mx) // 2
        diag_eq_max += (am == i and abs(dd) == mx)
        w(f'{i:>3} {dd:>9} {mx:>9} {am:>7} {str(am == i):>12} {ham:>16}')
    w(f'rows where best Walsh char is its OWN index: {diag_eq_max}/16')
    absvals = Counter(int(abs(x)) for x in D.flatten())
    w(f'|dot| value histogram over all 256 pairs: {dict(sorted(absvals.items()))}')
    w('')
    for N, SN, HN in ((32, S32, H32), (64, S64, H64)):
        DN = SN @ HN
        mx_per_row = [int(np.max(np.abs(DN[i]))) for i in range(N)]
        diag_per_row = [int(DN[i, i]) for i in range(N)]
        w(f'order {N}: per-row max|dot| histogram {dict(sorted(Counter(mx_per_row).items()))}; '
          f'diag dot histogram {dict(sorted(Counter(diag_per_row).items()))}; '
          f'argmax==row for {sum(1 for i in range(N) if int(np.argmax(np.abs(DN[i]))) == i)}/{N} rows')
    w('')
    # twin pairing = binary tree coordinates (since row 0 stays all-ones, no scrambling)
    w('Twin structure: row 0 of S_tower is all-ones at every order (verified above as')
    w('normalize no-op), so core T_{2^k-1} vertex v (1..2^k-1) carries the BINARY index')
    w('of the tower: twins at the last doubling are (v, v+2^(k-1)); vertex 2^(k-1) is the')
    w("border's twin.")
    w('VERDICT on "Walsh plus low-weight correction": REFUTED in the entrywise sense —')
    w('nnz(C_N) = N^2/2 exactly (HALF of all entries differ; the naive N log2 N recursion')
    w('argument fails because the -2I/+2I corrections cancel/seed on diagonals).')
    w('PARTIAL in the spectral sense: every row of S_tower has |<row, some Walsh char>| >= N/2;')
    w('per-row max|dot| histogram is the exact ladder {N/2: 4 rows, N/2+4: 8, ..., N-4: 8, N: 4};')
    w('diagonal dot spectrum = full arithmetic progression -N..N step 4, each value x2 (+-N x1).')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
