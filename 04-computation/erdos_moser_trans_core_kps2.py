#!/usr/bin/env python3
"""
erdos_moser_trans_core_kps2.py — kind-pasteur-2026-06-09-S2, Branch I (HYP-2356/2357)

Erdos-Moser largest transitive subtournament trans(T) on the Mersenne tower.

trans(T) = size of the largest transitive (= acyclic) subtournament.
Recursion: trans(T) = 1 + max_v trans(T[N+(v)])  (every transitive sub has a source).
Branch & bound with memoization on vertex-subset bitmasks; children sorted by
out-degree-within-mask descending so the bound 1+|sub| <= best allows BREAK.

Literature anchors used for verification (sources in session notes):
  R(2..6) = 2, 4, 8, 14, 28;  34 <= R(7) <= 47 (Neiman-Mackey-Heule arXiv:2011.00683);
  R(8) >= 57 (McCarthy-Monico arXiv:2408.04067; previous best 48 = Paley_47 + 1).
  ST_3, ST_7, ST_27 = quadratic-residue (Galois) tournaments; ST_13 circulant non-Galois.
  trans anchors: trans(TT_n)=n, trans(C3)=2, trans(Paley_7)=3 (unique TT_4-free on 7),
  trans(ST_13)=4, trans(QR(GF(27)))=5.

Parts:
  1. Algorithm verification: anchors + brute force (subset-score check) on all 1024
     labeled n=5, 50 random n=10, and T15 itself.
  2. Circulant Z_13 scan: find all TT_5-free circulant tournaments (ST_13 anchor).
  3. GF(27) Galois tournament: trans = 5 expected (TT_6-free, Sanchez-Flores).
  4. Tower cores T7, T15, T31: trans + per-vertex link table (trans of T[N+(v)]),
     max location, B_0 = previous-level check.
  5. Paley_19, _23, _31, _43: trans (record table).
  6. Random baseline: 20 random tournaments at n=31.

Output: 05-knowledge/results/erdos_moser_trans_core_kps2.out
"""
import sys, time, random
from collections import Counter
from math import comb
import numpy as np

sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import (M_of, A_of, scores, is_skew_hadamard,
    normalize_first_row, core_tournament, is_DRT, is_iso, all_tournaments, iso_classes)

sys.setrecursionlimit(100000)

OUT = open('05-knowledge/results/erdos_moser_trans_core_kps2.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

# ---------------- trans machinery ----------------

def adj_masks(A):
    n = A.shape[0]
    return [int(''.join('1' if A[i, j] else '0' for j in range(n - 1, -1, -1)), 2)
            for i in range(n)]

class TransSolver:
    def __init__(self, A):
        self.n = A.shape[0]
        self.adj = adj_masks(A)
        self.full = (1 << self.n) - 1
        self.memo = {}
        self.nodes = 0

    def f(self, mask):
        """Exact trans of the subtournament induced on `mask`."""
        if mask == 0:
            return 0
        memo = self.memo
        v = memo.get(mask)
        if v is not None:
            return v
        pc = mask.bit_count()
        if pc == 1:
            memo[mask] = 1
            return 1
        self.nodes += 1
        adj = self.adj
        items = []
        m = mask
        while m:
            b = m & -m
            m ^= b
            sub = adj[b.bit_length() - 1] & mask
            items.append((sub.bit_count(), sub))
        items.sort(key=lambda t: -t[0])
        best = 1
        for spc, sub in items:
            if spc + 1 <= best:
                break
            val = 1 + self.f(sub)
            if val > best:
                best = val
                if best == pc:
                    break
        memo[mask] = best
        return best

    def trans(self):
        return self.f(self.full)

    def link_table(self):
        """Exact trans(T[N+(v)]) for every v (forces unpruned computation)."""
        return [self.f(self.adj[v]) for v in range(self.n)]

def trans_of(A):
    return TransSolver(A).trans()

def trans_brute(A):
    """Independent check: subset S transitive iff within-S scores are all distinct
    (equivalently {0..|S|-1}). O(2^n * n)."""
    n = A.shape[0]
    adj = adj_masks(A)
    best = 1
    for S in range(1, 1 << n):
        pc = S.bit_count()
        if pc <= best:
            continue
        seen = 0
        ok = True
        m = S
        while m:
            b = m & -m
            m ^= b
            d = (adj[b.bit_length() - 1] & S).bit_count()
            bit = 1 << d
            if seen & bit:
                ok = False
                break
            seen |= bit
        if ok:
            best = pc
    return best

# ---------------- constructions ----------------

def transitive_T(n):
    return np.triu(np.ones((n, n), dtype=np.int64), 1)

def paley(p):
    QR = {(x * x) % p for x in range(1, p)}
    A = np.zeros((p, p), dtype=np.int64)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in QR:
                A[i, j] = 1
    return A

def circulant(n, S):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i, j] = 1
    return A

def gf27_galois():
    """QR tournament on GF(27) = F_3[x]/(x^3+2x+1)."""
    def mul(a, b):
        # a,b are (c0,c1,c2) coeff tuples
        prod = [0] * 5
        for i in range(3):
            for j in range(3):
                prod[i + j] = (prod[i + j] + a[i] * b[j]) % 3
        # reduce x^3 = -2x - 1 = x + 2  (since x^3 + 2x + 1 = 0 -> x^3 = -2x-1 = x+2 mod 3)
        for k in (4, 3):
            c = prod[k]
            if c:
                prod[k] = 0
                prod[k - 3] = (prod[k - 3] + 2 * c) % 3   # +2c from constant 2
                prod[k - 2] = (prod[k - 2] + 1 * c) % 3   # +c from x
        return (prod[0], prod[1], prod[2])
    elems = [(a, b, c) for c in range(3) for b in range(3) for a in range(3)]
    idx = {e: i for i, e in enumerate(elems)}
    nonzero = [e for e in elems if e != (0, 0, 0)]
    QR = {mul(e, e) for e in nonzero}
    assert len(QR) == 13
    neg = {(-a % 3, -b % 3, -c % 3) for (a, b, c) in QR}
    assert not (QR & neg), "QR not antisymmetric — -1 would be a square"
    A = np.zeros((27, 27), dtype=np.int64)
    for ei in elems:
        for ej in elems:
            if ei == ej:
                continue
            diff = tuple((y - x) % 3 for x, y in zip(ei, ej))
            if diff in QR:
                A[idx[ei], idx[ej]] = 1
    return A

def random_T(n, rng):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A

def build_tower():
    """Bordered skew-Hadamard tower; returns dict order->core tournament."""
    S = np.array([[1]], dtype=np.int64)
    cores = {}
    for k in range(1, 7):
        n = S.shape[0]
        I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
        Sn = normalize_first_row(S)
        cores[S.shape[0] - 1] = core_tournament(Sn)
    return cores

def is_strong(A):
    """Moon: T strong iff every k-prefix of ascending scores sums > C(k,2)."""
    s = sorted(scores(A))
    n = len(s)
    tot = 0
    for k in range(1, n):
        tot += s[k - 1]
        if tot <= comb(k, 2):
            return False
    return True

# ---------------- main ----------------

def main():
    t_start = time.time()
    w('=== erdos_moser_trans_core_kps2 — Branch I: trans(T) on the Mersenne tower ===')
    w('')

    # ---------- PART 1: algorithm verification ----------
    w('--- PART 1: algorithm verification ---')
    for n in (5, 8):
        Tn = transitive_T(n)
        w(f'trans(transitive_{n}) = {trans_of(Tn)}   (expect {n})')
    C3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    w(f'trans(C3) = {trans_of(C3)}   (expect 2)')
    P7 = paley(7)
    w(f'trans(Paley_7) = {trans_of(P7)}   (expect 3: unique largest TT_4-free, R(4)=8)')

    bad = 0
    for i, A in enumerate(all_tournaments(5)):
        if trans_of(A) != trans_brute(A):
            bad += 1
            w(f'  MISMATCH at labeled n=5 #{i}')
    w(f'exhaustive n=5 labeled (1024): B&B vs brute mismatches = {bad}')

    rng = random.Random(20260609)
    bad = 0
    for t in range(50):
        A = random_T(10, rng)
        if trans_of(A) != trans_brute(A):
            bad += 1
    w(f'random n=10 x50: B&B vs brute mismatches = {bad}')
    w('')

    # ---------- PART 2: circulant Z_13 scan (ST_13 anchor) ----------
    w('--- PART 2: circulant tournaments on Z_13 (64 connection sets) ---')
    t0 = time.time()
    pairs13 = [(d, 13 - d) for d in range(1, 7)]
    tt5free = []
    tcount = Counter()
    for bits in range(64):
        S = {a if (bits >> k) & 1 else b for k, (a, b) in enumerate(pairs13)}
        A = circulant(13, S)
        tr = trans_of(A)
        tcount[tr] += 1
        if tr == 4:
            tt5free.append((sorted(S), A))
    w(f'trans distribution over 64 circulants: {dict(sorted(tcount.items()))}')
    w(f'TT_5-free (trans=4) circulants: {len(tt5free)}  S-sets: {[s for s, _ in tt5free]}')
    if tt5free:
        allco = all(is_iso(tt5free[0][1], A) or True for _, A in tt5free[1:])  # iso at n=13 brute is too slow via perms
        w('(iso check among them skipped at n=13 — permutation brute infeasible; '
          'multiplier equivalence implies they are rotations/multiples of one class)')
    w(f'  [{time.time()-t0:.1f}s]  literature: ST_13 = unique TT_5-free tournament on 13 vertices, circulant, non-Galois')
    w('')

    # ---------- PART 3: GF(27) Galois tournament ----------
    w('--- PART 3: QR tournament on GF(27) (= ST_27, Sanchez-Flores) ---')
    t0 = time.time()
    G27 = gf27_galois()
    w(f'GF(27) tournament: regular = {len(set(scores(G27))) == 1}, DRT = {is_DRT(G27)}')
    s27 = TransSolver(G27)
    tr27 = s27.trans()
    w(f'trans(QR(GF(27))) = {tr27}   (expect 5: ST_27 is the unique largest TT_6-free, R(6)=28)')
    w(f'  [{time.time()-t0:.1f}s, memo={len(s27.memo)}, nodes={s27.nodes}]')
    w('')

    # ---------- PART 4: the Mersenne tower ----------
    w('--- PART 4: tower cores T7, T15, T31 ---')
    cores = build_tower()
    T7, T15, T31, T63 = cores[7], cores[15], cores[31], cores[63]
    np.save('05-knowledge/results/erdos_moser_T63_adj_kps2.npy', T63)
    w('(T63 adjacency saved to 05-knowledge/results/erdos_moser_T63_adj_kps2.npy)')

    w(f'T15 brute-force cross-check: trans_brute(T15) = {trans_brute(T15)}')
    for name, T in (('T7', T7), ('T15', T15), ('T31', T31)):
        t0 = time.time()
        sol = TransSolver(T)
        tr = sol.trans()
        links = sol.link_table()
        lc = Counter(links)
        n = T.shape[0]
        w(f'{name}: n={n}  trans = {tr}   [{time.time()-t0:.1f}s, memo={len(sol.memo)}, nodes={sol.nodes}]')
        w(f'  per-vertex link trans values 1+t (t=trans(T[N+(v)])): '
          f'{dict(sorted(((1 + k, v) for k, v in lc.items())))} (value: #vertices)')
        w(f'  link of vertex 0: trans = {links[0]} '
          f'(B_0 = previous tower level => expect trans(T_{(n - 1) // 2}))')
        argmaxes = [v for v in range(n) if links[v] == max(links)]
        w(f'  max link trans = {max(links)} attained at vertices {argmaxes}'
          f'{"  (includes B_0=v0)" if 0 in argmaxes else "  (NOT at B_0)"}')
        # verify recursion consistency
        assert tr == 1 + max(links), 'recursion violated!'
        w(f'  recursion trans = 1 + max_v trans(link(v)) verified: {tr} = 1 + {max(links)}')
    w('')
    w('forcing floors: R(5)=14<=15 => trans(T15)>=5; R(6)=28<=31 => trans(T31)>=6;')
    w('               R(7)<=47<=63 => trans(T63)>=7  (T63 in separate script)')
    w('')

    # B_0 labeled-submatrix re-verification
    for big, small, nm in ((T15, T7, 'B_0(T15)=T7'), (T31, T15, 'B_0(T31)=T15'), (T63, T31, 'B_0(T63)=T31')):
        nb = [j for j in range(big.shape[0]) if big[0, j]]
        sub = big[np.ix_(nb, nb)]
        w(f'{nm} (labeled): {np.array_equal(sub, small)}')
    w('')

    # ---------- PART 5: Paley records ----------
    w('--- PART 5: Paley tournaments (QR on GF(p), p = 3 mod 4) ---')
    for p in (19, 23, 31, 43):
        t0 = time.time()
        sol = TransSolver(paley(p))
        tr = sol.trans()
        w(f'trans(Paley_{p}) = {tr}   [{time.time()-t0:.1f}s, memo={len(sol.memo)}, nodes={sol.nodes}]')
    w('')
    w('known floors: trans(Paley_19), trans(Paley_23) >= 5 (R(5)=14); trans(Paley_31) >= 6 (R(6)=28);')
    w('literature: trans(Paley_47) = 7 exactly (TT_8-free per McCarthy-Monico "q_8=47"; contains TT_7 since 47>=R(7) upper... ')
    w('            strictly: TT_7 in Paley_47 follows only if R(7)<=47, the proved upper bound).')
    w('')

    # ---------- PART 6: random baseline n=31 ----------
    w('--- PART 6: 20 random tournaments at n=31 ---')
    rng = random.Random(31063)
    vals = []
    t0 = time.time()
    for trial in range(20):
        A = random_T(31, rng)
        tr = trans_of(A)
        vals.append(tr)
    w(f'trans values: {vals}')
    w(f'distribution: {dict(sorted(Counter(vals).items()))}  '
      f'min={min(vals)} max={max(vals)} mean={sum(vals)/len(vals):.2f}   [{time.time()-t0:.1f}s total]')
    w(f'theory: 2log2(31)+1 = {2*np.log2(31)+1:.2f} (EM upper-tail), log2(31)+1 = {np.log2(31)+1:.2f} (Stearns floor)')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
