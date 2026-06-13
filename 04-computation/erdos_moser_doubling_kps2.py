#!/usr/bin/env python3
"""
erdos_moser_doubling_kps2.py — kind-pasteur-2026-06-09-S2, Branch I task 4 (HYP-2357)

Doubling laws for trans(T) under the skew-Sylvester double
D(T): M' = [[M, M+I], [M-I, -M]] (THM-447).

Facts to establish:
  (a) trans(D(T)) >= trans(T) + 1 always (sink-twin append: chain v1..vk plus twin vk').
  (b) trans(D(T)) <= 2 trans(T) (copy1 part transitive in T, copy2 part transitive in T^op).
  (c) When is it exactly trans+1? Test predictor: strong connectivity of T.
  (d) A-priori limit (from directed Ramsey floors): the law trans(D(T)) = trans(T)+1
      cannot hold universally: D(Paley_7) has 14 vertices and R(5)=14 forces
      trans(D(Paley_7)) >= 5 = trans(Paley_7) + 2.
  (e) Iterated doubling D^k(C3): N = 3*2^k vertices; if "+1 per doubling" held, the
      family would have trans = k+2 = floor(log2 N)+1 (Stearns-tight), but R(7)<=47
      forces a break at N=48 (trans >= 7 > 6). Find the exact break schedule.

Parts:
  1. Exhaustive iso classes n=3..5: trans(T), trans(D(T)), strong(T); delta distribution.
  2. Exhaustive labeled n=6 (32768): delta distribution cross-tabbed vs strong(T).
  3. Samples n=7 (300), n=8 (150): same.
  4. Transitive line: trans(D(TT_n)) for n=1..10.
  5. Iterated doubling: D^k(C3) k=1..5 (N=6..96), D^k(Paley_7) k=1..2 (N=14,28).
  6. Tower doubles: trans(D(T7)), trans(D(T15)); D(T31) = T63 minus border-twin check.

Output: 05-knowledge/results/erdos_moser_doubling_kps2.out
"""
import sys, time, random
from collections import Counter
from math import comb
import numpy as np

sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import (M_of, scores, normalize_first_row,
    core_tournament, is_iso, all_tournaments, iso_classes, D_skew)

sys.setrecursionlimit(100000)

OUT = open('05-knowledge/results/erdos_moser_doubling_kps2.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

# ---------------- trans machinery (with optional deadline) ----------------

def adj_masks(A):
    n = A.shape[0]
    return [int(''.join('1' if A[i, j] else '0' for j in range(n - 1, -1, -1)), 2)
            for i in range(n)]

class Timeout(Exception):
    pass

class TransSolver:
    def __init__(self, A, deadline=None):
        self.n = A.shape[0]
        self.adj = adj_masks(A)
        self.full = (1 << self.n) - 1
        self.memo = {}
        self.nodes = 0
        self.deadline = deadline

    def f(self, mask):
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
        if self.deadline is not None and self.nodes % 65536 == 0 \
                and time.time() > self.deadline:
            raise Timeout()
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

def trans_of(A):
    return TransSolver(A).trans()

def is_strong(A):
    s = sorted(scores(A))
    n = len(s)
    tot = 0
    for k in range(1, n):
        tot += s[k - 1]
        if tot <= comb(k, 2):
            return False
    return True

def D_of(A):
    return D_skew(A)[0]

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
    S = np.array([[1]], dtype=np.int64)
    cores = {}
    for k in range(1, 7):
        n = S.shape[0]
        I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
        cores[S.shape[0] - 1] = core_tournament(normalize_first_row(S))
    return cores

def main():
    t_start = time.time()
    w('=== erdos_moser_doubling_kps2 — trans(D(T)) vs trans(T) (HYP-2357) ===')
    w('')

    # ---------- PART 1: exhaustive iso classes n=3..5 ----------
    w('--- PART 1: iso classes n=3..5 (exhaustive) ---')
    w(f'{"n":>2} {"idx":>3} {"scores":>16} {"strong":>6} {"tr(T)":>6} {"tr(D)":>6} {"delta":>6}')
    viol_lower = []
    deltas_by_strong = {True: Counter(), False: Counter()}
    for n in (3, 4, 5):
        for idx, A in enumerate(iso_classes(n)):
            trT = trans_of(A)
            trD = trans_of(D_of(A))
            st = is_strong(A)
            d = trD - trT
            deltas_by_strong[st][d] += 1
            if d < 1:
                viol_lower.append((n, idx))
            w(f'{n:>2} {idx:>3} {str(scores(A)):>16} {str(st):>6} {trT:>6} {trD:>6} {d:>6}')
    w(f'violations of trans(D) >= trans(T)+1: {len(viol_lower)} {viol_lower}')
    w(f'delta distribution | strong T:     {dict(sorted(deltas_by_strong[True].items()))}')
    w(f'delta distribution | non-strong T: {dict(sorted(deltas_by_strong[False].items()))}')
    w('')

    # ---------- PART 2: exhaustive labeled n=6 ----------
    w('--- PART 2: exhaustive labeled n=6 (32768 tournaments) ---')
    t0 = time.time()
    cross = Counter()   # (strong, delta) -> count
    ex2 = {}
    for i, A in enumerate(all_tournaments(6)):
        trT = trans_of(A)
        trD = trans_of(D_of(A))
        st = is_strong(A)
        d = trD - trT
        cross[(st, d)] += 1
        key = (st, d)
        if key not in ex2:
            ex2[key] = (i, trT, trD)
    w(f'cross-tab (strong, delta) -> count: {dict(sorted(cross.items()))}')
    w(f'first examples per (strong, delta): {ex2}')
    w(f'  [{time.time()-t0:.1f}s]')
    w('')

    # ---------- PART 3: samples n=7, n=8 ----------
    rng = random.Random(977)
    for n, cnt in ((7, 300), (8, 150)):
        t0 = time.time()
        cross = Counter()
        for _ in range(cnt):
            A = random_T(n, rng)
            d = trans_of(D_of(A)) - trans_of(A)
            cross[(is_strong(A), d)] += 1
        w(f'n={n} random x{cnt}: (strong, delta) -> {dict(sorted(cross.items()))}   [{time.time()-t0:.1f}s]')
    w('')

    # ---------- PART 4: transitive line ----------
    w('--- PART 4: trans(D(TT_n)), n=1..10 ---')
    for n in range(1, 11):
        Tn = transitive_T(n)
        trD = trans_of(D_of(Tn))
        w(f'n={n:>2}: trans(TT_n)={n}  trans(D)={trD}  delta={trD-n}  '
          f'(2*trans bound = {2*n})')
    w('')

    # ---------- PART 5: iterated doubling ----------
    w('--- PART 5: iterated doubling D^k ---')
    w('floors: R(2..6)=2,4,8,14,28; 34<=R(7)<=47; R(8)>=57.')
    C3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    A = C3
    tr = trans_of(A)
    w(f'C3: N=3 trans={tr}')
    for k in range(1, 6):
        A = D_of(A)
        N = A.shape[0]
        if N > 100:
            w(f'D^{k}(C3): N={N} — skipped (size)')
            break
        t0 = time.time()
        sol = TransSolver(A, deadline=time.time() + 420)
        try:
            tr = sol.trans()
            stearns = int(np.floor(np.log2(N))) + 1
            w(f'D^{k}(C3): N={N:>3} trans={tr}  strong={is_strong(A)}  '
              f'Stearns floor={stearns}  [{time.time()-t0:.1f}s, memo={len(sol.memo)}]')
        except Timeout:
            w(f'D^{k}(C3): N={N:>3} TIMEOUT after {time.time()-t0:.0f}s '
              f'(nodes={sol.nodes}, memo={len(sol.memo)})')
            break
    P7 = paley(7)
    A = P7
    tr = trans_of(A)
    w(f'Paley_7: N=7 trans={tr}')
    for k in range(1, 3):
        A = D_of(A)
        N = A.shape[0]
        t0 = time.time()
        sol = TransSolver(A, deadline=time.time() + 420)
        try:
            tr = sol.trans()
            w(f'D^{k}(P7): N={N:>3} trans={tr}  strong={is_strong(A)}  [{time.time()-t0:.1f}s]')
        except Timeout:
            w(f'D^{k}(P7): N={N:>3} TIMEOUT after {time.time()-t0:.0f}s')
            break
    w('')

    # ---------- PART 6: tower doubles ----------
    w('--- PART 6: tower doubles ---')
    cores = build_tower()
    T7, T15, T31, T63 = cores[7], cores[15], cores[31], cores[63]
    for name, T in (('D(T7) [N=14]', T7), ('D(T15) [N=30]', T15)):
        t0 = time.time()
        DT = D_of(T)
        sol = TransSolver(DT, deadline=time.time() + 420)
        try:
            tr = sol.trans()
            w(f'trans({name}) = {tr}  (trans(T)={trans_of(T)})  [{time.time()-t0:.1f}s, memo={len(sol.memo)}]')
        except Timeout:
            w(f'{name}: TIMEOUT after {time.time()-t0:.0f}s')
    # D(T31) vs T63 minus border-twin (core vertex 31 = S row 32)
    DT31 = D_of(T31)
    keep = [i for i in range(63) if i != 31]
    sub = T63[np.ix_(keep, keep)]
    w(f'D(T31) == T63 minus core-vertex 31 (border twin), as labeled matrices: '
      f'{np.array_equal(DT31, sub)}')
    t0 = time.time()
    sol = TransSolver(DT31, deadline=time.time() + 480)
    try:
        tr = sol.trans()
        w(f'trans(D(T31)) [N=62] = {tr}   [{time.time()-t0:.1f}s, memo={len(sol.memo)}, nodes={sol.nodes}]')
    except Timeout:
        w(f'trans(D(T31)): TIMEOUT after {time.time()-t0:.0f}s (nodes={sol.nodes})')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
