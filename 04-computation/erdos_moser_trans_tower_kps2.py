#!/usr/bin/env python3
"""
erdos_moser_trans_tower_kps2.py — kind-pasteur-2026-06-09-S2, Branch I (HYP-2356/2357)

Erdos-Moser largest transitive subtournament trans(T) on the Mersenne tower
T7, T15, T31, T63 (skew-Hadamard doubling cores, THM-448), vs Paley tournaments,
random baselines, and the known directed-Ramsey anchors:
  R(2)=2, R(3)=4, R(4)=8, R(5)=14, R(6)=28 [Sanchez-Flores 1994],
  34 <= R(7) <= 47 [Neiman-Mackey-Heule, arXiv:2011.00683, G&C 2022].
ST3, ST7, ST27 are the quadratic-residue (Galois) tournaments; ST13 is NOT Galois.
QR23 is TT6-free (order 23); Galois tournaments on 31, 43, 47 all contain TT7 (NMH).

Algorithm: trans(T) = 1 + max_v trans(T[N+(v)]) (every transitive sub has a source),
branch&bound on bitmask sets, memoized, descending-out-degree order with the
break-prune (1 + outdeg <= best so far at this node => no remaining candidate wins).

Parts:
 1. Anchor verification (transitive_n, C3, QR7, QR23, ST27/GF27, exhaustive n=4 and n=7).
 2. Mersenne tower: trans(T7), trans(T15), trans(T31), trans(T63).
 3. Link recursion: per-vertex trans(N+(v)) for T31 and T63; is max at B_0?
 4. Paley controls: QR11, QR19, QR23, QR31, QR43 (+47, 59, 67 time permitting).
 5. Random baselines: 20 random tournaments at n=31 and n=63.
 6. Circulant scans: Z13 (find ST13 among 64 circulants), Z31 and Z33 TT7-free scans.

Output: 05-knowledge/results/erdos_moser_trans_tower_kps2.out
"""
import sys, time, itertools
from collections import Counter
import numpy as np

sys.setrecursionlimit(100000)
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\04-computation')
from skew_doubling_core_kps1 import (M_of, A_of, scores, is_skew_hadamard,
    normalize_first_row, core_tournament, is_DRT)

OUT = None   # opened lazily in main() so that importing this module never truncates results
def w(s=''):
    global OUT
    if OUT is None:
        OUT = open(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\erdos_moser_trans_tower_kps2.out',
                   'w', encoding='utf-8')
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

# ---------------- trans engine ----------------

def masks_of(A):
    n = A.shape[0]
    return [sum(1 << j for j in range(n) if A[i, j]) for i in range(n)], n

class TransSolver:
    """Exact trans + decision has_TT on one tournament (bitmask masks)."""
    def __init__(self, masks):
        self.masks = masks
        self.memo = {}
        self.dmemo = {}
        self.nodes = 0

    def trans(self, S):
        if S == 0:
            return 0
        m = self.memo.get(S)
        if m is not None:
            return m
        self.nodes += 1
        items = []
        T = S
        mk = self.masks
        while T:
            b = T & -T
            u = b.bit_length() - 1
            T ^= b
            items.append(((mk[u] & S).bit_count(), u))
        items.sort(reverse=True)
        best = 0
        for sz, u in items:
            if sz + 1 <= best:
                break
            val = 1 + self.trans(mk[u] & S)
            if val > best:
                best = val
        self.memo[S] = best
        return best

    def has_TT(self, S, k):
        """Decision: does S contain a transitive subtournament on >= k vertices?"""
        if k <= 0:
            return True
        if S.bit_count() < k:
            return False
        key = (S, k)
        m = self.dmemo.get(key)
        if m is not None:
            return m
        self.nodes += 1
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
        self.dmemo[key] = res
        return res

def trans_of(A):
    mk, n = masks_of(A)
    sv = TransSolver(mk)
    return sv.trans((1 << n) - 1), sv

# ---------------- constructions ----------------

def transitive_T(n):
    return np.triu(np.ones((n, n), dtype=np.int64), 1)

def paley(p):
    qr = {(x * x) % p for x in range(1, p)}
    A = np.zeros((p, p), dtype=np.int64)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i, j] = 1
    return A

def circulant(n, Sset):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in Sset:
                A[i, j] = 1
    return A

def random_T(n, rng):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.integers(2):
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A

# GF(27) = GF(3)[x]/(x^3 - x - 1) (irreducible: no roots in GF(3))
def gf27_elements():
    return [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]

def gf27_mul(u, v):
    coeffs = [0] * 5
    for i, a in enumerate(u):
        for j, b in enumerate(v):
            coeffs[i + j] = (coeffs[i + j] + a * b) % 3
    for d in (4, 3):
        c = coeffs[d]
        if c:
            coeffs[d] = 0
            coeffs[d - 3] = (coeffs[d - 3] + c) % 3   # x^d -> x^{d-3} * (x + 1)
            coeffs[d - 2] = (coeffs[d - 2] + c) % 3
    return tuple(coeffs[:3])

def gf27_sub(u, v):
    return tuple((a - b) % 3 for a, b in zip(u, v))

def galois_27():
    els = gf27_elements()
    nz = [e for e in els if e != (0, 0, 0)]
    squares = {gf27_mul(g, g) for g in nz}
    assert len(squares) == 13
    neg = lambda u: tuple((-a) % 3 for a in u)
    assert all((d in squares) != (neg(d) in squares) for d in nz), 'not antisymmetric'
    idx = {e: i for i, e in enumerate(els)}
    A = np.zeros((27, 27), dtype=np.int64)
    for u in els:
        for v in els:
            if u != v and gf27_sub(v, u) in squares:
                A[idx[u], idx[v]] = 1
    return A

def double_S(S):
    n = S.shape[0]
    I = np.eye(n, dtype=np.int64)
    return np.block([[S, S], [S - 2 * I, 2 * I - S]])

def build_tower():
    S = np.array([[1]], dtype=np.int64)
    cores = {}
    for k in range(1, 7):
        S = double_S(S)
        assert is_skew_hadamard(S)
        T = core_tournament(normalize_first_row(S))
        cores[S.shape[0] - 1] = T
    return cores   # keys 3, 7, 15, 31, 63

def all_tournaments_n(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        yield A

# ---------------- main ----------------

def main():
    t_start = time.time()
    w('=== erdos_moser_trans_tower_kps2 — Branch I: Erdos-Moser trans on the Mersenne tower ===')
    w('Anchors (literature): R(3)=4, R(4)=8, R(5)=14, R(6)=28; 34 <= R(7) <= 47 (NMH 2022).')
    w('Forcing floor f(n) >= floor(log2 n) + 1 (Stearns): n=7:3  n=15:4  n=31:5  n=63:6.')
    w('Since 63 >= 47 >= R(7): trans(T63) >= 7 guaranteed. Since 31 >= 28 = R(6): trans(T31) >= 6.')
    w('')

    # ---------- PART 1: anchors ----------
    w('--- PART 1: algorithm anchors ---')
    for n in range(1, 11):
        tr, _ = trans_of(transitive_T(n))
        assert tr == n, (n, tr)
    w('trans(transitive_n) == n for n=1..10: PASS')
    C3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    trC3, _ = trans_of(C3)
    w(f'trans(C3) = {trC3} (expect 2): {"PASS" if trC3 == 2 else "FAIL"}')
    QR7 = paley(7)
    tr7, _ = trans_of(QR7)
    w(f'trans(Paley_7) = {tr7} (expect 3, unique TT4-free of order 7): {"PASS" if tr7 == 3 else "FAIL"}')

    # exhaustive n=4: every tournament contains TT3 (R(3)=4)
    bad = sum(1 for A in all_tournaments_n(4) if trans_of(A)[0] < 3)
    w(f'exhaustive n=4: tournaments with trans < 3: {bad}/64 (expect 0, R(3)=4): '
      f'{"PASS" if bad == 0 else "FAIL"}')

    # exhaustive n=7: count labeled tournaments with trans=3 (no TT4);
    # expect 7!/|Aut(QR7)| = 5040/21 = 240 (uniqueness of ST7 = Paley_7)
    t0 = time.time()
    pairs7 = [(i, j) for i in range(7) for j in range(i + 1, 7)]
    cnt3 = 0
    histo = Counter()
    for bits in range(1 << 21):
        mk = [0] * 7
        for k, (i, j) in enumerate(pairs7):
            if (bits >> k) & 1:
                mk[i] |= 1 << j
            else:
                mk[j] |= 1 << i
        sv = TransSolver(mk)
        if not sv.has_TT(127, 4):
            cnt3 += 1
        if bits < 50000:  # exact trans histogram on a deterministic 50k prefix
            histo[sv.trans(127)] += 1
    w(f'exhaustive n=7 (2^21 labeled): no-TT4 count = {cnt3} '
      f'(expect 240 = 5040/21, uniqueness of ST7): {"PASS" if cnt3 == 240 else "FAIL"}  '
      f'({time.time()-t0:.1f}s)')
    w(f'  exact trans histogram on first 50000 labeled: {dict(sorted(histo.items()))}')

    QR23 = paley(23)
    tr23, _ = trans_of(QR23)
    w(f'trans(Paley_23) = {tr23} (expect 5: QR23 is TT6-free, NMH sec. 3.1): '
      f'{"PASS" if tr23 == 5 else "FAIL"}')
    ST27 = galois_27()
    assert is_DRT(ST27)
    tr27, _ = trans_of(ST27)
    w(f'trans(Galois_27) = {tr27} (expect 5: ST27 unique TT6-free of order 27): '
      f'{"PASS" if tr27 == 5 else "FAIL"}')
    w('')

    # ---------- PART 2: the tower ----------
    w('--- PART 2: Mersenne tower cores ---')
    cores = build_tower()
    tower_solvers = {}
    tower_trans = {}
    for m in (3, 7, 15, 31, 63):
        t0 = time.time()
        tr, sv = trans_of(cores[m])
        tower_solvers[m] = sv
        tower_trans[m] = tr
        floor_f = m.bit_length()  # floor(log2 m) + 1 == bit_length for m = 2^k - 1
        w(f'trans(T{m}) = {tr}   (floor(log2 n)+1 = {floor_f}; nodes={sv.nodes}, '
          f'memo={len(sv.memo)}, {time.time()-t0:.1f}s)')
    w('')
    w(f'tower trans sequence (n=3,7,15,31,63): '
      f'{[tower_trans[m] for m in (3, 7, 15, 31, 63)]}')
    w(f'per-tower-step growth: '
      f'{[tower_trans[b] - tower_trans[a] for a, b in ((3,7),(7,15),(15,31),(31,63))]}')
    w('')

    # ---------- PART 3: link recursion ----------
    w('--- PART 3: link recursion trans(T) = 1 + max_v trans(N+(v)) ---')
    for m, prev in ((31, 15), (63, 31)):
        A = cores[m]
        sv = tower_solvers[m]
        mk, n = masks_of(A)
        link_tr = []
        for v in range(n):
            link_tr.append(sv.trans(mk[v]))
        c = Counter(link_tr)
        w(f'T{m}: link trans values (per out-neighborhood): {dict(sorted(c.items()))}')
        w(f'T{m}: trans = 1 + max(link) = {1 + max(link_tr)} '
          f'(direct: {tower_trans[m]})  consistent={1 + max(link_tr) == tower_trans[m]}')
        # B_0 check: out-neighborhood of vertex 0 should be the previous tower level (labeled)
        nbr0 = [j for j in range(n) if A[0, j]]
        sub0 = A[np.ix_(nbr0, nbr0)]
        same = np.array_equal(sub0, cores[prev])
        w(f'T{m}: B_0 (N+(0), size {len(nbr0)}) == T{prev} as labeled submatrix: {same}')
        w(f'T{m}: trans(B_0) = {link_tr[0]}  vs trans(T{prev}) = {tower_trans[prev]}  '
          f'(max link trans = {max(link_tr)}; attained at B_0: {link_tr[0] == max(link_tr)})')
        argmaxes = [v for v in range(n) if link_tr[v] == max(link_tr)]
        w(f'T{m}: vertices attaining max link trans: {len(argmaxes)} of {n} '
          f'(first few: {argmaxes[:10]})')
        w('')

    # ---------- PART 4: Paley controls ----------
    w('--- PART 4: Paley (quadratic residue) tournament controls ---')
    paley_results = {}
    for p in (7, 11, 19, 23, 31, 43, 47, 59, 67):
        if time.time() - t_start > 2400 and p > 47:
            w(f'Paley_{p}: SKIPPED (time budget)')
            continue
        t0 = time.time()
        tr, sv = trans_of(paley(p))
        paley_results[p] = tr
        w(f'trans(Paley_{p}) = {tr}   (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
    w('')
    w('NMH claim check: Galois tournaments on 31, 43, 47 all contain TT7 '
      f'=> trans >= 7: Paley31={paley_results.get(31)}, Paley43={paley_results.get(43)}, '
      f'Paley47={paley_results.get(47)}')
    w('')

    # ---------- PART 5: random baselines ----------
    w('--- PART 5: random baselines ---')
    rng = np.random.default_rng(20260609)
    for n in (31, 63):
        vals = []
        t0 = time.time()
        for trial in range(20):
            tr, _ = trans_of(random_T(n, rng))
            vals.append(tr)
        w(f'n={n}: 20 random tournaments trans = {sorted(vals)}  '
          f'min={min(vals)} max={max(vals)} mean={sum(vals)/len(vals):.2f}  '
          f'({time.time()-t0:.1f}s)  [2*log2(n)+1 = {2*np.log2(n)+1:.2f}]')
    w('')

    # ---------- PART 6: circulant scans ----------
    w('--- PART 6: circulant scans ---')
    # Z13: find ST13 (unique TT5-free of order 13, NOT Galois). Among 64 circulants?
    t0 = time.time()
    found13 = []
    pairs13 = [(d, 13 - d) for d in range(1, 7)]
    for bits in range(64):
        Sset = {a if (bits >> k) & 1 else b for k, (a, b) in enumerate(pairs13)}
        A = circulant(13, Sset)
        tr, _ = trans_of(A)
        if tr == 4:
            found13.append(sorted(Sset))
    w(f'Z13 scan (64 circulants): TT5-free (trans=4) symbol sets: {len(found13)} '
      f'({time.time()-t0:.1f}s)')
    for s in found13:
        w(f'  S = {s}')
    if found13:
        w('  => ST13 IS circulant (by Sanchez-Flores uniqueness, all of these are ST13).')
        ST13 = circulant(13, set(found13[0]))
        np.save(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\st13_adj_kps2.npy', ST13)
        w('  saved ST13 adjacency to 05-knowledge/results/st13_adj_kps2.npy')
    else:
        w('  => no circulant on Z13 is TT5-free; ST13 is NOT circulant.')
    w('')

    # Z31: TT7-free circulants (trans = 6). Budget-limited.
    for n_c, budget in ((31, 900), (33, 900)):
        t0 = time.time()
        half = (n_c - 1) // 2
        pairs_c = [(d, n_c - d) for d in range(1, half + 1)]
        total = 1 << half
        free_syms = []
        done = 0
        full_c = (1 << n_c) - 1
        for bits in range(total):
            if time.time() - t0 > budget:
                break
            Sset = [a if (bits >> k) & 1 else b for k, (a, b) in enumerate(pairs_c)]
            mk = [0] * n_c
            for i in range(n_c):
                acc = 0
                for d in Sset:
                    acc |= 1 << ((i + d) % n_c)
                mk[i] = acc
            sv = TransSolver(mk)
            if not sv.has_TT(full_c, 7):
                free_syms.append(frozenset(Sset))
            done += 1
        # dedupe by multiplier equivalence
        units = [u for u in range(1, n_c) if np.gcd(u, n_c) == 1]
        canon_syms = set()
        for fs in free_syms:
            best = None
            for u in units:
                img = frozenset((u * x) % n_c for x in fs)
                key = tuple(sorted(img))
                if best is None or key < best:
                    best = key
            canon_syms.add(best)
        w(f'Z{n_c} scan: {done}/{total} symbol sets checked in {time.time()-t0:.1f}s; '
          f'TT7-free circulants found: {len(free_syms)} '
          f'({len(canon_syms)} up to multiplier equivalence)')
        for cs in sorted(canon_syms):
            A = circulant(n_c, set(cs))
            tr, _ = trans_of(A)
            w(f'  S = {list(cs)}  trans = {tr}')
        w('')

    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
