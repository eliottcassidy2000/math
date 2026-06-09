#!/usr/bin/env python3
"""
skew_doubling_tiling_walsh_kps1.py -- kind-pasteur-2026-06-09-S1, BRANCH C
Tiling / Walsh decomposition of the skew-Sylvester double D(T)  (HYP-2335, "BLUE conjecture").

CONVENTIONS (match repo / tournament-tiling-explorer.html EXACTLY):
  * Vertex labels 1..N.  Base path: N -> N-1 -> ... -> 1 (arcs (k,k-1), k=N..2).
  * TILES = pairs (x,y), x - y >= 2, x = upper vertex, y = lower vertex.
  * Tile enumeration order: for y=1..N-2: for x=N down to y+2  (explorer order).
  * BIT CONVENTION (this script): bit = 1  iff  arc points x -> y (downward).
      NOTE: the explorer stores the complement (bits[i]==0 => arc x->y). Grid-symmetry
      is polarity-invariant (sigma compares bits pairwise), so BLUE verdicts agree.
  * Grid symmetry sigma: tile (x,y) -> (N-y+1, N-x+1); tiling BLUE iff invariant.

RELABELING (Task 1):
  D(T) on 2n vertices: copy-1 = u_1..u_n (indices 0..n-1), copy-2 = v_1..v_n (indices n..2n-1).
  Canonical Ham path: u_n,...,u_1, (twin), v_1,...,v_n.
  Relabel:  u_i -> label n+i ;  v_j -> label n+1-j.
  (path position t, 0-based, gets label 2n - t).  Copy-1 occupies labels n+1..2n,
  copy-2 occupies labels 1..n.

DERIVED BLOCK TRANSFORM (verified exhaustively below):
  Let x = tiling of T (frame labels 1..n), A = full adjacency of T (label space).
  Tiling t of D(T) in the 2n frame:
    copy-1 block:  t(n+i, n+j)     = x(i,j)                       [identity copy]
    copy-2 block:  t(X, Y)         = x(sigma_n(X,Y))              [grid transpose, NO negation]
    cross block:   t(n+i, n+1-j)   = 1 if i==j (twin tiles, the ANTI-DIAGONAL X+Y=2n+1)
                                     else A[i][j]  (T's FULL arc matrix, each pair twice,
                                                    complementary: A[i][j] and A[j][i])
  => sigma_{2n} maps copy-1 <-> copy-2 bits EQUAL (symmetric), but maps cross tile (i,j) to
     cross tile (j,i) with COMPLEMENTARY bits.  So t XOR sigma(t) = indicator of the
     n(n-1) off-twin cross tiles, a CONSTANT vector c_n independent of T.
  => D(T) is NEVER blue (n>=2); it lies in the fixed coset {t : t ^ sigma(t) = c_n}
     of the blue subspace.  dim(blue subspace at 2n) = n(n-1);  D-image is an AFFINE
     subspace of dim C(n-1,2) inside that coset; relative density 2^{-(n-1)(n+2)/2}.

Output -> 05-knowledge/results/skew_doubling_tiling_walsh_kps1.out
"""
import sys, time, itertools
import numpy as np

sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import (all_tournaments, iso_classes, H_count, M_of, A_of,
                                     scores, D_skew, is_skew_hadamard,
                                     normalize_first_row, core_tournament, is_DRT, ham_path)

OUT = open('05-knowledge/results/skew_doubling_tiling_walsh_kps1.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(str(s) + '\n'); OUT.flush(); print(s)

# ---------------- tiling machinery ----------------

def tiles_list(N):
    return [(x, y) for y in range(1, N - 1) for x in range(N, y + 1, -1)]

def sigma_map(N, tiles):
    idx = {t: i for i, t in enumerate(tiles)}
    return [idx[(N - y + 1, N - x + 1)] for (x, y) in tiles]

def tiling_to_A(bits, N, tiles):
    """bits in explorer order, OUR polarity (1 = arc x->y). Returns adjacency, idx = label-1."""
    A = np.zeros((N, N), dtype=np.int64)
    for k in range(N, 1, -1):
        A[k - 1, k - 2] = 1                      # base path k -> k-1
    for i, (x, y) in enumerate(tiles):
        if bits[i]:
            A[x - 1, y - 1] = 1
        else:
            A[y - 1, x - 1] = 1
    return A

def A_to_tiling(A, N, tiles):
    for k in range(N, 1, -1):
        assert A[k - 1, k - 2] == 1, "base path missing"
    return [int(A[x - 1, y - 1]) for (x, y) in tiles]

def is_gridsym(bits, sig):
    return all(bits[i] == bits[sig[i]] for i in range(len(bits)))

def relabel_D(Ad, n):
    """Relabel D(T) (indices: 0..n-1 = u_1..u_n, n..2n-1 = v_1..v_n) into the 2n tiling frame."""
    N = 2 * n
    newlab = [0] * N
    for idx in range(n):
        newlab[idx] = n + idx + 1            # u_{idx+1} -> n+idx+1
    for idx in range(n, N):
        newlab[idx] = 2 * n - idx            # v_{idx-n+1} -> n+1-(idx-n+1) = 2n-idx
    inv = [0] * N
    for i in range(N):
        inv[newlab[i] - 1] = i
    return Ad[np.ix_(inv, inv)]

def frame_by_path(A, p):
    """Relabel tournament A so Ham path p becomes base path N -> ... -> 1."""
    N = A.shape[0]
    inv = [0] * N            # inv[label-1] = old index; label of p[t] = N - t
    for t, v in enumerate(p):
        inv[N - t - 1] = v
    return A[np.ix_(inv, inv)]

# ---------------- S2: exhaustive block-transform + BLUE verdict ----------------

def analyze_one(Aframed, n, tiles_n, tiles_2n, sig2n, cross_offtwin_idx):
    """Aframed = adjacency of T containing base path. Returns (t_bits, blue?, formulas_ok, defect_ok, weight_ok)."""
    x = {(X, Y): int(Aframed[X - 1, Y - 1]) for (X, Y) in tiles_n}
    Ad, _ = D_skew(Aframed)
    ADl = relabel_D(Ad, n)
    t = A_to_tiling(ADl, 2 * n, tiles_2n)
    # block formulas
    ok = True
    for i, (X, Y) in enumerate(tiles_2n):
        if X > n and Y > n:
            pred = x[(X - n, Y - n)]
        elif X <= n and Y <= n:
            pred = x[(n + 1 - Y, n + 1 - X)]
        else:
            ii, jj = X - n, n + 1 - Y
            pred = 1 if ii == jj else int(Aframed[ii - 1, jj - 1])
        if t[i] != pred:
            ok = False
            break
    blue = is_gridsym(t, sig2n)
    defect = frozenset(i for i in range(len(t)) if t[i] != t[sig2n[i]])
    defect_ok = (defect == cross_offtwin_idx)
    wt_ok = (sum(t) == 2 * sum(x.values()) + (n - 1) + n * (n - 1) // 2)
    return t, blue, ok, defect_ok, wt_ok

w('=== BRANCH C: tiling / Walsh decomposition of D(T) (HYP-2335, BLUE conjecture) ===')
w('Bit convention: bit=1 iff arc x->y (explorer stores complement; grid-sym unaffected).')
w('Relabeling: u_i -> n+i, v_j -> n+1-j  (path position t gets label 2n-t).')
w('')

w('--- S2: exhaustive verification of block transform + BLUE verdict ---')
for n in (3, 4, 5, 6):
    tiles_n = tiles_list(n)
    m = len(tiles_n)
    tiles_2n = tiles_list(2 * n)
    sig2n = sigma_map(2 * n, tiles_2n)
    cross_offtwin = frozenset(i for i, (X, Y) in enumerate(tiles_2n)
                              if X > n and Y <= n and X + Y != 2 * n + 1)
    nb = nf = nd = nw = 0
    tot = 1 << m
    for mask in range(tot):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = tiling_to_A(bits, n, tiles_n)
        t, blue, fok, dok, wok = analyze_one(A, n, tiles_n, tiles_2n, sig2n, cross_offtwin)
        nb += blue; nf += fok; nd += dok; nw += wok
    w(f'n={n}: tilings={tot}  BLUE(D)={nb}/{tot}  block-formula exact={nf}/{tot}  '
      f'defect==cross-offtwin={nd}/{tot}  weight law 2|x|+(n-1)+C(n,2)={nw}/{tot}  '
      f'[|defect| = n(n-1) = {n*(n-1)} tiles = {n*(n-1)//2} sigma-pairs of {(len(tiles_2n)-(n-1))//2}]')
w('')

w('--- S2b: same over ALL labeled tournaments (insertion-path frame), n=3..5 ---')
for n in (3, 4, 5):
    tiles_n = tiles_list(n)
    tiles_2n = tiles_list(2 * n)
    sig2n = sigma_map(2 * n, tiles_2n)
    cross_offtwin = frozenset(i for i, (X, Y) in enumerate(tiles_2n)
                              if X > n and Y <= n and X + Y != 2 * n + 1)
    nb = nf = nd = tot = 0
    for A in all_tournaments(n):
        Af = frame_by_path(A, ham_path(A))
        t, blue, fok, dok, wok = analyze_one(Af, n, tiles_n, tiles_2n, sig2n, cross_offtwin)
        nb += blue; nf += fok; nd += dok; tot += 1
    w(f'n={n}: labeled={tot}  BLUE(D)={nb}/{tot}  formula={nf}/{tot}  defect-const={nd}/{tot}')
w('')

# ---------------- S3: what the sigma-image IS ----------------

w('--- S3: sigma-image identification ---')
for n in (3, 4, 5):
    tiles_n = tiles_list(n)
    m = len(tiles_n)
    tiles_2n = tiles_list(2 * n)
    sig2n = sigma_map(2 * n, tiles_2n)
    c_vec = [1 if (X > n and Y <= n and X + Y != 2 * n + 1) else 0 for (X, Y) in tiles_2n]
    ok_coset = ok_swapop = ok_neverD = 0
    tot = 1 << m
    for mask in range(tot):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = tiling_to_A(bits, n, tiles_n)
        Ad, Mp = D_skew(A)
        ADl = relabel_D(Ad, n)
        t = A_to_tiling(ADl, 2 * n, tiles_2n)
        s = [t[sig2n[i]] for i in range(len(t))]
        # (a) sigma(t) = t XOR c   (constant-coset law)
        ok_coset += all(s[i] == (t[i] ^ c_vec[i]) for i in range(len(t)))
        # (b) sigma(t) = tiling of swap-conjugated D(T)^op = A_of([[M,-M+I],[-M-I,-M]])
        M = M_of(A); I = np.eye(n, dtype=np.int64)
        Mpp = np.block([[M, -M + I], [-M - I, -M]])
        # check Mpp == P (-M') P with P = copy swap; dominance of D(T)^op is -M' (M' is skew):
        P = np.block([[np.zeros((n, n), dtype=np.int64), I], [I, np.zeros((n, n), dtype=np.int64)]])
        assert np.array_equal(Mpp, P @ (-Mp) @ P), 'swap-op identity failed'
        App = A_of(Mpp)
        ADl2 = relabel_D(App, n)   # same relabeling map
        s2 = A_to_tiling(ADl2, 2 * n, tiles_2n)
        ok_swapop += (s2 == s)
        # (c) sigma(t) is NEVER itself a D-image (cross block would need A[i][j]==A[j][i])
        ok_neverD += (s != t or n < 2)   # s==t would be required for s = t(D(x)) with same copy-1 block
    w(f'n={n}: sigma(t)=t^c (const coset) {ok_coset}/{tot};  '
      f'sigma(t) = tiling(swap-relabeled D(T)^op) {ok_swapop}/{tot};  '
      f'sigma(t) != t always (so sigma-image never a D-image) {ok_neverD}/{tot}')
w('sigma(t(D(T))) IS the canonical tiling of D(T)^op (copy-swap relabeled = path-reversal frame).')
w('')

# ---------------- S4: Walsh 2x2 tile schema (printed example, n=4) ----------------

w('--- S4: Walsh schema at tile level, example n=4 (frame size 8) ---')
n = 4
tiles_n = tiles_list(n)
ex_bits = [1, 0, 1]   # x: tile A=(4,1)->1, B=(3,1)->0, C=(4,2)->1
A = tiling_to_A(ex_bits, n, tiles_n)
Ad, _ = D_skew(A)
ADl = relabel_D(Ad, n)
tiles_8 = tiles_list(8)
t = A_to_tiling(ADl, 8, tiles_8)
tmap = {tiles_8[i]: t[i] for i in range(len(t))}
w('Region map (rows X=8..2, cols Y=1..7; 1/0 = tile bit, P = base path, . = none):')
w('   regions: X,Y in 5..8 -> COPY-1 = x ; X,Y in 1..4 -> COPY-2 = sigma_4(x) ;')
w('            X in 5..8, Y in 1..4 -> CROSS = A(T) full arc matrix, anti-diag X+Y=9 twins=1')
hdr = '     ' + ' '.join(f'Y={y}' for y in range(1, 8))
w(hdr)
for X in range(8, 1, -1):
    row = [f'X={X}:']
    for Y in range(1, 8):
        if X - Y >= 2:
            cell = str(tmap[(X, Y)])
        elif X - Y == 1:
            cell = 'P'
        else:
            cell = '.'
        row.append(f'  {cell}')
    w(' '.join(row))
w(f'x(T) = {ex_bits} (tiles A=(4,1),B=(3,1),C=(4,2));  A(T) arcs include base path 4->3->2->1')
w('Tile-level Walsh schema [[H,H],[H,-H]] reading:')
w('  +copy 1: copy-1 block = x ;  +copy 2: copy-2 block = sigma(x) (grid transpose, NO negation);')
w('  cross block = +copy 3 and -copy 4 INTERLEAVED: pair {i,j} appears as A[i][j] AND A[j][i]=NOT A[i][j];')
w('  the all-ones twin anti-diagonal (X+Y=2n+1) is the identity/border (the "+I" of the matrix schema).')
w('')

# ---------------- S5: affine-linearity of x -> t(D(x)) ----------------

w('--- S5: x -> t(D(x)) is an injective AFFINE-LINEAR map over GF(2) ---')
for n in (3, 4, 5):
    tiles_n = tiles_list(n)
    m = len(tiles_n)
    tiles_2n = tiles_list(2 * n)
    def tmap_of(mask):
        bits = [(mask >> k) & 1 for k in range(m)]
        return np.array(A_to_tiling(relabel_D(D_skew(tiling_to_A(bits, n, tiles_n))[0], n),
                                    2 * n, tiles_2n), dtype=np.int64)
    t0 = tmap_of(0)
    basis = [(tmap_of(1 << k) ^ t0) for k in range(m)]
    # exactness: t(x) == t0 ^ XOR of basis vectors
    rng = np.random.default_rng(1)
    ok = 0; trials = min(1 << m, 200)
    for mask in (range(1 << m) if (1 << m) <= 200 else rng.integers(0, 1 << m, 200)):
        mask = int(mask)
        pred = t0.copy()
        for k in range(m):
            if (mask >> k) & 1:
                pred ^= basis[k]
        ok += np.array_equal(pred, tmap_of(mask))
    Bmat = np.array(basis, dtype=np.int64) % 2
    # GF(2) rank
    Bm = Bmat.copy(); r = 0
    for col in range(Bm.shape[1]):
        piv = next((i for i in range(r, Bm.shape[0]) if Bm[i, col]), None)
        if piv is None:
            continue
        Bm[[r, piv]] = Bm[[piv, r]]
        for i in range(Bm.shape[0]):
            if i != r and Bm[i, col]:
                Bm[i] ^= Bm[r]
        r += 1
    w(f'n={n}: affine law holds {ok}/{trials} (sampled);  GF(2) rank of linear part = {r} '
      f'(= C(n-1,2) = {m} => injective);  blue-subspace dim at 2n = n(n-1) = {n*(n-1)};  '
      f'density of D-image in its coset = 2^-{(n-1)*(n+2)//2}')
w('')

# ---------------- S6: blue-FRAMABILITY of D(T) (other Ham paths) ----------------

w('--- S6: does D(T) admit ANY blue frame (some Ham path with grid-symmetric tiling)? ---')

def all_ham_paths(Alist):
    N = len(Alist)
    adj = [[j for j in range(N) if Alist[i][j]] for i in range(N)]
    res = []
    def dfs(path, used):
        if len(path) == N:
            res.append(path.copy()); return
        for j in adj[path[-1]]:
            if not used[j]:
                used[j] = True; path.append(j)
                dfs(path, used)
                path.pop(); used[j] = False
    for s0 in range(N):
        used = [False] * N; used[s0] = True
        dfs([s0], used)
    return res

def anti_involutions(A, cap=200000):
    """All involutive anti-automorphisms r (r^2=id, A[a][b]==A[r(b)][r(a)])."""
    N = A.shape[0]
    res = []
    r = [-1] * N
    nodes = [0]
    def ok_pair(a, b):
        for c in range(N):
            d = r[c]
            if d == -1 or c == a or c == b:
                continue
            if A[a][c] != A[d][b]: return False
            if A[a][d] != A[c][b]: return False
            if b != a:
                if A[b][c] != A[d][a]: return False
                if A[b][d] != A[c][a]: return False
        return True
    def bt():
        nodes[0] += 1
        if nodes[0] > cap:
            return
        a = next((v for v in range(N) if r[v] == -1), None)
        if a is None:
            res.append(r.copy()); return
        if ok_pair(a, a):
            r[a] = a; bt(); r[a] = -1
        for b in range(a + 1, N):
            if r[b] != -1:
                continue
            if ok_pair(a, b):
                r[a] = b; r[b] = a
                bt()
                r[a] = r[b] = -1
    bt()
    return res, nodes[0] <= cap

def reversed_ham_path(A, r):
    """Ham path p with r(p_i) = p_{N+1-i}, if one exists (works for odd and even N).
    A = adjacency as list-of-lists."""
    N = len(A)
    fixed = [v for v in range(N) if r[v] == v]
    pairs = sorted({tuple(sorted((v, r[v]))) for v in range(N) if r[v] != v})
    npairs = len(pairs)
    if N % 2 == 1:
        if len(fixed) != 1:
            return None
        f = fixed[0]
        right = []; used = set()
        def dfs1(cur):
            if len(right) == npairs:
                return True
            for pq in pairs:
                if pq in used:
                    continue
                for wv in pq:
                    if A[cur][wv]:
                        used.add(pq); right.append(wv)
                        if dfs1(wv):
                            return True
                        used.discard(pq); right.pop()
            return False
        if not dfs1(f):
            return None
        left = [r[v] for v in reversed(right)]
        p = left + [f] + right
    else:
        if fixed:
            return None
        p = None
        for pq in pairs:
            for (s0, t0) in (pq, pq[::-1]):
                if r[s0] != t0 or not A[s0][t0]:
                    continue
                right = [t0]; used = {pq}
                def dfs2(cur):
                    if len(right) == npairs:
                        return True
                    for qq in pairs:
                        if qq in used:
                            continue
                        for wv in qq:
                            if A[cur][wv]:
                                used.add(qq); right.append(wv)
                                if dfs2(wv):
                                    return True
                                used.discard(qq); right.pop()
                    return False
                if dfs2(t0):
                    left = [r[v] for v in reversed(right)]
                    p = left + right
                    break
            if p:
                break
        if not p:
            return None
    assert all(A[p[i]][p[i + 1]] for i in range(N - 1)), 'reversed path invalid'
    return p

for n in (3, 4, 5):
    tiles_2n = tiles_list(2 * n)
    sig2n = sigma_map(2 * n, tiles_2n)
    N = 2 * n
    w(f'  n={n} (D on {N} vertices), per iso class of T:')
    for idx, A in enumerate(iso_classes(n)):
        Ad = D_skew(A)[0]
        Alist = [[int(v) for v in row] for row in Ad]
        paths = all_ham_paths(Alist)
        nblue = 0
        for p in paths:
            bits = [Alist[p[N - X]][p[N - Y]] for (X, Y) in tiles_2n]
            if is_gridsym(bits, sig2n):
                nblue += 1
        # cross-validate with anti-involution method
        invs, complete = anti_involutions(Ad)
        method_blue = any(reversed_ham_path(Alist, r) for r in invs)
        agree = (nblue > 0) == method_blue
        w(f'    T#{idx} scores={scores(A)} H(T)={H_count(A)}: H(D)={len(paths)}, '
          f'blue Ham paths of D = {nblue}  ({100.0*nblue/len(paths):.2f}%)  '
          f'anti-involutions(D)={len(invs)}{"" if complete else "(capped)"}  '
          f'method-agrees={agree}')
w('')

# ---------------- S7: iterated doubling tower W_2^k ----------------

w('--- S7: tower W_{2^k} = D^k(point); Gray-code canonical frame; Walsh structure ---')

def gray(t):
    return t ^ (t >> 1)

def lowbit_pos(v):
    return (v & -v).bit_length() - 1

def pred_arc(a, b):
    """Closed form: arc a->b in W_{2^k} iff (b_L) XOR parity(popcount((a&b)>>(L+1))) == 1,
       L = lowest differing bit of a,b."""
    L = lowbit_pos(a ^ b)
    S = bin((a & b) >> (L + 1)).count('1') & 1
    return ((b >> L) & 1) ^ S

Wt = np.zeros((1, 1), dtype=np.int64)
towers = {1: Wt}
for k in range(1, 5):
    Wt = D_skew(Wt)[0]
    towers[1 << k] = Wt

for sz in (2, 4, 8, 16):
    A = towers[sz]
    # closed-form check
    bad = sum(1 for a in range(sz) for b in range(sz) if a != b and int(A[a, b]) != pred_arc(a, b))
    # recursive canonical path == Gray code
    p = [0]
    while len(p) < sz:
        p = p + [len(p) + v for v in reversed(p)]
    gray_ok = (p == [gray(t) for t in range(sz)])
    path_ok = all(A[p[i], p[i + 1]] for i in range(sz - 1))
    w(f'order {sz}: closed-form arc formula errors = {bad}/{sz*(sz-1)};  '
      f'canonical recursive path == Gray code: {gray_ok};  path valid: {path_ok}')

w('')
w('  W16 tiling in canonical (Gray) frame: label X <-> vertex gray(16-X)')
sz = 16
A16 = towers[16]
tiles16 = tiles_list(16)
sig16 = sigma_map(16, tiles16)
vert_of_label = [None] + [int(gray(16 - X)) for X in range(1, 17)]   # 1-indexed
for kk in range(16, 1, -1):
    assert A16[vert_of_label[kk], vert_of_label[kk - 1]] == 1, 'Gray base path missing'
t16 = [int(A16[vert_of_label[X], vert_of_label[Y]]) for (X, Y) in tiles16]
wt16 = sum(t16)
viol = [i for i in range(len(t16)) if t16[i] != t16[sig16[i]]]
fixedt = [i for i, (X, Y) in enumerate(tiles16) if X + Y == 17]
crosspred = set(i for i, (X, Y) in enumerate(tiles16) if X > 8 and Y <= 8 and X + Y != 17)
w(f'  m = {len(tiles16)} tiles; weight(t16) = {wt16}  '
  f'(predicted by recursion 2w+(n-1)+C(n,2) from 0: 0,0,2,13,61 -> {wt16==61})')
w(f'  grid-symmetric: {is_gridsym(t16, sig16)};  violated tiles = {len(viol)} '
  f'(predicted top-level cross off-twin = 56: {set(viol)==crosspred});  '
  f'sigma-fixed tiles (anti-diag) = {len(fixedt)}, all bits 1: {all(t16[i] for i in fixedt)}')
w('  105-bit vector (explorer order, bit=1 iff arc x->y):')
w('  ' + ''.join(map(str, t16)))
w('  ASCII grid (rows X=16..2, cols Y=1..15; P=base path, .=upper triangle):')
for X in range(16, 1, -1):
    row = []
    for Y in range(1, 16):
        if X - Y >= 2:
            row.append(str(int(A16[vert_of_label[X], vert_of_label[Y]])))
        elif X - Y == 1:
            row.append('P')
        else:
            row.append('.')
    w(f'   X={X:2d}: ' + ' '.join(row))

w('')
w('  Walsh battery on the 105 tiling bits (a=gray(16-X), b=gray(16-Y)):')
TM = lambda v: bin(v).count('1') & 1
cands = {
    'exact skew-Walsh ((b_L) ^ par(pc((a&b)>>(L+1))))': lambda X, Y, a, b: pred_arc(a, b),
    'pure Walsh par(pc(a&b))':                          lambda X, Y, a, b: TM(a & b),
    'frame Walsh par(pc((X-1)&(Y-1)))':                 lambda X, Y, a, b: TM((X - 1) & (Y - 1)),
    'frame Walsh par(pc((16-X)&(16-Y)))':               lambda X, Y, a, b: TM((16 - X) & (16 - Y)),
    'Thue-Morse TM(a)^TM(b)':                           lambda X, Y, a, b: TM(a) ^ TM(b),
    'TM(a|b)':                                          lambda X, Y, a, b: TM(a | b),
    'twin-direction only (b_L)':                        lambda X, Y, a, b: (b >> lowbit_pos(a ^ b)) & 1,
    'sign only par(pc((a&b)>>(L+1)))':                  lambda X, Y, a, b: TM((a & b) >> (lowbit_pos(a ^ b) + 1)),
}
for name, f in cands.items():
    match = sum(1 for i, (X, Y) in enumerate(tiles16)
                if t16[i] == f(X, Y, vert_of_label[X], vert_of_label[Y]))
    w(f'    {name:52s}: {match}/105')
w('')

# blue-framability of the tower W8, W16
w('  Blue-framability of tower tournaments (anti-involutions + reversed Ham path):')
for sz in (4, 8, 16):
    A = towers[sz]
    invs, complete = anti_involutions(A, cap=2000000)
    Alist = [[int(v) for v in row] for row in A]
    found = None
    for r in invs:
        p = reversed_ham_path(Alist, r)
        if p:
            found = (r, p)
            break
    if found:
        r, p = found
        tilesz = tiles_list(sz)
        sigz = sigma_map(sz, tilesz)
        bits = [Alist[p[sz - X]][p[sz - Y]] for (X, Y) in tilesz]
        ok = is_gridsym(bits, sigz)
        w(f'    W{sz}: anti-involutions={len(invs)}{"" if complete else "(capped)"};  BLUE FRAME EXISTS, '
          f'path={p}, tiling weight={sum(bits)}/{len(bits)}, gridsym verified={ok}')
    else:
        w(f'    W{sz}: anti-involutions={len(invs)}{"" if complete else "(capped)"};  no blue frame found '
          f'=> class PURE BLACK' if complete else f'    W{sz}: search capped, inconclusive')
w('')

# ---------------- S8: T15 = core of normalized order-16 tower ----------------

w('--- S8: T15 (core of normalized S16 from iterated S-doubling of [1]) ---')
S = np.array([[1]], dtype=np.int64)
for k in range(4):
    nS = S.shape[0]
    I = np.eye(nS, dtype=np.int64)
    S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
w(f'S16 skew-Hadamard: {is_skew_hadamard(S)}')
S16n = normalize_first_row(S)
w(f'normalized S16 skew-Hadamard: {is_skew_hadamard(S16n)}')
T15 = core_tournament(S16n)
w(f'T15: scores = {scores(T15)}  DRT = {is_DRT(T15)}')
t0 = time.time()
H15 = H_count(T15)
w(f'H(T15) = {H15}   ({time.time()-t0:.1f}s)')

# variant chain: normalize at 8, double, normalize, core
S4v = np.array([[1]], dtype=np.int64)
for k in range(3):
    nS = S4v.shape[0]; I = np.eye(nS, dtype=np.int64)
    S4v = np.block([[S4v, S4v], [S4v - 2 * I, 2 * I - S4v]])
S8v = normalize_first_row(S4v)
I8 = np.eye(8, dtype=np.int64)
S16v = normalize_first_row(np.block([[S8v, S8v], [S8v - 2 * I8, 2 * I8 - S8v]]))
T15v = core_tournament(S16v)
t0 = time.time()
H15v = H_count(T15v)
w(f'variant (normalize at 8 first): scores={sorted(set(scores(T15v)))} DRT={is_DRT(T15v)} '
  f'H={H15v}  ({time.time()-t0:.1f}s)  same H as T15: {H15v == H15}')

tiles15 = tiles_list(15)
sig15 = sigma_map(15, tiles15)
p15 = ham_path(T15)
T15l = [[int(v) for v in row] for row in T15]
bits15 = [T15l[p15[15 - X]][p15[15 - Y]] for (X, Y) in tiles15]
viol15 = [i for i in range(91) if bits15[i] != bits15[sig15[i]]]
w(f'insertion-frame Ham path p = {p15}')
w(f'T15 tiling (insertion frame), m=91: weight={sum(bits15)}  gridsym={is_gridsym(bits15, sig15)}  '
  f'violated tiles={len(viol15)} of {91 - 7} non-fixed (={ (91-7)//2 } sigma-pairs)')
w('  91-bit vector: ' + ''.join(map(str, bits15)))

w('  Walsh battery on insertion-frame T15 bits (frame labels only):')
cands15 = {
    'par(pc((X-1)&(Y-1)))': lambda X, Y: TM((X - 1) & (Y - 1)),
    'par(pc((15-X)&(15-Y)))': lambda X, Y: TM((15 - X) & (15 - Y)),
    'par(pc(X&Y))': lambda X, Y: TM(X & Y),
    'constant 1': lambda X, Y: 1,
}
for name, f in cands15.items():
    match = sum(1 for i, (X, Y) in enumerate(tiles15) if bits15[i] == f(X, Y))
    w(f'    {name:26s}: {match}/91')

t0 = time.time()
invs15, complete15 = anti_involutions(T15, cap=5000000)
w(f'anti-involutions of T15: {len(invs15)} {"(complete)" if complete15 else "(CAPPED)"}  '
  f'({time.time()-t0:.1f}s)')
blue15 = None
for r in invs15:
    p = reversed_ham_path(T15l, r)
    if p:
        blue15 = (r, p)
        break
if blue15:
    r, p = blue15
    bitsB = [T15l[p[15 - X]][p[15 - Y]] for (X, Y) in tiles15]
    okB = is_gridsym(bitsB, sig15)
    w(f'BLUE FRAME of T15 EXISTS: anti-involution r={r}')
    w(f'  reversed Ham path p={p}')
    w(f'  blue tiling: weight={sum(bitsB)}/91  gridsym verified={okB}')
    w('  blue 91-bit vector: ' + ''.join(map(str, bitsB)))
else:
    w('NO blue frame of T15 found' + ('' if complete15 else ' (search capped, inconclusive)') +
      (' => T15 class PURE BLACK (no involutive anti-automorphism reverses a Ham path)'
       if complete15 else ''))
w('')
w('=== done ===')
OUT.close()
