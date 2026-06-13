#!/usr/bin/env python3
"""
verify_C_tiling_kps1.py -- kind-pasteur-2026-06-09-S1
ADVERSARIAL VERIFICATION of branch C-tiling (skew_doubling_tiling_walsh_kps1).

Independent recomputation (fresh code, pure-python lists, own conventions derived
from definitions, NOT reusing the sibling's relabel/frame/tiling functions):

  V1  D(T) built from arc-level definition; cross-checked vs core D_skew.
  V2  For ALL (labeled T, Ham-path frame) pairs at n=3,4,5:
        - canonical Ham path of D valid
        - BLUE verdict of the canonical tiling of D  (claim: always 0)
        - defect t XOR sigma(t) == c_n (off-twin cross indicator)  (claim: always)
        - block transform formula (copy1=x, copy2=x.sigma, cross=A / twin=1)
        - weight law |t| = 2|x| + (n-1) + C(n,2)
      NOTE: this covers strictly more frames than the sibling's S2 (every Ham path
      of every labeled T, not just one insertion path / base-path-framed tilings).
  V3  Blue-Ham-path counts of D(T) for ALL iso classes n=3,4,5 by brute force
      (own path enumeration + own grid-sym check).  Claims: C3 9/45, regT5 185/15505,
      all 16 other classes 0.  Plus the score-multiset necessary condition
      (anti-aut of D needs T regular + n odd) as a structural cross-check.
  V4  Tower W_{2^k} via own doubling; closed-form arc law re-implemented and tested
      at orders 2..16 AND extended to 32; Gray-code path; canonical-frame tiling
      weights at orders 4,8,16 (claim 2,13,61); W16 105-bit vector compared to the
      sibling's published string; violated-tile set vs top-level cross-offtwin;
      fractal all-ones anti-diagonals at scales 16,8,4; Walsh battery recount.
  V5  T15: own S-doubling/normalization/core; skew-Hadamard, scores, DRT checks;
      H(T15) recomputed TWO ways (core DP + own inclusion-exclusion walk-count
      algorithm, a genuinely different algorithm); ALL isomorphisms T15 -> T15^op
      via own backtracking (superset of anti-involutions); published insertion-frame
      path validity + tiling weight 48 / 32 violations.

Output -> 05-knowledge/results/verify_C_tiling_kps1.out
"""
import sys, time
import numpy as np

sys.path.insert(0, '04-computation')
from skew_doubling_core_kps1 import D_skew, H_count, iso_classes, scores as core_scores

OUT = open('05-knowledge/results/verify_C_tiling_kps1.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(str(s) + '\n'); OUT.flush(); print(s)

t_start = time.time()

# ---------------- independent primitives ----------------

def my_all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for b in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (b >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def my_D(A):
    """D(T) from the arc-level definition of M' = [[M, M+I],[M-I, -M]]:
       copy1 = T;  copy2 = T^op;  cross: u_i->v_i (twin), u_i->v_j iff A[i][j] (i!=j)."""
    n = len(A); N = 2 * n
    Ad = [[0] * N for _ in range(N)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if A[i][j]:
                Ad[i][j] = 1            # copy 1: u_i -> u_j
                Ad[n + j][n + i] = 1    # copy 2 = T^op: v_j -> v_i
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j]:
                Ad[i][n + j] = 1        # (M+I)[i][j] > 0
            else:
                Ad[n + j][i] = 1        # (M-I)^T entry: v_j -> u_i iff A[j][i], j!=i
    return Ad

def my_ham_paths(A):
    N = len(A)
    out = [[j for j in range(N) if A[i][j]] for i in range(N)]
    res = []
    def go(p, used):
        if len(p) == N:
            res.append(p[:]); return
        for v in out[p[-1]]:
            if not (used >> v) & 1:
                p.append(v); go(p, used | (1 << v)); p.pop()
    for s in range(N):
        go([s], 1 << s)
    return res

def tiles_of(N):
    return [(x, y) for y in range(1, N - 1) for x in range(N, y + 1, -1)]

def sigma_of(N, tiles):
    pos = {t: i for i, t in enumerate(tiles)}
    return [pos[(N - y + 1, N - x + 1)] for (x, y) in tiles]

def frame_bits(Ad, P, tiles):
    """Frame digraph Ad by Ham path P: label(P[t]) = N - t; bit=1 iff arc x->y."""
    N = len(Ad)
    vert = [None] * (N + 1)
    for t, v in enumerate(P):
        vert[N - t] = v
    for k in range(N, 1, -1):
        assert Ad[vert[k]][vert[k - 1]] == 1, 'frame path not a path'
    return [Ad[vert[X]][vert[Y]] for (X, Y) in tiles], vert

# ---------------- V1: my_D vs core D_skew ----------------

w('=== ADVERSARIAL VERIFICATION of branch C-tiling ===')
w('')
w('--- V1: independent D(T) construction vs core D_skew ---')
for n in (3, 4):
    bad = 0; tot = 0
    for A in my_all_tournaments(n):
        tot += 1
        Ad1 = my_D(A)
        Ad2 = D_skew(np.array(A, dtype=np.int64))[0].tolist()
        if Ad1 != Ad2:
            bad += 1
    w(f'n={n}: my_D == D_skew in {tot - bad}/{tot} labeled tournaments')
# spot check n=5
bad = 0
for k, A in enumerate(my_all_tournaments(5)):
    if k % 37 == 0:
        if my_D(A) != D_skew(np.array(A, dtype=np.int64))[0].tolist():
            bad += 1
w(f'n=5 spot checks (every 37th): mismatches = {bad}')
w('')

# anchors
C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
regT5 = [[1 if (j - i) % 5 in (1, 2) else 0 for j in range(5)] for i in range(5)]
P7 = [[1 if (j - i) % 7 in (1, 2, 4) else 0 for j in range(7)] for i in range(7)]
w('--- anchors (own path enumeration) ---')
w(f'H(C3) = {len(my_ham_paths(C3))} (expect 3)')
w(f'H(regular T5) = {len(my_ham_paths(regT5))} (expect 15)')
w(f'H(Paley T7) = {len(my_ham_paths(P7))} (expect 189)')
w(f'H(D(C3)) = {len(my_ham_paths(my_D(C3)))} (expect 45)')
w('')

# ---------------- V2: ALL (T, frame) pairs, n=3..5 ----------------

w('--- V2: BLUE / defect / block formula / weight over ALL (labeled T, Ham path) pairs ---')
for n in (3, 4, 5):
    N = 2 * n
    tiles2 = tiles_of(N)
    sig2 = sigma_of(N, tiles2)
    tiles1 = tiles_of(n)
    cset = frozenset(i for i, (X, Y) in enumerate(tiles2)
                     if X > n and Y <= n and X + Y != N + 1)
    blue_cnt = total = formula_bad = defect_bad = weight_bad = path_bad = 0
    for A in my_all_tournaments(n):
        Ad = my_D(A)
        for q in my_ham_paths(A):
            total += 1
            P = q + [n + v for v in reversed(q)]
            if not all(Ad[P[t]][P[t + 1]] for t in range(N - 1)):
                path_bad += 1
                continue
            t_bits, _ = frame_bits(Ad, P, tiles2)
            if all(t_bits[i] == t_bits[sig2[i]] for i in range(len(t_bits))):
                blue_cnt += 1
            dft = frozenset(i for i in range(len(t_bits)) if t_bits[i] != t_bits[sig2[i]])
            if dft != cset:
                defect_bad += 1
            # framed T arc matrix AF (labels 1..n): label(q[t]) = n - t
            vT = [None] * (n + 1)
            for t, v in enumerate(q):
                vT[n - t] = v
            AF = [[0] * (n + 1) for _ in range(n + 1)]
            for a in range(1, n + 1):
                for b in range(1, n + 1):
                    if a != b:
                        AF[a][b] = A[vT[a]][vT[b]]
            ok = True
            for i, (X, Y) in enumerate(tiles2):
                if X > n and Y > n:
                    pred = AF[X - n][Y - n]
                elif X <= n and Y <= n:
                    pred = AF[n + 1 - Y][n + 1 - X]
                else:
                    ii, jj = X - n, n + 1 - Y
                    pred = 1 if ii == jj else AF[ii][jj]
                if t_bits[i] != pred:
                    ok = False
                    break
            if not ok:
                formula_bad += 1
            xw = sum(AF[x][y] for (x, y) in tiles1)
            if sum(t_bits) != 2 * xw + (n - 1) + n * (n - 1) // 2:
                weight_bad += 1
    w(f'n={n}: (T,frame) pairs = {total}  canonical-D-path failures = {path_bad}  '
      f'BLUE = {blue_cnt}  defect!=c_n = {defect_bad}  formula fails = {formula_bad}  '
      f'weight fails = {weight_bad}')
w('')

# ---------------- V3: blue Ham paths of D(T), all iso classes n=3..5 ----------------

w('--- V3: blue Ham paths of D(T), brute force over ALL Ham paths of D, per iso class ---')
w('    structural necessary condition: anti-aut of D needs score multiset of D invariant')
w('    under s -> 2n-1-s, which forces n odd AND T regular.')
for n in (3, 4, 5):
    N = 2 * n
    tiles2 = tiles_of(N)
    sig2 = sigma_of(N, tiles2)
    for idx, Anp in enumerate(iso_classes(n)):
        A = Anp.tolist()
        Ad = my_D(A)
        paths = my_ham_paths(Ad)
        nb = 0
        for p in paths:
            vert = [None] * (N + 1)
            for t, v in enumerate(p):
                vert[N - t] = v
            bits = [Ad[vert[X]][vert[Y]] for (X, Y) in tiles2]
            if all(bits[i] == bits[sig2[i]] for i in range(len(bits))):
                nb += 1
        sc = sorted(core_scores(Anp))
        regular = len(set(sc)) == 1
        struct_possible = regular and (n % 2 == 1)
        w(f'  n={n} T#{idx} scores={sc} H(T)={len(my_ham_paths(A))}: '
          f'H(D)={len(paths)}  blue={nb}  struct-possible={struct_possible}')
w('')

# ---------------- V4: tower ----------------

w('--- V4: tower W_{2^k}, own doubling; closed form at 2..32; Gray path; W16 frame ---')

def gray(t):
    return t ^ (t >> 1)

def my_pred(a, b):
    d = a ^ b
    L = (d & -d).bit_length() - 1
    s = bin((a & b) >> (L + 1)).count('1') & 1
    return ((b >> L) & 1) ^ s

W = [[0]]
orders = {1: W}
for k in range(1, 6):
    W = my_D(W)
    orders[1 << k] = W

# cross-check tower vs core-built
Wc = np.zeros((1, 1), dtype=np.int64)
for k in range(1, 6):
    Wc = D_skew(Wc)[0]
w(f'order-32 tower: my_D-built == D_skew-built: {orders[32] == Wc.tolist()}')

for sz in (2, 4, 8, 16, 32):
    A = orders[sz]
    bad = sum(1 for a in range(sz) for b in range(sz) if a != b and A[a][b] != my_pred(a, b))
    p = [gray(t) for t in range(sz)]
    pv = all(A[p[i]][p[i + 1]] for i in range(sz - 1))
    w(f'order {sz}: closed-form errors = {bad}/{sz * (sz - 1)}; Gray path valid = {pv}'
      + ('   <-- EXTENSION beyond sibling claim (they tested <=16)' if sz == 32 else ''))

w('')
w('  canonical (Gray) frame tilings of W4, W8, W16:')
sib_w16 = ('1000000011111111011001011001110001100111011010101010111010101110'
           '11001110101101000111100011110101110101101')
for sz in (4, 8, 16):
    A = orders[sz]
    tiles = tiles_of(sz)
    sig = sigma_of(sz, tiles)
    P = [gray(t) for t in range(sz)]
    bits, vert = frame_bits(A, P, tiles)
    h = sz // 2
    crosspred = frozenset(i for i, (X, Y) in enumerate(tiles)
                          if X > h and Y <= h and X + Y != sz + 1)
    viol = frozenset(i for i in range(len(bits)) if bits[i] != bits[sig[i]])
    gridsym = (len(viol) == 0)
    # fractal all-ones anti-diagonals at every scale
    fractal_ok = True
    s = sz
    while s >= 2:
        for blk in range(sz // s):
            lo = blk * s
            for i, (X, Y) in enumerate(tiles):
                if lo < X <= lo + s and lo < Y <= lo + s and X + Y == 2 * lo + s + 1:
                    if bits[i] != 1:
                        fractal_ok = False
        s //= 2
    msg = (f'  W{sz}: weight = {sum(bits)}  gridsym = {gridsym}  '
           f'violations = {len(viol)} (== top cross-offtwin: {viol == crosspred})  '
           f'fractal anti-diagonals all-ones at every scale: {fractal_ok}')
    if sz == 16:
        msg += f'  105-bit == sibling string: {"".join(map(str, bits)) == sib_w16}'
    w(msg)

w('')
w('  Walsh battery recount on W16 canonical-frame bits:')
tiles16 = tiles_of(16)
P16 = [gray(t) for t in range(16)]
bits16, vert16 = frame_bits(orders[16], P16, tiles16)
TM = lambda v: bin(v).count('1') & 1
batt = {
    'exact skew-Walsh': lambda a, b: my_pred(a, b),
    'pure Walsh par(pc(a&b))': lambda a, b: TM(a & b),
    'Thue-Morse TM(a)^TM(b)': lambda a, b: TM(a) ^ TM(b),
    'twin-direction only (b_L)': lambda a, b: (b >> ((a ^ b) & -(a ^ b)).bit_length() - 1) & 1,
    'sign only': lambda a, b: TM((a & b) >> (((a ^ b) & -(a ^ b)).bit_length())),
}
for name, f in batt.items():
    m = sum(1 for i, (X, Y) in enumerate(tiles16)
            if bits16[i] == f(vert16[X], vert16[Y]))
    w(f'    {name:28s}: {m}/105')
w('')

# ---------------- V5: T15 ----------------

w('--- V5: T15 (core of normalized order-16 S-tower) ---')

def sdouble(S):
    m = len(S)
    N = 2 * m
    R = [[0] * N for _ in range(N)]
    for i in range(m):
        for j in range(m):
            R[i][j] = S[i][j]
            R[i][m + j] = S[i][j]
            R[m + i][j] = S[i][j] - (2 if i == j else 0)
            R[m + i][m + j] = (2 if i == j else 0) - S[i][j]
    return R

def my_is_skh(S):
    m = len(S)
    for i in range(m):
        for j in range(m):
            if abs(S[i][j]) != 1:
                return False
            if S[i][j] + S[j][i] != (2 if i == j else 0):
                return False
    Snp = np.array(S, dtype=np.int64)
    return np.array_equal(Snp @ Snp.T, m * np.eye(m, dtype=np.int64))

S = [[1]]
for k in range(4):
    S = sdouble(S)
w(f'S16 skew-Hadamard (own check): {my_is_skh(S)}')
# normalize first row
for j in range(1, 16):
    if S[0][j] == -1:
        for i in range(16):
            S[i][j] *= -1
        for i in range(16):
            S[j][i] *= -1
w(f'normalized S16 skew-Hadamard: {my_is_skh(S)}  first row all +1: {all(v == 1 for v in S[0])}')
T15 = [[1 if (S[i + 1][j + 1] - (1 if i == j else 0)) > 0 else 0 for j in range(15)]
       for i in range(15)]
sc15 = [sum(row) for row in T15]
w(f'T15 scores: {sorted(set(sc15))} x {len(sc15)}  (expect all 7)')
lam = set()
for u in range(15):
    for v in range(u + 1, 15):
        lam.add(sum(T15[u][k] * T15[v][k] for k in range(15)))
w(f'common out-neighbour counts over pairs: {sorted(lam)}  DRT = {len(lam) == 1}')

t0 = time.time()
H15_dp = H_count(np.array(T15, dtype=np.int64))
t_dp = time.time() - t0
w(f'H(T15) core DP        = {H15_dp}   ({t_dp:.1f}s)')

# inclusion-exclusion (independent algorithm): H = sum_S (-1)^(15-|S|) #walks_14(T[S])
t0 = time.time()
A15 = np.array(T15, dtype=np.int64)
nv = 15
tot = 0
for Smask in range(1, 1 << nv):
    idx = [v for v in range(nv) if (Smask >> v) & 1]
    k = len(idx)
    B = A15[np.ix_(idx, idx)]
    vec = np.ones(k, dtype=np.int64)
    for _ in range(nv - 1):
        vec = B @ vec
    wk = int(vec.sum())
    if wk:
        tot += wk if ((nv - k) % 2 == 0) else -wk
t_ie = time.time() - t0
w(f'H(T15) inclusion-excl = {tot}   ({t_ie:.1f}s)   match = {tot == H15_dp}   '
  f'sibling value 198335025 match = {tot == 198335025}')

# all isomorphisms T15 -> T15^op (= all anti-automorphisms), own backtracking
w('')
t0 = time.time()
res_anti = []
nodes = [0]
f = [-1] * 15
used = [False] * 15
def bt(k):
    nodes[0] += 1
    if k == 15:
        res_anti.append(f[:]); return
    for u in range(15):
        if used[u]:
            continue
        ok = True
        for j in range(k):
            if T15[j][k] != T15[u][f[j]]:
                ok = False
                break
        if ok:
            f[k] = u; used[u] = True
            bt(k + 1)
            used[u] = False; f[k] = -1
bt(0)
n_anti = len(res_anti)
n_inv = sum(1 for r in res_anti if all(r[r[v]] == v for v in range(15)))
w(f'ALL anti-automorphisms (isos T15 -> T15^op): {n_anti}  '
  f'(backtrack nodes = {nodes[0]}, {time.time()-t0:.1f}s)')
w(f'  involutive among them: {n_inv}   => T15 self-converse: {n_anti > 0};  '
  f'anti-involutions = {n_inv} (sibling claims 0)')
w(f'  pure-black verdict (no blue frame in ANY Ham-path frame): {n_inv == 0}')

# published insertion-frame path and tiling
p15 = [14, 13, 10, 12, 8, 9, 11, 7, 6, 4, 5, 3, 2, 1, 0]
valid = all(T15[p15[i]][p15[i + 1]] for i in range(14))
w(f'sibling insertion path valid in T15: {valid}')
tiles15 = tiles_of(15)
sig15 = sigma_of(15, tiles15)
bits15, _ = frame_bits(T15, p15, tiles15)
viol15 = sum(1 for i in range(len(bits15)) if bits15[i] != bits15[sig15[i]])
w(f'insertion-frame tiling: weight = {sum(bits15)}/91 (claim 48)  '
  f'gridsym = {all(bits15[i] == bits15[sig15[i]] for i in range(91))} (claim False)  '
  f'violated tiles = {viol15} (claim 32)')

# anti-automorphism counts for D(C3), D(regT5) -- sibling claims exactly 1 anti-involution each
w('')
for name, T in (('D(C3)', C3), ('D(regT5)', regT5)):
    Ad = my_D(T)
    Nn = len(Ad)
    res2 = []
    g = [-1] * Nn
    used2 = [False] * Nn
    def bt2(k):
        if k == Nn:
            res2.append(g[:]); return
        for u in range(Nn):
            if used2[u]:
                continue
            ok = True
            for j in range(k):
                if Ad[j][k] != Ad[u][g[j]]:
                    ok = False
                    break
            if ok:
                g[k] = u; used2[u] = True
                bt2(k + 1)
                used2[u] = False; g[k] = -1
    bt2(0)
    ninv = sum(1 for r in res2 if all(r[r[v]] == v for v in range(Nn)))
    w(f'{name}: anti-automorphisms = {len(res2)}, involutive = {ninv} (sibling claims 1)')

w('')
w(f'=== done ({time.time()-t_start:.1f}s total) ===')
OUT.close()
