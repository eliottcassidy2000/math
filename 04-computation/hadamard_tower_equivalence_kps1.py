# hadamard_tower_equivalence_kps1.py
# kind-pasteur-2026-06-09-S1, Branch D: Hadamard equivalence of skew tower vs Sylvester/Walsh tower.
#
# Tasks:
#  1. Build skew tower S_k (S -> [[S,S],[S-2I,2I-S]] from [[1]]) and Sylvester tower H_k
#     at orders 2,4,8,16,32. Verify Hadamard + skew-type at each order.
#  2. Orders 4, 8: find explicit monomial (signed-permutation) equivalence skew_k ~ sylvester_k.
#  3. Order 16: invariants — quadruple ("Hall set") profile, Smith normal form over Z of H
#     and of (H+J)/2, GF(2) rank of (H+J)/2. Verdict: same class or different?
#  4. Order 32: FULL quadruple profile (C(32,4)=35960 — cheap, better than sampling).
#  5. Walsh-domain shape of S: N = W_k^T S (so C = N/2^k), k=3,4. Histogram, sparsity, structure.
#  6. Butterfly: S_{2m} x = [S_m(x1+x2); S_m(x1-x2) - 2(x1-x2)]. Verify + time vs dense at 1024.

import sys, os, time, itertools
import numpy as np

OUTPATH = os.path.join("05-knowledge", "results", "hadamard_tower_equivalence_kps1.out")

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
    def w(self, *args):
        s = " ".join(str(a) for a in args)
        print(s)
        self.f.write(s + "\n")
        self.f.flush()

T = Tee(OUTPATH)
W = T.w

# ---------------------------------------------------------------- towers
def skew_tower(order):
    S = np.array([[1]], dtype=np.int64)
    while S.shape[0] < order:
        n = S.shape[0]
        I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
    return S

def sylvester(order):
    H = np.array([[1]], dtype=np.int64)
    while H.shape[0] < order:
        H = np.block([[H, H], [H, -H]])
    return H

W("=" * 78)
W("PART 1: BUILD TOWERS, VERIFY HADAMARD + SKEW-TYPE")
W("=" * 78)
for order in (2, 4, 8, 16, 32):
    S = skew_tower(order)
    Hs = sylvester(order)
    n = order
    I = np.eye(n, dtype=np.int64)
    had_S = np.array_equal(S @ S.T, n * I)
    skew_S = np.array_equal(S + S.T, 2 * I)
    had_H = np.array_equal(Hs @ Hs.T, n * I)
    sym_H = np.array_equal(Hs, Hs.T)
    W(f"order {n:3d}: skew tower Hadamard={had_S} skew-type(S+S^T=2I)={skew_S} | "
      f"Sylvester Hadamard={had_H} symmetric={sym_H}")
    assert had_S and skew_S and had_H

def matstr(M):
    return "\n".join("".join("+" if v > 0 else "-" for v in row) for row in M)

W("\nskew tower S_16 (+/- form):")
W(matstr(skew_tower(16)))
W("\nskew tower S_32 (+/- form):")
W(matstr(skew_tower(32)))

# ---------------------------------------------------------------- monomial equivalence search
def colsig(rows_list, n):
    # canonical multiset of partial columns up to column sign
    cols = []
    for c in range(n):
        v = tuple(r[c] for r in rows_list)
        nv = tuple(-x for x in v)
        cols.append(v if v <= nv else nv)
    cols.sort()
    return cols

def monomial_equiv(A, B, node_cap=5_000_000, time_cap=300.0):
    """Find signed row perm + signed col perm with (R A) C = B.
    Returns (status, witness) where status in {'FOUND','NONE','CAP'} and witness =
    (rowperm sigma, rowsigns eps, colperm tau, colsigns delta) meaning
    B[i,j] = eps[i] * delta[j] * A[sigma[i], tau[j]]."""
    n = A.shape[0]
    Arows = [tuple(int(v) for v in A[i]) for i in range(n)]
    Brows = [tuple(int(v) for v in B[i]) for i in range(n)]
    sigB = [None] * (n + 1)
    for k in range(1, n + 1):
        sigB[k] = colsig(Brows[:k], n)
    nodes = [0]
    t0 = time.time()
    used = [False] * n
    chosen = []     # signed A-row tuples, in target row order
    assign = []     # (source row idx, sign)
    capped = [False]

    def bt(k):
        if k == n:
            return True
        for i in range(n):
            if used[i]:
                continue
            for s in (1, -1):
                nodes[0] += 1
                if nodes[0] > node_cap or time.time() - t0 > time_cap:
                    capped[0] = True
                    return False
                row = Arows[i] if s == 1 else tuple(-x for x in Arows[i])
                chosen.append(row)
                if colsig(chosen, n) == sigB[k + 1]:
                    used[i] = True
                    assign.append((i, s))
                    if bt(k + 1):
                        return True
                    used[i] = False
                    assign.pop()
                chosen.pop()
                if capped[0]:
                    return False
        return False

    ok = bt(0)
    if not ok:
        return ("CAP" if capped[0] else "NONE"), None, nodes[0], time.time() - t0
    # reconstruct column matching: P = R A (rows permuted+signed); find C with P C = B
    P = np.array(chosen, dtype=np.int64)
    usedc = [False] * n
    tau = [None] * n
    delta = [None] * n
    for c in range(n):
        target = B[:, c]
        found = False
        for c2 in range(n):
            if usedc[c2]:
                continue
            if np.array_equal(P[:, c2], target):
                tau[c] = c2; delta[c] = 1; usedc[c2] = True; found = True; break
            if np.array_equal(P[:, c2], -target):
                tau[c] = c2; delta[c] = -1; usedc[c2] = True; found = True; break
        if not found:
            return "NONE", None, nodes[0], time.time() - t0
    sigma = [a[0] for a in assign]
    eps = [a[1] for a in assign]
    # verify
    chk = np.empty((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            chk[i, j] = eps[i] * delta[j] * A[sigma[i], tau[j]]
    assert np.array_equal(chk, B), "witness verification FAILED"
    return "FOUND", (sigma, eps, tau, delta), nodes[0], time.time() - t0

W("\n" + "=" * 78)
W("PART 2: EXPLICIT MONOMIAL EQUIVALENCE AT ORDERS 2, 4, 8 (skew_k ~ sylvester_k)")
W("=" * 78)
for order in (2, 4, 8):
    A = skew_tower(order)
    B = sylvester(order)
    status, wit, nodes, dt = monomial_equiv(A, B)
    W(f"\norder {order}: status={status}  nodes={nodes}  time={dt:.3f}s")
    if status == "FOUND":
        sigma, eps, tau, delta = wit
        W(f"  row perm sigma (target i <- source sigma[i]): {sigma}")
        W(f"  row signs eps:  {eps}")
        W(f"  col perm tau  (target j <- source tau[j]):   {tau}")
        W(f"  col signs delta:{delta}")
        W(f"  verified: sylvester[i,j] = eps[i]*delta[j]*skew[sigma[i],tau[j]]  EXACT")

# ---------------------------------------------------------------- invariants
def profile(H):
    """Full quadruple profile: counts of |sum_m H_im H_jm H_km H_lm| over all row 4-subsets."""
    n = H.shape[0]
    pairs = list(itertools.combinations(range(n), 2))
    PV = {p: H[p[0]] * H[p[1]] for p in pairs}
    counts = {}
    for (i, j, k, l) in itertools.combinations(range(n), 4):
        s = abs(int(np.dot(PV[(i, j)], PV[(k, l)])))
        counts[s] = counts.get(s, 0) + 1
    return dict(sorted(counts.items()))

def smith_normal_form_z(Ain):
    """SNF over Z, returns list of nonzero invariant factors (divisibility chain)."""
    A = [[int(x) for x in row] for row in Ain]
    m, n = len(A), len(A[0])
    t = 0
    while t < min(m, n):
        # find pivot = min abs nonzero in A[t:,t:]
        piv = None; best = None
        for i in range(t, m):
            for j in range(t, n):
                a = A[i][j]
                if a != 0 and (best is None or abs(a) < best):
                    best = abs(a); piv = (i, j)
        if piv is None:
            break
        i0, j0 = piv
        A[t], A[i0] = A[i0], A[t]
        for r in range(m):
            A[r][t], A[r][j0] = A[r][j0], A[r][t]
        while True:
            # clear column t
            again = False
            for i in range(t + 1, m):
                if A[i][t] != 0:
                    q = A[i][t] // A[t][t]
                    for j in range(t, n):
                        A[i][j] -= q * A[t][j]
                    if A[i][t] != 0:
                        A[t], A[i] = A[i], A[t]
                        again = True
            # clear row t
            for j in range(t + 1, n):
                if A[t][j] != 0:
                    q = A[t][j] // A[t][t]
                    for i in range(t, m):
                        A[i][j] -= q * A[i][t]
                    if A[t][j] != 0:
                        for i in range(t, m):
                            A[i][t], A[i][j] = A[i][j], A[i][t]
                        again = True
            if again:
                continue
            colc = all(A[i][t] == 0 for i in range(t + 1, m))
            rowc = all(A[t][j] == 0 for j in range(t + 1, n))
            if colc and rowc:
                # divisibility check
                d = A[t][t]
                bad = None
                for i in range(t + 1, m):
                    for j in range(t + 1, n):
                        if A[i][j] % d != 0:
                            bad = (i, j); break
                    if bad:
                        break
                if bad:
                    i, _ = bad
                    for jj in range(t, n):
                        A[t][jj] += A[i][jj]
                    continue
                break
        t += 1
    # collect
    diag = []
    for i in range(min(m, n)):
        if A[i][i] != 0:
            diag.append(abs(A[i][i]))
    return diag

def snf_compact(d):
    out = []
    for v in d:
        if out and out[-1][0] == v:
            out[-1][1] += 1
        else:
            out.append([v, 1])
    return " ".join(f"{v}^{c}" if c > 1 else f"{v}" for v, c in out)

def gf2_rank(M):
    rows = []
    for r in M:
        x = 0
        for v in r:
            x = (x << 1) | (int(v) & 1)
        rows.append(x)
    rank = 0
    for b in range(M.shape[1] - 1, -1, -1):
        piv = None
        for idx, r in enumerate(rows):
            if (r >> b) & 1:
                piv = idx; break
        if piv is None:
            continue
        pr = rows.pop(piv)
        rows = [r ^ pr if (r >> b) & 1 else r for r in rows]
        rank += 1
    return rank

# sympy cross-check for SNF
def snf_sympy(Ain):
    try:
        from sympy import Matrix, ZZ
        from sympy.matrices.normalforms import smith_normal_form
        Msnf = smith_normal_form(Matrix(Ain.tolist()), domain=ZZ)
        d = [abs(int(Msnf[i, i])) for i in range(min(Msnf.shape)) if Msnf[i, i] != 0]
        return d
    except Exception as e:
        return None

W("\n" + "=" * 78)
W("PART 3: ORDER-16 INVARIANTS — skew_16 vs sylvester_16")
W("=" * 78)
S16 = skew_tower(16)
H16 = sylvester(16)
J16 = np.ones((16, 16), dtype=np.int64)

for name, M in (("sylvester_16", H16), ("skew_16", S16)):
    prof = profile(M)
    proft = profile(M.T)
    W(f"\n{name}:")
    W(f"  row-quadruple profile (|sum| -> count, total C(16,4)=1820): {prof}")
    W(f"  col-quadruple profile (transpose):                          {proft}")
    snf_h = smith_normal_form_z(M)
    snf_h2 = snf_sympy(M)
    W(f"  SNF(H)        own: {snf_compact(snf_h)}   sympy: {snf_compact(snf_h2) if snf_h2 else 'n/a'}  agree={snf_h == snf_h2}")
    B = (M + J16) // 2
    snf_b = smith_normal_form_z(B)
    snf_b2 = snf_sympy(B)
    W(f"  SNF((H+J)/2)  own: {snf_compact(snf_b)}   sympy: {snf_compact(snf_b2) if snf_b2 else 'n/a'}  agree={snf_b == snf_b2}")
    W(f"  GF(2) rank of (H+J)/2: {gf2_rank(B % 2)}")

W("\nHall-set counts (quadruples with |sum|=16):")
p_s = profile(S16).get(16, 0)
p_h = profile(H16).get(16, 0)
W(f"  sylvester_16: {p_h}   skew_16: {p_s}")
W(f"VERDICT order 16: {'SAME profile' if profile(S16)==profile(H16) else 'DIFFERENT profile -> INEQUIVALENT Hadamard matrices'}")

# If invariants happen to agree, try witness search; else skip (provably inequivalent).
if profile(S16) == profile(H16):
    st, wit, nd, dt = monomial_equiv(S16, H16, node_cap=5_000_000, time_cap=300)
    W(f"  witness search at 16: {st} nodes={nd} time={dt:.1f}s")

W("\n" + "=" * 78)
W("PART 4: ORDER-32 FULL QUADRUPLE PROFILE (C(32,4)=35960 quadruples, exact, no sampling)")
W("=" * 78)
S32 = skew_tower(32)
H32 = sylvester(32)
t0 = time.time()
p32s = profile(S32)
p32h = profile(H32)
W(f"sylvester_32 profile: {p32h}")
W(f"skew_32      profile: {p32s}")
W(f"profile equal: {p32s == p32h}   (computed in {time.time()-t0:.1f}s)")
B32s = (S32 + np.ones((32, 32), dtype=np.int64)) // 2
B32h = (H32 + np.ones((32, 32), dtype=np.int64)) // 2
W(f"GF(2) rank of (H+J)/2: sylvester_32 = {gf2_rank(B32h % 2)}, skew_32 = {gf2_rank(B32s % 2)}")
# own pure-Python SNF suffers coefficient blowup at 32x32; sympy (instant) is used here.
# (own vs sympy cross-checked and agreeing at order 16 above.)
snf32s = snf_sympy(S32)
snf32h = snf_sympy(H32)
W(f"SNF(H) order 32 (sympy): sylvester: {snf_compact(snf32h)}")
W(f"                         skew:      {snf_compact(snf32s)}")
snf32bs = snf_sympy(B32s)
snf32bh = snf_sympy(B32h)
W(f"SNF((H+J)/2) order 32: sylvester: {snf_compact(snf32bh)}")
W(f"                       skew:      {snf_compact(snf32bs)}")

W("\n" + "=" * 78)
W("PART 5: WALSH-DOMAIN SHAPE OF S — N = W_k^T S (C = N / 2^k), k = 3, 4")
W("=" * 78)
for k in (3, 4, 5):
    m = 2 ** k
    Wm = sylvester(m)
    Sm = skew_tower(m)
    N = Wm.T @ Sm
    nnz = int(np.count_nonzero(N))
    vals, cts = np.unique(N, return_counts=True)
    hist = {int(v): int(c) for v, c in zip(vals, cts)}
    W(f"\nk={k} (m={m}): N = W^T S, C = N/{m}")
    W(f"  nonzeros: {nnz}/{m*m} ({100*nnz/(m*m):.1f}%)   entry histogram of N: {hist}")
    W(f"  diag(N): {list(np.diag(N))}")
    # triangularity checks in natural order
    lower = np.allclose(np.triu(N, 1), 0)
    upper = np.allclose(np.tril(N, -1), 0)
    W(f"  N lower-triangular: {lower}   upper-triangular: {upper}")
    # monomial?
    monomial = all(np.count_nonzero(N[i]) == 1 for i in range(m))
    W(f"  N monomial (1 nonzero/row): {monomial}")
    if k <= 4:
        W("  N matrix:")
        for row in N:
            W("   " + " ".join(f"{int(v):4d}" for v in row))
    # nonzero counts per row/col
    W(f"  nnz per row: {[int(np.count_nonzero(N[i])) for i in range(m)]}")
    W(f"  nnz per col: {[int(np.count_nonzero(N[:,j])) for j in range(m)]}")
    # recursive structure check: N_2m vs blocks of N_m?
    # also check conjugation C2 = W^T S W / m^2 ... print histogram of W^T S W / m
    N2 = Wm.T @ Sm @ Wm
    assert np.all(N2 % m == 0)
    N2 = N2 // m
    vals2, cts2 = np.unique(N2, return_counts=True)
    hist2 = {int(v): int(c) for v, c in zip(vals2, cts2)}
    nnz2 = int(np.count_nonzero(N2))
    W(f"  conjugation (W^T S W)/{m}: nnz {nnz2}/{m*m}, histogram {hist2}, "
      f"lower-tri={np.allclose(np.triu(N2,1),0)}, upper-tri={np.allclose(np.tril(N2,-1),0)}")

# Structural identity: S = W + Corr, Corr_{2m} = [[C_m, C_m],[C_m-2I, -C_m+2I]], Corr_1 = 0
W("\nStructural check: S_m = W_m + Corr_m with Corr_{2m}=[[C,C],[C-2I,2I-C]], C_1=0:")
def corr_tower(order):
    Cm = np.array([[0]], dtype=np.int64)
    while Cm.shape[0] < order:
        n = Cm.shape[0]
        I = np.eye(n, dtype=np.int64)
        Cm = np.block([[Cm, Cm], [Cm - 2 * I, 2 * I - Cm]])
    return Cm
for order in (2, 4, 8, 16, 32):
    ok = np.array_equal(skew_tower(order), sylvester(order) + corr_tower(order))
    W(f"  order {order}: S = W + Corr exact: {ok}")

W("\n" + "=" * 78)
W("PART 6: SKEW-WALSH BUTTERFLY — verify + time at order 1024 (and 4096)")
W("=" * 78)

def skew_rec(x):
    m = len(x)
    if m == 1:
        return x.copy()
    h = m // 2
    a = x[:h]; b = x[h:]
    u = a + b; v = a - b
    return np.concatenate([skew_rec(u), skew_rec(v) - 2 * v])

def skew_iter(x):
    y = x.copy()
    C = np.zeros_like(y)
    m = len(y)
    s = m
    while s >= 2:
        Y = y.reshape(-1, s)
        h = s // 2
        a = Y[:, :h].copy()
        b = Y[:, h:].copy()
        Y[:, :h] = a + b
        v = a - b
        Y[:, h:] = v
        C.reshape(-1, s)[:, h:] -= 2 * v
        s = h
    return y + C

# exact correctness on integer vectors
rng = np.random.default_rng(20260609)
for m in (2, 4, 8, 16, 64, 256, 1024):
    Sm = skew_tower(m)
    x = rng.integers(-9, 10, size=m).astype(np.int64)
    yd = Sm @ x
    yr = skew_rec(x)
    yi = skew_iter(x)
    ok = np.array_equal(yd, yr) and np.array_equal(yd, yi)
    if not ok:
        W(f"  m={m}: MISMATCH! dense={yd[:8]} rec={yr[:8]} iter={yi[:8]}")
    assert ok
W("  butterfly EXACT match vs dense matvec at m = 2,4,8,16,64,256,1024 (int64 vectors): True")

# block identity check: S_{2m} x == [S_m(x1+x2); S_m(x1-x2)-2(x1-x2)]
for m in (8, 64, 512):
    S2m = skew_tower(2 * m)
    Sm = skew_tower(m)
    x = rng.integers(-9, 10, size=2 * m).astype(np.int64)
    a, b = x[:m], x[m:]
    lhs = S2m @ x
    rhs = np.concatenate([Sm @ (a + b), Sm @ (a - b) - 2 * (a - b)])
    assert np.array_equal(lhs, rhs)
W("  block identity S_2m x = [S_m(x1+x2); S_m(x1-x2)-2(x1-x2)] EXACT at m=8,64,512: True")

def bench(fn, x, reps):
    best = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        fn(x)
        dt = time.perf_counter() - t0
        if dt < best:
            best = dt
    return best

for m in (1024, 4096):
    Sm = skew_tower(m).astype(np.float64)
    x = rng.standard_normal(m)
    tb = time.perf_counter()
    _ = skew_tower(m)
    tbuild = time.perf_counter() - tb
    t_dense = bench(lambda v: Sm @ v, x, 50)
    t_rec = bench(skew_rec, x, 50)
    t_iter = bench(skew_iter, x, 50)
    W(f"\n  m={m}: dense matvec (BLAS, prebuilt): {t_dense*1e6:9.1f} us | "
      f"recursive butterfly: {t_rec*1e6:9.1f} us | iterative butterfly: {t_iter*1e6:9.1f} us")
    W(f"         matrix build time: {tbuild*1e3:.1f} ms, memory: {m*m*8/1e6:.0f} MB (float64)"
      f" | speedup iter vs dense: {t_dense/t_iter:.2f}x, rec vs dense: {t_dense/t_rec:.2f}x")

# big orders, butterfly only (dense matrix would be 2GB at 16384)
for m in (16384, 2 ** 20):
    x = rng.standard_normal(m)
    t_rec = bench(skew_rec, x, 5)
    t_iter = bench(skew_iter, x, 5)
    W(f"  m={m}: recursive butterfly {t_rec*1e3:8.2f} ms | iterative butterfly {t_iter*1e3:8.2f} ms"
      f"  (dense matvec impossible/impractical: {m*m*8/1e9:.0f} GB matrix)")

# pure-python naive single-shot at 1024 for the algorithmic reference point
m = 1024
Sm = skew_tower(m)
Slist = [[int(v) for v in row] for row in Sm]
xl = [float(v) for v in rng.standard_normal(m)]
t0 = time.perf_counter()
yn = [sum(Slist[i][j] * xl[j] for j in range(m)) for i in range(m)]
t_py = time.perf_counter() - t0
W(f"\n  pure-Python naive O(m^2) matvec at m=1024 (reference): {t_py*1e3:.1f} ms")

W("\nDONE.")
T.f.close()
