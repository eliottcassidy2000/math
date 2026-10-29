"""Independent audit library (audit2, 2026-10-29).  Written from scratch; does NOT import the producer's runs.py.

Path counts: H[s, t, U] = number of directed Hamiltonian paths of T[U] from s to t (vertex set exactly U).
Computed by a PREPEND-vertex recursion organised by path length (the producer appends at the end, by subset):
    H_1[s, s, {s}] = 1,   H_l[s, t, U] = sum_r A[s, r] H_{l-1}[r, t, U - {s}]   (s in U).
Run observers: M_k[(s_1..s_k), (t_1..t_k)] = sum over ordered set partitions (U_1..U_k) of V of prod H[s_i, t_i, U_i].
"""
import itertools
import subprocess
import numpy as np
import flint


def gentourng(n, extra=()):
    out = subprocess.run(['gentourng', '-q'] + list(extra) + [str(n)], capture_output=True, text=True,
                         check=True).stdout
    return [ln.strip() for ln in out.split('\n') if ln.strip()]


def adj_from_ascii(n, s):
    """gentourng default output: upper triangle row by row; '1' at (i,j), i<j, read as the arc i -> j.
    (The opposite reading gives the converse tournament; every statistic audited here is converse-invariant:
    M_k(T^op) = M_k(T)^T.)"""
    A = np.zeros((n, n), dtype=np.int64)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1':
                A[i, j] = 1
            else:
                A[j, i] = 1
            k += 1
    assert k == len(s)
    return A


def circulant(n, S):
    return np.array([[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)], dtype=np.int64)


def path_table(A):
    """H[s, t, U] (int64 numpy array of shape (n, n, 2^n))."""
    n = A.shape[0]
    N = 1 << n
    H = np.zeros((n, n, N), dtype=np.int64)
    cur = np.zeros((n, n, N), dtype=np.int64)
    for s in range(n):
        cur[s, s, 1 << s] = 1
    H += cur
    without = [np.array([U for U in range(N) if not (U >> s) & 1], dtype=np.int64) for s in range(n)]
    for _ in range(2, n + 1):
        Q = np.tensordot(A, cur, axes=([1], [0]))      # Q[s, t, U'] = sum_r A[s, r] cur[r, t, U']
        nxt = np.zeros_like(cur)
        for s in range(n):
            idx = without[s]
            row = nxt[s]                     # view, shape (n, N)
            row[:, idx | (1 << s)] = Q[s][:, idx]
        cur = nxt
        H += cur
    return H


def ham_cycles(A, H):
    n = A.shape[0]
    full = (1 << n) - 1
    return int(sum(H[0, t, full] for t in range(1, n) if A[t, 0]))


def is_strong(A):
    n = A.shape[0]
    def reach(M):
        seen = {0}
        stack = [0]
        while stack:
            x = stack.pop()
            for y in range(n):
                if M[x, y] and y not in seen:
                    seen.add(y)
                    stack.append(y)
        return len(seen) == n
    return reach(A) and reach(A.T)


def pairs_ordered(n):
    """Ordered pairs grouped by base 2-set: (b0,b1), (b1,b0)."""
    out = []
    for b in itertools.combinations(range(n), 2):
        out.append((b[0], b[1]))
        out.append((b[1], b[0]))
    return out


def m2_matrix(A, H=None):
    n = A.shape[0]
    if H is None:
        H = path_table(A)
    N = 1 << n
    full = N - 1
    Hf = H.reshape(n * n, N)
    Hc = H[:, :, full ^ np.arange(N)].reshape(n * n, N)
    T = Hf @ Hc.T                      # T[(a1,b1),(a2,b2)] = sum_U H[a1,b1,U] H[a2,b2,V-U]
    pr = pairs_ordered(n)
    M = np.array([[T[a[0] * n + b[0], a[1] * n + b[1]] for b in pr] for a in pr], dtype=np.int64)
    return M, pr


def m2_sub_table(A, H):
    """M2W[W][a2, b2, a3, b3] = sum_{U subset W} H[a2,b2,U] H[a3,b3,W-U] for every subset W."""
    n = A.shape[0]
    N = 1 << n
    out = np.zeros((N, n, n, n, n), dtype=np.int64)
    for W in range(N):
        subs = []
        U = W
        while True:
            subs.append(U)
            if U == 0:
                break
            U = (U - 1) & W
        subs = np.array(subs, dtype=np.int64)
        X = H[:, :, subs]            # (n, n, m)
        Y = H[:, :, W ^ subs]
        out[W] = np.einsum('abm,cdm->abcd', X, Y)
    return out


S3 = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]   # image tuples, lexicographic


def triples_by_orbit(n):
    bases = list(itertools.combinations(range(n), 3))
    tup = [tuple(b[g[k]] for k in range(3)) for b in bases for g in S3]
    return bases, tup


def m3_matrix(A, H=None):
    n = A.shape[0]
    if H is None:
        H = path_table(A)
    N = 1 << n
    full = N - 1
    M2W = m2_sub_table(A, H)                      # (N, n, n, n, n): [W, a2, b2, a3, b3]
    Hf = H.reshape(n * n, N)                      # [(a1,b1), U1]
    R = M2W[full ^ np.arange(N)].reshape(N, n ** 4)   # [U1, (a2,b2,a3,b3)]
    T = (Hf @ R).reshape(n, n, n, n, n, n)        # [a1,b1,a2,b2,a3,b3]
    bases, tup = triples_by_orbit(n)
    idx = np.array(tup, dtype=np.int64)
    a1, a2, a3 = idx[:, 0][:, None], idx[:, 1][:, None], idx[:, 2][:, None]
    b1, b2, b3 = idx[:, 0][None, :], idx[:, 1][None, :], idx[:, 2][None, :]
    M = T[a1, b1, a2, b2, a3, b3]
    return M.astype(np.int64), bases, tup


def brute_mk(A, k):
    """Brute force: enumerate ordered set partitions into k nonempty blocks and Hamiltonian paths inside blocks."""
    n = A.shape[0]
    from collections import defaultdict
    cnt = defaultdict(int)
    def paths(block):
        res = []
        for p in itertools.permutations(block):
            if all(A[p[i], p[i + 1]] for i in range(len(p) - 1)):
                res.append((p[0], p[-1]))
        return res
    for lab in itertools.product(range(k), repeat=n):
        blocks = [[v for v in range(n) if lab[v] == i] for i in range(k)]
        if any(not b for b in blocks):
            continue
        plist = [paths(b) for b in blocks]
        for combo in itertools.product(*plist):
            s = tuple(c[0] for c in combo)
            t = tuple(c[1] for c in combo)
            cnt[(s, t)] += 1
    return cnt


def rank_mod_p(M, p):
    M = np.array(M, dtype=np.int64) % p
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + nz[0]
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        col = M[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % p
        r += 1
    return r


def nullspace_mod_p(M, p):
    """Basis (rows) of the right kernel of M over F_p."""
    M = np.array(M, dtype=np.int64) % p
    rows, cols = M.shape
    r = 0
    pivcols = []
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + nz[0]
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        col = M[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % p
        pivcols.append(c)
        r += 1
    free = [c for c in range(cols) if c not in pivcols]
    basis = []
    for f in free:
        v = np.zeros(cols, dtype=np.int64)
        v[f] = 1
        for i, pc in enumerate(pivcols):
            v[pc] = (-M[i, f]) % p
        basis.append(v)
    return np.array(basis, dtype=np.int64).reshape(len(basis), cols)


def inv_mod_p(M, p):
    M = np.array(M, dtype=np.int64) % p
    n = M.shape[0]
    aug = np.concatenate([M, np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        nz = np.nonzero(aug[c:, c])[0]
        if nz.size == 0:
            raise ValueError('singular mod p')
        piv = c + nz[0]
        if piv != c:
            aug[[c, piv]] = aug[[piv, c]]
        inv = pow(int(aug[c, c]), p - 2, p)
        aug[c] = (aug[c] * inv) % p
        col = aug[:, c].copy()
        col[c] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            aug[nzr] = (aug[nzr] - np.outer(col[nzr], aug[c])) % p
    return aug[:, n:]


def snf_diag(M):
    """Smith normal form diagonal via FLINT (independent of PARI); returns list of nonnegative ints."""
    M = np.asarray(M)
    fm = flint.fmpz_mat([[int(x) for x in row] for row in M])
    S = fm.snf()
    k = min(S.nrows(), S.ncols())
    d = [abs(int(S[i, i])) for i in range(k)]
    # pad: for a square matrix k = n; zeros included
    return d


def snf_summary(d):
    out = {}
    for x in d:
        out[x] = out.get(x, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: (kv[0] == 0, kv[0])))


def fmt_snf(summ):
    return ' '.join('%d^%d' % (k, v) for k, v in summ.items())


def det_exact(M):
    fm = flint.fmpz_mat([[int(x) for x in row] for row in np.asarray(M)])
    return int(fm.det())
