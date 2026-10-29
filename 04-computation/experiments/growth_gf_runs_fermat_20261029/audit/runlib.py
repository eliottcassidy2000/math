"""Audit library for THM-4846: run observers by subset dynamic programming (no permutation cutting),
Smith forms by FLINT (python-flint fmpz_mat.snf), independent of the producer's PARI pipeline.

hp[S][(s,t)] = number of Hamiltonian paths of T[S] from s to t, built by PREPENDING the start vertex:
    hp[S][(s,t)] = sum_{u in S-s, s->u} hp[S-s][(u,t)].
M_k[(s_1..s_k),(t_1..t_k)] = sum over ORDERED set partitions (S_1..S_k) of V of prod_i hp[S_i][(s_i,t_i)],
enumerated as surjections V -> {0..k-1}.
"""
import itertools
import subprocess
from collections import Counter
import flint


def tournaments(n):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True, check=True).stdout
    return [ln.strip() for ln in out.split() if ln.strip()]


def adj_from_bits(n, bits):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[idx] == '1':
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def circulant(n, S):
    return [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]


def hp_table(A):
    n = len(A)
    hp = {}
    for v in range(n):
        hp[1 << v] = {(v, v): 1}
    for size in range(2, n + 1):
        for comb in itertools.combinations(range(n), size):
            S = 0
            for v in comb:
                S |= 1 << v
            d = {}
            for s in comb:
                rest = S & ~(1 << s)
                for (u, t), c in hp[rest].items():
                    if A[s][u]:
                        d[(s, t)] = d.get((s, t), 0) + c
            hp[S] = d
    return hp


def ham_cycles(A, hp):
    n = len(A)
    full = (1 << n) - 1
    return sum(c for (s, t), c in hp[full].items() if s == 0 and A[t][0])


def ham_paths(hp, n):
    return sum(hp[(1 << n) - 1].values())


def is_strong(A):
    n = len(A)
    for src in range(n):
        seen = {src}
        stack = [src]
        while stack:
            x = stack.pop()
            for y in range(n):
                if A[x][y] and y not in seen:
                    seen.add(y)
                    stack.append(y)
        if len(seen) < n:
            return False
    return True


def ordered_tuples(n, k):
    """Ordered k-tuples of distinct vertices grouped by sorted base; within a base the S_k elements in lexicographic
    order of image tuples (for k = 3 this is 012,021,102,120,201,210, the order of THREE_RUN_GLUE (R2))."""
    perms = list(itertools.permutations(range(k)))
    out = []
    for b in itertools.combinations(range(n), k):
        for g in perms:
            out.append(tuple(b[g[i]] for i in range(k)))
    return out


def run_matrix(A, k, hp=None, tuples=None):
    n = len(A)
    if hp is None:
        hp = hp_table(A)
    if tuples is None:
        tuples = ordered_tuples(n, k)
    index = {t: i for i, t in enumerate(tuples)}
    N = len(tuples)
    M = [[0] * N for _ in range(N)]
    for f in itertools.product(range(k), repeat=n):
        masks = [0] * k
        for v, b in enumerate(f):
            masks[b] |= 1 << v
        if any(m == 0 for m in masks):
            continue
        lists = [list(hp[m].items()) for m in masks]
        for combo in itertools.product(*lists):
            s = tuple(c[0][0] for c in combo)
            t = tuple(c[0][1] for c in combo)
            val = 1
            for c in combo:
                val *= c[1]
            M[index[s]][index[t]] += val
    return tuples, M


def plus_sector(A, k, hp=None):
    """Unordered covers by k paths: start SET x end SET (k-subsets lexicographic)."""
    n = len(A)
    if hp is None:
        hp = hp_table(A)
    sets = list(itertools.combinations(range(n), k))
    idx = {s: i for i, s in enumerate(sets)}
    P = [[0] * len(sets) for _ in sets]
    # unordered set partitions: canonical ordering of blocks by their minimum vertex
    for f in itertools.product(range(k), repeat=n):
        masks = [0] * k
        for v, b in enumerate(f):
            masks[b] |= 1 << v
        if any(m == 0 for m in masks):
            continue
        mins = [(m & -m).bit_length() for m in masks]
        if mins != sorted(mins):
            continue
        lists = [list(hp[m].items()) for m in masks]
        for combo in itertools.product(*lists):
            s = tuple(sorted(c[0][0] for c in combo))
            t = tuple(sorted(c[0][1] for c in combo))
            val = 1
            for c in combo:
                val *= c[1]
            P[idx[s]][idx[t]] += val
    return sets, P


def snf_diag(M):
    """Invariant factors (with zeros) of an integer matrix via FLINT."""
    if not M:
        return []
    X = flint.fmpz_mat(M)
    D = X.snf()
    r = min(D.nrows(), D.ncols())
    return [abs(int(D[i, i])) for i in range(r)]


def snf_summary(d):
    c = Counter(d)
    keys = sorted(c, key=lambda x: (x == 0, x))
    return [(x, c[x]) for x in keys]


def fmt(summary):
    return ' '.join('%d^%d' % kv for kv in summary)


def rank_mod_p(M, p):
    """Rank over F_p (FLINT nmod_mat); rank_mod_p_py is the pure-Python cross-check."""
    if not M or not M[0]:
        return 0
    return flint.nmod_mat([[x % p for x in row] for row in M], p).rank()


def rank_mod_p_py(M, p):
    rows = [[x % p for x in row] for row in M]
    r = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c], p - 2, p)
        rows[r] = [x * inv % p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(x - f * y) % p for x, y in zip(rows[i], rows[r])]
        r += 1
    return r


def rank_Q(M):
    """Exact rational rank via FLINT."""
    if not M:
        return 0
    return flint.fmpz_mat(M).rank()


def matmul(X, Y):
    return (flint.fmpz_mat(X) * flint.fmpz_mat(Y)).tolist()


# THREE_RUN_GLUE (R2) chart and H = 6 F^-1 (R4)
F3 = [[1, 1, 1, 1, 1, 1],
      [1, -1, -1, 1, 1, -1],
      [1, 1, 0, -1, 0, -1],
      [0, 0, 1, -1, 1, -1],
      [0, -1, 1, 1, -1, 0],
      [1, -1, 0, 0, -1, 1]]
H3 = [[1, 1, 2, 0, 0, 2],
      [1, -1, 2, -2, 0, -2],
      [1, -1, 0, 2, 2, 0],
      [1, 1, 0, -2, 2, -2],
      [1, 1, -2, 2, -2, 0],
      [1, -1, -2, 0, -2, 2]]
F2 = [[1, 1], [1, -1]]
H2 = [[1, 1], [1, -1]]          # 2 F2^-1


def block_diag(B, r):
    s = len(B)
    out = [[0] * (s * r) for _ in range(s * r)]
    for o in range(r):
        for i in range(s):
            for j in range(s):
                out[o * s + i][o * s + j] = B[i][j]
    return out


def unglued(M, k):
    """N = F_big M F_big^{-1} for k = 2 (S_2 chart) or k = 3 (index-54 chart); checks integrality."""
    if k == 3:
        F, H, den, s = F3, H3, 6, 6
    else:
        F, H, den, s = F2, H2, 2, 2
    r = len(M) // s
    X = matmul(matmul(block_diag(F, r), M), block_diag(H, r))
    N = []
    for row in X:
        nr = []
        for v in row:
            v = int(v)
            if v % den:
                raise ValueError('not integral')
            nr.append(v // den)
        N.append(nr)
    return N


def sector(N, k, coords):
    s = 6 if k == 3 else 2
    r = len(N) // s
    idx = [o * s + c for o in range(r) for c in coords]
    return [[N[i][j] for j in idx] for i in idx]


def ndiv(summary, p):
    return sum(m for d, m in summary if d == 0 or d % p == 0)
