"""Ordered k-run response matrices of tournaments, their S_3 sectors, and Smith forms.

M_k[s,t] counts ordered spanning covers of T by k vertex-disjoint directed paths, path i running from
s_i to t_i (singletons allowed).  This is the observer of
05-knowledge/results/overnight_rank_return_20261029/THREE_RUN.md (T1) and
05-knowledge/results/overnight_research_20261029/THREE_RUN_GLUE.md (R1), recomputed here independently
by cutting permutations.

For k = 3 the ordered triples are arranged as (sorted base b, g in S_3) -> (b[g0], b[g1], b[g2]), with
S_3 listed by image tuples 012, 021, 102, 120, 201, 210, so that permuting run positions by h sends
(b, g) to (b, g o h).  The Fourier chart F (index 54) of THREE_RUN_GLUE (R2) is applied orbitwise.

Smith forms and ranks are computed by PARI/gp (matsnf, matrank, matker).
"""
import itertools
import subprocess
import tempfile
import os

S3 = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
F3 = [[1, 1, 1, 1, 1, 1],
      [1, -1, -1, 1, 1, -1],
      [1, 1, 0, -1, 0, -1],
      [0, 0, 1, -1, 1, -1],
      [0, -1, 1, 1, -1, 0],
      [1, -1, 0, 0, -1, 1]]
H3 = [[1, 1, 2, 0, 0, 2],          # H = 6 F^-1 (THREE_RUN_GLUE (R4))
      [1, -1, 2, -2, 0, -2],
      [1, -1, 0, 2, 2, 0],
      [1, 1, 0, -2, 2, -2],
      [1, 1, -2, 2, -2, 0],
      [1, -1, -2, 0, -2, 2]]


def tour_from_bits(n, bits):
    """gentourng upper-triangle bit string -> adjacency matrix (1 means i -> j for i < j)."""
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


def tournaments(n):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True, check=True).stdout
    return [line.strip() for line in out.split() if line.strip()]


def circulant(n, S):
    return [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]


def ordered_tuples(n, k):
    """Ordered k-tuples of distinct vertices; for k = 3 grouped in S_3 orbits as (base, g)."""
    if k == 3:
        out = []
        for b in itertools.combinations(range(n), 3):
            for g in S3:
                out.append((b[g[0]], b[g[1]], b[g[2]]))
        return out
    if k == 2:
        out = []
        for b in itertools.combinations(range(n), 2):
            out.append((b[0], b[1]))
            out.append((b[1], b[0]))
        return out
    return list(itertools.permutations(range(n), k))


def run_matrix(A, k, tuples=None):
    n = len(A)
    if tuples is None:
        tuples = ordered_tuples(n, k)
    index = {t: i for i, t in enumerate(tuples)}
    N = len(tuples)
    M = [[0] * N for _ in range(N)]
    for w in itertools.permutations(range(n)):
        breaks = [i for i in range(n - 1) if not A[w[i]][w[i + 1]]]
        if len(breaks) > k - 1:
            continue
        free = [i for i in range(n - 1) if A[w[i]][w[i + 1]]]
        for extra in itertools.combinations(free, k - 1 - len(breaks)):
            cuts = sorted(breaks + list(extra))
            starts = tuple([w[0]] + [w[c + 1] for c in cuts])
            ends = tuple([w[c] for c in cuts] + [w[n - 1]])
            M[index[starts]][index[ends]] += 1
    return tuples, M


def matmul(X, Y):
    m, p, q = len(X), len(Y), len(Y[0])
    out = [[0] * q for _ in range(m)]
    for i in range(m):
        Xi = X[i]
        oi = out[i]
        for t in range(p):
            a = Xi[t]
            if a:
                Yt = Y[t]
                for j in range(q):
                    if Yt[j]:
                        oi[j] += a * Yt[j]
    return out


def block_diag(B, r):
    s = len(B)
    out = [[0] * (s * r) for _ in range(s * r)]
    for o in range(r):
        for i in range(s):
            for j in range(s):
                out[o * s + i][o * s + j] = B[i][j]
    return out


def unglued_k3(M, norbits):
    """N = F_big M F_big^-1 (must be integral); returned with the sector blocks P, S, Lrow0, Lrow1."""
    Fb = block_diag(F3, norbits)
    Hb = block_diag(H3, norbits)
    X = matmul(matmul(Fb, M), Hb)
    N = []
    for row in X:
        nr = []
        for v in row:
            if v % 6:
                raise ValueError('F M F^-1 not integral')
            nr.append(v // 6)
        N.append(nr)
    # coordinates within an orbit: 0 = u (plus), 1 = v (minus), W = [[2,3],[4,5]] row-major.
    # M acts by LEFT multiplication of each W by L_ab = sum_x K_ab(x) rho(x^-1), so the standard
    # sector splits into the two columns (2,4) and (3,5) of W.
    def pick(coords):
        idx = [o * 6 + c for o in range(norbits) for c in coords]
        return [[N[i][j] for j in idx] for i in idx]
    P = pick([0])
    S = pick([1])
    L0 = pick([2, 4])
    L1 = pick([3, 5])
    # sanity: N is block diagonal for this split
    sets = [[0], [1], [2, 4], [3, 5]]
    for a in range(4):
        for b in range(4):
            if a == b:
                continue
            for o in range(norbits):
                for o2 in range(norbits):
                    for c in sets[a]:
                        for c2 in sets[b]:
                            if N[o * 6 + c][o2 * 6 + c2]:
                                raise ValueError('sector split not block diagonal')
    return N, P, S, L0, L1


def gp_matrix(M):
    # Mat() keeps a single row a t_MAT (a bare one-row literal is a t_VEC, which matrank rejects)
    return 'Mat([' + ';'.join(','.join(str(x) for x in row) for row in M) + '])'


def gp_batch(commands, timeout=3600):
    """Run gp on a list of expressions, one print per expression; returns list of output lines."""
    with tempfile.NamedTemporaryFile('w', suffix='.gp', delete=False) as f:
        # parisizemax on its own line: a default(parisize) on the same line discards the rest of it
        f.write('default(parisizemax, 4000000000);\n')
        for c in commands:
            f.write('print(' + c + ');\n')
        f.write('quit;\n')
        path = f.name
    try:
        res = subprocess.run(['gp', '-q', '-s', '800000000', path], capture_output=True, text=True,
                             timeout=timeout, check=True)
        err = [ln for ln in res.stderr.split('\n') if ln.strip() and 'stack size' not in ln]
        if err:
            raise RuntimeError('gp error: ' + ' | '.join(err)[:500])
        out = res.stdout
    finally:
        os.unlink(path)
    lines = [ln for ln in out.split('\n') if ln.strip()]
    if len(lines) != len(commands):
        raise RuntimeError('gp output mismatch: %d lines for %d commands' % (len(lines), len(commands)))
    return lines


def parse_vec(s):
    s = s.strip().strip('[]')
    return [int(x) for x in s.split(',')] if s else []


def snf_summary(d):
    """Elementary divisors -> dict value:multiplicity (0 = kernel)."""
    out = {}
    for x in d:
        x = abs(x)
        out[x] = out.get(x, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: (kv[0] == 0, kv[0])))


def fmt_snf(summary):
    return ' '.join('%d^%d' % (k, v) for k, v in summary.items())
