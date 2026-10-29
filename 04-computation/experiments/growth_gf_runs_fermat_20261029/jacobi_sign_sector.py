"""Sign sector of the two-run observer as a multilinear Jacobi complementary minor (check for every class n <= 6).

Work in R = Z[x_v]/(x_v^2), X = diag(x), B = I - A X (det B = 1 + nilpotent, a unit), P = X B^-1 (the path matrix of
THREE_RUN_GLUE (R1): P_st = sum over directed s-t paths of the product of their vertex variables).  Then
    S_2[alpha, beta] = M_2[(a1,a2),(b1,b2)] - M_2[(a1,a2),(b2,b1)] = [x_V] det P[alpha, beta]
                     = [x_V] x_a1 x_a2 (-1)^(a1+a2+b1+b2) det B[V-beta, V-alpha] / det B          (Jacobi),
and det B[.,.] expands over permutations whose non-fixed cycles are directed cycles of T (A has zero diagonal and
no 2-cycles), so the sign sector is a ratio of signed path-cycle-cover counts.  Plus sector likewise = [x_V] per P[alpha,beta].
"""
import itertools
from runs import tournaments, tour_from_bits
from fast_runs import m2_fast


def rmul(f, g):
    out = {}
    for m1, c1 in f.items():
        for m2, c2 in g.items():
            if m1 & m2:
                continue
            k = m1 | m2
            out[k] = out.get(k, 0) + c1 * c2
    return {k: v for k, v in out.items() if v}


def radd(f, g, s=1):
    out = dict(f)
    for k, v in g.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


def rdet(Mat):
    """Determinant over R by Laplace expansion along the first row (small sizes)."""
    k = len(Mat)
    if k == 0:
        return {0: 1}
    if k == 1:
        return Mat[0][0]
    tot = {}
    for j in range(k):
        if not Mat[0][j]:
            continue
        minor = [row[:j] + row[j + 1:] for row in Mat[1:]]
        term = rmul(Mat[0][j], rdet(minor))
        tot = radd(tot, term, -1 if j % 2 else 1)
    return tot


def rinv_unit(f, n):
    """(1 + N)^-1 for f = 1 + N with N nilpotent (no constant term besides 1)."""
    assert f.get(0, 0) == 1
    N = {k: v for k, v in f.items() if k}
    out = {0: 1}
    term = {0: 1}
    for _ in range(n + 1):
        term = rmul(term, {k: -v for k, v in N.items()})
        if not term:
            break
        out = radd(out, term)
    return out


checked = 0
for n in (3, 4, 5, 6):
    full = (1 << n) - 1
    for bits in tournaments(n):
        A = tour_from_bits(n, bits)
        T, _ = m2_fast(A)
        B = [[({0: 1} if i == j else {}) if not A[i][j] else radd({0: 1} if i == j else {}, {1 << j: 1}, -1) for j in range(n)] for i in range(n)]
        detB = rdet(B)
        inv = rinv_unit(detB, n)
        for al in itertools.combinations(range(n), 2):
            for be in itertools.combinations(range(n), 2):
                rows = [v for v in range(n) if v not in be]
                cols = [v for v in range(n) if v not in al]
                sub = [[B[r][c] for c in cols] for r in rows]
                val = rmul(rmul({(1 << al[0]) | (1 << al[1]): 1}, rdet(sub)), inv)
                sign = (-1) ** (al[0] + al[1] + be[0] + be[1])
                jac = sign * val.get(full, 0)
                direct = int(T[al[0], al[1], be[0], be[1]] - T[al[0], al[1], be[1], be[0]])
                assert jac == direct, (n, bits, al, be, jac, direct)
                checked += 1
    print('n=%d: Jacobi form of the sign sector verified on every class and every (alpha, beta)' % n)
print('entries checked:', checked)
