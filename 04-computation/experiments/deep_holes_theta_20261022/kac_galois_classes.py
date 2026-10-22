"""Galois action on the level-2 Kac classes of the shell-1 elements, exactly (session opus-2026-10-22-S53; THM-4760 (4')).

For a simply-laced component of Coxeter number h and m = h + 2, the shell-1 element has Kac class
    P = rho/m (Kac (3; 1, ..., 1), the level-2 vacuum) on non-carriers,
    C = (rho + theta)/m (the level-2 adjoint) on the carrier.
For each unit d mod m, dx is reduced into the fundamental alcove by the affine Weyl group W x Q (simple reflections
a_j -> a_j - a_i A_ij, and the affine reflection a_j -> a_j - (theta(x) - 1)(theta, alpha_j)), with exact fractions.
The resulting Kac coordinates (s_0; s_1, ..., s_n) are compared with those of P and C, and with their images under
the symmetries of the extended Dynkin diagram (which the glue realises up to the glue code). This upgrades the
root-spectrum check of root_multisets.py to conjugacy in W x Q (exact match) or W x P (match up to a diagram
symmetry). Usage: python kac_galois_classes.py"""
from fractions import Fraction as Fr
import math, itertools


def cartan(t, n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    if t == 'A':
        edges = [(i, i + 1) for i in range(n - 1)]
    elif t == 'D':
        edges = [(i, i + 1) for i in range(n - 2)] + [(n - 3, n - 1)]
    elif t == 'E':
        edges = [(0, 2), (2, 3), (1, 3)] + [(i, i + 1) for i in range(3, n - 1)]
    for a, b in edges:
        A[a][b] = A[b][a] = -1
    return A


def highest_root(A):
    n = len(A)
    roots = {tuple(1 if j == i else 0 for j in range(n)) for i in range(n)}
    frontier = list(roots)
    while frontier:
        new = []
        for r in frontier:
            for i in range(n):
                ip = sum(r[j] * A[j][i] for j in range(n))
                s = list(r); s[i] -= ip; s = tuple(s)
                if all(c >= 0 for c in s) and s not in roots:
                    roots.add(s); new.append(s)
        frontier = new
    return max(roots, key=sum)


def affine_cartan(A, theta):
    """Extended Cartan matrix, node 0 = -theta."""
    n = len(A)
    th_a = [sum(theta[j] * A[j][i] for j in range(n)) for i in range(n)]       # (theta, alpha_i)
    B = [[0] * (n + 1) for _ in range(n + 1)]
    B[0][0] = 2
    for i in range(n):
        B[0][i + 1] = B[i + 1][0] = -th_a[i]
        for j in range(n):
            B[i + 1][j + 1] = A[i][j]
    return B


def diagram_symmetries(B):
    """Permutations of the extended nodes preserving the extended Cartan matrix (brute force on small diagrams via
    refinement by degree)."""
    N = len(B)
    deg = [sum(1 for j in range(N) if j != i and B[i][j] != 0) for i in range(N)]
    syms = []

    def extend(perm, used):
        k = len(perm)
        if k == N:
            syms.append(tuple(perm)); return
        for c in range(N):
            if c in used or deg[c] != deg[k]:
                continue
            if all(B[k][j] == B[c][perm[j]] for j in range(k)):
                perm.append(c); used.add(c)
                extend(perm, used)
                perm.pop(); used.discard(c)

    extend([], set())
    return syms


def reduce_alcove(a, A, theta):
    n = len(A)
    th_a = [sum(theta[j] * A[j][i] for j in range(n)) for i in range(n)]
    a = list(a)
    for _ in range(100000):
        neg = [i for i in range(n) if a[i] < 0]
        if neg:
            i = neg[0]
            ai = a[i]
            a = [a[j] - ai * A[i][j] for j in range(n)]
            continue
        th = sum(theta[j] * a[j] for j in range(n))
        if th > 1:
            c = th - 1
            a = [a[j] - c * th_a[j] for j in range(n)]
            continue
        return a
    raise RuntimeError('no convergence')


def kac(a, theta, m):
    s0 = m * (1 - sum(t * x for t, x in zip(theta, a)))
    return tuple(int(v) for v in [s0] + [m * x for x in a])


for (t, n) in [('A', 1), ('D', 4), ('D', 6), ('D', 8), ('D', 10), ('D', 12), ('D', 16), ('D', 24), ('E', 6), ('E', 7), ('E', 8)]:
    A = cartan(t, n)
    theta = highest_root(A)
    h = sum(theta) + 1
    m = h + 2
    th_a = [sum(theta[j] * A[j][i] for j in range(n)) for i in range(n)]
    P = [Fr(1, m)] * n                                    # alpha_i(rho/m) = 1/m
    C = [Fr(1 + th_a[i], m) for i in range(n)]            # alpha_i((rho + theta)/m)
    kP, kC = kac(P, theta, m), kac(C, theta, m)
    B = affine_cartan(A, theta)
    syms = diagram_symmetries(B)
    orbit = lambda k: {tuple(k[p] for p in s) for s in syms}
    oP, oC = orbit(kP), orbit(kC)
    units = [d for d in range(1, m) if math.gcd(d, m) == 1]
    rowP, rowC = [], []
    for d in units:
        for base, row in ((P, rowP), (C, rowC)):
            kd = kac(reduce_alcove([d * x for x in base], A, theta), theta, m)
            row.append('P' if kd == kP else 'C' if kd == kC else 'p' if kd in oP else 'c' if kd in oC else '?')
    print(f'{t}{n:<2d} (h = {h}, m = {m}): P = {kP}, C = {kC}; |Omega-type symmetries| = {len(syms)}')
    print(f'      d:        {" ".join(f"{d:>2d}" for d in units)}')
    print(f'      d * P ->  {" ".join(f"{c:>2s}" for c in rowP)}')
    print(f'      d * C ->  {" ".join(f"{c:>2s}" for c in rowC)}')
print('legend: P, C exact (conjugate under W x Q); p, c up to an extended-diagram symmetry; ? another class')
