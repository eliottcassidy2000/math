#!/usr/bin/env python3
"""Lane twoapex (S93 Chan-Pak wave): core objects of the two-apex family G(2,0), re-derived.

Tournament T(w), w = l_1 ... l_n, l_i = a_i + 2 b_i in {0,1,2,3}:
  spine s_1 -> s_2 -> ... -> s_n transitive (s_i -> s_j for i < j), apexes u, v with u -> v,
  a_i = 1 iff u -> s_i (else s_i -> u), b_i = 1 iff v -> s_i (else s_i -> v).
A Hamiltonian path is [R1] x [R2] y [R3] with {x, y} = {u, v} and R1, R2, R3 increasing runs partitioning the spine:
  last(R1) -> x, x -> first(R2), last(R2) -> y, y -> first(R3); R2 empty needs x -> y, i.e. (x, y) = (u, v).
Arrangement 0: (x, y) = (u, v), bits (alpha, beta) = (a, b).  Arrangement 1: (x, y) = (v, u), (alpha, beta) = (b, a).

36-state transfer (the count of 3-colourings of the prefix by the status of each run):
  r1 in {0 empty, 1 last ok (alpha = 0), 2 last bad}, r2 in {0 empty, 1 last ok (beta = 0), 2 last bad},
  r3 in {0 empty, 1 nonempty}; a run can only be opened by a vertex satisfying its first-vertex condition.
24-state aggregation (what a suffix can see): c1 = [r1 == 2] (R1 needs a later vertex), r2, r3.
"""
import itertools
from fractions import Fraction


def brute(n, a, b):
    """Hamiltonian paths by subset dynamic programming on the n + 2 vertices."""
    V = n + 2
    u, v = n, n + 1
    adj = [[0] * V for _ in range(V)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = 1
        if a[i]:
            adj[u][i] = 1
        else:
            adj[i][u] = 1
        if b[i]:
            adj[v][i] = 1
        else:
            adj[i][v] = 1
    adj[u][v] = 1
    full = (1 << V) - 1
    dp = [[0] * V for _ in range(1 << V)]
    for x in range(V):
        dp[1 << x][x] = 1
    for S in range(1, 1 << V):
        row = dp[S]
        for x in range(V):
            c = row[x]
            if c:
                ax = adj[x]
                for y in range(V):
                    if ax[y] and not (S >> y) & 1:
                        dp[S | (1 << y)][y] += c
    return sum(dp[full])


def brute_perm(n, a, b):
    """Literal enumeration of vertex orders (independent of the subset DP); for n <= 5."""
    V = n + 2
    u, v = n, n + 1
    beats = [[False] * V for _ in range(V)]
    for i in range(n):
        for j in range(i + 1, n):
            beats[i][j] = True
        if a[i]:
            beats[u][i] = True
        else:
            beats[i][u] = True
        if b[i]:
            beats[v][i] = True
        else:
            beats[i][v] = True
    beats[u][v] = True
    tot = 0
    for p in itertools.permutations(range(V)):
        ok = True
        for k in range(V - 1):
            if not beats[p[k]][p[k + 1]]:
                ok = False
                break
        if ok:
            tot += 1
    return tot


# ---------------------------------------------------------------------------------------------------------------
# 36-state transfer, as explicit integer matrices (row-vector convention: state' = state @ M36[l]).
def idx36(arr, r1, r2, r3):
    return arr * 18 + r1 * 6 + r2 * 2 + r3


def build36():
    M = [[[0] * 36 for _ in range(36)] for _ in range(4)]
    for l in range(4):
        a, b = l & 1, (l >> 1) & 1
        for arr in range(2):
            al, be = (a, b) if arr == 0 else (b, a)
            for r1 in range(3):
                for r2 in range(3):
                    for r3 in range(2):
                        s = idx36(arr, r1, r2, r3)
                        # vertex joins R1
                        M[l][s][idx36(arr, 1 if al == 0 else 2, r2, r3)] += 1
                        # vertex joins R2 (opening it needs x -> vertex)
                        if r2 != 0 or al == 1:
                            M[l][s][idx36(arr, r1, 1 if be == 0 else 2, r3)] += 1
                        # vertex joins R3 (opening it needs y -> vertex)
                        if r3 != 0 or be == 1:
                            M[l][s][idx36(arr, r1, r2, 1)] += 1
    init = [0] * 36
    init[idx36(0, 0, 0, 0)] = 1
    init[idx36(1, 0, 0, 0)] = 1
    fin = [0] * 36
    for arr in range(2):
        for r1 in (0, 1):
            for r2 in (0, 1):
                if r2 == 0 and arr == 1:
                    continue
                for r3 in range(2):
                    fin[idx36(arr, r1, r2, r3)] = 1
    return M, init, fin


M36, INIT36, FIN36 = build36()


def vecmat(v, M):
    n = len(M[0])
    out = [0] * n
    for i, c in enumerate(v):
        if c:
            Mi = M[i]
            for j in range(n):
                if Mi[j]:
                    out[j] += c * Mi[j]
    return out


def matvec(M, v):
    return [sum(Mi[j] * v[j] for j in range(len(v)) if Mi[j]) for Mi in M]


def dot(x, y):
    return sum(p * q for p, q in zip(x, y))


def H36(word):
    v = INIT36
    for l in word:
        v = vecmat(v, M36[l])
    return dot(v, FIN36)


# ---------------------------------------------------------------------------------------------------------------
# Direct formula as a sum over the markers (no transfer): used as a third independent count.
def H_direct(word):
    """sum over arrangements and 3-colourings, literally (3^n colourings); for n <= 9."""
    n = len(word)
    a = [l & 1 for l in word]
    b = [(l >> 1) & 1 for l in word]
    tot = 0
    for arr in range(2):
        al, be = (a, b) if arr == 0 else (b, a)
        for col in itertools.product((0, 1, 2), repeat=n):
            R = [[i for i in range(n) if col[i] == k] for k in range(3)]
            if R[0] and al[R[0][-1]] != 0:
                continue
            if R[1]:
                if al[R[1][0]] != 1 or be[R[1][-1]] != 0:
                    continue
            elif arr == 1:
                continue
            if R[2] and be[R[2][0]] != 1:
                continue
            tot += 1
    return tot


# ---------------------------------------------------------------------------------------------------------------
# exact linear algebra over Q
def rank_basis(rows):
    """Row-reduce a list of integer/Fraction vectors; return (rank, list of pivot row indices in input order)."""
    basis = []          # list of (pivot_col, reduced_row)
    picked = []
    for idx, r in enumerate(rows):
        r = [Fraction(x) for x in r]
        for pc, br in basis:
            if r[pc] != 0:
                f = r[pc]
                r = [x - f * y for x, y in zip(r, br)]
        pc = next((j for j, x in enumerate(r) if x != 0), None)
        if pc is None:
            continue
        inv = 1 / r[pc]
        r = [x * inv for x in r]
        basis.append((pc, r))
        picked.append(idx)
    return len(basis), picked


def solve_in_basis(basis_rows, target):
    """Express target as a rational combination of the linearly independent basis_rows (exact)."""
    k = len(basis_rows)
    n = len(target)
    # solve sum_i c_i basis_rows[i] = target : (n equations, k unknowns)
    A = [[Fraction(basis_rows[i][j]) for i in range(k)] + [Fraction(target[j])] for j in range(n)]
    row = 0
    piv = []
    for col in range(k):
        p = next((r for r in range(row, n) if A[r][col] != 0), None)
        if p is None:
            raise ValueError("basis rows dependent")
        A[row], A[p] = A[p], A[row]
        inv = 1 / A[row][col]
        A[row] = [x * inv for x in A[row]]
        for r in range(n):
            if r != row and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[row])]
        piv.append(col)
        row += 1
    for r in range(row, n):
        if A[r][k] != 0:
            raise ValueError("target not in span")
    return [A[i][k] for i in range(k)]


def words(n):
    return itertools.product(range(4), repeat=n)


def letters_ab(word):
    return [l & 1 for l in word], [(l >> 1) & 1 for l in word]
