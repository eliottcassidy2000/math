#!/usr/bin/env python3
"""S91 helper (g20_local): exact integer lattice tools (Python integers only).

  Lattice        incremental echelon Z-basis of the lattice spanned by inserted integer vectors
  smith_diagonal elementary divisors of an integer matrix
  rank_mod       rank of an integer matrix modulo a prime
  solve_q        one rational solution of a linear system (Fractions)
"""
from fractions import Fraction


def xgcd(a, b):
    """returns (g, s, t) with s a + t b = g = gcd(a, b) >= 0"""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


class Lattice:
    """Z-span of the inserted vectors, kept as an echelon basis {pivot column: row}."""

    def __init__(self, dim):
        self.dim = dim
        self.basis = {}

    def add(self, vec):
        """inserts vec; returns True iff the lattice grew."""
        v = [int(x) for x in vec]
        changed = False
        for j in range(self.dim):
            if v[j] == 0:
                continue
            b = self.basis.get(j)
            if b is None:
                if v[j] < 0:
                    v = [-x for x in v]
                self.basis[j] = v
                return True
            if v[j] % b[j] == 0:
                q = v[j] // b[j]
                v = [x - q * y for x, y in zip(v, b)]
            else:
                g, s, t = xgcd(b[j], v[j])
                nb = [s * x + t * y for x, y in zip(b, v)]
                bj, vj = b[j] // g, v[j] // g
                v = [bj * y - vj * x for x, y in zip(b, v)]
                self.basis[j] = nb
                changed = True
        return changed

    def contains(self, vec):
        v = [int(x) for x in vec]
        for j in range(self.dim):
            if v[j] == 0:
                continue
            b = self.basis.get(j)
            if b is None or v[j] % b[j] != 0:
                return False
            q = v[j] // b[j]
            v = [x - q * y for x, y in zip(v, b)]
        return True

    def coords(self, vec):
        """integer coordinates of vec in the echelon basis (rows ordered by pivot); raises if not a member."""
        v = [int(x) for x in vec]
        out = []
        for j in sorted(self.basis):
            b = self.basis[j]
            if v[j] % b[j] != 0:
                raise ValueError("vector not in lattice")
            q = v[j] // b[j]
            out.append(q)
            if q:
                v = [x - q * y for x, y in zip(v, b)]
        if any(v):
            raise ValueError("vector not in lattice")
        return out

    def rows(self):
        return [self.basis[j] for j in sorted(self.basis)]

    @property
    def rank(self):
        return len(self.basis)


def smith_diagonal(mat):
    """elementary divisors (nonzero diagonal of the Smith normal form) of an integer matrix."""
    A = [[int(x) for x in row] for row in mat]
    m = len(A)
    n = len(A[0]) if m else 0
    diag = []
    r = 0
    c = 0
    while r < m and c < n:
        # find a nonzero entry of smallest absolute value in the remaining block
        best = None
        for i in range(r, m):
            for j in range(c, n):
                if A[i][j] != 0 and (best is None or abs(A[i][j]) < abs(A[best[0]][best[1]])):
                    best = (i, j)
        if best is None:
            break
        i, j = best
        A[r], A[i] = A[i], A[r]
        for row in A:
            row[c], row[j] = row[j], row[c]
        while True:
            piv = A[r][c]
            dirty = False
            for i in range(r + 1, m):
                if A[i][c] != 0:
                    q = A[i][c] // piv
                    if q:
                        A[i] = [x - q * y for x, y in zip(A[i], A[r])]
                    if A[i][c] != 0:
                        A[r], A[i] = A[i], A[r]
                        dirty = True
                        break
            if dirty:
                continue
            for j in range(c + 1, n):
                if A[r][j] != 0:
                    q = A[r][j] // piv
                    if q:
                        for row in A:
                            row[j] -= q * row[c]
                    if A[r][j] != 0:
                        for row in A:
                            row[c], row[j] = row[j], row[c]
                        dirty = True
                        break
            if dirty:
                continue
            # the pivot must divide every remaining entry
            bad = None
            for i in range(r + 1, m):
                for j in range(c + 1, n):
                    if A[i][j] % piv != 0:
                        bad = i
                        break
                if bad is not None:
                    break
            if bad is None:
                break
            A[r] = [x + y for x, y in zip(A[r], A[bad])]
        diag.append(abs(A[r][c]))
        r += 1
        c += 1
    diag.sort()
    # normalise to a divisibility chain
    out = list(diag)
    changed = True
    while changed:
        changed = False
        for i in range(len(out) - 1):
            a, b = out[i], out[i + 1]
            if b % a != 0:
                g = xgcd(a, b)[0]
                out[i], out[i + 1] = g, a * b // g
                changed = True
        out.sort()
    return out


def rank_mod(mat, p):
    A = [[int(x) % p for x in row] for row in mat]
    m = len(A)
    n = len(A[0]) if m else 0
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        r += 1
        if r == m:
            break
    return r


def rank_q(mat):
    A = [[Fraction(x) for x in row] for row in mat]
    m = len(A)
    n = len(A[0]) if m else 0
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == m:
            break
    return r


def solve_q(mat, rhs):
    """one rational solution x of mat x = rhs (free variables set to 0), or None if inconsistent."""
    m = len(mat)
    n = len(mat[0]) if m else 0
    A = [[Fraction(x) for x in row] + [Fraction(b)] for row, b in zip(mat, rhs)]
    piv_cols = []
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv_cols.append(c)
        r += 1
        if r == m:
            break
    for i in range(r, m):
        if A[i][n] != 0:
            return None
    x = [Fraction(0)] * n
    for i, c in enumerate(piv_cols):
        x[c] = A[i][n]
    return x
