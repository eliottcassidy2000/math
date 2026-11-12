#!/usr/bin/env python3
"""Exact integer lattice helpers (row echelon bases over Z) for the minimal representation of the two-apex count."""
from fractions import Fraction


def ext_gcd(a, b):
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0


class Lattice:
    """A sublattice of Z^n kept as an echelon Z-basis (pivot column -> row)."""

    def __init__(self, n):
        self.n = n
        self.rows = {}          # pivot column -> row (list of ints), zero before the pivot

    def rank(self):
        return len(self.rows)

    def add(self, v):
        """Add a generator; return True if the lattice changed."""
        v = list(v)
        changed = False
        for col in range(self.n):
            if v[col] == 0:
                continue
            if col not in self.rows:
                if v[col] < 0:
                    v = [-x for x in v]
                self.rows[col] = v
                return True
            r = self.rows[col]
            if v[col] % r[col] == 0:
                f = v[col] // r[col]
                v = [x - f * y for x, y in zip(v, r)]
                continue
            g, x, y = ext_gcd(r[col], v[col])
            nr = [x * p + y * q for p, q in zip(r, v)]
            nv = [(v[col] // g) * p - (r[col] // g) * q for p, q in zip(r, v)]
            self.rows[col] = nr
            v = nv
            changed = True
        return changed

    def basis(self):
        return [self.rows[c] for c in sorted(self.rows)]

    def coords(self, v):
        """Integer coordinates of v in basis(); raises if v is not in the lattice."""
        v = list(v)
        cols = sorted(self.rows)
        out = []
        for c in cols:
            r = self.rows[c]
            if v[c] % r[c] != 0:
                raise ValueError("not in lattice")
            f = v[c] // r[c]
            out.append(f)
            if f:
                v = [x - f * y for x, y in zip(v, r)]
        if any(v):
            raise ValueError("not in lattice (outside the span)")
        return out

    def contains(self, v):
        try:
            self.coords(v)
            return True
        except ValueError:
            return False


def det_int(M):
    """Exact determinant of an integer square matrix (fraction-free Bareiss)."""
    n = len(M)
    A = [list(r) for r in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            p = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if p is None:
                return 0
            A[k], A[p] = A[p], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def inv_frac(M):
    n = len(M)
    A = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def matmul(A, B):
    n, m, k = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(m) if A[i][t]) for j in range(k)] for i in range(n)]


def rank_frac(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    rank = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        p = next((r for r in range(rank, len(rows)) if rows[r][c] != 0), None)
        if p is None:
            continue
        rows[rank], rows[p] = rows[p], rows[rank]
        inv = 1 / rows[rank][c]
        rows[rank] = [x * inv for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][c] != 0:
                f = rows[r][c]
                rows[r] = [x - f * y for x, y in zip(rows[r], rows[rank])]
        rank += 1
    return rank


def nullspace_frac(M):
    """Basis of {x : M x = 0} for a rational matrix M (list of rows)."""
    rows = [[Fraction(x) for x in r] for r in M]
    ncol = len(rows[0])
    piv = []
    rank = 0
    for c in range(ncol):
        p = next((r for r in range(rank, len(rows)) if rows[r][c] != 0), None)
        if p is None:
            continue
        rows[rank], rows[p] = rows[p], rows[rank]
        inv = 1 / rows[rank][c]
        rows[rank] = [x * inv for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][c] != 0:
                f = rows[r][c]
                rows[r] = [x - f * y for x, y in zip(rows[r], rows[rank])]
        piv.append(c)
        rank += 1
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for fc in free:
        x = [Fraction(0)] * ncol
        x[fc] = Fraction(1)
        for i, pc in enumerate(piv):
            x[pc] = -rows[i][fc]
        basis.append(x)
    return basis
