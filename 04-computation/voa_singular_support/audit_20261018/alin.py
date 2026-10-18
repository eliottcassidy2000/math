"""Audit helper (independent of the author's code): sparse vectors as dicts and
exact (fmpq_mat) / modular (nmod_mat) row reduction.

A 'vector' is a dict  key -> coefficient.  A 'basis' is a list of keys fixing
column order.  All exact arithmetic uses flint.fmpq.
"""
import flint

Q = flint.fmpq


def q(x, y=1):
    return Q(x, y)


def vadd(acc, vec, s=1):
    """acc += s*vec (in place)."""
    for k, c in vec.items():
        nv = acc.get(k, 0) + s * c
        if nv == 0:
            acc.pop(k, None)
        else:
            acc[k] = nv
    return acc


def vscale(vec, s):
    return {k: s * c for k, c in vec.items() if s * c != 0}


def to_mat(vectors, colindex, ncols, modulus=None):
    """Rows = vectors (dicts), columns via colindex: key -> column."""
    if modulus is None:
        M = flint.fmpq_mat(len(vectors), ncols)
        for i, v in enumerate(vectors):
            for k, c in v.items():
                M[i, colindex[k]] = c
        return M
    ent = [0] * (len(vectors) * ncols)
    for i, v in enumerate(vectors):
        for k, c in v.items():
            if isinstance(c, flint.fmpq):
                cc = (int(c.p) * pow(int(c.q), -1, modulus)) % modulus
            else:  # int or nmod
                cc = int(c) % modulus
            ent[i * ncols + colindex[k]] = cc
    return flint.nmod_mat(len(vectors), ncols, ent, modulus)


def rref_vectors(vectors, basis, modulus=None):
    """Row-reduce the vectors with columns in the order of `basis`.
    Returns list of (pivot_col, rowdict) for the nonzero rows (reduced echelon form)."""
    if not vectors:
        return []
    colindex = {k: i for i, k in enumerate(basis)}
    M = to_mat(vectors, colindex, len(basis), modulus)
    R, r = M.rref()
    out = []
    n = len(basis)
    ent = R.entries()
    for i in range(r):
        row = {}
        piv = None
        base = i * n
        for j in range(n):
            c = ent[base + j]
            if c != 0:
                if piv is None:
                    piv = j
                row[basis[j]] = c
        out.append((piv, row))
    return out


def rank(vectors, basis, modulus=None):
    if not vectors:
        return 0
    colindex = {k: i for i, k in enumerate(basis)}
    return to_mat(vectors, colindex, len(basis), modulus).rank()


def nullspace(vectors, basis, modulus=None):
    """Kernel of the linear map whose matrix has COLUMNS indexed by basis and rows = vectors,
    i.e. returns all x (dict over basis) with  sum_k row[k] x[k] = 0  for every row."""
    rows = rref_vectors(vectors, basis, modulus)
    pivcols = [p for p, _ in rows]
    pivset = set(pivcols)
    free = [j for j in range(len(basis)) if j not in pivset]
    one = 1 if modulus is not None else Q(1)
    ker = []
    for fj in free:
        x = {basis[fj]: one}
        for p, row in rows:
            c = row.get(basis[fj], 0)
            if c != 0:
                x[basis[p]] = -c if modulus is None else (-int(c)) % modulus
        ker.append(x)
    return ker


def in_span(vec, vectors, basis, modulus=None):
    r0 = rank(vectors, basis, modulus)
    r1 = rank(vectors + [vec], basis, modulus)
    return r1 == r0


def leading_forms(vectors, basis, degree, modulus=None):
    """Associated graded of span(vectors) for the filtration by `degree` (key -> int),
    where HIGHER degree = more leading.  Returns list of leading-form dicts (a basis of gr)."""
    order = sorted(basis, key=lambda k: -degree(k))
    rows = rref_vectors(vectors, order, modulus)
    out = []
    for p, row in rows:
        d0 = degree(order[p])
        out.append({k: c for k, c in row.items() if degree(k) == d0})
    return out
