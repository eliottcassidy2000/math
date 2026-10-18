#!/usr/bin/env python3
"""Independent computation of dim L_D (the simple quotient of a universal vacuum module), with no singular vectors,
no character formula and no contravariant form.

M = maximal proper submodule = { v : (U(n_+^) v) has zero vacuum component }.  With X ranging over a generating
set of n_+^ (positive modes, plus positive-root zero modes for affine/sl2 parts):
    v in M_(D,ch)  iff  X v in M_(D - w(X), ch + chg(X)) for all X   (D > 0).
We keep, for each block, a matrix Q (rows = coordinates on L = V/M).  Blocks are processed by increasing weight and,
within a weight, in an order where X v lands in an already-processed block (zero modes raise a height function).
"""
import flint

def simple_quotient_dims(alg, Dmax, nplus, height=None, only=None, verbose=False):
    """alg: ModeAlgebra; nplus: list of modes (w <= 0) generating n_+^; height(ch): zero modes X must raise it.
    Returns dict (D, ch) -> dim L.  'only' = optional set of blocks to restrict to (with their cone)."""
    P = alg.P
    basis = alg.pbw_basis(Dmax)
    idx = {key: {mo: i for i, mo in enumerate(ms)} for key, ms in basis.items()}
    Q = {}
    dims = {}
    zero = tuple(0 for _ in alg.chg[0])
    hfun = height if height is not None else (lambda ch: 0)
    keys = sorted(basis.keys(), key=lambda k: (k[0], -hfun(k[1])))
    for key in keys:
        D, ch = key
        if only is not None and key not in only:
            continue
        ms = basis[key]
        if D == 0:
            Q[key] = flint.nmod_mat(1, 1, [1], P) if ch == zero else flint.nmod_mat(0, len(ms), [], P)
            dims[key] = 1 if ch == zero else 0
            continue
        cols = []          # one column (list of rows) per basis vector, stacked over X
        blocks = []
        for X in nplus:
            w, t = X
            tgt = (D + w, tuple(a + b for a, b in zip(ch, alg.chg[t])))
            if tgt not in Q:
                continue
            Qt = Q[tgt]
            if Qt.nrows() == 0:
                continue
            blocks.append((X, tgt, Qt))
        rows_total = sum(Qt.nrows() for _, _, Qt in blocks)
        if rows_total == 0:
            Q[key] = flint.nmod_mat(0, len(ms), [], P)
            dims[key] = 0
            continue
        big = [[0] * len(ms) for _ in range(rows_total)]
        r0 = 0
        for X, tgt, Qt in blocks:
            it = idx[tgt]
            Ql = [[int(x) for x in row] for row in Qt.tolist()]
            for j, mo in enumerate(ms):
                img = alg.apply(X, mo)
                if not img:
                    continue
                for rr in range(len(Ql)):
                    row = Ql[rr]
                    s = 0
                    for mo2, co in img.items():
                        s += row[it[mo2]] * co
                    big[r0 + rr][j] = s % P
            r0 += len(Ql)
        M = flint.nmod_mat(rows_total, len(ms), [x for row in big for x in row], P)
        R, rk = M.rref()
        Rl = R.tolist()[:rk]
        Q[key] = flint.nmod_mat(rk, len(ms), [int(x) for row in Rl for x in row], P) if rk else flint.nmod_mat(0, len(ms), [], P)
        dims[key] = rk
        if verbose:
            print(key, len(ms), rk, flush=True)
    return dims

def per_weight(dims):
    out = {}
    for (D, ch), d in dims.items():
        out[D] = out.get(D, 0) + d
    return dict(sorted(out.items()))
