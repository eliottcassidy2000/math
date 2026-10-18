#!/usr/bin/env python3
"""Independent exact audit of THM-4671 (does not import the author's modules).

For each (p,q):
  * singular vector v of V^c at weight 2N from L_1 v = L_2 v = 0 (exact, fmpq), check the
    solution space is 1-dimensional and the L_{-2}^N coefficient is nonzero;
  * Lemma 1: alpha, beta (coefficients of L_{-4}L_{-2}^{N-2}, L_{-3}^2 L_{-2}^{N-3});
  * Lemma 2 / part 1: u = L_{-2}L_{-1}v - N L_{-3}v, its top length, and its length-N symbol
    against A b_N + B d^3 f;
  * maximal submodule M_D of V^c computed from scratch (w in M_D iff L_1 w in M_{D-1} and
    L_2 w in M_{D-2}), compared with I_D = U(Vir_-) v  (no character formula used);
  * grI_D (leading forms for the PBW-length filtration), arc ideal J_D = (x_0^N)_d,
    dim ker(pi)_D = dim grI_D - dim J_D for D <= 2N+3 (+ extra weights), and b_N in grI \\ J.
Usage: python vir_audit.py [extra_weights]
"""
import sys
import time
import flint
from alin import rref_vectors, rank, nullspace, leading_forms, in_span
from virmod import VirVac, partitions_min2, central_charge, poly_mul, poly_d, symbol

Q = flint.fmpq


def pnum(k):
    # number of partitions of k
    from sympy import partition
    return int(partition(k))


def op_matrix_rows(vir, n, basis_src, basis_tgt):
    """Return fmpq_mat A with A[i,j] = coeff of basis_tgt[i] in L_n basis_src[j]."""
    idx = {m: i for i, m in enumerate(basis_tgt)}
    A = flint.fmpq_mat(len(basis_tgt), len(basis_src))
    for j, b in enumerate(basis_src):
        for m, c in vir.act(n, (b)).items():
            A[idx[m], j] = c
    return A


def maximal_submodule_functionals(vir, Dmax):
    """K[D] = fmpq_mat whose rows span the annihilator of M_D (dim L_D = rank)."""
    bases = {D: partitions_min2(D) for D in range(Dmax + 1)}
    K = {0: flint.fmpq_mat(1, 1, [1])}
    for D in range(1, Dmax + 1):
        B = bases[D]
        if not B:
            K[D] = flint.fmpq_mat(0, 0)
            continue
        blocks = []
        for n in (1, 2):
            if D - n < 0:
                continue
            Bt = bases[D - n]
            if not Bt or K[D - n].nrows() == 0:
                continue
            A = op_matrix_rows(vir, n, B, Bt)
            blocks.append(K[D - n] * A)
        rows = []
        for Mb in blocks:
            ent = Mb.entries()
            nc = Mb.ncols()
            for i in range(Mb.nrows()):
                rows.append({B[j]: ent[i * nc + j] for j in range(nc) if ent[i * nc + j] != 0})
        rr = rref_vectors(rows, B)
        Km = flint.fmpq_mat(len(rr), len(B))
        col = {m: j for j, m in enumerate(B)}
        for i, (_, row) in enumerate(rr):
            for m, c in row.items():
                Km[i, col[m]] = c
        K[D] = Km
    return K, bases


def apply_functionals(K, basis, vec):
    col = {m: j for j, m in enumerate(basis)}
    x = flint.fmpq_mat(len(basis), 1)
    for m, c in vec.items():
        x[col[m], 0] = c
    y = K * x
    return [y[i, 0] for i in range(y.nrows())]


def audit_model(p, q_, extra=0):
    t0 = time.time()
    N = (p - 1) * (q_ - 1) // 2
    c = central_charge(p, q_)
    vir = VirVac(c, Q(1))
    W = 2 * N
    B2N = partitions_min2(W)
    # --- singular vector -------------------------------------------------------------
    rows = []
    for n in (1, 2):
        # equations: for each target monomial m, sum_j (L_n b_j)[m] x_j = 0
        eq = {}
        for b in B2N:
            for m, cc in vir.act(n, b).items():
                eq.setdefault(m, {})[b] = cc
        rows.extend(eq.values())
    ker = nullspace(rows, B2N)
    res = {'p': p, 'q': q_, 'N': N, 'c': c, 'dim_sing_space': len(ker)}
    assert len(ker) == 1, (p, q_, len(ker))
    v = ker[0]
    top = (2,) * N
    lead = v.get(top, 0)
    res['L2^N coeff nonzero'] = (lead != 0)
    v = {m: cc / lead for m, cc in v.items()}
    # verify singular: L_1 v = L_2 v = 0
    assert not vir.actv(1, v) and not vir.actv(2, v)
    # --- Lemma 1 -------------------------------------------------------------------
    alpha = v.get((4,) + (2,) * (N - 2), Q(0)) if N >= 2 else None
    beta = v.get((3, 3) + (2,) * (N - 3), Q(0)) if N >= 3 else None
    alpha_pred = -Q(N) * (8 * (N - 1) + c) / 12
    res['alpha'] = alpha
    res['alpha_ok'] = (alpha == alpha_pred)
    if beta is not None:
        res['beta'] = beta
        res['5a+8b_ok'] = (5 * alpha + 8 * beta == -Q(3 * N * (N - 1), 2))
    # length profile of v
    res['v_max_len'] = max(len(m) for m in v)
    # --- u and its symbol --------------------------------------------------------------
    u = vir.actv(-2, vir.actv(-1, v))
    for m, cc in vir.actv(-3, v).items():
        nv = u.get(m, 0) - N * cc
        if nv == 0:
            u.pop(m, None)
        else:
            u[m] = nv
    res['u_max_len'] = max(len(m) for m in u) if u else None
    sN = symbol(u, N)
    f = {top: Q(1)}
    d3f = poly_d(poly_d(poly_d(f)))
    A = Q(N * (2 * N + 1), 8) * (c - (1 - Q(3 * (2 * N - 1) ** 2, 2 * N + 1)))
    B = -(6 * N + c - 10) / 24
    A_pq = Q(3 * (p - 1) * (p - 2) * (q_ - 1) * (q_ - 2) * (p * p * q_ * q_ - (p - q_) ** 2), 16 * p * q_)
    res['A'] = A
    res['A_pq_matches'] = (A == A_pq)
    res['B'] = B
    if N >= 3:
        bN = {(4, 3) + (2,) * (N - 2): Q(1), (3, 3, 3) + (2,) * (N - 3): Q(2 * N - 7, 12)}
    else:
        bN = None
    pred = {}
    if bN is not None:
        for m, cc in bN.items():
            pred[m] = pred.get(m, 0) + A * cc
    for m, cc in d3f.items():
        pred[m] = pred.get(m, 0) + B * cc
    pred = {m: cc for m, cc in pred.items() if cc != 0}
    res['sigma_N(u) == A b_N + B d3f'] = (sN == pred)
    res['sigma_N(u)'] = sN
    # p = 2: is u - B L_{-1}^3 v zero in V^c ?
    L13v = vir.actv(-1, vir.actv(-1, vir.actv(-1, v)))
    diff = dict(u)
    for m, cc in L13v.items():
        nv = diff.get(m, 0) - B * cc
        if nv == 0:
            diff.pop(m, None)
        else:
            diff[m] = nv
    res['u - B L_{-1}^3 v == 0'] = (not diff)
    # --- maximal submodule vs I, grI, J -------------------------------------------------
    Dmax = W + 3 + extra
    K, bases = maximal_submodule_functionals(vir, Dmax)
    I = {W: [v]}
    for D in range(W + 1, Dmax + 1):
        gens = [vir.actv(-1, w) for w in I[D - 1]]
        if D - 2 >= W:
            gens += [vir.actv(-2, w) for w in I[D - 2]]
        rr = rref_vectors(gens, bases[D])
        I[D] = [row for _, row in rr]
    table = []
    for D in range(W, Dmax + 1):
        Bd = bases[D]
        dimV = len(Bd)
        dimL = K[D].rank() if K[D].nrows() else 0
        dimI = len(I[D])
        # I inside M: all functionals vanish on I
        inM = all(all(x == 0 for x in apply_functionals(K[D], Bd, w)) for w in I[D])
        gr = leading_forms(I[D], Bd, degree=len)
        # J_D: span of m * d^s f
        Jg = []
        for s in range(0, D - W + 1):
            ds = dict(f)
            for _ in range(s):
                ds = poly_d(ds)
            for mm in partitions_min2(D - W - s):
                Jg.append(poly_mul({mm: Q(1)}, ds))
        dimJ = rank(Jg, Bd) if Jg else 0
        J_in_gr = all(in_span(g, gr, Bd) for g in Jg) if Jg else True
        table.append((D, dimV, dimL, dimV - dimL, dimI, inM, len(gr), dimJ, J_in_gr, len(gr) - dimJ))
        if D == W + 3 and bN is not None:
            res['b_N in grI'] = in_span(bN, gr, Bd)
            res['b_N in J'] = in_span(bN, Jg, Bd)
    res['table'] = table
    # Rocha-Caridi count check: dim I_{2N+k} = p(k) for k < 2p-1
    res['dimI_vs_p(k)'] = [(row[0] - W, row[4], pnum(row[0] - W)) for row in table]
    res['time'] = round(time.time() - t0, 1)
    return res


if __name__ == '__main__':
    extra = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    models = [(3, 4), (3, 5), (4, 5), (2, 7), (3, 7), (2, 9), (3, 8), (2, 11), (4, 7), (5, 6), (3, 10), (5, 7), (2, 13)]
    if len(sys.argv) > 2:
        models = [tuple(map(int, s.split(','))) for s in sys.argv[2:]]
    allok = True
    for (p, q_) in models:
        r = audit_model(p, q_, extra)
        print(f"=== Vir({p},{q_})  N={r['N']}  c={r['c']}  [{r['time']}s]")
        for k in ['dim_sing_space', 'L2^N coeff nonzero', 'alpha', 'alpha_ok', 'beta', '5a+8b_ok', 'v_max_len',
                  'u_max_len', 'A', 'A_pq_matches', 'B', 'sigma_N(u) == A b_N + B d3f',
                  'u - B L_{-1}^3 v == 0', 'b_N in grI', 'b_N in J']:
            if k in r:
                print(f"   {k}: {r[k]}")
        print("   D  dimV  dimL  dimM  dimI  I<=M  dim grI  dim J  J<=grI  dim ker(pi)")
        for row in r['table']:
            print("   " + "  ".join(str(x) for x in row))
        print("   dim I_{2N+k} vs p(k):", r['dimI_vs_p(k)'])
        ok = r['alpha_ok'] and r['A_pq_matches'] and r['sigma_N(u) == A b_N + B d3f'] and r['u_max_len'] == r['N']
        if r['N'] >= 3:
            ok = ok and r['5a+8b_ok']
        for row in r['table']:
            D, dimV, dimL, dimM, dimI, inM, dgr, dJ, Jin, dker = row
            ok = ok and (dimM == dimI) and inM and Jin
            if D < 2 * r['N'] + 3:
                ok = ok and dker == 0
            if D == 2 * r['N'] + 3:
                ok = ok and dker == (1 if p >= 3 else 0)
        print("   ALL CHECKS OK:", ok)
        allok = allok and ok
    print("GLOBAL:", allok)
