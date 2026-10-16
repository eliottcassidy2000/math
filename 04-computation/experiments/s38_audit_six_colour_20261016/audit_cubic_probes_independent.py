#!/usr/bin/env python3
"""S38 audit (independent referee): independent reconstruction of the four cubic-carrier probes of
cubic_carrier_three_colour_probes_20261016.py, with my own field arithmetic, my own edge finder and my own
colourings.

* Field: L' = Q(zeta55)[w]/(f), w = 2 theta, f(w) = w^3 + 2 w^2 - 4 w - 4 (monic; theta-basis is NOT used).
  Elements are 3 x 40 integer arrays (coefficients of zeta^i w^k) over a common denominator.
* v is built from its DEFINING FORMULA v = (theta - theta xi + xi^2) / (theta - theta xi^-1 + xi^-2), xi = zeta55^5,
  by solving the 120 x 120 multiplication system over Q (flint fmpq_mat), not from the JSON receipt; it is then
  compared with the receipt coefficients used by the S38 script.
* Probes: Cayley balls over the same direction sets.  Edge finder: EVERY pair of points is tested at three primes
  p = 1 mod 55 at which f has a root (ring maps zeta -> r, w -> s and the conjugate map zeta -> r^-1; a true unit
  pair maps to (A_a - A_b)(B_a - B_b) = 1 mod p, so no unit pair is lost), and every survivor is confirmed by the
  exact norm (x - y)(bar x - bar y) = 1 computed in Z[zeta][w].
* Cartesian structure: every point is tracked with its decomposition (cyclotomic part, v-coordinates); consistency
  on every edge is checked, every exact edge is classified, and an explicit 3-colouring is built from the product
  structure (residue maps zeta5 -> 1 mod 5, zeta11 -> 1 mod 11; SAT only for the zeta55 projection of probe 3).
* Q-ranks by elimination modulo two primes (exact: mod-p rank <= Q-rank <= the trivial upper bound).
* LLL on the SATURATED integer relation lattice (kernel-LLL on [N*B | I]); compared with the S38 numbers.
Output: audit_cubic_probes_independent.out (stdout) and .json in this folder.
"""
import itertools, json, math, os, sys, time
import numpy as np
import sympy as sp
from flint import fmpz_poly, fmpq_mat, fmpq, fmpz_mat, fmpz
from pysat.solvers import Solver

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
OUTJ = os.path.join(HERE, 'audit_cubic_probes_independent.json')
RES = {}
NZ = 40


def log(*a):
    print(*a, flush=True)


zs = sp.symbols('zs')
PHI55 = fmpz_poly([int(c) for c in reversed(sp.Poly(sp.cyclotomic_poly(55, zs), zs).all_coeffs())])
assert PHI55.degree() == 40


# ---------------------------------------------------------------- exact arithmetic in Z[zeta55][w]/(f)
def zred(p):
    return p % PHI55


def zmono(i):
    cs = [0] * (i % 55 + 1)
    cs[i % 55] = 1
    return zred(fmpz_poly(cs))


def zconj(p):
    cs = [0] * 56
    for i in range(p.degree() + 1):
        cs[(-i) % 55] += int(p[i])
    return zred(fmpz_poly(cs))


def emul(a, b):
    """a, b: lists of 3 fmpz_poly (coefficients of w^0, w^1, w^2). f: w^3 = -2w^2 + 4w + 4, w^4 = 8w^2 - 4w - 8."""
    c = [fmpz_poly(0)] * 5
    for i in range(3):
        for k in range(3):
            c[i + k] = c[i + k] + a[i] * b[k]
    c = [zred(x) for x in c]
    return [zred(c[0] + 4 * c[3] - 8 * c[4]), zred(c[1] + 4 * c[3] - 4 * c[4]), zred(c[2] - 2 * c[3] + 8 * c[4])]


def econj(a):
    return [zconj(x) for x in a]


def evec_to_elem(vec):
    """integer vector of length 120 (blocks k=0,1,2 of 40 zeta-coefficients) -> element (numerator)."""
    return [fmpz_poly([int(x) for x in vec[k * NZ:(k + 1) * NZ]]) for k in range(3)]


def elem_to_vec(a):
    out = np.zeros(3 * NZ, dtype=object)
    for k in range(3):
        for i in range(a[k].degree() + 1):
            out[k * NZ + i] = int(a[k][i])
    return out


def exact_unit(vec, den):
    a = evec_to_elem(vec)
    n = emul(a, econj(a))
    return n[1] == 0 and n[2] == 0 and n[0] == fmpz_poly([den * den])


# ---------------------------------------------------------------- v from its defining formula (rational solve)
def rational_elem_from_sympy(expr_z, expr_w=None):
    raise NotImplementedError


def mult_matrix(a_q):
    """a_q: rational element as 3 lists of 40 fmpq; returns the 120 x 120 fmpq_mat of multiplication by a."""
    # work with a common denominator
    den = 1
    for blk in a_q:
        for x in blk:
            den = den * int(x.q) // math.gcd(den, int(x.q))
    a_int = [fmpz_poly([int(x * den) for x in blk]) for blk in a_q]
    cols = []
    for k in range(3):
        for i in range(NZ):
            e = [fmpz_poly(0)] * 3
            e[k] = zmono(i)
            prod = emul(a_int, e)
            cols.append(elem_to_vec(prod))
    M = fmpq_mat(3 * NZ, 3 * NZ)
    for c, col in enumerate(cols):
        for r in range(3 * NZ):
            M[r, c] = fmpq(int(col[r]), den)
    return M


def build_v():
    # theta = w/2.  numerator N = theta - theta xi + xi^2, denominator Nb = theta - theta xi^-1 + xi^-2, xi = zeta^5.
    half = fmpq(1, 2)
    N = [[fmpq(0)] * NZ for _ in range(3)]
    Nb = [[fmpq(0)] * NZ for _ in range(3)]

    def add(E, k, i, c):
        p = zmono(i)
        for t in range(p.degree() + 1):
            E[k][t] += c * int(p[t])

    add(N, 1, 0, half); add(N, 1, 5, -half); add(N, 0, 10, fmpq(1))
    add(Nb, 1, 0, half); add(Nb, 1, -5, -half); add(Nb, 0, -10, fmpq(1))
    M = mult_matrix(Nb)
    rhs = fmpq_mat(3 * NZ, 1, [N[k][i] for k in range(3) for i in range(NZ)])
    x = M.solve(rhs)
    vq = [[x[k * NZ + i, 0] for i in range(NZ)] for k in range(3)]
    return vq, N, Nb


def q_to_intvec(eq, den):
    out = np.zeros(3 * NZ, dtype=np.int64)
    for k in range(3):
        for i in range(NZ):
            y = eq[k][i] * den
            assert y.q == 1
            out[k * NZ + i] = int(y.p)
    return out


def receipt_v_in_w_basis():
    """The S38 script's v (from the G7 JSON receipt), converted to the w-basis for comparison."""
    z, th, j = sp.symbols('z th j')
    coef_str = json.load(open(os.path.join(ROOT, '05-knowledge/results/elliptic_cubic_direction_20261016c.json')))['unit_coefficients']
    coef = [sp.sympify(s.replace('theta', 'th')) for s in coef_str]
    gauss = sum(sp.legendre_symbol(i, 11) * z ** (5 * i) for i in range(1, 11))
    vexpr = sp.expand(sum(c.subs(j, gauss) * z ** (5 * h) for h, c in enumerate(coef)))
    PH = sp.Poly(sp.cyclotomic_poly(55, z), z)
    out = [[fmpq(0)] * NZ for _ in range(3)]
    ww = sp.symbols('ww')
    vexpr_w = sp.expand(vexpr.subs(th, ww / 2))
    for k in range(3):
        ck = sp.Poly(sp.expand(vexpr_w).coeff(ww, k), z)
        ck = ck.rem(PH) if ck.degree() >= NZ else ck
        cs = ck.all_coeffs()[::-1]
        for i, c in enumerate(cs):
            c = sp.Rational(c)
            out[k][i] = fmpq(int(c.p), int(c.q))
    return out


# ---------------------------------------------------------------- Galois action, roots of unity
def sigma_q(eq, m):
    out = []
    for k in range(3):
        acc = [fmpq(0)] * NZ
        for i in range(NZ):
            if eq[k][i] == 0:
                continue
            p = zmono(i * m)
            for t in range(p.degree() + 1):
                acc[t] += eq[k][i] * int(p[t])
        out.append(acc)
    return out


def rou_q(i):
    out = [[fmpq(0)] * NZ for _ in range(3)]
    p = zmono(i)
    for t in range(p.degree() + 1):
        out[0][t] = fmpq(int(p[t]))
    return out


def qmul(a, b):
    den = 1
    for E in (a, b):
        for blk in E:
            for x in blk:
                den = den * int(x.q) // math.gcd(den, int(x.q))
    ai = [fmpz_poly([int(x * den) for x in blk]) for blk in a]
    bi = [fmpz_poly([int(x * den) for x in blk]) for blk in b]
    pr = emul(ai, bi)
    out = []
    for k in range(3):
        blk = [fmpq(0)] * NZ
        for t in range(pr[k].degree() + 1):
            blk[t] = fmpq(int(pr[k][t]), den * den)
        out.append(blk)
    return out


# ---------------------------------------------------------------- modular maps
def find_primes(cnt, avoid):
    ps = []
    p = 56
    while len(ps) < cnt:
        p += 55
        if not sp.isprime(p) or any(p % a == 0 for a in avoid):
            continue
        roots = [s for s in range(p) if (s ** 3 + 2 * s * s - 4 * s - 4) % p == 0]
        if not roots:
            continue
        g = sp.primitive_root(p)
        r = pow(g, (p - 1) // 55, p)
        ps.append((p, r, roots[0]))
    return ps


def modmaps(P, den, p, r, s):
    rinv = pow(r, -1, p)
    bas = np.zeros(3 * NZ, dtype=np.int64)
    basb = np.zeros(3 * NZ, dtype=np.int64)
    for k in range(3):
        for i in range(NZ):
            bas[k * NZ + i] = pow(r, i, p) * pow(s, k, p) % p
            basb[k * NZ + i] = pow(rinv, i, p) * pow(s, k, p) % p
    Pm = np.mod(P, p)
    dinv = pow(den, -1, p)
    A = (Pm @ bas) % p * dinv % p
    B = (Pm @ basb) % p * dinv % p
    return A.astype(np.int64), B.astype(np.int64)


def modular_unit_pairs(P, den, primes):
    maps = [modmaps(P, den, p, r, s) for (p, r, s) in primes]
    n = len(P)
    cand = []
    for a in range(n - 1):
        mask = None
        for (A, B), (p, _, _) in zip(maps, primes):
            m = ((A[a] - A[a + 1:]) % p) * ((B[a] - B[a + 1:]) % p) % p == 1
            mask = m if mask is None else (mask & m)
        idx = np.nonzero(mask)[0]
        for b in idx:
            cand.append((a, a + 1 + int(b)))
    return cand


# ---------------------------------------------------------------- probes
def main():
    T = time.time()
    vq, N, Nb = build_v()
    # checks on v
    # (1) v * Nb == N exactly
    prod = qmul(vq, Nb)
    assert all(prod[k][i] == N[k][i] for k in range(3) for i in range(NZ)), 'v * Nb != N'
    # (2) norm one
    nv = qmul(vq, sigma_conj(vq))
    assert nv[0][0] == 1 and all(nv[k][i] == 0 for k in range(3) for i in range(NZ) if (k, i) != (0, 0)), 'norm'
    rv = receipt_v_in_w_basis()
    same = all(rv[k][i] == vq[k][i] for k in range(3) for i in range(NZ))
    dens = sorted(set(int(x.q) for blk in vq for x in blk))
    log(f'v from the defining formula: v*Nb = N exactly, norm one exactly; equals the S38 receipt v: {same}; '
        f'denominators in the w-basis: {dens}')
    RES['v_checks'] = dict(defining_formula_ok=True, norm_one=True, equals_receipt=same, w_basis_denominators=dens)
    assert same
    # conjugates sigma_m, m = 1 mod 5
    ms = [m for m in range(1, 55) if math.gcd(m, 55) == 1 and m % 5 == 1]
    vconj = {m % 11: sigma_q(vq, m) for m in ms}
    # sigma_(-1 on xi)(v) = conj(v) = 1/v
    for r in vconj:
        if (11 - r) in vconj:
            pr = qmul(vconj[r], vconj[11 - r])
            assert pr[0][0] == 1 and all(pr[k][i] == 0 for k in range(3) for i in range(NZ) if (k, i) != (0, 0))
    log('v_r * v_(11-r) = 1 for all r (conjugate pairs)')
    DEN = 1
    for r in vconj:
        for blk in vconj[r]:
            for x in blk:
                DEN = DEN * int(x.q) // math.gcd(DEN, int(x.q))
    log(f'common denominator of the conjugates in the w-basis: {DEN}')
    RES['DEN_w_basis'] = DEN

    def dirs_of(spec):
        D = []  # (vector, tag) tag = ('c', i) for root zeta^i (sign s) or ('v', r, k5) etc.
        if 'mu10' in spec:
            D += [(rou_q(11 * k), ('c', 11 * k)) for k in range(5)]
        if 'mu22' in spec:
            D += [(rou_q(5 * k), ('c', 5 * k)) for k in range(11)]
        if 'mu110' in spec:
            D += [(rou_q(k), ('c', k)) for k in range(55)]
        if 'v' in spec:
            D += [(vconj[r], ('v', r, 0)) for r in sorted(vconj)]
        if 'v5' in spec:
            D += [(qmul(vconj[r], rou_q(11 * k)), ('v', r, k)) for r in sorted(vconj) for k in range(1, 5)]
        out, seen = [], set()
        for E, tag in D:
            for sgn in (1, -1):
                vec = q_to_intvec([[x * sgn for x in blk] for blk in E], DEN)
                t = tuple(vec.tolist())
                if t not in seen:
                    seen.add(t)
                    out.append((vec, tag, sgn))
        return out

    primes = find_primes(3, avoid=[DEN])
    log(f'filter primes (p, zeta55 -> r, w -> s): {primes}')
    RES['filter_primes'] = primes
    PROBES = [('mu10+v', ['mu10', 'v'], 3), ('mu22+v', ['mu22', 'v'], 3),
              ('mu110+v', ['mu110', 'v'], 2), ('mu10+v+v5', ['mu10', 'v', 'v5'], 2)]
    S38 = {'mu10+v': (30, 4971, 14220), 'mu22+v': (42, 13287, 37086), 'mu110+v': (130, 8581, 17230),
           'mu10+v+v5': (110, 6161, 12430)}
    RES['probes'] = []
    for name, spec, rad in PROBES:
        t0 = time.time()
        D = dirs_of(spec)
        # decomposition: probe 4 uses (alpha_0..alpha_10) in Z[zeta5]^11 (4 coords each, basis 1,z5,z5^2,z5^3);
        # probes 1-3 use (cyclotomic exponent-count vector in Z^55 reduced later, n in Z^10)
        rs = sorted(vconj)
        rindex = {r: i for i, r in enumerate(rs)}
        is4 = 'v5' in spec

        def z5vec(k):  # zeta5^k in basis 1, z5, z5^2, z5^3
            k %= 5
            v = np.zeros(4, dtype=np.int64)
            if k < 4:
                v[k] = 1
            else:
                v[:] = -1
            return v

        def zvec55(i):  # reduced coefficient vector of zeta55^i
            p_ = zmono(i)
            return np.array([int(p_[t]) if t <= p_.degree() else 0 for t in range(NZ)], dtype=np.int64)

        def dec_of(tag, sgn):
            """probe 4: (alpha_0..alpha_10) in Z[zeta5]^11, basis 1,z5,z5^2,z5^3 (canonical), flattened to 44;
            probes 1-3: [ell (40, canonical coefficients of the cyclotomic part in Z[zeta55]), n (10), rsum (1)]."""
            if is4:
                d = np.zeros((11, 4), dtype=np.int64)
                if tag[0] == 'c':
                    d[0] = sgn * z5vec(tag[1] // 11)
                else:
                    d[1 + rindex[tag[1]]] = sgn * z5vec(tag[2])
                return d.reshape(-1)
            d = np.zeros(NZ + 10 + 1, dtype=np.int64)
            if tag[0] == 'c':
                d[:NZ] = sgn * zvec55(tag[1])
                d[NZ + 10] = sgn
            else:
                d[NZ + rindex[tag[1]]] = sgn
            return d

        CANON = 44 if is4 else NZ + 10          # canonical part of the decomposition
        gens = [(vec, dec_of(tag, sgn)) for vec, tag, sgn in D]
        pts = {tuple([0] * (3 * NZ)): 0}
        order = [np.zeros(3 * NZ, dtype=np.int64)]
        decs = [gens[0][1] * 0]
        frontier = [0]
        rediscovery_conflicts = 0
        for _ in range(rad):
            new = []
            for pi in frontier:
                for vec, dd in gens:
                    q = order[pi] + vec
                    t = tuple(q.tolist())
                    if t not in pts:
                        pts[t] = len(order)
                        order.append(q)
                        decs.append(decs[pi] + dd)
                        new.append(pts[t])
                    else:
                        if not np.array_equal((decs[pi] + dd)[:CANON], decs[pts[t]][:CANON]):
                            rediscovery_conflicts += 1
            frontier = new
        P = np.array(order)
        n = len(P)
        # exact edge finding (all pairs, modular filter at three primes, exact norm on survivors)
        cand = modular_unit_pairs(P, DEN, primes)
        E = [(a, b) for a, b in cand if exact_unit(P[a] - P[b], DEN)]
        gset = {tuple(vec.tolist()): i for i, (vec, _) in enumerate(gens)}
        nongen = 0
        dec_consistent = True
        for a, b in E:
            t = tuple((P[b] - P[a]).tolist())
            if t not in gset:
                nongen += 1
                continue
            if not np.array_equal((decs[b] - decs[a])[:CANON], gens[gset[t]][1][:CANON]):
                dec_consistent = False
        # every pair at a generator difference is an edge (the converse inclusion)
        gen_pairs = 0
        for a in range(n):
            for vec, _ in gens:
                t = tuple((P[a] + vec).tolist())
                if t in pts and pts[t] > a:
                    gen_pairs += 1
        adj = [set() for _ in range(n)]
        for a, b in E:
            adj[a].add(b); adj[b].add(a)
        tri = sum(len(adj[a] & adj[b]) for a, b in E) // 3
        # explicit 3-colouring from the product structure
        if is4:
            val = [int(d.sum()) % 5 for d in decs]       # zeta5 -> 1 on each of the eleven Z[zeta5] coordinates
            cyc = [0, 1, 0, 1, 2]
            col = [cyc[x] for x in val]
            method = 'sum of residues zeta5 -> 1 (mod 5) over the eleven pentagonal factors, then C5 -> 3 colours'
        elif 'mu110' not in spec:
            mod = 5 if 'mu10' in spec else 11
            val = [(int(d[NZ + 10]) + int(d[NZ:NZ + 10].sum())) % mod for d in decs]
            cyc = [0, 1, 0, 1, 2] if mod == 5 else [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2]
            col = [cyc[x] for x in val]
            method = f'residue map zeta{mod} -> 1 (mod {mod}) plus sum of v-coordinates, then C{mod} -> 3 colours'
        else:
            proj = {}
            pid = []
            for d in decs:
                key = tuple(int(x) for x in d[:NZ])
                if key not in proj:
                    proj[key] = len(proj)
                pid.append(proj[key])
            keys = list(proj)
            kv = np.array(keys, dtype=np.int64)
            roots = []
            for i in range(55):
                rv_ = zvec55(i)
                roots += [rv_, -rv_]
            kd = {k: i for i, k in enumerate(keys)}
            pE = set()
            for i, k in enumerate(kv):
                for rv_ in roots:
                    t = tuple((k + rv_).tolist())
                    if t in kd:
                        j_ = kd[t]
                        pE.add((min(i, j_), max(i, j_)))
            m = len(keys)
            var = lambda x, c: 3 * x + c + 1
            with Solver(name='glucose4') as S:
                for x in range(m):
                    S.add_clause([var(x, c) for c in range(3)])
                for a_, b_ in pE:
                    for c in range(3):
                        S.add_clause([-var(a_, c), -var(b_, c)])
                assert S.solve()
                mm = set(l for l in S.get_model() if l > 0)
            pc = [next(c for c in range(3) if var(x, c) in mm) for x in range(m)]
            col = [(pc[pid[i]] + int(decs[i][NZ:NZ + 10].sum())) % 3 for i in range(n)]
            method = (f'SAT 3-colouring of the cyclotomic projection ({m} points, {len(pE)} root-of-unity edges) '
                      f'plus sum of v-coordinates mod 3')
        col_ok = all(col[a] != col[b] for a, b in E)
        avgdeg = 2 * len(E) / n
        rec = dict(probe=name, directions=len(D), radius=rad, vertices=n, modular_candidates=len(cand),
                   exact_unit_edges=len(E), generator_difference_pairs=gen_pairs, non_generator_edges=nongen,
                   decomposition_consistent_on_edges=dec_consistent, bfs_rediscovery_conflicts=rediscovery_conflicts,
                   triangles=tri, average_degree=round(avgdeg, 3), product_colouring_valid=col_ok,
                   colouring_method=method, s38_reported=S38[name],
                   matches_s38=(len(D), n, len(E)) == S38[name], seconds=round(time.time() - t0, 1))
        RES['probes'].append(rec)
        log(json.dumps(rec))

    # ---------------- Q-ranks modulo primes
    def rank_mod(rows, p):
        A = np.array([[int(x) % p for x in r] for r in rows], dtype=np.int64)
        rk = 0
        rows_, cols_ = A.shape
        for c in range(cols_):
            piv = None
            for r in range(rk, rows_):
                if A[r, c] % p:
                    piv = r
                    break
            if piv is None:
                continue
            A[[rk, piv]] = A[[piv, rk]]
            inv = pow(int(A[rk, c]), -1, p)
            A[rk] = A[rk] * inv % p
            for r in range(rows_):
                if r != rk and A[r, c]:
                    A[r] = (A[r] - A[r, c] * A[rk]) % p
            rk += 1
        return rk

    Z = [q_to_intvec(rou_q(k), DEN) for k in range(55)]
    V = [q_to_intvec(vconj[r], DEN) for r in sorted(vconj)]
    VZ = [q_to_intvec(qmul(vconj[r], rou_q(k)), DEN) for r in sorted(vconj) for k in range(55)]
    M5 = [q_to_intvec(rou_q(11 * k), DEN) for k in range(5)]
    VM5 = [q_to_intvec(qmul(vconj[r], rou_q(11 * k)), DEN) for r in sorted(vconj) for k in range(5)]
    ranks = {}
    for p in (1000003, 998244353):
        ranks[p] = dict(roots=rank_mod(Z, p), roots_plus_v=rank_mod(Z + V, p), roots_plus_all_zv=rank_mod(Z + VZ, p),
                        mu5_plus_vmu5=rank_mod(M5 + VM5, p))
    log('ranks mod p (lower bounds for Q-ranks; trivial upper bounds 40, 50, 120, 44):', ranks)
    RES['ranks_mod_p'] = ranks

    # ---------------- saturated relation lattice + LLL for the four triples
    RES['lll'] = []
    one = rou_q(0)
    for trip in [(1, 3, 4), (2, 9, 3), (1, 10, 2), (1, 2, 3)]:
        t0 = time.time()
        rows = []
        for w in [one] + [vconj[r] for r in trip]:
            rows += [q_to_intvec(qmul(w, rou_q(k)), DEN) for k in range(40)]
        nr = len(rows)
        Nbig = 10 ** 12
        M = fmpz_mat(nr, 3 * NZ + nr)
        for i, r in enumerate(rows):
            for c in range(3 * NZ):
                M[i, c] = int(r[c]) * Nbig
            M[i, 3 * NZ + i] = 1
        R = M.lll()
        kern = []
        for i in range(nr):
            if all(R[i, c] == 0 for c in range(3 * NZ)):
                kern.append([int(R[i, 3 * NZ + c]) for c in range(nr)])
        # verify kernel vectors are relations and that the lattice is saturated (det of kernel basis Gram vs
        # nothing: saturation is automatic for the kernel-LLL construction)
        for kv_ in kern[:3]:
            s = np.zeros(3 * NZ, dtype=object)
            for c, cf in enumerate(kv_):
                if cf:
                    s = s + np.array(rows[c], dtype=object) * cf
            assert all(int(x) == 0 for x in s)
        L1 = sorted(sum(abs(x) for x in kv_) for kv_ in kern)
        Linf = [max(abs(x) for x in kv_) for kv_ in kern]
        rec = dict(triple=list(trip), kernel_rank=len(kern), min_L1=L1[0], L1_first5=L1[:5],
                   min_Linf=min(Linf), seconds=round(time.time() - t0, 1))
        RES['lll'].append(rec)
        log('saturated kernel LLL:', json.dumps(rec))
    RES['seconds'] = round(time.time() - T, 1)
    json.dump(RES, open(OUTJ, 'w'), indent=1, default=str)
    log('wrote', OUTJ)


def sigma_conj(eq):
    """complex conjugation zeta -> zeta^-1, w fixed."""
    return sigma_q(eq, -1)


if __name__ == '__main__':
    main()
