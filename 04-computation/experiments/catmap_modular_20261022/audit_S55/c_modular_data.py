"""Audit C (THM-4762 (2)-(4)): independent modular data, modular-relation checks, orders of C = T S T^-1 S^-1.
Auditor for opus-2026-10-22-S55.  Written independently of catmap_mtc.py.

Builders (own code):
  * Fibonacci; SU(2)_k; Virasoro M(p,p') (Di Francesco 10.134, own Kac-table enumeration);
  * SU(3)_k from the Kac-Peterson formula in orthonormal R^3 coordinates (Weyl group = S_3 acting by permutations);
  * pointed categories of even lattices (any signature) from the discriminant form, own enumeration of L*/L;
  * D(Z_N); D(G) for G = S_3 from the Coste-Gannon-Ruelle group-theoretic formula (no hard-coded S-matrix);
  * Deligne products.
Every modular datum is checked: S symmetric unitary, S^2 a permutation matrix, (S T)^3 = Theta S^2 with
Theta = Gauss sum / D = exp(2 pi i c/8), Verlinde coefficients non-negative integers (rank <= 100).
Orders: eigenvalues of the unitary C are roots of unity; order = lcm of denominators, projective order = lcm of
denominators of ratios; both confirmed by explicit powers (C^o = I, C^po = lambda I, and minimality by prime divisors).
High precision (mpmath, 50 digits) for rank <= 16, numpy double otherwise.
Usage: python c_modular_data.py [all]   ('all' also rebuilds the full 135-entry census and compares with the JSON)"""
import sys, math, cmath, itertools, json
from fractions import Fraction
import numpy as np
import mpmath as mp

mp.mp.dps = 50
TWO_PI_I = 2j * mp.pi


def e(x):                                         # exp(2 pi i x), x a Fraction or mpf
    return mp.expjpi(2 * mp.mpf(x.numerator) / x.denominator) if isinstance(x, Fraction) else mp.expjpi(2 * x)


class MD:
    def __init__(self, name, S, h, c, unit=0):
        self.name, self.h, self.c, self.unit = name, [Fraction(x) for x in h], Fraction(c), unit
        self.r = len(self.h)
        self.S = S if isinstance(S, mp.matrix) else mp.matrix(S)
        self.T = mp.diag([e(x) for x in self.h])
        hu = self.h[unit]
        self.N = math.lcm(*[(x - hu).denominator for x in self.h])
        self.Nlin = math.lcm(*[(x - self.c / 24).denominator for x in self.h])


def fib():
    phi = (1 + mp.sqrt(5)) / 2
    D = mp.sqrt(2 + phi)
    return MD('Fib', [[1 / D, phi / D], [phi / D, -1 / D]], [0, Fraction(2, 5)], Fraction(14, 5))


def su2(k):
    n = k + 2
    S = [[mp.sqrt(mp.mpf(2) / n) * mp.sinpi(mp.mpf((a + 1) * (b + 1)) / n) for b in range(k + 1)] for a in range(k + 1)]
    return MD(f'SU(2)_{k}', S, [Fraction(a * (a + 2), 4 * n) for a in range(k + 1)], Fraction(3 * k, k + 2))


def minimal(p, pp):
    assert math.gcd(p, pp) == 1 and 2 <= p < pp
    reps = []
    for r in range(1, pp):
        for s in range(1, p):
            if (r, s) < (pp - r, p - s):
                reps.append((r, s))
    reps.sort(key=lambda rs: (rs != (1, 1), rs))
    h = [Fraction((p * r - pp * s) ** 2 - (p - pp) ** 2, 4 * p * pp) for r, s in reps]
    S = [[2 * mp.sqrt(mp.mpf(2) / (p * pp)) * (-1) ** (1 + s * rho + r * sig) * mp.sinpi(mp.mpf(p * r * rho) / pp)
          * mp.sinpi(mp.mpf(pp * s * sig) / p) for (rho, sig) in reps] for (r, s) in reps]
    md = MD(f'M({p},{pp})', S, h, 1 - Fraction(6 * (p - pp) ** 2, p * pp))
    md.reps = reps
    return md


def su3(k):
    """Kac-Peterson: S_{lm} ~ sum_{w in S_3} sgn(w) exp(-2 pi i <w(l+rho), m+rho>/(k+3)), orthonormal R^3 coordinates."""
    n = k + 3
    w1 = (Fraction(2, 3), Fraction(-1, 3), Fraction(-1, 3)); w2 = (Fraction(1, 3), Fraction(1, 3), Fraction(-2, 3))
    wts = [(a, b) for a in range(k + 1) for b in range(k + 1 - a)]

    def vec(a, b):
        return tuple(a * x + b * y for x, y in zip(w1, w2))
    rho = vec(1, 1)
    perms = list(itertools.permutations(range(3)))

    def sgn(pm):
        return (-1) ** sum(1 for i in range(3) for j in range(i + 1, 3) if pm[i] > pm[j])
    S = mp.matrix(len(wts), len(wts))
    for i, (a, b) in enumerate(wts):
        lr = tuple(x + y for x, y in zip(vec(a, b), rho))
        for j, (c, d) in enumerate(wts):
            mr = tuple(x + y for x, y in zip(vec(c, d), rho))
            S[i, j] = mp.fsum(sgn(pm) * e(-sum(lr[pm[t]] * mr[t] for t in range(3)) / n) for pm in perms)
    nrm = mp.sqrt(mp.fsum(abs(S[0, j]) ** 2 for j in range(len(wts))))
    ph = S[0, 0] / abs(S[0, 0])
    S = S / (nrm * ph)
    h = []
    for (a, b) in wts:
        v = vec(a, b); v2 = tuple(x + 2 * y for x, y in zip(v, rho))
        h.append(sum(x * y for x, y in zip(v, v2)) / (2 * n))
    return MD(f'SU(3)_{k}', S, h, Fraction(8 * k, k + 3))


def disc_form(G):
    G = [[Fraction(x) for x in row] for row in G]
    import sympy as sp
    M = sp.Matrix(G)
    Mi = M.inv()
    d = abs(int(M.det()))
    n = len(G)
    seen = {}
    box = range(d) if d > 1 else range(1)
    for v in itertools.product(box, repeat=n):
        x = Mi * sp.Matrix(v)
        key = tuple(Fraction(int(sp.numer(c)), int(sp.denom(c))) % 1 for c in x)
        seen.setdefault(key, None)
        if len(seen) == d:
            break
    reps = sorted(seen)
    def bil(x, y):
        return sum(x[i] * G[i][j] * y[j] for i in range(n) for j in range(n))
    q = [(bil(x, x) / 2) % 1 for x in reps]
    B = [[bil(x, y) % 1 for y in reps] for x in reps]
    return reps, q, B


def pointed(name, Gram):
    reps, q, B = disc_form(Gram)
    n = len(reps)
    S = [[e(-B[i][j]) / mp.sqrt(n) for j in range(n)] for i in range(n)]
    sig = int(sum(np.sign(np.linalg.eigvalsh(np.array(Gram, dtype=float)))))
    md = MD(name, S, q, Fraction(sig))
    md.reps, md.q = reps, q
    return md


def dzn(N):
    an = [(a, b) for a in range(N) for b in range(N)]
    S = [[e(Fraction(-(a * d + c * b), N)) / N for (c, d) in an] for (a, b) in an]
    return MD(f'D(Z_{N})', S, [Fraction(a * b, N) for a, b in an], 0)


# ---------------- D(G) from the group-theoretic formula (G = S_3 as permutations of {0,1,2})
def pmul(x, y):                                    # (x y)(i) = x(y(i))
    return tuple(x[y[i]] for i in range(len(y)))


def pinv(x):
    out = [0] * len(x)
    for i, xi in enumerate(x):
        out[xi] = i
    return tuple(out)


def double_S3():
    G = list(itertools.permutations(range(3)))
    eid = (0, 1, 2)
    classes, done = [], set()
    for g in G:
        if g in done:
            continue
        cl = sorted({pmul(pmul(x, g), pinv(x)) for x in G})
        done |= set(cl); classes.append(cl)
    classes.sort(key=lambda cl: (len(cl), cl))     # [e], [3-cycles], [transpositions] after sorting by size? sizes 1,2,3
    reps = [cl[0] for cl in classes]

    def cyc_type(g):
        fixed = sum(1 for i in range(3) if g[i] == i)
        return {3: 'e', 0: '3', 1: '2'}[fixed]

    def centralizer(a):
        return [x for x in G if pmul(x, a) == pmul(a, x)]

    def irreps(a):
        """list of (label, character function) for the centralizer of a."""
        Ca = centralizer(a)
        if len(Ca) == 6:
            sign = lambda g: (-1) ** sum(1 for i in range(3) for j in range(i + 1, 3) if g[i] > g[j])
            fixed = lambda g: sum(1 for i in range(3) if g[i] == i)
            return [('1', lambda g: 1), ('sgn', sign), ('2', lambda g: fixed(g) - 1)]
        # cyclic centralizer generated by a
        order = len(Ca)
        powers = {}
        x = eid
        for m in range(order):
            powers[x] = m; x = pmul(x, a)
        assert len(powers) == order
        z = mp.expjpi(mp.mpf(2) / order)
        return [(f'chi{j}', (lambda j: (lambda g: z ** (j * powers[g])))(j)) for j in range(order)]

    conj_x = {}
    for cl, a in zip(classes, reps):
        for g in cl:
            conj_x[g] = next(x for x in G if pmul(pmul(x, a), pinv(x)) == g)
    objs = []
    for cl, a in zip(classes, reps):
        for lab, ch in irreps(a):
            objs.append((a, cl, lab, ch))
    r = len(objs)
    S = mp.matrix(r, r)
    for i, (a, Ka, la, ca) in enumerate(objs):
        for j, (b, Kb, lb, cb) in enumerate(objs):
            tot = 0
            for g in Ka:
                for hh in Kb:
                    if pmul(g, hh) != pmul(hh, g):
                        continue
                    xg, xh = conj_x[g], conj_x[hh]
                    tot += mp.conj(mp.mpmathify(ca(pmul(pmul(pinv(xg), hh), xg)))) * mp.conj(mp.mpmathify(cb(pmul(pmul(pinv(xh), g), xh))))
            S[i, j] = mp.mpmathify(tot) / 6
    h = []
    for (a, Ka, la, ca) in objs:
        th = complex(ca(a)) / complex(ca(eid))
        f = Fraction(cmath.phase(th) / (2 * math.pi)).limit_denominator(12) % 1
        h.append(f)
    md = MD('D(S_3)', S, h, 0)
    md.labels = [(cyc_type(a), lab) for (a, Ka, lab, ca) in objs]
    return md


def deligne(*mds):
    S = mds[0].S; h = mds[0].h; c = mds[0].c; name = mds[0].name
    for m in mds[1:]:
        Snew = mp.matrix(S.rows * m.S.rows, S.cols * m.S.cols)
        for i in range(S.rows):
            for j in range(S.cols):
                for k in range(m.S.rows):
                    for l in range(m.S.cols):
                        Snew[i * m.S.rows + k, j * m.S.cols + l] = S[i, j] * m.S[k, l]
        S = Snew; h = [x + y for x in h for y in m.h]; c = c + m.c; name = name + ' x ' + m.name
    return MD(name, S, h, c)


# ---------------- checks and orders
def to_np(M):
    return np.array(M.tolist(), dtype=complex)


def check_md(md, verlinde_max=100):
    S = to_np(md.S); T = np.diag([complex(e(x)) for x in md.h]); r = md.r
    out = {}
    out['sym'] = np.abs(S - S.T).max()
    out['unit'] = np.abs(S @ S.conj().T - np.eye(r)).max()
    S2 = S @ S
    perm_ok = np.allclose(np.abs(S2), np.round(np.abs(S2)), atol=1e-9) and np.allclose(np.abs(S2).sum(0), 1, atol=1e-9)
    out['S2_perm'] = bool(perm_ok and np.allclose(np.linalg.matrix_power(S2, 2), np.eye(r), atol=1e-9))
    ST3 = np.linalg.matrix_power(S @ T, 3)
    i0 = md.unit
    j0 = int(np.argmax(np.abs(S2[i0])))
    Theta = ST3[i0, j0] / S2[i0, j0]
    out['ST3'] = np.abs(ST3 - Theta * S2).max()
    d = S[i0] / S[i0, i0]
    gauss = np.sum(d ** 2 * np.diag(T)) / math.sqrt(np.sum(np.abs(d) ** 2).real) if True else None
    # for non-pseudo-unitary data the categorical D may carry a sign: compare up to sign
    th_c = cmath.exp(2j * math.pi * float(md.c) / 8)
    out['Theta_vs_c'] = min(abs(Theta - th_c), abs(Theta + th_c))
    out['gauss_vs_Theta'] = min(abs(gauss - Theta), abs(gauss + Theta))
    if r <= verlinde_max:
        Nijk = np.einsum('im,jm,km,m->ijk', S, S, S.conj(), 1 / S[i0])
        out['verlinde_int'] = np.abs(Nijk - np.round(Nijk.real)).max()
        out['verlinde_nonneg'] = bool((np.round(Nijk.real) >= 0).all())
    return out


def rat_phase(z, cap=20000):
    a = (cmath.phase(complex(z)) / (2 * math.pi)) % 1.0
    f = Fraction(a).limit_denominator(cap)
    if abs(float(f) - a) > 1e-9:
        raise ValueError((z, a, f))
    return f % 1


def prime_factors(n):
    ps, d = [], 2
    while d * d <= n:
        if n % d == 0:
            ps.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        ps.append(n)
    return ps


def catmap(md):
    S, T = md.S, md.T
    if md.r <= 16:
        Si = mp.inverse(S); Ti = mp.diag([1 / T[i, i] for i in range(md.r)])
        C = T * S * Ti * Si
        if md.r == 1:
            return C, [complex(C[0, 0])], 'mp'
        ev = mp.eig(C, left=False, right=False)
        return C, [complex(z) for z in ev], 'mp'
    Snp = to_np(S); Tnp = np.diag([complex(e(x)) for x in md.h])
    C = Tnp @ Snp @ Tnp.conj() @ np.linalg.inv(Snp)
    return C, list(np.linalg.eigvals(C)), 'np'


def mpow(C, m, kind):
    if kind == 'mp':
        X = mp.eye(C.rows); B = C; k = m
        while k:
            if k & 1:
                X = X * B
            B = B * B; k >>= 1
        return X
    return np.linalg.matrix_power(C, m)


def dist_scalar(X, kind):
    if kind == 'mp':
        lam = X[0, 0]
        return float(max(abs(X[i, j] - (lam if i == j else 0)) for i in range(X.rows) for j in range(X.cols))), complex(lam)
    lam = X[0, 0]
    return float(np.abs(X - lam * np.eye(len(X))).max()), complex(lam)


def orders(md):
    C, ev, kind = catmap(md)
    ph = [rat_phase(z) for z in ev]
    assert all(abs(abs(z) - 1) < 1e-9 for z in ev)
    o = math.lcm(*[f.denominator for f in ph])
    po = math.lcm(*[((f - ph[0]) % 1).denominator for f in ph])
    tol = 1e-25 if kind == 'mp' else 1e-8
    dI, lam = dist_scalar(mpow(C, o, kind), kind)
    assert dI < tol and abs(lam - 1) < tol, (md.name, 'C^o != I', dI, lam)
    for pr in prime_factors(o):
        dd, lam2 = dist_scalar(mpow(C, o // pr, kind), kind)
        assert not (dd < 1e-6 and abs(lam2 - 1) < 1e-6), (md.name, 'order not minimal')
    dP, _ = dist_scalar(mpow(C, po, kind), kind)
    assert dP < tol, (md.name, 'C^po not scalar', dP)
    for pr in prime_factors(po):
        dd, _ = dist_scalar(mpow(C, po // pr, kind), kind)
        assert dd > 1e-6, (md.name, 'proj order not minimal')
    return o, po, C, kind


def ordA(N):
    if N == 1:
        return 1
    X = (1, 0, 0, 1); k = 0
    while True:
        X = ((2 * X[0] + X[1]) % N, (X[0] + X[1]) % N, (2 * X[2] + X[3]) % N, (X[2] + X[3]) % N); k += 1
        if X == (1 % N, 0, 0, 1 % N):
            return k


def report(md, census=None):
    chk = check_md(md)
    o, po, C, kind = orders(md)
    mN, mL = ordA(md.N), ordA(md.Nlin)
    row = dict(name=md.name, rank=md.r, N=md.N, mN=mN, po=po, o=o, Nl=md.Nlin, mL=mL)
    ok = (chk['sym'] < 1e-9 and chk['unit'] < 1e-9 and chk['S2_perm'] and chk['ST3'] < 1e-9 and chk['Theta_vs_c'] < 1e-9
          and chk.get('verlinde_int', 0) < 1e-6 and chk.get('verlinde_nonneg', True))
    cmp = ''
    if census is not None:
        cmp = 'census MATCH' if all(census.get(k) == row[k] for k in ('rank', 'N', 'mN', 'po', 'o', 'Nl', 'mL')) else f'census MISMATCH {census}'
    print(f'{md.name:34s} rank {md.r:3d} N {md.N:4d} ord_N(A) {mN:4d} proj-ord C {po:4d} ord C {o:4d} N\' {md.Nlin:5d} '
          f'ord_N\'(A) {mL:4d} | relations {"OK" if ok else "FAIL " + str(chk)} [{kind}] po|mN {mN % po == 0} o|mL {mL % o == 0} {cmp}',
          flush=True)
    return row, C, kind


if __name__ == '__main__':
    census = {r['name']: r for r in json.load(open('../catmap_census.json'))}
    alias = {'Fib': 'Fib = (G2)_1'}
    print('=== required list (audit item C) ===')
    req = [fib(), minimal(3, 4), minimal(2, 5), minimal(2, 7), minimal(4, 5)] + [su2(k) for k in range(1, 13)] + \
          [su3(1), su3(2)] + [dzn(N) for N in range(2, 11)]
    for md in req:
        report(md, census.get(alias.get(md.name, md.name)))
    d6 = double_S3()
    report(d6, census.get('D(S_3)'))
    print('D(S_3) labels (class type, centralizer irrep):', d6.labels, '; h =', [str(x) for x in d6.h])
    lat = {'A1 (definite)': [[2]], 'A1 + U (Lorentzian)': [[2, 0, 0], [0, 0, 1], [0, 1, 0]],
           'A1(-1) (negative definite)': [[-2]], 'A2 (definite)': [[2, -1], [-1, 2]],
           'A2 + U (Lorentzian)': [[2, -1, 0, 0], [-1, 2, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
           'golden [[2,1],[1,-2]] (indef.)': [[2, 1], [1, -2]], 'U = II_1,1 (Lorentzian)': [[0, 1], [1, 0]]}
    Cs = {}
    for nm, G in lat.items():
        md = pointed(nm, G)
        row, C, kind = report(md, census.get(nm))
        Cs[nm] = (md, C)
        print('     discriminant form: q =', [str(x) for x in md.q], '; c = signature =', md.c)
    # ---- compare the C matrices of lattices with (non-)isomorphic discriminant forms
    print('\n=== C for definite vs indefinite lattices (cyclic discriminant groups, ordered by a generator) ===')

    def ordered_C(md, C, gen_q):
        """reorder the discriminant group as 0, g, 2g, ... with q(g) = gen_q; return C in that basis (numpy)."""
        reps = md.reps
        n = len(reps)
        idx = {x: i for i, x in enumerate(reps)}
        for gi, g in enumerate(reps):
            if md.q[gi] != gen_q:
                continue
            orb = [tuple((k * c) % 1 for c in g) for k in range(n)]
            if len(set(orb)) == n:
                perm = [idx[x] for x in orb]
                Cn = to_np(C) if not isinstance(C, np.ndarray) else C
                return Cn[np.ix_(perm, perm)]
        return None
    CA1 = ordered_C(*Cs['A1 (definite)'], Fraction(1, 4))
    CA1U = ordered_C(*Cs['A1 + U (Lorentzian)'], Fraction(1, 4))
    CA1m = ordered_C(*Cs['A1(-1) (negative definite)'], Fraction(3, 4))
    CA2 = ordered_C(*Cs['A2 (definite)'], Fraction(1, 3))
    CA2U = ordered_C(*Cs['A2 + U (Lorentzian)'], Fraction(1, 3))
    np.set_printoptions(precision=4, suppress=True)
    print('C(A1) =\n', CA1, '\nmax|C(A1) - C(A1+U)| =', np.abs(CA1 - CA1U).max(), '; max|C(A1) - C(A1(-1))| =',
          np.abs(CA1 - CA1m).max(), '; max|conj C(A1) - C(A1(-1))| =', np.abs(CA1.conj() - CA1m).max())
    print('max|C(A2) - C(A2+U)| =', np.abs(CA2 - CA2U).max())
    # ---- Fibonacci x D(Z_2) x D(Z_3)
    print()
    report(deligne(fib(), dzn(2), dzn(3)), census.get('Fib = (G2)_1 x D(Z_2) x D(Z_3)'))
    # ---- check the session's hard-coded D(S_3) S-matrix against the group-theoretic one
    w = cmath.exp(2j * math.pi / 3)
    S_sess = (1 / 6) * np.array([[1, 1, 2, 2, 2, 2, 3, 3], [1, 1, 2, 2, 2, 2, -3, -3], [2, 2, 4, -2, -2, -2, 0, 0],
                                 [2, 2, -2, 4, -2, -2, 0, 0], [2, 2, -2, -2, -2, 4, 0, 0], [2, 2, -2, -2, 4, -2, 0, 0],
                                 [3, -3, 0, 0, 0, 0, 3, -3], [3, -3, 0, 0, 0, 0, -3, 3]], dtype=complex)
    h_sess = [Fraction(0)] * 4 + [Fraction(1, 3), Fraction(2, 3), Fraction(0), Fraction(1, 2)]
    S_mine = to_np(d6.S)
    found = False
    for perm in itertools.permutations(range(8)):
        if [d6.h[p] for p in perm] != h_sess:
            continue
        if np.allclose(S_mine[np.ix_(perm, perm)], S_sess, atol=1e-12):
            print('session D(S_3) S-matrix = group-theoretic S-matrix under relabelling', perm, [d6.labels[p] for p in perm])
            found = True
            break
    print('session D(S_3) hard-coded data reproduced:', found)
    if len(sys.argv) > 1 and sys.argv[1] == 'all':
        print('\n=== full census rebuild ===')
        allm = [fib()] + [minimal(p, pp) for p in range(2, 7) for pp in range(p + 1, 13) if math.gcd(p, pp) == 1 and p * pp <= 66]
        allm += [su2(k) for k in range(1, 61)] + [su3(k) for k in range(1, 11)]
        E8 = [[2, -1, 0, 0, 0, 0, 0, 0], [-1, 2, -1, 0, 0, 0, 0, 0], [0, -1, 2, -1, 0, 0, 0, -1], [0, 0, -1, 2, -1, 0, 0, 0],
              [0, 0, 0, -1, 2, -1, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0], [0, 0, 0, 0, 0, -1, 2, 0], [0, 0, -1, 0, 0, 0, 0, 2]]
        lat2 = dict(lat)
        lat2['A4 (definite)'] = [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2]]
        lat2['D4 (definite)'] = [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]]
        lat2['E8 (definite, unimodular)'] = E8
        allm += [pointed(nm, G) for nm, G in lat2.items()]
        allm += [dzn(N) for N in range(2, 31)] + [d6, deligne(fib(), dzn(2), dzn(3))]
        rows = []
        for md in allm:
            nm = {'Fib': 'Fib = (G2)_1', 'Fib x D(Z_2) x D(Z_3)': 'Fib = (G2)_1 x D(Z_2) x D(Z_3)'}.get(md.name, md.name)
            row, C, kind = report(md, census.get(nm))
            row['name'] = nm
            rows.append(row)
        json.dump(rows, open('c_census_audit.json', 'w'), indent=1)
        mism = [r['name'] for r in rows if any(census[r['name']].get(k) != r[k] for k in ('rank', 'N', 'mN', 'po', 'o', 'Nl', 'mL'))]
        sat = [r for r in rows if r['po'] == r['mN']]
        satN3 = [r for r in sat if r['N'] >= 3]
        print(f'\n{len(rows)} modular data rebuilt; census mismatches: {mism}; po | ord_N(A) in all: '
              f'{all(r["mN"] % r["po"] == 0 for r in rows)}; o | ord_N\'(A) in all: {all(r["mL"] % r["o"] == 0 for r in rows)}')
        print(f'saturating po = ord_N(A): {len(sat)} (of which N >= 3, i.e. po = pi(N)/2 literally: {len(satN3)}; '
              f'N <= 2: {[r["name"] for r in sat if r["N"] < 3]})')
        defs = {}
        for r in rows:
            defs.setdefault(r['mN'] // r['po'], []).append(r['name'])
        print('deficiency r = ord_N(A)/po value counts:', {k: len(v) for k, v in defs.items()}, '; r = 4:', defs.get(4))
        print('projective order 15:', [r['name'] for r in rows if r['po'] == 15])
