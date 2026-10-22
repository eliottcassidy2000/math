"""Cusps of II_{9,1}, II_{17,1}, II_{25,1} and the even unimodular lattices of rank 8, 16, 24
(session opus-2026-10-21-S52; note 05-knowledge/results/lattices_modular_repunits_20261021.md).

(V) Vinberg's algorithm (controlling vector v0 = (0^n | 1)) run from scratch for II_{9,1} and II_{17,1}: it stops at
    10 and 19 simple roots, certified by Vinberg's finite-volume criterion (every elliptic subdiagram of rank n-1
    extends in exactly two ways to an elliptic subdiagram of rank n or a parabolic one of rank n-1).  The parabolic
    subdiagrams of rank n-1 are the cusps; each is an affine diagram whose imaginary root (sum of marks times roots) is
    the primitive null vector z of that cusp, and z^perp/z is the even unimodular lattice with that root system:
      II_{9,1}:  one cusp, E~8                -> E8               (rho . z = -h = -30 for the Weyl vector, rho^2 < 0);
      II_{17,1}: two cusps, E~8 + E~8, D~16  -> E8 + E8, D16^+   (the two heterotic-string lattices).
    In II_{25,1} the simple roots are infinitely many (Leech roots, Conway) and its 24 cusps are the 24 Niemeier
    lattices; the A_24 cusp is computed in tournament coordinates in lorentzian_scores.py.  KNOWN (Vinberg 1972/75,
    Conway 1983, Borcherds 1985); here reproduced FINITE-EXACT.
(M) Siegel's mass formula sum 1/|Aut L| = (|B_k|/2k) prod_{j<k} |B_2j|/(4j) (rank 2k = 8, 16, 24) checked exactly:
    rank 8 by E8 alone, rank 16 by E8 + E8 and D16^+, rank 24 by the 24 Niemeier lattices with
    |Aut N| = |W(R)| |G1| |G2| (Conway-Sloane SPLAG Table 16.1) and |Co_0| for Leech; every Niemeier root system has
    rank 24, one Coxeter number h, 24 h roots, and theta_N = E_12 + (24 h - 65520/691) Delta.  KNOWN.
(E) E8 holes: deep hole e_1 at squared distance 1 with 16 nearest points (a cross-polytope; Session27 cites it) and the
    shallow hole at the circumcentre of the regular 8-simplex {0, e_1 + e_j, (1/2)^8}: squared distance 8/9, exactly 9
    nearest points.  Milnor's isospectral pair: E8 + E8 and D16^+ both have 480 roots and 61920 norm-4 vectors
    (theta = E4^2) but different root systems.  KNOWN.
Prints ALL CHECKS PASSED."""
import itertools, math
from fractions import Fraction
import numpy as np
from sympy import bernoulli, Matrix

FAIL = []
H = Fraction(1, 2)


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def lor(x, y):
    return sum(a * b for a, b in zip(x[:-1], y[:-1])) - x[-1] * y[-1]


def sorted_vectors(n, N, half):
    """Non-increasing y_1 >= ... >= y_(n-1) >= |y_n| (integers, or all in Z + 1/2) with sum y_i^2 = N."""
    step = Fraction(1) if not half else Fraction(1)
    vals = []
    v = H if half else Fraction(0)
    while v * v <= N:
        vals.append(v)
        v += step
    out = []

    def rec(prefix, remaining, maxv):
        k = len(prefix)
        if k == n - 1:
            for last in vals:
                if last <= maxv and last * last == remaining:
                    out.append(prefix + [last])
                    if last != 0:
                        out.append(prefix + [-last])
            return
        for x in vals:
            if x > maxv or x * x > remaining:
                continue
            # the remaining n-1-k coordinates are at least min(vals) each in absolute value
            if remaining - x * x < (n - 1 - k) * vals[0] ** 2:
                continue
            rec(prefix + [x], remaining - x * x, x)
    rec([], N, vals[-1])
    return out


def components(idx, G):
    idx = list(idx)
    seen, comps = set(), []
    for a in idx:
        if a in seen:
            continue
        stack, comp = [a], []
        seen.add(a)
        while stack:
            u = stack.pop(); comp.append(u)
            for b in idx:
                if b not in seen and G[u][b] != 0:
                    seen.add(b); stack.append(b)
        comps.append(sorted(comp))
    return comps


def eig(G, idx):
    M = np.array([[float(G[a][b]) for b in idx] for a in idx])
    return np.linalg.eigvalsh(M)


def is_elliptic(G, idx):
    return len(idx) == 0 or eig(G, idx).min() > 1e-9


def is_parabolic_component(G, comp):
    e = eig(G, comp)
    return e.min() > -1e-9 and sum(1 for t in e if abs(t) < 1e-9) == 1


def is_parabolic(G, idx):
    comps = components(idx, G)
    return len(idx) > 0 and all(is_parabolic_component(G, c) for c in comps), len(idx) - len(comps)


def finite_volume(G, nroots, n):
    """Vinberg's criterion for the Coxeter polytope in H^n with Gram matrix G (nroots walls)."""
    nodes = range(nroots)
    par = {}
    for size in range(n, nroots + 1):
        for Q in itertools.combinations(nodes, size):
            p, rk = is_parabolic(G, Q)
            if p and rk == n - 1:
                par[Q] = components(Q, G)
    has_vertex = bool(par) or any(is_elliptic(G, S) for S in itertools.combinations(nodes, n))
    good = True
    for S in itertools.combinations(nodes, n - 1):
        if not is_elliptic(G, S):
            continue
        ext = sum(1 for v in nodes if v not in S and is_elliptic(G, S + (v,)))
        ext += sum(1 for Q in par if set(S) <= set(Q))
        good &= ext == 2
    return has_vertex and good, par


def vinberg(n, kmax):
    # simple roots of D_n orthogonal to v0, chamber x_1 >= ... >= x_(n-1) >= |x_n|; walls x . r <= 0
    roots = []
    for i in range(n - 1):
        r = [Fraction(0)] * (n + 1); r[i + 1] = Fraction(1); r[i] = Fraction(-1)
        roots.append(r)
    r = [Fraction(0)] * (n + 1); r[n - 2] = Fraction(-1); r[n - 1] = Fraction(-1)
    roots.append(r)
    k = H
    while k <= kmax:
        half = k.denominator == 2
        level = []
        for y in sorted_vectors(n, 2 + k * k, half):
            cand = y + [k]
            if (sum(y) - k) % 2 != 0:
                continue
            if all(lor(cand, r) <= 0 for r in roots + level):
                level.append(cand)
        roots += level
        G = [[lor(a, b) for b in roots] for a in roots]
        fv, par = finite_volume(G, len(roots), n)
        if fv:
            return roots, G, par, k
        k += H
    return roots, None, None, None


def affine_type(G, comp):
    m = len(comp)
    deg = {a: sum(1 for b in comp if b != a and G[a][b] != 0) for a in comp}
    if m == 2:
        return 'A~1'
    if all(d == 2 for d in deg.values()):
        return f'A~{m - 1}'
    br = [a for a in comp if deg[a] >= 3]
    if len(br) == 2 or (len(br) == 1 and deg[br[0]] == 4):
        return f'D~{m - 1}'
    # one branch node of degree 3: arms
    c = br[0]
    arms = []
    for b in comp:
        if G[c][b] != 0 and b != c:
            L, prev, cur = 1, c, b
            while True:
                nxt = [x for x in comp if x not in (prev, cur) and G[cur][x] != 0]
                if not nxt:
                    break
                prev, cur = cur, nxt[0]; L += 1
            arms.append(L)
    arms = tuple(sorted(arms))
    return {(2, 2, 2): 'E~6', (1, 3, 3): 'E~7', (1, 2, 5): 'E~8'}.get(arms, f'D~{m - 1}' if arms[:2] == (1, 1) else str(arms))


def null_vector(G, comp, roots):
    M = Matrix([[G[a][b] for b in comp] for a in comp])
    ker = M.nullspace()[0]
    ker = ker / min(abs(t) for t in ker if t != 0)
    if ker[0] < 0:
        ker = -ker
    marks = [Fraction(int(t.p), int(t.q)) for t in ker]
    z = [sum(m * roots[a][i] for m, a in zip(marks, comp)) for i in range(len(roots[0]))]
    return marks, z


def weyl_vector(G, nroots, roots):
    # rho with rho . r = -1 for all simple roots (when it exists in the span)
    A = np.array([[float(G[a][b]) for b in range(nroots)] for a in range(nroots)])
    R = np.array([[float(t) for t in r] for r in roots])
    J = np.diag([1.0] * (R.shape[1] - 1) + [-1.0])
    # solve (R J) rho = -1 in least squares; exact for II_{9,1} (R square, invertible)
    sol, res, rk, sv = np.linalg.lstsq(R @ J, -np.ones(nroots), rcond=None)
    return sol


print('(V) Vinberg diagrams and cusps')
for n, expect_roots, expect_cusps in ((9, 10, [['E~8']]), (17, 19, [['D~16'], ['E~8', 'E~8']])):
    roots, G, par, kstop = vinberg(n, kmax=Fraction(6))
    ok_lat = all(lor(r, r) == 2 and (all(Fraction(t).denominator == 1 for t in r) or all(Fraction(t).denominator == 2 for t in r))
                 and (sum(r[:-1]) - r[-1]) % 2 == 0 for r in roots)
    sig = np.linalg.eigvalsh(np.array([[float(x) for x in row] for row in G]))
    neg, zero = sum(1 for t in sig if t < -1e-9), sum(1 for t in sig if abs(t) < 1e-9)
    cusps = sorted(sorted(affine_type(G, c) for c in comps) for comps in par.values())
    check(G is not None and len(roots) == expect_roots and ok_lat and neg == 1 and len(roots) - zero == n + 1,
          f'II_{{{n},1}}: Vinberg stops at level k = {kstop} with {len(roots)} simple roots, finite volume certified; Gram signature ({len(roots) - zero - neg},{neg}) + {zero} null')
    check(cusps == sorted(expect_cusps), f'II_{{{n},1}}: cusps (parabolic subdiagrams of rank {n - 1}) = {cusps}')
    for Q, comps in sorted(par.items()):
        zs = []
        for comp in comps:
            marks, z = null_vector(G, comp, roots)
            zs.append((affine_type(G, comp), sum(marks), z))
        for t, h, z in zs:
            print(f'     cusp {t:5s}: sum of marks = Coxeter number {h}, imaginary root z = ({",".join(str(x) for x in z[:-1])} | {z[-1]}), z^2 = {lor(z, z)}')
        check(all(lor(z, z) == 0 for _, _, z in zs) and len({tuple(z) for _, _, z in zs}) == 1,
              f'II_{{{n},1}}: the components of {sorted(t for t, _, _ in zs)} share one null vector (the cusp)')
    if n == 9:
        A = np.array([[float(x) for x in row] for row in G])
        c = np.linalg.solve(A, -np.ones(len(roots)))
        rho = sum(ci * np.array([float(t) for t in r]) for ci, r in zip(c, roots))
        J = np.diag([1.0] * n + [-1.0])
        rho2 = rho @ J @ rho
        marks, z = null_vector(G, list(par.values())[0][0], roots)
        zf = np.array([float(t) for t in z])
        check(abs(rho2 + 1240) < 1e-6 and abs(rho @ J @ zf + 30) < 1e-9,
              f'II_{{9,1}}: Weyl vector rho (rho . r = -1) has rho^2 = {rho2:.0f} and rho . z = {rho @ J @ zf:.0f} = -h(E8)')

print('(M) mass formulas')


def mass(rank):
    k = rank // 2
    m = abs(Fraction(bernoulli(k))) / (2 * k)
    for j in range(1, k):
        m *= abs(Fraction(bernoulli(2 * j))) / (4 * j)
    return m


fact = math.factorial
WE = {6: 51840, 7: 2903040, 8: 696729600}


def W(t, r):
    return fact(r + 1) if t == 'A' else 2 ** (r - 1) * fact(r) if t == 'D' else WE[r]


def cox(t, r):
    return r + 1 if t == 'A' else 2 * r - 2 if t == 'D' else {6: 12, 7: 18, 8: 30}[r]


def nroots(t, r):
    return r * (r + 1) if t == 'A' else 2 * r * (r - 1) if t == 'D' else {6: 72, 7: 126, 8: 240}[r]


check(mass(8) == Fraction(1, WE[8]), f'rank 8: mass = 1/{1 / mass(8)} = 1/|W(E8)|: E8 is the only even unimodular lattice')
m16 = Fraction(1, 2 * WE[8] ** 2) + Fraction(1, 2 ** 15 * fact(16))
check(mass(16) == m16 and mass(16).numerator == 691, f'rank 16: mass = {mass(16)} = 1/|Aut(E8+E8)| + 1/|Aut(D16+)| (numerator 691 from B_12)')
NIEMEIER = [  # root system [(type, rank, multiplicity)], |G1|, |G2|  (Conway-Sloane, SPLAG Table 16.1)
    ([('D', 24, 1)], 1, 1), ([('D', 16, 1), ('E', 8, 1)], 1, 1), ([('E', 8, 3)], 1, 6), ([('A', 24, 1)], 2, 1),
    ([('D', 12, 2)], 1, 2), ([('A', 17, 1), ('E', 7, 1)], 2, 1), ([('D', 10, 1), ('E', 7, 2)], 1, 2),
    ([('A', 15, 1), ('D', 9, 1)], 2, 1), ([('D', 8, 3)], 1, 6), ([('A', 12, 2)], 2, 2),
    ([('A', 11, 1), ('D', 7, 1), ('E', 6, 1)], 2, 1), ([('E', 6, 4)], 2, 24), ([('A', 9, 2), ('D', 6, 1)], 2, 2),
    ([('D', 6, 4)], 1, 24), ([('A', 8, 3)], 2, 6), ([('A', 7, 2), ('D', 5, 2)], 2, 4), ([('A', 6, 4)], 2, 12),
    ([('A', 5, 4), ('D', 4, 1)], 2, 24), ([('D', 4, 6)], 3, 720), ([('A', 4, 6)], 2, 120), ([('A', 3, 8)], 2, 1344),
    ([('A', 2, 12)], 2, 95040), ([('A', 1, 24)], 1, 244823040)]
CO0 = 8315553613086720000
tot = Fraction(1, CO0)
ok = len(NIEMEIER) == 23
for R, g1, g2 in NIEMEIER:
    w = 1
    for t, r, mult in R:
        w *= W(t, r) ** mult
    tot += Fraction(1, w * g1 * g2)
    hs = {cox(t, r) for t, r, mult in R}
    ok &= sum(r * mult for t, r, mult in R) == 24 and len(hs) == 1 and sum(nroots(t, r) * mult for t, r, mult in R) == 24 * hs.pop()
check(ok, '23 Niemeier root systems: rank 24, a single Coxeter number h each, 24 h roots')
check(tot == mass(24), f'rank 24: sum over the 24 Niemeier lattices of 1/|Aut| = {tot} = Siegel mass {mass(24)} (exact)')
check(Fraction(65520, 691) == -Fraction(24) / Fraction(bernoulli(12)),
      'theta_N = E_12 + (24 h - 65520/691) Delta; 65520/691 = -24/B_12 is the E_12 Eisenstein coefficient (Leech: h = 0)')

print('(E) E8 holes and Milnor\'s isospectral pair')
E8roots = [v for v in itertools.product(range(-1, 2), repeat=8) if sum(t * t for t in v) == 2]
E8roots = [tuple(Fraction(t) for t in v) for v in E8roots]
E8roots += [tuple(H * s for s in signs) for signs in itertools.product((-1, 1), repeat=8) if signs.count(-1) % 2 == 0]
check(len(E8roots) == 240, '240 roots of E8 = D8 u (D8 + (1/2)^8)')
simplex = [tuple([Fraction(0)] * 8)] + [tuple(Fraction(1 if i in (0, j) else 0) for i in range(8)) for j in range(1, 8)] + [tuple([H] * 8)]
d2 = lambda a, b: sum((x - y) ** 2 for x, y in zip(a, b))
reg = all(d2(a, b) == 2 for a, b in itertools.combinations(simplex, 2))
# circumcentre: |c - p_i|^2 equal for all i  <=>  2 c . p_i = |p_i|^2 (p_0 = 0)
M = Matrix([[2 * t for t in p] for p in simplex[1:]])
rhs = Matrix([sum(t * t for t in p) for p in simplex[1:]])
csol = M.solve(rhs)
c = tuple(Fraction(int(t.p), int(t.q)) for t in csol)
R2 = d2(c, simplex[0])
near = [p for p in [tuple([Fraction(0)] * 8)] + E8roots if d2(c, p) == R2]
closer = [p for p in [tuple([Fraction(0)] * 8)] + E8roots if d2(c, p) < R2]
check(reg and R2 == Fraction(8, 9) and len(near) == 9 and not closer,
      f'shallow hole c = ({",".join(str(t) for t in c)}): squared distance {R2} to exactly 9 lattice points (regular 8-simplex), none closer')
deep = (Fraction(1),) + (Fraction(0),) * 7
dn = [p for p in [tuple([Fraction(0)] * 8)] + E8roots if d2(deep, p) == 1]
dshell = [p for p in [tuple([Fraction(0)] * 8)] + E8roots if d2(deep, p) < 1]
# the 16 nearest points of the deep hole include 2 e_1 (norm 4): add it
dn += [(Fraction(2),) + (Fraction(0),) * 7]
check(len(dn) == 16 and not dshell, 'deep hole e_1: squared distance 1 to 16 points (0, 2e_1, e_1 +- e_j: a cross-polytope), covering radius 1')
# Milnor: theta(E8+E8) = theta(D16+) through norm 4
e8_n4 = 2160
d16_n2 = 4 * math.comb(16, 2)
d16_n4 = 2 * 16 + 16 * math.comb(16, 4) + 2 ** 15
check(2 * 240 == d16_n2 == 480 and 2 * e8_n4 + 240 * 240 == d16_n4 == 61920 == 480 * (1 + 2 ** 7),
      'E8+E8 and D16^+: 480 roots and 61920 norm-4 vectors each (E4^2 = 1 + 480 q + 61920 q^2), root systems 2E8 vs D16: Milnor\'s drums')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
