# -*- coding: utf-8 -*-
# G1 groundwork: THE DEGREE ARITHMETIC OF 2-JET KELLER MAPS
# kind-pasteur 2026-07-26 (S133 scratch) -- scratch only, no canon writes.
#
# Legs:
#  1. Known 2-jet Keller maps: exact polynomial inverses (field degree 1), det J,
#     graded D-system; THM-1310 wild control (field degree 3).
#  2. Frame/projected-system lemma: Res_z(E_i,E_j) = a_k^2 - u_k b_k, kernel identities.
#  3. Counting engines (space Groebner mod p vs projected plane system) cross-validation.
#  4. TABLE A: generic degree per box vs BKK formula N_max = s^3-(s-2)^3 (staircase).
#  5. TABLE C/E: shape strata of box (1,2,3) (+ (2,3,4) conic cap), D-staircase refinements.
#  6. Reachability in minimal box (0,1,2): witnesses for N=1..8; quartic S4 witness + stats.
import random, time, itertools, sys
from collections import Counter
from sympy import (symbols, Matrix, Poly, groebner, expand, factor_list, resultant,
                   Integer, Rational, simplify, factor)
from sympy.polys.orderings import grevlex

x, y, z, t = symbols('x y z t')
SEED = 20260726
random.seed(SEED)
PR  = 10007
PR2 = 31013

def say(*a):
    print(*a); sys.stdout.flush()

def br(u, v, w):
    return Matrix.hstack(u, v, w).det()

def dsystem(A, B, C):
    Ax, Ay = A.diff(x), A.diff(y)
    Bx, By = B.diff(x), B.diff(y)
    Cx, Cy = C.diff(x), C.diff(y)
    D5 = 2*br(Ax, Ay, A)
    D4 = br(Ax, Ay, B) + 2*br(Ax, By, A) + 2*br(Bx, Ay, A)
    D3 = br(Ax, By, B) + br(Bx, Ay, B) + 2*(br(Ax, Cy, A) + br(Cx, Ay, A) + br(Bx, By, A))
    D2 = br(Bx, By, B) + br(Ax, Cy, B) + br(Cx, Ay, B) + 2*(br(Bx, Cy, A) + br(Cx, By, A))
    D1 = br(Bx, Cy, B) + br(Cx, By, B) + 2*br(Cx, Cy, A)
    D0 = br(Cx, Cy, B)
    return [expand(D) for D in (D5, D4, D3, D2, D1, D0)]

def detJ_of(A, B, C):
    F = A*z**2 + B*z + C
    return expand(F.jacobian([x, y, z]).det())

# ---------- counting machinery ----------
def qdim(G, nvars):
    polys = G.polys
    if len(polys) == 1 and polys[0].is_one:
        return 0
    lms = [max(pl.monoms(), key=grevlex) for pl in polys]
    bounds = []
    for i in range(nvars):
        b = None
        for m in lms:
            if all(m[j] == 0 for j in range(nvars) if j != i):
                b = m[i] if b is None else min(b, m[i])
        if b is None:
            return None  # not zero-dimensional
        bounds.append(b)
    cnt = 0
    for mono in itertools.product(*[range(b) for b in bounds]):
        if not any(all(mono[j] >= lm[j] for j in range(nvars)) for lm in lms):
            cnt += 1
    return cnt

def space_count(A, B, C, p, target=None):
    if target is None:
        target = [random.randrange(p) for _ in range(3)]
    F = A*z**2 + B*z + C
    polys = [expand(F[i] - target[i]) for i in range(3)]
    t0 = time.time()
    G = groebner(polys, x, y, z, modulus=p, order='grevlex')
    n = qdim(G, 3)
    return n, time.time() - t0

def plane_counts(A, B, C, p, target=None, k=None):
    """returns (raw, sat, k, secs): affine count of {g1=0, g2_k=0} and its u_k-saturation."""
    if target is None:
        target = [random.randrange(p) for _ in range(3)]
    PM = Matrix(3, 1, [Integer(v) for v in target])
    Ct = C - PM
    u = A.cross(B)
    a = A.cross(PM - C)
    b = B.cross(Ct)
    if k is None:
        for kk in (2, 0, 1):
            if not Poly(expand(u[kk]), x, y, modulus=p).is_zero:
                k = kk; break
        if k is None:
            return None, None, None, 0.0
    g1 = expand(Matrix.hstack(A, B, Ct).det())
    g2 = expand(a[k]**2 - u[k]*b[k])
    t0 = time.time()
    Graw = groebner([g1, g2], x, y, modulus=p, order='grevlex')
    raw = qdim(Graw, 2)
    Gsat = groebner([g1, g2, expand(1 - t*u[k])], x, y, t, modulus=p, order='grevlex')
    sat = qdim(Gsat, 3)
    return raw, sat, k, time.time() - t0

# ---------- random sampling ----------
def rpoly2(deg, p):
    return sum(random.randrange(p)*x**i*y**j for i in range(deg + 1) for j in range(deg + 1 - i))

def rvec(deg, p):
    return Matrix(3, 1, [rpoly2(deg, p) for _ in range(3)])

# ---------- linear algebra mod p ----------
def lin_solve_modp(rows, rhs, nunk, p):
    m = len(rows)
    M = [[rows[i].get(j, 0) % p for j in range(nunk)] + [rhs[i] % p] for i in range(m)]
    piv_cols = []
    r = 0
    for c in range(nunk):
        pr_ = None
        for rr in range(r, m):
            if M[rr][c] % p:
                pr_ = rr; break
        if pr_ is None:
            continue
        M[r], M[pr_] = M[pr_], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(v*inv) % p for v in M[r]]
        for rr in range(m):
            if rr != r and M[rr][c]:
                f = M[rr][c]
                M[rr] = [(M[rr][j] - f*M[r][j]) % p for j in range(nunk + 1)]
        piv_cols.append(c)
        r += 1
        if r == m:
            break
    for rr in range(m):
        if all(M[rr][j] % p == 0 for j in range(nunk)) and M[rr][nunk] % p:
            return None
    free = [c for c in range(nunk) if c not in piv_cols]
    part = [0]*nunk
    for i, c in enumerate(piv_cols):
        part[c] = M[i][nunk] % p
    basis = []
    for fc in free:
        v = [0]*nunk
        v[fc] = 1
        for i, c in enumerate(piv_cols):
            v[c] = (-M[i][fc]) % p
        basis.append(v)
    return part, basis, len(piv_cols)

def linear_conditions(expr, usyms, p):
    """expr must be affine-linear in usyms with integer (x,y)-poly coefficients.
       returns rows, rhs encoding [expr == 0 as poly in x,y]."""
    e = expand(expr)
    if e == 0:
        return [], []
    pe = Poly(e, x, y)
    rows, rhs = [], []
    for coefexpr in pe.coeffs():
        if not coefexpr.free_symbols & set(usyms):
            rows.append({}); rhs.append((-int(coefexpr)) % p)
            continue
        cp = Poly(coefexpr, *usyms)
        row = {}; const = 0
        for mono, cf in zip(cp.monoms(), cp.coeffs()):
            d = sum(mono)
            if d == 0:
                const = int(cf)
            elif d == 1:
                row[mono.index(1)] = (row.get(mono.index(1), 0) + int(cf)) % p
            else:
                raise ValueError('condition not linear in unknowns')
        rows.append(row); rhs.append((-const) % p)
    return rows, rhs

def generic_vec_syms(prefix, deg):
    n = (deg + 1)*(deg + 2)//2
    syms = symbols(prefix + '0:%d' % (3*n))
    monos = [x**i*y**j for i in range(deg + 1) for j in range(deg + 1 - i)]
    V = Matrix(3, 1, [sum(syms[c*n + m]*monos[m] for m in range(n)) for c in range(3)])
    return V, list(syms)

def random_solution(part, basis, p):
    vec = part[:]
    for bv in basis:
        r = random.randrange(p)
        vec = [(vec[i] + r*bv[i]) % p for i in range(len(vec))]
    return vec

def solve_D4_for_B(A, dB, p, forbid_const_B3=False, tries=6):
    Bm, bsyms = generic_vec_syms('b', dB)
    Z3 = Matrix(3, 1, [0, 0, 0])
    D = dsystem(A, Bm, Z3)
    rows, rhs = linear_conditions(D[1], bsyms, p)     # D4
    sol = lin_solve_modp(rows, rhs, len(bsyms), p)
    if sol is None:
        return None, None
    part, basis, rank = sol
    for _ in range(tries):
        vec = random_solution(part, basis, p)
        subsd = {bsyms[i]: vec[i] for i in range(len(bsyms))}
        Bs = Bm.subs(subsd)
        if forbid_const_B3 and Poly(expand(Bs[2].diff(x)), x, y, modulus=p).is_zero \
           and Poly(expand(Bs[2].diff(y)), x, y, modulus=p).is_zero:
            continue
        return Bs, rank
    return Bs, rank

def solve_lin_for_C(A, B, dC, p, which):
    """which: list of D-indices among [2]=D3,[3]=D2,[4]=D1 that are affine-linear in C."""
    Cm, csyms = generic_vec_syms('c', dC)
    D = dsystem(A, B, Cm)
    rows, rhs = [], []
    for w in which:
        r_, h_ = linear_conditions(D[w], csyms, p)
        rows += r_; rhs += h_
    sol = lin_solve_modp(rows, rhs, len(csyms), p)
    if sol is None:
        return None, None, None
    part, basis, rank = sol
    vec = random_solution(part, basis, p)
    subsd = {csyms[i]: vec[i] for i in range(len(csyms))}
    return Cm.subs(subsd), rank, len(csyms)

def assert_zero_modp(expr, p, label):
    assert Poly(expand(expr), x, y, modulus=p).is_zero, 'NOT ZERO mod p: ' + label

def polydeg_modp(expr, p):
    pp = Poly(expand(expr), x, y, modulus=p)
    return None if pp.is_zero else pp.total_degree()

# =====================================================================
say('='*78)
say('LEG 1: known 2-jet Keller maps -- exact field degrees and graded system')
say('='*78)

P1s, P2s, P3s = symbols('P1 P2 P3')

# T1 = (x + z^2, y, z)
A1 = Matrix(3, 1, [1, 0, 0]); B1 = Matrix(3, 1, [0, 0, 1]); C1 = Matrix(3, 1, [x, y, 0])
F1 = A1*z**2 + B1*z + C1
inv1 = [P1s - P3s**2, P2s, P3s]
comp = [expand(F1[i].subs({x: inv1[0], y: inv1[1], z: inv1[2]})) for i in range(3)]
assert comp == [P1s, P2s, P3s]
comp2 = [expand(e.subs({P1s: F1[0], P2s: F1[1], P3s: F1[2]})) for e in inv1]
assert comp2 == [x, y, z]
D_T1 = dsystem(A1, B1, C1)
assert D_T1[:5] == [0, 0, 0, 0, 0] and D_T1[5] == 1
assert expand(detJ_of(A1, B1, C1)) == 1
say('T1 = (x+z^2, y, z): polynomial inverse (P1-P3^2, P2, P3) verified both ways.')
say('    det J = 1, D5..D1 = 0, D0 = 1, box (dA,dB,dC) = (0,0,1). FIELD DEGREE = 1.')

# T2 = P2-refutation tame Keller map (THM-2446 S5)
A2 = Matrix(3, 1, [x, 1, 0]); B2 = Matrix(3, 1, [y, 0, 1])
C2 = Matrix(3, 1, [x*y - x**3 + x, 3*x**2 - 2*y, x])
F2 = A2*z**2 + B2*z + C2
assert expand(detJ_of(A2, B2, C2)) == -2
D_T2 = dsystem(A2, B2, C2)
assert D_T2[:5] == [0, 0, 0, 0, 0] and D_T2[5] == -2
Xi = P1s - (P3s**3 - P3s*P2s)/2
Zi = P3s - Xi
Yi = ((P3s - Xi)**2 + 3*Xi**2 - P2s)/2
inv2 = [Xi, Yi, Zi]
comp = [expand(F2[i].subs({x: inv2[0], y: inv2[1], z: inv2[2]}, simultaneous=True)) for i in range(3)]
assert comp == [P1s, P2s, P3s], comp
comp2 = [expand(e.subs({P1s: F2[0], P2s: F2[1], P3s: F2[2]}, simultaneous=True)) for e in inv2]
assert comp2 == [x, y, z], comp2
say('T2 = (xz^2+yz+xy-x^3+x, z^2+3x^2-2y, x+z): NEW EXPLICIT POLYNOMIAL INVERSE')
say('    x = P1 - (P3^3 - P2*P3)/2 ;  z = P3 - x ;  y = ((P3-x)^2 + 3x^2 - P2)/2')
say('    verified both compositions = identity. det J = -2, D5..D1 = 0, D0 = -2.')
say('    box (dA,dB,dC) = (1,1,3). FIELD DEGREE = 1.')

# W = THM-1310 wild z-affine map
u_ = 1 + x*y
AW = Matrix(3, 1, [0, 0, 0])
BW = Matrix(3, 1, [u_**3, 3*x*u_**2, -x**3])
CW = Matrix(3, 1, [y**2*u_*(4 + 3*x*y), y + 3*x*y**2*(4 + 3*x*y), 2*x - 3*x**2*y])
assert expand(detJ_of(AW, BW, CW)) == -2
D_W = dsystem(AW, BW, CW)
assert D_W[:5] == [0, 0, 0, 0, 0] and D_W[5] == -2
assert expand(BW[1]**3 + 27*BW[0]**2*BW[2]) == 0   # cuspidal-cubic cone (THM-2446 / THM-1340)
n1, s1 = space_count(AW, BW, CW, PR)
n2, s2 = space_count(AW, BW, CW, PR2)
say('W = THM-1310 map: det J = -2, D-system (0,0,0,0,0,-2), B on cusp cone 27B1^2B3+B2^3=0.')
say('    fiber count mod %d: %s (%.1fs), mod %d: %s (%.1fs)  -- expected 3' % (PR, n1, s1, PR2, n2, s2))
assert n1 == 3 and n2 == 3
say('    D0(T2) = D0(W) = -2 with field degrees 1 vs 3: D0 carries NO field-degree info.')

# B nowhere zero / C immersion consequences of D0 = const (Nullstellensatz certs)
GB = groebner([BW[0], BW[1], BW[2]], x, y, order='grevlex')
assert GB.polys[0].is_one or any(pl.is_one for pl in GB.polys)
C2x, C2y = C2.diff(x), C2.diff(y)
w2 = C2x.cross(C2y)
GC = groebner([w2[0], w2[1], w2[2]], x, y, order='grevlex')
assert any(pl.is_one for pl in GC.polys)
say('D0-consequences verified: ideal(B_W) = (1) [B nowhere 0]; ideal(C2_x x C2_y) = (1)')
say('    [C is an everywhere-immersion for T2]. Both are forced by D0 = const != 0.')

say('')
say('='*78)
say('LEG 2: the projected plane fiber system (frame identities, proof-grade)')
say('='*78)
a1s, a2s, a3s, b1s, b2s, b3s, c1s, c2s, c3s = symbols('a1 a2 a3 b1 b2 b3 c1 c2 c3')
Av = Matrix(3, 1, [a1s, a2s, a3s]); Bv = Matrix(3, 1, [b1s, b2s, b3s]); Cv = Matrix(3, 1, [c1s, c2s, c3s])
E1 = Av[0]*z**2 + Bv[0]*z + Cv[0]
E2 = Av[1]*z**2 + Bv[1]*z + Cv[1]
E3 = Av[2]*z**2 + Bv[2]*z + Cv[2]
uv = Av.cross(Bv); av = -Av.cross(Cv); bv = Bv.cross(Cv)   # Cv plays C~ = C - P
r12 = expand(resultant(E1, E2, z))
assert expand(r12 - (av[2]**2 - uv[2]*bv[2])) == 0
r13 = expand(resultant(E1, E3, z))
assert expand(r13 - (av[1]**2 - uv[1]*bv[1])) == 0
r23 = expand(resultant(E2, E3, z))
assert expand(r23 - (av[0]**2 - uv[0]*bv[0])) == 0
say('Res_z(E_i, E_j) = a_k^2 - u_k b_k  for (i,j,k) cyclic  -- EXACT, 9 generic symbols.')
assert expand(uv.dot(Av)) == 0 and expand(uv.dot(Bv)) == 0
assert expand(uv.dot(Cv) - Matrix.hstack(Av, Bv, Cv).det()) == 0
say('u.A = u.B = 0 and u.C~ = det[A|B|C~] = g1: u = AxB is the LEFT KERNEL of the')
say('fiber matrix on {g1=0}; u_k != 0 makes E_k redundant, so the fiber bijects with')
say('{g1 = 0, g2_k = 0, u_k != 0} via z = a_k/u_k.  (Lemma proof in report.)')

say('')
say('='*78)
say('LEG 3: engine cross-validation (space vs projected-plane counting, mod p)')
say('='*78)
for box in [(0, 1, 2), (1, 2, 3)]:
    dA, dB, dC = box
    for inst in range(2):
        A = rvec(dA, PR); B = rvec(dB, PR); C = rvec(dC, PR)
        target = [random.randrange(PR) for _ in range(3)]
        ns, ts = space_count(A, B, C, PR, target)
        raw, sat, k, tp = plane_counts(A, B, C, PR, target)
        say('box %s inst %d: space=%s (%.1fs)  plane raw=%s sat=%s (chart k=%d, %.1fs)'
            % (box, inst, ns, ts, raw, sat, k, tp))
        assert ns == sat, 'plane saturated count must equal space count'
n_t2, _ = space_count(A2, B2, C2, PR)
assert n_t2 == 1
say('T2 control through engine: fiber count = 1 OK')

say('')
say('='*78)
say('LEG 4: TABLE A -- generic (non-Keller) field degree per box vs BKK')
say('='*78)
say('BKK/hull formula (proved in report): N_max(dA,dB,dC) = (dC^2+dC*b+b^2)+(b^2+b*dA+dA^2),')
say('b = max(dB,(dA+dC)/2); staircase (d,d+1,d+2): N_max = (d+2)^3 - d^3 = 6d^2+12d+8.')
predictions = {(0, 1, 2): 8, (1, 2, 3): 26, (2, 3, 4): 56, (1, 1, 3): 26}
for box in [(0, 1, 2), (1, 2, 3), (1, 1, 3)]:
    dA, dB, dC = box
    vals = []
    for inst in range(3):
        A = rvec(dA, PR); B = rvec(dB, PR); C = rvec(dC, PR)
        n_, t_ = space_count(A, B, C, PR)
        vals.append(n_)
    # one cross-prime check on the last instance
    n2_, t2_ = space_count(A, B, C, PR2)
    say('box %s: measured %s (+ mod %d: %s), predicted N_max = %d'
        % (box, vals, PR2, n2_, predictions[box]))
    assert max(vals) == predictions[box]
# (2,3,4): plane-saturated engine (validated in LEG 3); one space attempt, time-boxed by fate
box = (2, 3, 4)
vals = []
for inst in range(2):
    A = rvec(2, PR); B = rvec(3, PR); C = rvec(4, PR)
    raw, sat, k, tp = plane_counts(A, B, C, PR)
    say('box (2,3,4) inst %d: plane raw=%s sat=%s (chart k=%s, %.1fs), predicted 56'
        % (inst, raw, sat, k, tp))
    vals.append(sat)
say('box (2,3,4): measured %s, predicted N_max = 56' % (vals,))
assert max(vals) == predictions[box]
say('Projected-system Bezout bookkeeping per staircase box: deg g1 = 3d+3, deg g2 = 4d+4,')
say('product 12(d+1)^2 = 12/48/108 vs true 8/26/56: excess 6d^2+12d+4 sits entirely at')
say('plane-infinity (raw affine counts above already equal N on random instances).')

say('')
say('='*78)
say('LEG 6: reachability in the MINIMAL box (0,1,2) (all total-degree-2 maps)')
say('='*78)
witnesses = {
    1: (Matrix(3,1,[1,0,0]), Matrix(3,1,[0,0,1]), Matrix(3,1,[x,y,0]),      '(x+z^2, y, z)'),
    2: (Matrix(3,1,[0,0,0]), Matrix(3,1,[0,0,1]), Matrix(3,1,[x**2,y,0]),   '(x^2, y, z)'),
    3: (Matrix(3,1,[0,1,0]), Matrix(3,1,[0,0,y]), Matrix(3,1,[x,y,0]),      '(x, z^2+y, yz)'),
    4: (Matrix(3,1,[0,0,1]), Matrix(3,1,[0,0,0]), Matrix(3,1,[x**2,y,0]),   '(x^2, y, z^2)'),
}
for N, (A, B, C, name) in sorted(witnesses.items()):
    n_, _ = space_count(A, B, C, PR)
    n2_, _ = space_count(A, B, C, PR)
    say('  N = %d witness %-22s measured %s,%s' % (N, name, n_, n2_))
    assert n_ == N and n2_ == N
say('  N = 8: generic (Table A).  Searching coefficient range {-1,0,1} for 5,6,7 ...')
found = {}
tries = 0
while len(found) < 8 and tries < 220:
    tries += 1
    A = Matrix(3, 1, [random.randint(-1, 1) for _ in range(3)])
    B = Matrix(3, 1, [sum(random.randint(-1, 1)*m for m in (1, x, y)) for _ in range(3)])
    C = Matrix(3, 1, [sum(random.randint(-1, 1)*m for m in (1, x, y, x**2, x*y, y**2)) for _ in range(3)])
    n_, _ = space_count(A, B, C, PR)
    if n_ is None or n_ == 0 or n_ in found:
        continue
    n2_, _ = space_count(A, B, C, PR)
    if n2_ != n_:
        continue
    found[n_] = (list(A), list(B), list(C))
say('  observed degree spectrum in %d small-coefficient samples: %s' % (tries, sorted(found)))
for N in sorted(found):
    if N in (5, 6, 7):
        A_, B_, C_ = found[N]
        say('    N=%d witness: A=%s B=%s C=%s' % (N, A_, B_, C_))

say('')
say('--- the quartic S4 witness Q4 = (x, z^2+y^2+y, yz+x), box (0,1,2), NON-Keller ---')
AQ = Matrix(3, 1, [0, 1, 0]); BQ = Matrix(3, 1, [0, 0, y]); CQ = Matrix(3, 1, [x, y**2 + y, x])
nQ, _ = space_count(AQ, BQ, CQ, PR)
say('  fiber count = %s (expected 4)' % nQ)
assert nQ == 4
q_, s_ = symbols('q_ s_')
elim = expand(resultant(z**2 + y**2 + y - q_, y*z + s_, y))
say('  z-eliminant (x=P1, s=P3-P1): %s' % elim)
assert expand(elim - (z**4 - q_*z**2 + s_*z + s_**2)) == 0
say('  = z^4 - q z^2 + s z + s^2 : DEPRESSED quartic (no z^3 -- structural for 2-jet route).')
pdep, qdep, rdep = -q_, s_, s_**2
Iq = pdep**2 + 12*rdep
Jq = 2*pdep**3 + 27*qdep**2 - 72*pdep*rdep
Delta4 = expand((4*Iq**3 - Jq**2)/27)
say('  I = %s,  J = %s' % (expand(Iq), expand(Jq)))
say('  Delta4 = %s' % factor(Delta4))
cyc = Counter(); tot = 0
t0 = time.time()
while tot < 2500:
    qq = random.randrange(PR); ss = random.randrange(1, PR)
    f = Poly(z**4 - qq*z**2 + ss*z + ss**2, z, modulus=PR)
    if f.discriminant() % PR == 0:
        continue
    fl = f.factor_list()[1]
    part = tuple(sorted(sum([[g.degree()]*m for g, m in fl], [])))
    cyc[part] += 1; tot += 1
say('  cycle statistics over %d unramified specializations mod %d (%.0fs):' % (tot, PR, time.time() - t0))
th_S4 = {(1,1,1,1): 1/24, (1,1,2): 6/24, (2,2): 3/24, (1,3): 8/24, (4,): 6/24}
for part in sorted(cyc, key=lambda s: (len(s), s)):
    say('    %-12s observed %.4f   S4-Chebotarev %.4f' % (part, cyc[part]/tot, th_S4.get(part, 0.0)))
say('  A4 predicts no (1,1,2)/(4,); V4 predicts only (1,1,1,1)/(2,2); D4/C4 have no (1,3)')
say('  at rate 1/3.  Verdict: monodromy S4.  ==> degree 4 with FULL S4 is reachable by')
say('  the 2-jet SHAPE already in the minimal box; only KELLER kills it there (Wang).')

say('')
say('='*78)
say('LEG 5: TABLE C/E -- Keller-shape strata of box (1,2,3) and the D-staircase')
say('='*78)
say('Classification (proved in report): D5=0 in box (1,2,3) <=> rank{A0,A1,A2} <= 2.')
say('r=2 (unique lam, lam.A=0): T (lam.B=const!=0, KELLER=>N=1 via Moh), H1 (lam.B nonconst),')
say('H0 (lam.B=0); r=1: A=ell*v; r=0: z-affine (B<=2, C<=3) with plane- or conic-valued B.')

def stage_measure(label, A, B, C, p=PR, ntar=2):
    vals = []
    for _ in range(ntar):
        n_, t_ = space_count(A, B, C, p)
        vals.append(n_)
    d1 = polydeg_modp(dsystem(A, B, C)[4], p)
    d0 = polydeg_modp(dsystem(A, B, C)[5], p)
    say('  %-34s N = %-12s deg D1 = %-4s deg D0 = %s' % (label, vals, d1, d0))
    return vals

# --- T-shape ceiling (Keller forces 1 there; ceiling is the non-Keller MV)
A = Matrix(3, 1, [rpoly2(1, PR), rpoly2(1, PR), Integer(0)])
B = Matrix(3, 1, [rpoly2(2, PR), rpoly2(2, PR), Integer(random.randrange(1, PR))])
C = rvec(3, PR)
stage_measure('T  (A3=0, B3 = const)', A, B, C)

# --- H1 shape and staircase
say('H1 (A3=0 affine A; B3 nonconstant):')
for inst in range(2):
    A = Matrix(3, 1, [rpoly2(1, PR), rpoly2(1, PR), Integer(0)])
    B = rvec(2, PR); C = rvec(3, PR)
    assert_zero_modp(dsystem(A, B, C)[0], PR, 'D5 H1')
    stage_measure('H1 stage D5=0 (inst %d)' % inst, A, B, C)
A = Matrix(3, 1, [rpoly2(1, PR), rpoly2(1, PR), Integer(0)])
B4, rankB = solve_D4_for_B(A, 2, PR, forbid_const_B3=True)
say('  [D4 linear system on 18 B-coeffs: rank %s]' % rankB)
C = rvec(3, PR)
assert_zero_modp(dsystem(A, B4, C)[1], PR, 'D4')
stage_measure('H1 + D4=0', A, B4, C)
C32, rankC, ncoef = solve_lin_for_C(A, B4, 3, PR, which=[2, 3])
if C32 is None:
    say('  H1 + D3=D2=0: affine system INCONSISTENT for this (A,B) -- B must be special')
else:
    say('  [D3,D2 affine-linear system on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(A, B4, C32)[2], PR, 'D3')
    assert_zero_modp(dsystem(A, B4, C32)[3], PR, 'D2')
    stage_measure('H1 + D4=D3=D2=0', A, B4, C32)

# --- H0 shape
say('H0 (A3=0, B3=0: third row z-free):')
for inst in range(2):
    A = Matrix(3, 1, [rpoly2(1, PR), rpoly2(1, PR), Integer(0)])
    B = Matrix(3, 1, [rpoly2(2, PR), rpoly2(2, PR), Integer(0)])
    C = rvec(3, PR)
    stage_measure('H0 stage D5=0 (inst %d)' % inst, A, B, C)

# --- r=1
say('r=1 (A = ell*v, ell affine):')
v = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
ell = rpoly2(1, PR)
A = v*ell
B = rvec(2, PR); C = rvec(3, PR)
assert_zero_modp(dsystem(A, B, C)[0], PR, 'D5 r1')
stage_measure('r1 stage D5=0', A, B, C)

# --- r=0 z-affine strata
say('r=0 (A=0, z-affine, dB<=2, dC<=3):')
A0 = Matrix(3, 1, [0, 0, 0])
B = rvec(2, PR); C = rvec(3, PR)
stage_measure('r0 generic (no conditions)', A0, B, C)
# plane-valued B: [Bx,By,B] = 0 automatic
v1 = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
v2 = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
Bpl = v1*rpoly2(2, PR) + v2*rpoly2(2, PR)
assert_zero_modp(dsystem(A0, Bpl, C)[3], PR, 'D2 zaff plane')  # D2 = [Bx,By,B]
stage_measure('r0 + plane-valued B (D2=0)', A0, Bpl, C)
Cd1, rankC, ncoef = solve_lin_for_C(A0, Bpl, 3, PR, which=[4])   # D1 linear in C when A=0
if Cd1 is not None:
    say('  [z-affine D1 linear system on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(A0, Bpl, Cd1)[4], PR, 'D1 zaff plane')
    stage_measure('r0 plane-B + D1=0', A0, Bpl, Cd1)
# conic-valued B (Veronese cap): [Bx,By,B] = 0 by Euler
f_ = rpoly2(1, PR); g_ = rpoly2(1, PR)
M = Matrix(3, 3, [random.randrange(PR) for _ in range(9)])
Bcon = M*Matrix(3, 1, [f_**2, f_*g_, g_**2])
assert_zero_modp(dsystem(A0, Bcon, C)[3], PR, 'D2 zaff conic')
stage_measure('r0 + conic-valued B (D2=0)', A0, Bcon, C)
Cd1, rankC, ncoef = solve_lin_for_C(A0, Bcon, 3, PR, which=[4])
if Cd1 is not None:
    say('  [z-affine D1 linear system on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(A0, Bcon, Cd1)[4], PR, 'D1 zaff conic')
    stage_measure('r0 conic-B + D1=0', A0, Bcon, Cd1)

say('')
say('--- box (2,3,4): the conic-cap door (THM-2446 S3) ---')
f_ = rpoly2(1, PR); g_ = rpoly2(1, PR)
M = Matrix(3, 3, [random.randrange(PR) for _ in range(9)])
Acon = M*Matrix(3, 1, [f_**2, f_*g_, g_**2])
Bfree = rvec(3, PR); Cfree = rvec(4, PR)
assert_zero_modp(dsystem(Acon, Bfree, Cfree)[0], PR, 'D5 conic cap 234')
raw, sat, k, tp = plane_counts(Acon, Bfree, Cfree, PR)
say('  conic-cap A, free B,C: plane raw=%s sat=%s (%.1fs)' % (raw, sat, tp))
B4c, rankB = solve_D4_for_B(Acon, 3, PR)
say('  [D4 linear system on 30 B-coeffs: rank %s]' % rankB)
assert_zero_modp(dsystem(Acon, B4c, Cfree)[1], PR, 'D4 conic cap')
raw2, sat2, k2, tp2 = plane_counts(Acon, B4c, Cfree, PR)
say('  conic-cap A + D4=0 B: plane raw=%s sat=%s (%.1fs)' % (raw2, sat2, tp2))
C32c, rankC, ncoef = solve_lin_for_C(Acon, B4c, 4, PR, which=[2, 3])
if C32c is None:
    say('  + D3=D2=0 on C: affine system INCONSISTENT for random D4-kernel B')
    say('    ==> in box (2,3,4) the staircase is NOT layer-by-layer solvable: B and C')
    say('        must be solved JOINTLY below D4 (hunt-design constraint).')
else:
    say('  [D3,D2 affine-linear system on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(Acon, B4c, C32c)[2], PR, 'D3 conic cap')
    assert_zero_modp(dsystem(Acon, B4c, C32c)[3], PR, 'D2 conic cap')
    raw3, sat3, k3, tp3 = plane_counts(Acon, B4c, C32c, PR)
    say('  conic-cap A + D4=D3=D2=0: plane raw=%s sat=%s (%.1fs)' % (raw3, sat3, tp3))

say('')
say('ALL CHECKS PASSED')
