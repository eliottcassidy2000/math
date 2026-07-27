# -*- coding: utf-8 -*-
# G1 groundwork PART 2 (parallel engine): strata tables with SPACE Groebner counts.
import random, time, itertools, sys
from collections import Counter
from sympy import symbols, Matrix, Poly, groebner, expand, factor, resultant, Integer
from sympy.polys.orderings import grevlex

x, y, z, t = symbols('x y z t')
random.seed(917331)
PR  = 10007
PR2 = 31013

def say(*a):
    print(*a); sys.stdout.flush()

def br(u, v, w): return Matrix.hstack(u, v, w).det()

def dsystem(A, B, C):
    Ax, Ay = A.diff(x), A.diff(y); Bx, By = B.diff(x), B.diff(y); Cx, Cy = C.diff(x), C.diff(y)
    D5 = 2*br(Ax, Ay, A)
    D4 = br(Ax, Ay, B) + 2*br(Ax, By, A) + 2*br(Bx, Ay, A)
    D3 = br(Ax, By, B) + br(Bx, Ay, B) + 2*(br(Ax, Cy, A) + br(Cx, Ay, A) + br(Bx, By, A))
    D2 = br(Bx, By, B) + br(Ax, Cy, B) + br(Cx, Ay, B) + 2*(br(Bx, Cy, A) + br(Cx, By, A))
    D1 = br(Bx, Cy, B) + br(Cx, By, B) + 2*br(Cx, Cy, A)
    D0 = br(Cx, Cy, B)
    return [expand(D) for D in (D5, D4, D3, D2, D1, D0)]

def qdim(G, nvars):
    polys = G.polys
    if any(p.is_one for p in polys): return 0
    lms = [max(pl.monoms(), key=grevlex) for pl in polys]
    bounds = []
    for i in range(nvars):
        b = None
        for m in lms:
            if all(m[j] == 0 for j in range(nvars) if j != i):
                b = m[i] if b is None else min(b, m[i])
        if b is None: return None
        bounds.append(b)
    cnt = 0
    for mono in itertools.product(*[range(b) for b in bounds]):
        if not any(all(mono[j] >= lm[j] for j in range(nvars)) for lm in lms):
            cnt += 1
    return cnt

def space_count(A, B, C, p, target=None):
    if target is None: target = [random.randrange(p) for _ in range(3)]
    F = A*z**2 + B*z + C
    t0 = time.time()
    G = groebner([expand(F[i] - target[i]) for i in range(3)], x, y, z, modulus=p, order='grevlex')
    return qdim(G, 3), time.time() - t0

def rpoly2(deg, p): return sum(random.randrange(p)*x**i*y**j for i in range(deg+1) for j in range(deg+1-i))
def rvec(deg, p):  return Matrix(3, 1, [rpoly2(deg, p) for _ in range(3)])

def lin_solve_modp(rows, rhs, nunk, p):
    m = len(rows)
    M = [[rows[i].get(j, 0) % p for j in range(nunk)] + [rhs[i] % p] for i in range(m)]
    piv_cols = []; r = 0
    for c in range(nunk):
        pr_ = None
        for rr in range(r, m):
            if M[rr][c] % p: pr_ = rr; break
        if pr_ is None: continue
        M[r], M[pr_] = M[pr_], M[r]
        inv = pow(M[r][c], p-2, p)
        M[r] = [(v*inv) % p for v in M[r]]
        for rr in range(m):
            if rr != r and M[rr][c]:
                f = M[rr][c]
                M[rr] = [(M[rr][j]-f*M[r][j]) % p for j in range(nunk+1)]
        piv_cols.append(c); r += 1
        if r == m: break
    for rr in range(m):
        if all(M[rr][j] % p == 0 for j in range(nunk)) and M[rr][nunk] % p:
            return None
    free = [c for c in range(nunk) if c not in piv_cols]
    part = [0]*nunk
    for i, c in enumerate(piv_cols): part[c] = M[i][nunk] % p
    basis = []
    for fc in free:
        v = [0]*nunk; v[fc] = 1
        for i, c in enumerate(piv_cols): v[c] = (-M[i][fc]) % p
        basis.append(v)
    return part, basis, len(piv_cols)

def linear_conditions(expr, usyms, p):
    e = expand(expr)
    if e == 0: return [], []
    pe = Poly(e, x, y)
    rows, rhs = [], []
    for coefexpr in pe.coeffs():
        if not getattr(coefexpr, 'free_symbols', set()) & set(usyms):
            rows.append({}); rhs.append((-int(coefexpr)) % p); continue
        cp = Poly(coefexpr, *usyms)
        row = {}; const = 0
        for mono, cf in zip(cp.monoms(), cp.coeffs()):
            d = sum(mono)
            if d == 0: const = int(cf)
            elif d == 1: row[mono.index(1)] = (row.get(mono.index(1), 0)+int(cf)) % p
            else: raise ValueError('not linear')
        rows.append(row); rhs.append((-const) % p)
    return rows, rhs

def generic_vec_syms(prefix, deg):
    n = (deg+1)*(deg+2)//2
    syms = symbols(prefix + '0:%d' % (3*n))
    monos = [x**i*y**j for i in range(deg+1) for j in range(deg+1-i)]
    V = Matrix(3, 1, [sum(syms[c*n+m]*monos[m] for m in range(n)) for c in range(3)])
    return V, list(syms)

def random_solution(part, basis, p):
    vec = part[:]
    for bv in basis:
        r = random.randrange(p)
        vec = [(vec[i]+r*bv[i]) % p for i in range(len(vec))]
    return vec

def solve_D4_for_B(A, dB, p, forbid_const_B3=False, tries=6):
    Bm, bsyms = generic_vec_syms('b', dB)
    D = dsystem(A, Bm, Matrix(3, 1, [0, 0, 0]))
    rows, rhs = linear_conditions(D[1], bsyms, p)
    sol = lin_solve_modp(rows, rhs, len(bsyms), p)
    if sol is None: return None, None
    part, basis, rank = sol
    Bs = None
    for _ in range(tries):
        vec = random_solution(part, basis, p)
        subsd = {bsyms[i]: vec[i] for i in range(len(bsyms))}
        Bs = Bm.subs(subsd)
        if forbid_const_B3 and Poly(expand(Bs[2].diff(x)), x, y, modulus=p).is_zero \
           and Poly(expand(Bs[2].diff(y)), x, y, modulus=p).is_zero:
            continue
        break
    return Bs, rank

def solve_lin_for_C(A, B, dC, p, which):
    Cm, csyms = generic_vec_syms('c', dC)
    D = dsystem(A, B, Cm)
    rows, rhs = [], []
    for w in which:
        r_, h_ = linear_conditions(D[w], csyms, p)
        rows += r_; rhs += h_
    sol = lin_solve_modp(rows, rhs, len(csyms), p)
    if sol is None: return None, None, None
    part, basis, rank = sol
    vec = random_solution(part, basis, p)
    return Cm.subs({csyms[i]: vec[i] for i in range(len(csyms))}), rank, len(csyms)

def assert_zero_modp(expr, p, label):
    assert Poly(expand(expr), x, y, modulus=p).is_zero, 'NOT ZERO mod p: '+label

def polydeg_modp(expr, p):
    pp = Poly(expand(expr), x, y, modulus=p)
    return None if pp.is_zero else pp.total_degree()

def stage(label, A, B, C, p=PR, ntar=2):
    vals = []; secs = 0.
    for _ in range(ntar):
        n_, t_ = space_count(A, B, C, p); vals.append(n_); secs += t_
    D = dsystem(A, B, C)
    d1 = polydeg_modp(D[4], p); d0 = polydeg_modp(D[5], p)
    say('  %-36s N = %-10s deg(D1,D0) = (%s,%s)  [%.0fs]' % (label, vals, d1, d0, secs))
    return vals

say('#'*78)
say('PART 2 (space engine): box (2,3,4) generic + strata tables + witnesses')
say('#'*78)

say('--- box (2,3,4) generic (space Groebner) ---')
for inst in range(2):
    A = rvec(2, PR); B = rvec(3, PR); C = rvec(4, PR)
    n_, t_ = space_count(A, B, C, PR)
    say('  inst %d: N = %s (%.0fs), predicted 56' % (inst, n_, t_))
A = rvec(2, PR2); B = rvec(3, PR2); C = rvec(4, PR2)
n_, t_ = space_count(A, B, C, PR2)
say('  cross-prime %d: N = %s (%.0fs)' % (PR2, n_, t_))

say('')
say('--- box (1,2,3) shape strata (D5=0 classification) ---')
# T-shape
A = Matrix(3, 1, [rpoly2(1, PR), rpoly2(1, PR), Integer(0)])
B = Matrix(3, 1, [rpoly2(2, PR), rpoly2(2, PR), Integer(random.randrange(1, PR))])
C = rvec(3, PR)
stage('T  (A3=0, B3 const; Keller=>N=1)', A, B, C)
# H1
for inst in range(2):
    A = Matrix(3, 1, [rpoly2(1, PR), rpoly2(1, PR), Integer(0)])
    B = rvec(2, PR); C = rvec(3, PR)
    assert_zero_modp(dsystem(A, B, C)[0], PR, 'D5 H1')
    stage('H1 (A3=0, B3 nonconst) inst %d' % inst, A, B, C)
A = Matrix(3, 1, [rpoly2(1, PR), rpoly2(1, PR), Integer(0)])
B4, rankB = solve_D4_for_B(A, 2, PR, forbid_const_B3=True)
say('  [H1 D4-linear system on 18 B-coeffs: rank %s]' % rankB)
C = rvec(3, PR)
assert_zero_modp(dsystem(A, B4, C)[1], PR, 'D4 H1')
stage('H1 + D4=0', A, B4, C)
C32, rankC, ncoef = solve_lin_for_C(A, B4, 3, PR, which=[2, 3])
if C32 is None:
    say('  H1 + D3=D2=0: affine system INCONSISTENT (B must be special)')
else:
    say('  [H1 D3&D2 affine-linear on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(A, B4, C32)[2], PR, 'D3 H1')
    assert_zero_modp(dsystem(A, B4, C32)[3], PR, 'D2 H1')
    stage('H1 + D4=D3=D2=0', A, B4, C32)
# H0
for inst in range(2):
    A = Matrix(3, 1, [rpoly2(1, PR), rpoly2(1, PR), Integer(0)])
    B = Matrix(3, 1, [rpoly2(2, PR), rpoly2(2, PR), Integer(0)])
    C = rvec(3, PR)
    stage('H0 (A3=B3=0, row 3 z-free) inst %d' % inst, A, B, C)
# r=1
v = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
A = v*rpoly2(1, PR)
B = rvec(2, PR); C = rvec(3, PR)
assert_zero_modp(dsystem(A, B, C)[0], PR, 'D5 r1')
stage('r1 (A = ell*v)', A, B, C)
# r=0 z-affine
A0 = Matrix(3, 1, [0, 0, 0])
B = rvec(2, PR); C = rvec(3, PR)
stage('r0 z-affine generic (B<=2,C<=3)', A0, B, C)
v1 = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
v2 = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
Bpl = v1*rpoly2(2, PR) + v2*rpoly2(2, PR)
assert_zero_modp(dsystem(A0, Bpl, C)[3], PR, 'D2 r0 plane')
stage('r0 + plane-valued B ([BxByB]=0)', A0, Bpl, C)
Cd1, rankC, ncoef = solve_lin_for_C(A0, Bpl, 3, PR, which=[4])
if Cd1 is not None:
    say('  [r0-plane D1-linear on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(A0, Bpl, Cd1)[4], PR, 'D1 r0 plane')
    stage('r0 plane-B + D1=0', A0, Bpl, Cd1)
f_ = rpoly2(1, PR); g_ = rpoly2(1, PR)
M = Matrix(3, 3, [random.randrange(PR) for _ in range(9)])
Bcon = M*Matrix(3, 1, [f_**2, f_*g_, g_**2])
assert_zero_modp(dsystem(A0, Bcon, C)[3], PR, 'D2 r0 conic')
stage('r0 + conic(Veronese)-valued B', A0, Bcon, C)
Cd1, rankC, ncoef = solve_lin_for_C(A0, Bcon, 3, PR, which=[4])
if Cd1 is not None:
    say('  [r0-conic D1-linear on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(A0, Bcon, Cd1)[4], PR, 'D1 r0 conic')
    stage('r0 conic-B + D1=0', A0, Bcon, Cd1)

say('')
say('--- box (2,3,4) doors ---')
# conic-cap A (THM-2446 d=4 door)
f_ = rpoly2(1, PR); g_ = rpoly2(1, PR)
M = Matrix(3, 3, [random.randrange(PR) for _ in range(9)])
Acon = M*Matrix(3, 1, [f_**2, f_*g_, g_**2])
Bfree = rvec(3, PR); Cfree = rvec(4, PR)
assert_zero_modp(dsystem(Acon, Bfree, Cfree)[0], PR, 'D5 conic cap')
stage('conic-cap A + free B,C', Acon, Bfree, Cfree, ntar=1)
B4c, rankB = solve_D4_for_B(Acon, 3, PR)
say('  [conic-cap D4-linear on 30 B-coeffs: rank %s]' % rankB)
assert_zero_modp(dsystem(Acon, B4c, Cfree)[1], PR, 'D4 conic cap')
stage('conic-cap A + D4=0 B', Acon, B4c, Cfree, ntar=1)
C32c, rankC, ncoef = solve_lin_for_C(Acon, B4c, 4, PR, which=[2, 3])
if C32c is None:
    say('  conic-cap + D3=D2=0 on C: INCONSISTENT for random D4-kernel B')
    say('  ==> B and C must be solved JOINTLY below D4 in box (2,3,4).')
else:
    say('  [conic-cap D3&D2 affine-linear on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(Acon, B4c, C32c)[2], PR, 'D3 cc')
    assert_zero_modp(dsystem(Acon, B4c, C32c)[3], PR, 'D2 cc')
    stage('conic-cap + D4=D3=D2=0', Acon, B4c, C32c, ntar=1)
# z-affine (dB,dC)=(3,4): generic, then cusp engine with AFFINE lambda,mu
A0 = Matrix(3, 1, [0, 0, 0])
B = rvec(3, PR); C = rvec(4, PR)
stage('r0(2,3,4): z-affine generic (3,4)', A0, B, C, ntar=1)
lam = rpoly2(1, PR); mu = rpoly2(1, PR)
Bcusp = Matrix(3, 1, [lam**3, 3*lam**2*mu, -mu**3])
assert expand(Bcusp[1]**3 + 27*Bcusp[0]**2*Bcusp[2]) == 0
assert_zero_modp(dsystem(A0, Bcusp, C)[3], PR, 'D2 cusp')
stage('r0(2,3,4): CUSP-ENGINE B, free C', A0, Bcusp, C, ntar=2)
Cd1, rankC, ncoef = solve_lin_for_C(A0, Bcusp, 4, PR, which=[4])
if Cd1 is not None:
    say('  [cusp-engine D1-linear on %d C-coeffs: rank %s]' % (ncoef, rankC))
    assert_zero_modp(dsystem(A0, Bcusp, Cd1)[4], PR, 'D1 cusp')
    stage('r0(2,3,4) cusp-B + D1=0', A0, Bcusp, Cd1, ntar=2)

say('')
say('--- box (0,1,2) witnesses + quartic S4 stats (early copy) ---')
wit = [
    (1, Matrix(3,1,[1,0,0]), Matrix(3,1,[0,0,1]), Matrix(3,1,[x,y,0]),    '(x+z^2, y, z)'),
    (2, Matrix(3,1,[0,0,0]), Matrix(3,1,[0,0,1]), Matrix(3,1,[x**2,y,0]), '(x^2, y, z)'),
    (3, Matrix(3,1,[0,1,0]), Matrix(3,1,[0,0,y]), Matrix(3,1,[x,y,0]),    '(x, z^2+y, yz)'),
    (4, Matrix(3,1,[0,0,1]), Matrix(3,1,[0,0,0]), Matrix(3,1,[x**2,y,0]), '(x^2, y, z^2)'),
    (4, Matrix(3,1,[0,1,0]), Matrix(3,1,[0,0,y]), Matrix(3,1,[x,y**2+y,x]),'Q4=(x, z^2+y^2+y, yz+x)'),
]
for N, A, B, C, name in wit:
    n_, _ = space_count(A, B, C, PR)
    say('  N=%d  %-26s measured %s' % (N, name, n_))
    assert n_ == N
found = {}
tries = 0
while len(found) < 8 and tries < 300:
    tries += 1
    A = Matrix(3, 1, [random.randint(-1, 1) for _ in range(3)])
    B = Matrix(3, 1, [sum(random.randint(-1, 1)*m for m in (1, x, y)) for _ in range(3)])
    C = Matrix(3, 1, [sum(random.randint(-1, 1)*m for m in (1, x, y, x**2, x*y, y**2)) for _ in range(3)])
    n_, _ = space_count(A, B, C, PR)
    if n_ is None or n_ == 0 or n_ in found: continue
    n2_, _ = space_count(A, B, C, PR)
    if n2_ != n_: continue
    found[n_] = (list(A), list(B), list(C))
say('  observed spectrum in %d samples (coeffs in {-1,0,1}): %s' % (tries, sorted(found)))
for N in sorted(found):
    A_, B_, C_ = found[N]
    say('    N=%d witness: A=%s B=%s C=%s' % (N, A_, B_, C_))

cyc = Counter(); tot = 0
while tot < 2500:
    qq = random.randrange(PR); ss = random.randrange(1, PR)
    f = Poly(z**4 - qq*z**2 + ss*z + ss**2, z, modulus=PR)
    if f.discriminant() % PR == 0: continue
    fl = f.factor_list()[1]
    part = tuple(sorted(sum([[g.degree()]*m for g, m in fl], [])))
    cyc[part] += 1; tot += 1
th = {(1,1,1,1): 1/24, (1,1,2): 6/24, (2,2): 3/24, (1,3): 8/24, (4,): 6/24}
say('  Q4 eliminant z^4 - q z^2 + s z + s^2 cycle stats over %d samples:' % tot)
for part in sorted(cyc, key=lambda s: (len(s), s)):
    say('    %-12s obs %.4f   S4 %.4f' % (part, cyc[part]/tot, th.get(part, 0.)))
q_, s_ = symbols('q_ s_')
pdep, qdep, rdep = -q_, s_, s_**2
Iq = pdep**2 + 12*rdep; Jq = 2*pdep**3 + 27*qdep**2 - 72*pdep*rdep
say('  Delta4(Q4 family) = %s' % factor(expand((4*Iq**3 - Jq**2)/27)))
say('')
say('PART2 ALL CHECKS PASSED')
