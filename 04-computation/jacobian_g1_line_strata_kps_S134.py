# -*- coding: utf-8 -*-
# PART 4: (i) exact Euler-kill verification for r2a (H1/H0 dead);
#         (ii) the LIVE line-A strata of box (1,2,3): r2b (A = A0 + ell*w) and
#              r1 (A = ellt*v, B in span(v,w)); THM-2451's family is r2b's skeleton.
#         (iii) box (2,3,4): plane-valued quadratic A -- does D4 still kill B3?
import random, time, itertools, sys
from sympy import symbols, Matrix, Poly, groebner, expand, Integer, Rational, simplify
from sympy.polys.orderings import grevlex

x, y, z = symbols('x y z')
random.seed(31415926)
PR = 10007

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

# ---------- exact symbolic proofs of the Euler kills ----------
say('EXACT (symbolic over QQ): the r2a Euler kills, A = (x, y, 0)')
Apl = Matrix(3, 1, [x, y, 0])
b1, b2, b3 = symbols('b1 b2 b3', cls=lambda n, **kw: symbols(n))
from sympy import Function
B3 = sum(symbols('t%d%d' % (i, j))*x**i*y**j for i in range(4) for j in range(4-i))
Bgen = Matrix(3, 1, [sum(symbols('u%d%d' % (i, j))*x**i*y**j for i in range(3) for j in range(3-i)),
                     sum(symbols('v%d%d' % (i, j))*x**i*y**j for i in range(3) for j in range(3-i)),
                     B3])
D = dsystem(Apl, Bgen, Matrix(3, 1, [0, 0, 0]))
lhs = expand(D[1])
rhs = expand(B3 - 2*x*B3.diff(x) - 2*y*B3.diff(y))
assert expand(lhs - rhs) == 0
say('  D4(A=(x,y,0); any B with this B3) == (1 - 2x d/dx - 2y d/dy) B3   [EXACT, deg<=3 B3]')
say('  eigenvalue on degree-d monomials: 1-2d != 0  ==> D4=0 <=> B3=0. H1 IS EMPTY (r2a).')
Cgen = Matrix(3, 1, [sum(symbols('p%d%d' % (i, j))*x**i*y**j for i in range(4) for j in range(4-i)),
                     sum(symbols('q%d%d' % (i, j))*x**i*y**j for i in range(4) for j in range(4-i)),
                     sum(symbols('r%d%d' % (i, j))*x**i*y**j for i in range(4) for j in range(4-i))])
Bpl = Matrix(3, 1, [Bgen[0], Bgen[1], 0])
D = dsystem(Apl, Bpl, Cgen)
C3 = Cgen[2]
assert expand(D[2] + 2*(x*C3.diff(x) + y*C3.diff(y))) == 0
say('  D3(A=(x,y,0), B3=0) == -2*(x d/dx + y d/dy) C3  [EXACT] ==> D3=0 <=> C3 = const.')
C3c = Cgen[2].subs({symbols('r%d%d' % (i, j)): 0 for i in range(4) for j in range(4-i) if (i, j) != (0, 0)})
Cc = Matrix(3, 1, [Cgen[0], Cgen[1], C3c])
Dc = dsystem(Apl, Bpl, Cc)
assert expand(Dc[3]) == 0 and expand(Dc[4]) == 0 and expand(Dc[5]) == 0
say('  and then D2 == D1 == D0 == 0 identically [EXACT] ==> det J == 0: H0 IS EMPTY (r2a).')
say('  ==> r2a Keller maps are exactly the T-shape ==> field degree 1 (Moh reduction).')
say('')

# ---------- machinery (mod p) ----------
def qdim(G, nvars):
    polys = G.polys
    if any(pl.is_one for pl in polys): return 0
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
    for mono in itertools.product(*[range(bb) for bb in bounds]):
        if not any(all(mono[j] >= lm[j] for j in range(nvars)) for lm in lms):
            cnt += 1
    return cnt

def space_count(A, B, C, p):
    target = [random.randrange(p) for _ in range(3)]
    F = A*z**2 + B*z + C
    G = groebner([expand(F[i]-target[i]) for i in range(3)], x, y, z, modulus=p, order='grevlex')
    return qdim(G, 3)

def rpoly2(deg, p): return sum(random.randrange(p)*x**i*y**j for i in range(deg+1) for j in range(deg+1-i))

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

def polydeg_modp(expr, p):
    pp = Poly(expand(expr), x, y, modulus=p)
    return None if pp.is_zero else pp.total_degree()

def report(label, A, B, C, p=PR, ntar=2):
    Ns = [space_count(A, B, C, p) for _ in range(ntar)]
    D = dsystem(A, B, C)
    degs = [polydeg_modp(D[i], p) for i in range(6)]
    say('  %-40s N=%s  deg(D5..D0)=%s' % (label, Ns, degs))
    return Ns

say('--- LIVE STRATUM r2b: A = A0 + ell*w (affine-line A, THM-2451 skeleton) ---')
for trial in range(3):
    A0v = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
    wv = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
    ell = rpoly2(1, PR)
    A = A0v + wv*ell
    # D5 = 2[A_x,A_y,A]: A_x = ell_x w, A_y = ell_y w parallel ==> D5 == 0 automatically
    Bm, bsyms = generic_vec_syms('b', 2)
    D = dsystem(A, Bm, Matrix(3, 1, [0, 0, 0]))
    rows, rhs = linear_conditions(D[1], bsyms, PR)
    sol = lin_solve_modp(rows, rhs, len(bsyms), PR)
    if sol is None:
        say('  trial %d: D4 system inconsistent (unexpected)' % trial); continue
    partB, basisB, rankB = sol
    B = Bm.subs({bsyms[i]: v for i, v in enumerate(random_solution(partB, basisB, PR))})
    Cm, csyms = generic_vec_syms('c', 3)
    D = dsystem(A, B, Cm)
    rows, rhs = [], []
    for w_ in (2, 3):
        r_, h_ = linear_conditions(D[w_], csyms, PR)
        rows += r_; rhs += h_
    solC = lin_solve_modp(rows, rhs, len(csyms), PR)
    if solC is None:
        say('  trial %d: rankB=%d; D3&D2 INCONSISTENT in C (B-source obstructed)' % (trial, rankB))
        continue
    partC, basisC, rankC = solC
    C = Cm.subs({csyms[i]: v for i, v in enumerate(random_solution(partC, basisC, PR))})
    for idx in (0, 1, 2, 3):
        assert Poly(expand(dsystem(A, B, C)[idx]), x, y, modulus=PR).is_zero
    say('  trial %d: rankB(D4)=%d kerB=%d; rankC(D3&D2)=%d solC=%d' %
        (trial, rankB, 18-rankB, rankC, 30-rankC))
    report('r2b + D4=D3=D2=0', A, B, C)

say('')
say('--- LIVE STRATUM r1: A = ellt*v, B in span(v,w) (D3 auto-zero), D2 linear in C ---')
for trial in range(3):
    v = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
    w = Matrix(3, 1, [random.randrange(PR) for _ in range(3)])
    ellt = rpoly2(1, PR)
    A = v*ellt
    B = v*rpoly2(2, PR) + w*rpoly2(2, PR)
    Dpre = dsystem(A, B, Matrix(3, 1, [0, 0, 0]))
    assert Poly(expand(Dpre[0]), x, y, modulus=PR).is_zero  # D5
    assert Poly(expand(Dpre[1]), x, y, modulus=PR).is_zero  # D4 (A line-valued)
    assert Poly(expand(Dpre[2]), x, y, modulus=PR).is_zero  # D3 (B in span(v,w) with A || v)
    Cm, csyms = generic_vec_syms('c', 3)
    D = dsystem(A, B, Cm)
    rows, rhs = linear_conditions(D[3], csyms, PR)   # D2 linear in C
    solC = lin_solve_modp(rows, rhs, len(csyms), PR)
    if solC is None:
        say('  trial %d: D2 INCONSISTENT in C' % trial); continue
    partC, basisC, rankC = solC
    C = Cm.subs({csyms[i]: vv for i, vv in enumerate(random_solution(partC, basisC, PR))})
    assert Poly(expand(dsystem(A, B, C)[3]), x, y, modulus=PR).is_zero
    say('  trial %d: rankC(D2)=%d solC=%d' % (trial, rankC, 30-rankC))
    report('r1(span) + D5=D4=D3=D2=0', A, B, C)

say('')
say('--- box (2,3,4): plane-valued QUADRATIC A: does D4 still kill B3? ---')
for trial in range(2):
    A = Matrix(3, 1, [rpoly2(2, PR), rpoly2(2, PR), Integer(0)])
    Bm, bsyms = generic_vec_syms('b', 3)
    D = dsystem(A, Bm, Matrix(3, 1, [0, 0, 0]))
    rows, rhs = linear_conditions(D[1], bsyms, PR)
    sol = lin_solve_modp(rows, rhs, len(bsyms), PR)
    partB, basisB, rankB = sol
    # is B3 forced to 0 on the kernel? check whether any kernel basis vector has B3-part
    n1 = 10  # coeffs per component for deg 3
    b3_free = any(any(bv[2*n1 + j] % PR for j in range(n1)) for bv in basisB)
    say('  trial %d: rank(D4) = %d on 30 unknowns; kernel dim %d; B3-direction free: %s'
        % (trial, rankB, 30-rankB, b3_free))
    B = Bm.subs({bsyms[i]: v for i, v in enumerate(random_solution(partB, basisB, PR))})
    d4chk = Poly(expand(dsystem(A, B, Matrix(3, 1, [0, 0, 0]))[1]), x, y, modulus=PR).is_zero
    b3deg = polydeg_modp(B[2], PR)
    say('           sample kernel B: D4==0: %s, deg B3 = %s' % (d4chk, b3deg))

say('')
say('PART4 DONE')
