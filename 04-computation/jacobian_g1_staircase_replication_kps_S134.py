# -*- coding: utf-8 -*-
# PART 3: replication of the staircase-degeneracy signal in box (1,2,3).
# Question: does {D5=D4=D3=D2=0} on shapes H1/H0/r1 force det J == 0 identically
# (i.e. D1 == D0 == 0) for RANDOM staircase choices, or was PART 2's N=0 a branch fluke?
import random, time, itertools, sys
from sympy import symbols, Matrix, Poly, groebner, expand, Integer
from sympy.polys.orderings import grevlex

x, y, z, t = symbols('x y z t')
random.seed(555001)
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

def generic_vec_syms(prefix, deg, mask=(1, 1, 1)):
    n = (deg+1)*(deg+2)//2
    syms = symbols(prefix + '0:%d' % (3*n))
    monos = [x**i*y**j for i in range(deg+1) for j in range(deg+1-i)]
    used = []
    comps = []
    for c in range(3):
        if mask[c]:
            comps.append(sum(syms[c*n+m]*monos[m] for m in range(n)))
            used += [syms[c*n+m] for m in range(n)]
        else:
            comps.append(Integer(0))
    return Matrix(3, 1, comps), used

def random_solution(part, basis, p):
    vec = part[:]
    for bv in basis:
        r = random.randrange(p)
        vec = [(vec[i]+r*bv[i]) % p for i in range(len(vec))]
    return vec

def polydeg_modp(expr, p):
    pp = Poly(expand(expr), x, y, modulus=p)
    return None if pp.is_zero else pp.total_degree()

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

def staircase_trial(shape, p=PR):
    """returns (rankB, kerB_dim, rankC, solC_dim, degD1, degD0, N)"""
    if shape == 'H1':
        A = Matrix(3, 1, [rpoly2(1, p), rpoly2(1, p), Integer(0)])
        Bmask = (1, 1, 1)
    elif shape == 'H0':
        A = Matrix(3, 1, [rpoly2(1, p), rpoly2(1, p), Integer(0)])
        Bmask = (1, 1, 0)
    elif shape == 'r1':
        v = Matrix(3, 1, [random.randrange(p) for _ in range(3)])
        A = v*rpoly2(1, p)
        Bmask = (1, 1, 1)
    else:
        raise ValueError(shape)
    Bm, bsyms = generic_vec_syms('b', 2, Bmask)
    D = dsystem(A, Bm, Matrix(3, 1, [0, 0, 0]))
    rows, rhs = linear_conditions(D[1], bsyms, p)
    sol = lin_solve_modp(rows, rhs, len(bsyms), p)
    if sol is None: return None
    partB, basisB, rankB = sol
    B = Bm.subs({bsyms[i]: v for i, v in enumerate(random_solution(partB, basisB, p))})
    Cm, csyms = generic_vec_syms('c', 3)
    D = dsystem(A, B, Cm)
    rows, rhs = [], []
    for w in (2, 3):
        r_, h_ = linear_conditions(D[w], csyms, p)
        rows += r_; rhs += h_
    solC = lin_solve_modp(rows, rhs, len(csyms), p)
    if solC is None:
        return (rankB, len(bsyms)-rankB, None, None, None, None, None)
    partC, basisC, rankC = solC
    C = Cm.subs({csyms[i]: v for i, v in enumerate(random_solution(partC, basisC, p))})
    Dfin = dsystem(A, B, C)
    for lbl, idx in (('D5', 0), ('D4', 1), ('D3', 2), ('D2', 3)):
        assert Poly(expand(Dfin[idx]), x, y, modulus=p).is_zero, lbl + ' not zero'
    d1 = polydeg_modp(Dfin[4], p); d0 = polydeg_modp(Dfin[5], p)
    N = space_count(A, B, C, p)
    return (rankB, len(bsyms)-rankB, rankC, len(csyms)-rankC, d1, d0, N)

say('PART 3: staircase replication, box (1,2,3); entries =')
say('(rankB, dim ker D4, rankC, dim sol D3&D2, deg D1, deg D0, N) -- deg None = identically 0')
for shape in ('H1', 'H0', 'r1'):
    for trial in range(4):
        res = staircase_trial(shape)
        say('  %s trial %d: %s' % (shape, trial, res))

say('')
say('Interpretation: if deg D1 = deg D0 = None across trials, the affine-linear')
say('staircase branch of {D5=D4=D3=D2=0} lies inside {det J == 0} for that shape;')
say('Keller solutions (D0 = const != 0) then require B outside the sampled kernel')
say('branch or joint (B,C) solving -- a hunt-design fact, not a box emptiness proof.')
say('PART3 DONE')
