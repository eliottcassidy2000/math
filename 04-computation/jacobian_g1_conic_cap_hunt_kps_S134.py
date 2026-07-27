"""
G1 conic-cap stratified construction hunt (kind-pasteur, 2026-07-26).

Strata (deg A <= 2 box, A conical i.e. image of [A] a full conic):
  C1: A = (x^2, xy, y^2)   (p,q independent affine -> source-affine normal form)
  C2: A = (x^2, x, 1)      (degenerate parametrization q ~ const; GL2-on-(p,q)
                            + source affine reduce every remaining case here)
Staircase per stratum:
  D4 = 0 linear in B  ->  kernel;  D3 = D2 = 0 linear in C given B;
  D1 = 0, D0 = const quadratic in the leftover C-parameters (Rabinowitsch).
Field degree of every consistent solution via mod-p fiber Groebner counts.

Usage: python g1_conic_cap_hunt_kps.py <stage> [args]
Stages: sanity | kernels | joint1 | joint2 | batchC1 i0 i1 [degC] | batchC2 | enlarge i0 i1
"""
import sys, time, random, itertools, traceback
import sympy as sp

x, y, z = sp.symbols('x y z')
PMOD1, PMOD2 = 10007, 31013
T0 = time.time()
RNG_SEED = 20260726

def log(*a):
    print('[%7.1fs]' % (time.time() - T0), *a, flush=True)

# ---------------- core algebra ----------------

def bracket(u, v, w):
    """[u,v,w] = det(u|v|w), columns."""
    return sp.expand(
        u[0]*(v[1]*w[2] - v[2]*w[1])
        - v[0]*(u[1]*w[2] - u[2]*w[1])
        + w[0]*(u[1]*v[2] - u[2]*v[1]))

def dvec(V, s):
    return [sp.diff(c, s) for c in V]

def D_system(A, B, C):
    Ax, Ay = dvec(A, x), dvec(A, y)
    Bx, By = dvec(B, x), dvec(B, y)
    Cx, Cy = dvec(C, x), dvec(C, y)
    D5 = 2*bracket(Ax, Ay, A)
    D4 = bracket(Ax, Ay, B) + 2*bracket(Ax, By, A) + 2*bracket(Bx, Ay, A)
    D3 = (bracket(Ax, By, B) + bracket(Bx, Ay, B)
          + 2*(bracket(Ax, Cy, A) + bracket(Cx, Ay, A) + bracket(Bx, By, A)))
    D2 = (bracket(Bx, By, B) + bracket(Ax, Cy, B) + bracket(Cx, Ay, B)
          + 2*(bracket(Bx, Cy, A) + bracket(Cx, By, A)))
    D1 = bracket(Bx, Cy, B) + bracket(Cx, By, B) + 2*bracket(Cx, Cy, A)
    D0 = bracket(Cx, Cy, B)
    return [sp.expand(t) for t in (D5, D4, D3, D2, D1, D0)]

def build_F(A, B, C):
    return [sp.expand(A[i]*z**2 + B[i]*z + C[i]) for i in range(3)]

def detJ(F):
    J = sp.Matrix([[sp.diff(F[i], v) for v in (x, y, z)] for i in range(3)])
    return sp.expand(J.det())

def monoms_upto(d):
    return [(i, t - i) for t in range(d + 1) for i in range(t + 1)]

def make_vec(prefix, d, skip_const=False):
    syms, vec = [], []
    for comp in range(3):
        e = sp.Integer(0)
        for (i, j) in monoms_upto(d):
            if skip_const and i == 0 and j == 0:
                continue
            s = sp.Symbol('%s%d_%d_%d' % (prefix, comp, i, j))
            syms.append(s)
            e += s * x**i * y**j
        vec.append(e)
    return vec, syms

def xy_coeffs(expr):
    expr = sp.expand(expr)
    if expr == 0:
        return []
    return list(sp.Poly(expr, x, y).coeffs())

def xy_split_const(expr):
    expr = sp.expand(expr)
    if expr == 0:
        return [], sp.Integer(0)
    p = sp.Poly(expr, x, y)
    nonconst, const = [], sp.Integer(0)
    for mon, coef in p.terms():
        if mon == (0, 0):
            const = coef
        else:
            nonconst.append(coef)
    return nonconst, const

# ---------------- groebner helpers ----------------

def clear_den(e, gens):
    try:
        pl = sp.Poly(e, *gens)
        _, p2 = pl.clear_denoms(convert=True)
        return p2.as_expr()
    except Exception:
        return sp.expand(e * sp.together(e).as_numer_denom()[1])

def gb(eqs, gens, p=None):
    eqs = [sp.expand(e) for e in eqs]
    eqs = [e for e in eqs if e != 0]
    if not eqs:
        return None
    if p is None:
        return sp.groebner(eqs, *gens, order='grevlex')
    eqs = [clear_den(e, gens) for e in eqs]
    return sp.groebner(eqs, *gens, modulus=p, order='grevlex')

def gb_is_one(G):
    if G is None:
        return False
    ex = list(G.exprs)
    return len(ex) == 1 and ex[0].is_number and ex[0] != 0

def quotient_dim(G, gens):
    """#standard monomials if zero-dimensional else None."""
    lms = [pg.monoms(order='grevlex')[0] for pg in G.polys]
    n = len(gens)
    bounds = []
    for i in range(n):
        pures = [m[i] for m in lms if all(m[j] == 0 for j in range(n) if j != i)]
        if not pures:
            return None
        bounds.append(min(pures))
    cnt = 0
    for e in itertools.product(*[range(b) for b in bounds]):
        if not any(all(e[k] >= m[k] for k in range(n)) for m in lms):
            cnt += 1
    return cnt

def fiber_counts(F, p, trials=4, seed=7):
    rng = random.Random(seed)
    out = []
    for _ in range(trials):
        Pt = [rng.randrange(1, p) for _ in range(3)]
        eqs = [sp.expand(F[i] - Pt[i]) for i in range(3)]
        G = gb(eqs, (x, y, z), p=p)
        if gb_is_one(G):
            out.append(0)
            continue
        q = quotient_dim(G, (x, y, z))
        out.append(q if q is not None else -1)
    return out

# ---------------- strata ----------------

A_C1 = [x**2, x*y, y**2]
A_C2 = [x**2, x, sp.Integer(1)]
B_SEED = [-2*x**3, x**2 - 2*x**2*y, x + 2*x*y - 2*x*y**2]   # THM-2446 S5 conic B

def kernelB(A, degB):
    Bv, bs = make_vec('b', degB)
    D5, D4, _, _, _, _ = D_system(A, Bv, [sp.Integer(0)]*3)
    assert sp.expand(D5) == 0, 'stratum A must kill D5'
    eqs = xy_coeffs(D4)
    M, rhs = sp.linear_eq_to_matrix(eqs, bs)
    assert rhs == sp.zeros(len(eqs), 1)
    ns = M.nullspace()
    basis = []
    for v in ns:
        Bcand = [sp.expand(sum(v[k]*x**i*y**j for k, (comp, (i, j)) in enumerate(
            [(c, m) for c in range(3) for m in monoms_upto(degB)]) if comp == c0))
            for c0 in range(3)]
        basis.append(Bcand)
    return basis, bs, M

def vec_from_coords(coords, degB):
    idx = [(c, m) for c in range(3) for m in monoms_upto(degB)]
    out = [sp.Integer(0)]*3
    for k, (c, (i, j)) in enumerate(idx):
        out[c] += coords[k]*x**i*y**j
    return [sp.expand(e) for e in out]

def in_span(vecs_coords, target_coords):
    if not vecs_coords:
        return all(t == 0 for t in target_coords)
    M = sp.Matrix(vecs_coords).T
    try:
        M.gauss_jordan_solve(sp.Matrix(target_coords))
        return True
    except ValueError:
        return False

def coords_of_B(B, degB):
    idx = [(c, m) for c in range(3) for m in monoms_upto(degB)]
    out = []
    for (c, (i, j)) in idx:
        out.append(sp.expand(B[c]).coeff(x, i).coeff(y, j))
    return out

# ---------------- stages ----------------

def stage_sanity():
    rng = random.Random(3)
    Av = [sum(rng.randint(-3, 3)*x**i*y**j for (i, j) in monoms_upto(2)) for _ in range(3)]
    Bv = [sum(rng.randint(-3, 3)*x**i*y**j for (i, j) in monoms_upto(2)) for _ in range(3)]
    Cv = [sum(rng.randint(-3, 3)*x**i*y**j for (i, j) in monoms_upto(2)) for _ in range(3)]
    Ds = D_system(Av, Bv, Cv)
    F = build_F(Av, Bv, Cv)
    lhs = detJ(F)
    rhs = sp.expand(sum(Ds[5-k]*z**k for k in range(6)))
    assert sp.expand(lhs - rhs) == 0, 'D-system mismatch'
    log('SANITY OK: D-system == det J on random quadratic instance')
    # seed B in D4-kernel of C1?
    D5, D4, D3, D2, D1, D0 = D_system(A_C1, B_SEED, [sp.Integer(0)]*3)
    log('seed B: D5 =', D5, ' D4 =', sp.expand(D4))
    # B == 0 lemma
    log('B==0 forces D0=[C_x,C_y,0]=0: no 2-jet Keller map with B identically 0 (trivial).')

def stage_kernels():
    for tag, A in (('C1', A_C1), ('C2', A_C2)):
        for degB in (1, 2, 3, 4):
            basis, bs, M = kernelB(A, degB)
            dim = len(basis)
            # gauge subspace 2*A*phi, phi affine, inside box?
            gauge = []
            for phi in (sp.Integer(1), x, y):
                Bg = [sp.expand(2*a*phi) for a in A]
                if all(sp.degree(sp.Poly(c, x, y)) <= degB if c != 0 else True for c in Bg):
                    Bgc = coords_of_B(Bg, degB)
                    gauge.append((phi, Bgc))
            allc = [coords_of_B(b, degB) for b in basis]
            ok = [str(phi) for phi, gc in gauge if in_span(allc, gc)]
            log('stratum %s degB=%d: dim ker D4 = %d (of %d), gauge 2*A*phi in kernel for phi in %s'
                % (tag, degB, dim, 3*len(monoms_upto(degB)), ok))
            if degB == 3:
                for k, b in enumerate(basis):
                    log('   %s degB3 kernel basis[%d] = %s' % (tag, k, b))

def solve_C_stage(A, Bv, degC, tag, do_qq=True):
    """Given A, concrete B in ker D4: impose D3=D2=0 (linear in C), then D1=0,
    D0=const (Rabinowitsch).  Returns verdict string."""
    Cv, cs = make_vec('c', degC, skip_const=True)
    D5, D4, D3, D2, D1, D0 = D_system(A, Bv, Cv)
    if sp.expand(D4) != 0:
        log(tag, 'B NOT in D4-kernel, skip')
        return 'not-in-kernel'
    eqs = xy_coeffs(D3) + xy_coeffs(D2)
    M, rhs = sp.linear_eq_to_matrix(eqs, cs)
    try:
        sol, params = M.gauss_jordan_solve(rhs)
    except ValueError:
        log(tag, 'KILL exact at D3/D2: linear system in C inconsistent (QQ)')
        return 'kill-D3D2'
    taus = sorted(list(params), key=str) if params else []
    nfree = len(taus)
    sub = dict(zip(cs, list(sol)))
    Cpar = [sp.expand(Cv[i].subs(sub)) for i in range(3)]
    D5b, D4b, D3b, D2b, D1p, D0p = D_system(A, Bv, Cpar)
    assert sp.expand(D3b) == 0 and sp.expand(D2b) == 0, 'staircase substitution broke'
    eqs1 = xy_coeffs(D1p)
    ncD0, d0c = xy_split_const(D0p)
    eqsF = eqs1 + ncD0
    log(tag, 'C-affine space after D3=D2=0: dim %d (of %d); D1 eqs %d, D0 nonconst eqs %d'
        % (nfree, len(cs), len(eqs1), len(ncD0)))
    w = sp.Symbol('w_rab')
    if not eqsF:
        if d0c == 0:
            log(tag, 'D1,D0-nonconst vacuous but D0 const == 0: KILL')
            return 'kill-D0zero'
        return finish_witness(A, Bv, Cpar, taus, [], d0c, tag)
    gens = taus + [w]
    G1 = gb(eqsF + [w*d0c - 1], gens, p=PMOD1)
    if gb_is_one(G1):
        G2 = gb(eqsF + [w*d0c - 1], gens, p=PMOD2)
        if gb_is_one(G2):
            verdict = 'kill-D1D0-modp'
            if do_qq and nfree <= 24:
                try:
                    GQ = gb(eqsF + [w*d0c - 1], gens)
                    if gb_is_one(GQ):
                        log(tag, 'KILL exact (QQ Groebner = <1>): no Keller completion, any C in box')
                        return 'kill-D1D0-QQ'
                except Exception as ex:
                    log(tag, 'QQ GB failed:', repr(ex))
            log(tag, 'KILL mod both primes (Rabinowitsch GB=<1> mod 10007 & 31013)')
            return verdict
        log(tag, 'prime discrepancy!? mod', PMOD1, 'empty, mod', PMOD2, 'not')
    qd = quotient_dim(G1, gens) if G1 is not None else None
    log(tag, 'CONSISTENT mod %d: Keller variety nonempty; quotient dim mod p = %s'
        % (PMOD1, ('%d (zero-dim)' % qd) if qd is not None else 'positive-dimensional'))
    return attempt_witness(A, Bv, Cpar, taus, eqsF, d0c, tag)

def attempt_witness(A, Bv, Cpar, taus, eqsF, d0c, tag):
    # try QQ solve on the (D1, D0-nonconst) system, then filter d0c != 0
    try:
        GQ = gb(eqsF, taus)
        if gb_is_one(GQ):
            log(tag, 'D1/D0-nonconst empty over QQ but nonempty mod p?? (bad prime) -- treat as kill')
            return 'kill-D1D0-QQ'
        sols = sp.solve(list(GQ.exprs), taus, dict=True)
    except Exception as ex:
        log(tag, 'witness solve failed:', repr(ex))
        return 'consistent-no-witness'
    rng = random.Random(11)
    for so in sols:
        for attempt in range(6):
            full = dict(so)
            freev = [t for t in taus if t not in so or (so.get(t) is not None and so[t].free_symbols & set(taus))]
            valmap = {}
            for t in taus:
                if t in so:
                    continue
                valmap[t] = sp.Integer(rng.randint(-2, 2)) if attempt else sp.Integer(0)
            trial = {}
            okk = True
            for t in taus:
                v = so.get(t, t)
                v = sp.simplify(sp.sympify(v).subs(valmap))
                if v.free_symbols:
                    okk = False
                    break
                trial[t] = v
            if not okk:
                continue
            if any(sp.simplify(sp.sympify(e).subs(trial)) != 0 for e in eqsF):
                continue
            dc = sp.simplify(sp.sympify(d0c).subs(trial))
            if dc != 0:
                return finish_witness(A, Bv, Cpar, taus, trial, dc, tag)
        # also try the raw parametric solution at zero
    log(tag, 'solutions exist but all sampled points have D0const = 0 (degenerate, det J == 0): KILL-on-slice')
    return 'consistent-only-degenerate-sampled'

def finish_witness(A, Bv, Cpar, taus, trial, dc, tag):
    if isinstance(trial, dict):
        Cw = [sp.expand(c.subs(trial)) for c in Cpar]
    else:
        Cw = Cpar
    if any(sp.sympify(c).free_symbols - {x, y} for c in Cw):
        # free parameters remain (vacuous case): set them 0
        rest = set().union(*[sp.sympify(c).free_symbols - {x, y} for c in Cw])
        Cw = [sp.expand(c.subs({r: 0 for r in rest})) for c in Cw]
        dc2 = None
    F = build_F(A, Bv, Cw)
    dj = detJ(F)
    if dj.free_symbols or dj == 0:
        log(tag, 'witness FAILED verification: det J =', dj)
        return 'witness-failed'
    fc1 = fiber_counts(F, PMOD1)
    fdeg = max(fc1)
    fc2 = fiber_counts(F, PMOD2) if fdeg >= 2 else None
    log(tag, 'WITNESS: det J =', dj, ' fiber counts mod %d: %s' % (PMOD1, fc1),
        ('mod %d: %s' % (PMOD2, fc2)) if fc2 else '')
    log(tag, '  B =', Bv)
    log(tag, '  C =', Cw)
    if fdeg >= 4:
        log(tag, '*** JACKPOT CANDIDATE: field degree >= 4 -- run exact certificate ***')
        return 'JACKPOT?fdeg=%d' % fdeg
    return 'witness-fdeg-%d' % fdeg

def candidates_C1(degB=3):
    basis, bs, M = kernelB(A_C1, degB)
    allc = [coords_of_B(b, degB) for b in basis]
    # gauge coords
    gaugec = []
    for phi in (sp.Integer(1), x, y):
        Bg = [sp.expand(2*a*phi) for a in A_C1]
        gaugec.append(coords_of_B(Bg, degB))
    # complement of gauge inside kernel: row-reduce kernel coords mod gauge
    Mg = sp.Matrix(gaugec)
    comp = []
    cur = [r[:] for r in gaugec]
    for c in allc:
        if not in_span(cur, c):
            comp.append(c)
            cur.append(c)
    cands = []
    seedc = coords_of_B(B_SEED, degB)
    cands.append(('seed-THM2446S5', B_SEED))
    for k, c in enumerate(comp):
        cands.append(('kerbasis[%d]' % k, vec_from_coords(c, degB)))
    rng = random.Random(RNG_SEED)
    for t in range(10):
        while True:
            coeffs = [rng.choice([-1, 0, 1]) for _ in comp]
            if any(coeffs):
                break
        c = [sum(coeffs[i]*comp[i][k] for i in range(len(comp))) for k in range(len(comp[0]))]
        cands.append(('rand[%d]' % t, vec_from_coords(c, degB)))
    return cands, len(basis), len(comp), in_span(allc, seedc)

def stage_batchC1(i0, i1, degC=4):
    cands, kdim, cdim, seed_in = candidates_C1(3)
    log('C1 kernel dim (degB<=3) = %d, complement-of-gauge dim = %d, seed in kernel: %s'
        % (kdim, cdim, seed_in))
    results = []
    for k in range(i0, min(i1, len(cands))):
        name, Bv = cands[k]
        tag = 'C1[%s]' % name
        log('--- instance %d: %s, B = %s' % (k, name, Bv))
        try:
            v = solve_C_stage(A_C1, Bv, degC, tag)
        except Exception:
            traceback.print_exc()
            v = 'error'
        results.append((name, v))
        log('verdict:', name, '->', v)
    log('BATCH C1 SUMMARY [%d:%d] degC=%d:' % (i0, i1, degC))
    for n, v in results:
        log('   ', n, '->', v)

def stage_batchC2(degC=4, i0=0, i1=999):
    basis, bs, M = kernelB(A_C2, 3)
    allc = [coords_of_B(b, 3) for b in basis]
    gaugec = []
    for phi in (sp.Integer(1), x, y):
        Bg = [sp.expand(2*a*phi) for a in A_C2]
        gaugec.append(coords_of_B(Bg, 3))
    comp, cur = [], [r[:] for r in gaugec]
    for c in allc:
        if not in_span(cur, c):
            comp.append(c)
            cur.append(c)
    log('C2 kernel dim (degB<=3) = %d, complement-of-gauge dim = %d' % (len(basis), len(comp)))
    cands = [('kerbasis[%d]' % k, vec_from_coords(c, 3)) for k, c in enumerate(comp)]
    rng = random.Random(RNG_SEED + 1)
    for t in range(8):
        while True:
            coeffs = [rng.choice([-1, 0, 1]) for _ in comp]
            if any(coeffs):
                break
        c = [sum(coeffs[i]*comp[i][k] for i in range(len(comp))) for k in range(len(comp[0]))]
        cands.append(('rand[%d]' % t, vec_from_coords(c, 3)))
    results = []
    for k, (name, Bv) in enumerate(cands):
        if not (i0 <= k < i1):
            continue
        tag = 'C2[%s]' % name
        log('--- instance %d: %s, B = %s' % (k, name, Bv))
        try:
            v = solve_C_stage(A_C2, Bv, degC, tag)
        except Exception:
            traceback.print_exc()
            v = 'error'
        results.append((name, v))
        log('verdict:', name, '->', v)
    log('BATCH C2 SUMMARY degC=%d:' % degC)
    for n, v in results:
        log('   ', n, '->', v)

def stage_joint(which, degB=1, degC=2):
    """Exhaustive box: joint Rabinowitsch Groebner over B and C coefficients."""
    A = A_C1 if which == 1 else A_C2
    Bv, bs = make_vec('b', degB)
    Cv, cs = make_vec('c', degC, skip_const=True)
    D5, D4, D3, D2, D1, D0 = D_system(A, Bv, Cv)
    ncD0, d0c = xy_split_const(D0)
    w = sp.Symbol('w_rab')
    eqs = (xy_coeffs(D4) + xy_coeffs(D3) + xy_coeffs(D2) + xy_coeffs(D1)
           + ncD0 + [w*d0c - 1])
    gens = bs + cs + [w]
    log('joint stratum C%d degB<=%d degC<=%d: %d equations, %d unknowns (C const dropped)'
        % (which, degB, degC, len([e for e in eqs if e != 0]), len(gens)))
    G1 = gb(eqs, gens, p=PMOD1)
    if gb_is_one(G1):
        log('joint K1 C%d: EMPTY mod %d' % (which, PMOD1))
        G2 = gb(eqs, gens, p=PMOD2)
        log('joint K1 C%d: mod %d ->' % (which, PMOD2), 'EMPTY' if gb_is_one(G2) else 'NONEMPTY?!')
        try:
            GQ = gb(eqs, gens)
            log('joint K1 C%d: QQ Groebner ->' % which, 'EMPTY (PROOF-GRADE)' if gb_is_one(GQ)
                else 'NONEMPTY over QQ !!')
        except Exception as ex:
            log('joint K1 C%d: QQ GB failed:' % which, repr(ex))
        return
    qd = quotient_dim(G1, gens)
    log('joint K1 C%d: NONEMPTY mod %d, quotient dim = %s' % (which, PMOD1, qd))

def stage_enlarge(i0, i1):
    """degB <= 4 kernel directions (new ones beyond degB<=3 + gauge), degC = 4."""
    basis4, bs4, M4 = kernelB(A_C1, 4)
    all4 = [coords_of_B(b, 4) for b in basis4]
    # old directions embedded
    basis3, _, _ = kernelB(A_C1, 3)
    old = [coords_of_B(b, 4) for b in basis3]
    for phi in (sp.Integer(1), x, y, x*x, x*y, y*y):
        Bg = [sp.expand(2*a*phi) for a in A_C1]
        if max(sp.total_degree(sp.Poly(c, x, y)) if c != 0 else 0 for c in Bg) <= 4:
            old.append(coords_of_B(Bg, 4))
    comp, cur = [], [r[:] for r in old]
    for c in all4:
        if not in_span(cur, c):
            comp.append(c)
            cur.append(c)
    log('C1 kernel dim (degB<=4) = %d; new-beyond-(degB3+quadratic gauge) directions = %d'
        % (len(basis4), len(comp)))
    cands = [('ker4basis[%d]' % k, vec_from_coords(c, 4)) for k, c in enumerate(comp)]
    rng = random.Random(RNG_SEED + 2)
    for t in range(6):
        while True:
            coeffs = [rng.choice([-1, 0, 1]) for _ in comp]
            if any(coeffs):
                break
        c = [sum(coeffs[i]*comp[i][k] for i in range(len(comp))) for k in range(len(comp[0]))]
        cands.append(('rand4[%d]' % t, vec_from_coords(c, 4)))
    results = []
    for k in range(i0, min(i1, len(cands))):
        name, Bv = cands[k]
        tag = 'C1e[%s]' % name
        log('--- enlarge instance %d: %s, B = %s' % (k, name, Bv))
        try:
            v = solve_C_stage(A_C1, Bv, 4, tag)
        except Exception:
            traceback.print_exc()
            v = 'error'
        results.append((name, v))
        log('verdict:', name, '->', v)
    log('ENLARGE SUMMARY [%d:%d]:' % (i0, i1))
    for n, v in results:
        log('   ', n, '->', v)

if __name__ == '__main__':
    st = sys.argv[1] if len(sys.argv) > 1 else 'sanity'
    if st == 'sanity':
        stage_sanity()
    elif st == 'kernels':
        stage_kernels()
    elif st == 'joint1':
        stage_joint(1)
    elif st == 'joint2':
        stage_joint(2)
    elif st == 'jointX':
        stage_joint(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    elif st == 'batchC1':
        i0, i1 = int(sys.argv[2]), int(sys.argv[3])
        degC = int(sys.argv[4]) if len(sys.argv) > 4 else 4
        stage_batchC1(i0, i1, degC)
    elif st == 'batchC2':
        if len(sys.argv) > 3:
            stage_batchC2(4, int(sys.argv[2]), int(sys.argv[3]))
        else:
            stage_batchC2()
    elif st == 'enlarge':
        stage_enlarge(int(sys.argv[2]), int(sys.argv[3]))
    log('STAGE %s DONE' % st)
