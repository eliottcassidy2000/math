#!/usr/bin/env python3
"""
lane J1 -- JC(2) leading-form tower: rigorous local law, Euclid chain, search order.

All arithmetic exact (sympy over QQ / Fraction).  Seven parts, each prints a
PASS/FAIL line.

PART 1   Lemma J   : local Jacobian valuation law  v_L(Jac(F,G)) = v_L F+v_L G-1
                    unless v_L F/deg F = v_L G/deg G.
PART 2   Tower     : graded identities, L0/L1, TTL, T1, Euclid chain, equal-slope
                    law, hull dichotomy -- all on genuine automorphisms.
PART 3   Obstruction: dim ker / dim coker of Phi_j for prescribed (H,a,b) with
                    a,b >= 2 -- the counterexample regime.
PART 3b  T2        : P_(n-1) determined by Q_(m-1) up to one H0-power scalar.
PART 3c  Reduced   : the explicit next orders (R_j), j <= floor((a-1)/b)+1.
PART 4   Feasibility: exact hull-matching arithmetic (single slope) -> t<=g-1.
PART 5   Table     : the search order over coprime (a,b), a,b>=2, by CF depth.

Written to be re-runnable: python3 jc2_tower_depth_laneJ1.py
"""
import sympy as sp
from sympy import symbols, Rational, Integer
from fractions import Fraction
import random, itertools, sys

x, y = symbols('x y')
oo = sp.oo
RESULTS = []


def jac(P, Q):
    return sp.expand(sp.diff(P, x)*sp.diff(Q, y) - sp.diff(P, y)*sp.diff(Q, x))


def hparts(F):
    """dict: total degree -> homogeneous component (nonzero ones only)."""
    F = sp.expand(F)
    if F == 0:
        return {}
    p = sp.Poly(F, x, y)
    d = {}
    for mon, c in zip(p.monoms(), p.coeffs()):
        deg = mon[0] + mon[1]
        d[deg] = d.get(deg, 0) + c*x**mon[0]*y**mon[1]
    return {k: sp.expand(v) for k, v in d.items() if sp.expand(v) != 0}


def vL(F, L):
    """L-adic valuation of a binary form F w.r.t. a linear form L; oo if F=0."""
    F = sp.expand(F)
    if F == 0:
        return oo
    v = 0
    PL = sp.Poly(L, x, y)
    while True:
        q, r = sp.div(sp.Poly(F, x, y), PL)
        if sp.expand(r.as_expr()) == 0:
            F = sp.expand(q.as_expr())
            v += 1
            if F == 0:
                return oo
        else:
            return v


def tdeg(F):
    F = sp.expand(F)
    return None if F == 0 else sp.Poly(F, x, y).total_degree()


# ----------------------------------------------------------------------------
# PART 1 -- Lemma J
# ----------------------------------------------------------------------------
def part1(ntrials=260, seed=17):
    rng = random.Random(seed)
    bad = []
    tested_generic = 0
    tested_exc = 0
    for _ in range(ntrials):
        c = rng.randint(-3, 3)
        L = x + c*y if rng.random() < 0.8 else y
        e1, e2 = rng.randint(0, 4), rng.randint(0, 4)
        d1, d2 = rng.randint(0, 4), rng.randint(0, 4)
        Fp = sum(rng.randint(-4, 4)*x**(d1-i)*y**i for i in range(d1+1))
        Gp = sum(rng.randint(-4, 4)*x**(d2-i)*y**i for i in range(d2+1))
        F = sp.expand(L**e1*Fp)
        G = sp.expand(L**e2*Gp)
        if F == 0 or G == 0:
            continue
        vF, vG = vL(F, L), vL(G, L)
        dF, dG = tdeg(F), tdeg(G)
        exc = (vF*dG == vG*dF)
        J = jac(F, G)
        vJ = vL(J, L)
        if exc:
            tested_exc += 1
            ok = (vJ == oo) or (vJ >= vF + vG)
        else:
            tested_generic += 1
            ok = (vJ == vF + vG - 1)
        if not ok:
            bad.append((L, F, G, vF, vG, dF, dG, exc, vJ))
    ok = (len(bad) == 0)
    print("PART 1  Lemma J (local Jacobian valuation law)")
    print("   generic cases tested : %d   exceptional cases tested : %d" %
          (tested_generic, tested_exc))
    print("   violations           : %d" % len(bad))
    for b in bad[:3]:
        print("   FAIL", b)
    print("   PART1 PASS =", ok)
    RESULTS.append(("PART1 LemmaJ", ok))
    return ok


# ----------------------------------------------------------------------------
# automorphism generator
# ----------------------------------------------------------------------------
def build_auto(rng, nfac=2, maxdeg=3, degcap=40):
    P, Q = x, y
    for _ in range(nfac):
        if max(tdeg(P) or 0, tdeg(Q) or 0)*maxdeg > degcap:
            break
        while True:
            a, b, c, d = [rng.randint(-2, 2) for _ in range(4)]
            if a*d - b*c != 0:
                break
        P, Q = sp.expand(a*P + b*Q), sp.expand(c*P + d*Q)
        dd = rng.randint(2, maxdeg)
        poly = sum(rng.randint(-2, 2)*x**i for i in range(dd)) + \
            rng.choice([1, -1, 2])*x**dd
        if rng.random() < 0.5:
            Q = sp.expand(Q + poly.subs(x, P))
        else:
            P = sp.expand(P + poly.subs(x, Q))
    return sp.expand(P), sp.expand(Q)


def linear_factors(F):
    """binary form -> list of (linear form, multiplicity), over QQ(i) if needed.
    We only use automorphism leading forms, which split over QQ in our sample;
    fall back to sympy factor_list and assert linearity."""
    fl = sp.factor_list(F)
    out = []
    for f, mult in fl[1]:
        if sp.Poly(f, x, y).total_degree() == 1:
            out.append((sp.expand(f), mult))
        else:
            return None            # non-split leading form: skip this sample
    return out


def hull_slopes(seq):
    """seq: list of values v_0..v_N in Z>=0 u {oo}. Returns the decreasing-slope
    list of the LOWER convex hull of the finite points, as Fractions, together
    with the vertex index list.  Slopes are returned as positive numbers
    sigma = -(slope), i.e. the descent rates."""
    pts = [(i, int(v)) for i, v in enumerate(seq) if v != oo]
    if not pts:
        return [], []
    # Andrew monotone chain, lower hull
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            # keep if turn is counter-clockwise (convex from below)
            if (x2-x1)*(p[1]-y1) - (y2-y1)*(p[0]-x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    slopes = []
    for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
        slopes.append(Fraction(y1-y2, x2-x1))
    return slopes, [h[0] for h in hull]


def sigma_first(seq):
    """sigma = max_{i>=1} (v_0 - v_i)/i  -- the first (steepest) descent rate."""
    v0 = seq[0]
    best = None
    for i, v in enumerate(seq):
        if i == 0 or v == oo:
            continue
        r = Fraction(int(v0) - int(v), i)
        if best is None or r > best:
            best = r
    return best


def crossing_vertex(slopes, verts, s):
    """B-hull vertex index kappa where the descent rate crosses the value s:
    rates strictly before kappa are > s, rates after are < s.  Returns None if
    s occurs as a rate of this hull."""
    for r in slopes:
        if r == s:
            return None
    kappa = verts[0]
    for r, v in zip(slopes, verts[1:]):
        if r > s:
            kappa = v
        else:
            break
    return kappa


def analyse_pair(P, Q, checks, verbose=False, tag=""):
    """Run every tower check on a Jacobian pair (Jac = nonzero constant)."""
    hp, hq = hparts(P), hparts(Q)
    n, m = max(hp), max(hq)
    Jc = jac(P, Q)
    if not (Jc.is_number and Jc != 0):
        return None
    if n + m < 3:
        return None

    # --- graded identities
    grad_ok = True
    for j in range(0, n+m-1):
        s = 0
        for i in range(0, j+1):
            k = j - i
            A = hp.get(n-i, 0)
            B = hq.get(m-k, 0)
            if A != 0 and B != 0:
                s += jac(A, B)
        want = Jc if j == n+m-2 else 0
        if sp.expand(s - want) != 0:
            grad_ok = False
    checks['graded'].append(grad_ok)

    g = sp.igcd(n, m)
    a, b = n//g, m//g
    Pn, Qm = hp[n], hq[m]

    # --- L0 : Pn = c H^a , Qm = c' H^b with deg H = g
    facs = linear_factors(Pn)
    if facs is None:
        return None
    K = len(facs)
    mult_ok = all(mu % a == 0 for _, mu in facs)
    Hfacs = [(L, mu//a) for L, mu in facs]
    H = 1
    for L, e in Hfacs:
        H = sp.expand(H*L**e)
    l0_ok = mult_ok and tdeg(H) == g
    if l0_ok:
        cP = sp.simplify(sp.expand(Pn/H**a))
        cQ = sp.simplify(sp.expand(Qm/H**b))
        l0_ok = cP.is_number and cQ.is_number and cP != 0 and cQ != 0
    checks['L0'].append(l0_ok)
    checks['K1'].append(K == 1)
    checks['minab1'].append(min(a, b) == 1)
    if not l0_ok:
        return None

    P1 = hp.get(1, 0)
    Q1 = hq.get(1, 0)

    for (L, e) in Hfacs:
        us = [vL(hp.get(n-i, 0), L) for i in range(0, n+1)]
        ws = [vL(hq.get(m-k, 0), L) for k in range(0, m+1)]
        assert us[0] == a*e and ws[0] == b*e, (us[0], a*e, ws[0], b*e)

        # --- TTL : for every order j, min over i+k=j attained >=2 times, or the
        #     unique minimiser is exceptional  u_i (m-k) = w_k (n-i)
        ttl_ok = True
        for j in range(1, n+m-2):
            vals = []
            for i in range(max(0, j-m), min(n, j)+1):
                k = j-i
                if us[i] == oo or ws[k] == oo:
                    continue
                vals.append((us[i]+ws[k], i, k))
            if not vals:
                continue
            mn = min(v[0] for v in vals)
            arg = [v for v in vals if v[0] == mn]
            if len(arg) >= 2:
                continue
            _, i, k = arg[0]
            if us[i]*(m-k) == ws[k]*(n-i):
                continue
            ttl_ok = False
        checks['TTL'].append(ttl_ok)

        # NOTE ON SCOPE.  Every claim below is derived from the graded identity
        # at some order j, which exists (and reads "=0") only for
        # 1 <= j <= n+m-3.  jtop is that bound.
        jtop = n+m-3

        # --- T1 : u_1 >= (a-b) e  and (non/doubly exceptional) u_1 - w_1 = (a-b)e
        if jtop >= 1:
            u1, w1 = us[1], ws[1]
            d = a-b
            t1_div = (u1 == oo) or (u1 >= d*e)
            checks['T1div'].append(t1_div)
            if not t1_div:
                print("   T1div FAIL", tag, L, e, a, b, n, m, us[:4], ws[:4])
            exc01 = (w1 != oo) and (w1*n == (a*e)*(m-1))   # (0,1) exceptional
            exc10 = (u1 != oo) and (u1*m == (b*e)*(n-1))   # (1,0) exceptional
            if (not exc01 and not exc10) or (exc01 and exc10):
                t1_exact = (u1 == oo or w1 == oo) or (u1 - w1 == d*e)
                checks['T1exact'].append(t1_exact)
                if not t1_exact:
                    print("   T1exact FAIL", tag, L, e, a, b, n, m, us[:4], ws[:4])

        # --- Euclid chain : u_i >= (a - i b) e  for 0 <= i <= floor(a/b),
        #     each step i needing order j=i to be inside the tower.
        ec = True
        if a >= b:
            for i in range(0, min(a//b, jtop) + 1):
                if i < len(us) and us[i] != oo and us[i] < (a-i*b)*e:
                    ec = False
        else:
            for k in range(0, min(b//a, jtop) + 1):
                if k < len(ws) and ws[k] != oo and ws[k] < (b-k*a)*e:
                    ec = False
        checks['euclid'].append(ec)
        if not ec:
            print("   euclid FAIL", tag, L, e, a, b, n, m, us, ws)

        # --- equal-slope law, in the scope where it is proved:
        #     needs n,m >= 2 and L nmid P_1, L nmid Q_1
        sP, sQ = sigma_first(us), sigma_first(ws)
        LdivP1 = (P1 != 0 and vL(P1, L) >= 1)
        LdivQ1 = (Q1 != 0 and vL(Q1, L) >= 1)
        in_scope = (n >= 2 and m >= 2 and not LdivP1 and not LdivQ1
                    and sP is not None and sQ is not None)
        if in_scope:
            checks['sigma_scoped'].append(sP == sQ)
            if sP != sQ:
                print("   sigma FAIL", tag, L, e, a, b, n, m, us, ws, sP, sQ)
        checks['sigma_all'].append(sP == sQ)

        # --- hull slope DICHOTOMY: every descent rate of one hull that is not a
        #     descent rate of the other must sit over an exceptional corner.
        slA, vertA = hull_slopes(us)
        slB, vertB = hull_slopes(ws)
        setA = set(s for s in slA if s > 0)
        setB = set(s for s in slB if s > 0)
        checks['hullsets_equal'].append(setA == setB)
        if n >= 2 and m >= 2:
            dich = True
            for (sl, vt, other_sl, other_vt, side) in [
                    (slA, vertA, slB, vertB, 'A'), (slB, vertB, slA, vertA, 'B')]:
                for idx, s in enumerate(sl):
                    if s <= 0 or s in set(other_sl):
                        continue
                    kap = crossing_vertex(other_sl, other_vt, s)
                    if kap is None:
                        continue
                    p, pa = vt[idx], vt[idx+1]
                    for pp in (p, pa):
                        if side == 'A':
                            i, k = pp, kap
                        else:
                            i, k = kap, pp
                        if i+k < 1 or i+k > jtop:
                            continue
                        if us[i] == oo or ws[k] == oo:
                            continue
                        if us[i]*(m-k) != ws[k]*(n-i):
                            dich = False
                            print("   dichotomy FAIL", tag, L, e, a, b, n, m,
                                  us, ws, s, i, k)
            checks['hull_dichotomy'].append(dich)

        if verbose:
            print("   %s L=%s e=%d  a,b=(%d,%d) g=%d" % (tag, L, e, a, b, g))
            print("      u =", us)
            print("      w =", ws)
            print("      sigma_P=%s sigma_Q=%s   hullA=%s hullB=%s" %
                  (sP, sQ, slA, slB))
    return dict(n=n, m=m, g=g, a=a, b=b, K=K, H=H)


def part2(verbose_first=3):
    print()
    print("PART 2  the tower on genuine Jacobian pairs (polynomial automorphisms)")
    checks = {k: [] for k in ['graded', 'L0', 'K1', 'minab1', 'TTL', 'T1div',
                              'T1exact', 'euclid', 'sigma_scoped', 'sigma_all',
                              'hullsets_equal', 'hull_dichotomy']}
    rng = random.Random(20260731)
    samples = []
    tries = 0
    while len(samples) < 30 and tries < 400:
        tries += 1
        nf = rng.choice([2, 2, 3, 3])
        md = rng.choice([2, 3, 3, 4])
        P, Q = build_auto(rng, nfac=nf, maxdeg=md)
        hp, hq = hparts(P), hparts(Q)
        if not hp or not hq:
            continue
        n, m = max(hp), max(hq)
        if n + m < 8 or n + m > 36 or min(n, m) < 2:
            continue
        r = analyse_pair(P, Q, checks, verbose=(len(samples) < verbose_first),
                         tag="[%02d]" % len(samples))
        if r is None:
            continue
        samples.append(r)

    print("   samples analysed          :", len(samples))
    degs = sorted(set((s['n'], s['m'], s['g'], s['a'], s['b'], s['K'])
                      for s in samples))
    print("   (n,m,g,a,b,K) multiset    :", degs)
    allok = True
    for k in ['graded', 'L0', 'TTL', 'T1div', 'T1exact', 'euclid',
              'sigma_scoped', 'hull_dichotomy']:
        v = checks[k]
        ok = all(v)
        allok &= ok
        print("   %-16s : %4d checks, all True = %s" % (k, len(v), ok))
    for k in ['K1', 'minab1', 'sigma_all', 'hullsets_equal']:
        v = checks[k]
        print("   %-16s : %4d checks, all True = %s   (observation, not a claim)"
              % (k, len(v), all(v)))
    print("   PART2 PASS =", allok)
    RESULTS.append(("PART2 tower on automorphisms", allok))
    return allok


# ----------------------------------------------------------------------------
# PART 3 -- deformation tower for prescribed (H,a,b) with a,b >= 2
# ----------------------------------------------------------------------------
def basis(d):
    return [x**(d-i)*y**i for i in range(d+1)] if d >= 0 else []


def coeff_vec(F, d):
    F = sp.expand(F)
    p = sp.Poly(F, x, y) if F != 0 else None
    out = []
    for i in range(d+1):
        out.append(0 if p is None else p.coeff_monomial(x**(d-i)*y**i))
    return out


def solve_tower(H, a, b, jmax=None, rng=None, trials=1):
    """P_n = H^a, Q_m = H^b.  Solve  Jac(P_n,Q_{m-j}) + Jac(P_{n-j},Q_m)
       = -sum_{i+k=j,i,k>=1} Jac(P_{n-i},Q_{m-k})  order by order over QQ,
       picking a random rational point of the solution set at each order.
       Returns (jbreak, log) where jbreak is the first order that is
       inconsistent (None if we reach jmax)."""
    g = tdeg(H)
    n, m = g*a, g*b
    if jmax is None:
        jmax = n+m-3
    Pn = sp.expand(H**a)
    Qm = sp.expand(H**b)
    best = -1
    bestlog = None
    for _t in range(trials):
        Pc = {0: Pn}
        Qc = {0: Qm}
        log = []
        broke = None
        for j in range(1, jmax+1):
            dP, dQ = n-j, m-j
            rhs = 0
            for i in range(1, j):
                k = j-i
                A = Pc.get(i, 0)
                B = Qc.get(k, 0)
                if A != 0 and B != 0:
                    rhs += jac(A, B)
            rhs = sp.expand(-rhs)
            tgt = n+m-2-j
            if tgt < 0:
                break
            cols = []
            for mo in basis(dP):
                cols.append(coeff_vec(jac(mo, Qm), tgt))
            for mo in basis(dQ):
                cols.append(coeff_vec(jac(Pn, mo), tgt))
            if not cols:
                if sp.expand(rhs) != 0:
                    broke = j
                    break
                continue
            M = sp.Matrix(cols).T          # (tgt+1) x (#unknowns)
            rv = sp.Matrix(coeff_vec(rhs, tgt))
            aug = M.row_join(rv)
            if M.rank() != aug.rank():
                broke = j
                log.append((j, 'inconsistent', M.rows, M.cols, M.rank()))
                break
            sol = M.gauss_jordan_solve(rv)
            part, params = sol
            if params:
                sub = {p: Rational(rng.randint(-3, 3) if rng else 1, 1)
                       for p in params}
                part = part.subs(sub)
            vals = list(part)
            A = sum(vals[i]*mo for i, mo in enumerate(basis(dP)))
            B = sum(vals[len(basis(dP))+i]*mo for i, mo in enumerate(basis(dQ)))
            Pc[j] = sp.expand(A)
            Qc[j] = sp.expand(B)
            log.append((j, 'ok', M.rows, M.cols, M.rank(), len(params)))
        reached = (broke-1) if broke else jmax
        if reached > best:
            best, bestlog, bb = reached, log, broke
    return bb if 'bb' in dir() else broke, best, bestlog


def primitive_part(H):
    """H = c * H0^d with H0 primitive (not a proper power).  Return (H0, d)."""
    fl = sp.factor_list(H)
    mults = [mu for _, mu in fl[1]]
    d = 0
    for mu in mults:
        d = sp.igcd(d, mu)
    H0 = 1
    for f, mu in fl[1]:
        H0 = sp.expand(H0*f**(mu//d))
    return sp.expand(H0), int(d)


def Phi_matrix(H, a, b, j):
    """matrix of  Phi_j(A,B) = Jac(A,H^b) + Jac(H^a,B)  on monomial bases."""
    g = tdeg(H)
    n, m = g*a, g*b
    Pn, Qm = sp.expand(H**a), sp.expand(H**b)
    dP, dQ, tgt = n-j, m-j, n+m-2-j
    cols = []
    for mo in basis(dP):
        cols.append(coeff_vec(jac(mo, Qm), tgt))
    for mo in basis(dQ):
        cols.append(coeff_vec(jac(Pn, mo), tgt))
    if not cols:
        return None
    return sp.Matrix(cols).T


def part3():
    print()
    print("PART 3  the graded obstruction map  Phi_j(A,B) = Jac(A,Q_m)+Jac(P_n,B)")
    print("   PREDICTION (proved in the artifact, a>=b):")
    print("       dim ker  Phi_j = max(0, m-j+1) + [g0 | j]")
    print("       dim coker Phi_j = (m-2) + [g0 | j]        for 1 <= j <= m")
    print("   where H = c*H0^d, g0 = deg H0.  The '+[g0|j]' is the ONLY free")
    print("   parameter beyond 'P_(n-j) is forced by Q_(m-j)'.")
    Hs = {
        'y^2'            : sp.expand(y**2),
        'y^3'            : sp.expand(y**3),
        'x^2*y'          : sp.expand(x**2*y),
        'x*y'            : sp.expand(x*y),
        'x*(x+y)'        : sp.expand(x*(x+y)),
        'x*y*(x+y)'      : sp.expand(x*y*(x+y)),
        '(x*y)^2'        : sp.expand((x*y)**2),
        'x*y*(x+y)*(x-y)': sp.expand(x*y*(x+y)*(x-y)),
    }
    ab_list = [(3, 2), (5, 2), (4, 3), (5, 3), (2, 1), (3, 1)]
    ok = True
    nchk = 0
    print()
    print("   %-16s %2s %2s %2s %2s %3s %3s | %s" %
          ("H", "g", "g0", "a", "b", "n", "m", "j: (dimker,dimcoker) pred/obs"))
    for hn, H in Hs.items():
        g = tdeg(H)
        H0, d = primitive_part(H)
        g0 = tdeg(H0)
        for (a, b) in ab_list:
            if sp.igcd(a, b) != 1 or a < b:
                continue
            n, m = g*a, g*b
            if n+m > 24:
                continue
            obs = []
            for j in range(1, min(m, n+m-3)+1):
                M = Phi_matrix(H, a, b, j)
                if M is None:
                    continue
                # structural check: the described set really lies in ker Phi_j
                for mo in basis(m-j):
                    B = mo
                    A = sp.expand(sp.Rational(a, b)*H**(a-b)*B)
                    if sp.expand(jac(A, H**b) + jac(H**a, B)) != 0:
                        ok = False
                        print("   KER-CONTAINMENT FAIL", hn, a, b, j, mo)
                if (n-j) % g0 == 0:
                    A = sp.expand(sp.Rational(1, b)*H0**((n-j)//g0))
                    if sp.expand(jac(A, H**b)) != 0:
                        ok = False
                        print("   KER-LAMBDA FAIL", hn, a, b, j)
                r = M.rank()
                dker = M.cols - r
                dcok = M.rows - r
                pker = max(0, m-j+1) + (1 if j % g0 == 0 else 0)
                pcok = (m-2) + (1 if j % g0 == 0 else 0)
                nchk += 1
                good = (dker == pker and dcok == pcok)
                if not good:
                    ok = False
                obs.append("%d:(%d,%d)%s" % (j, dker, dcok, "" if good else "!!"))
            print("   %-16s %2d %2d %2d %2d %3d %3d | %s" %
                  (hn, g, g0, a, b, n, m, " ".join(obs)))
    print()
    print("   orders checked : %d      dimension law holds everywhere : %s"
          % (nchk, ok))
    print("   PART3 PASS =", ok)
    RESULTS.append(("PART3 obstruction dimension law", ok))
    return ok


def part3b():
    """T2 / T3: the exact determination of P_(n-1), P_(n-2) by Q's tail."""
    print()
    print("PART 3b  T2 and T3, the exact determination relations")
    print("   T2 (a>=b): P_(n-1) = H^(a-b) A_1 and  c a Q_(m-1) - c' b A_1 =")
    print("              lambda * H0^((m-1)/g0), with lambda = 0 unless g0 = 1.")
    checks = []
    rng = random.Random(3)
    samples = []
    tries = 0
    while len(samples) < 22 and tries < 300:
        tries += 1
        P, Q = build_auto(rng, nfac=rng.choice([2, 3]), maxdeg=rng.choice([2, 3, 4]))
        hp, hq = hparts(P), hparts(Q)
        if not hp or not hq:
            continue
        n, m = max(hp), max(hq)
        if n+m < 8 or n+m > 34 or min(n, m) < 2:
            continue
        Jc = jac(P, Q)
        if not (Jc.is_number and Jc != 0):
            continue
        g = sp.igcd(n, m)
        a, b = n//g, m//g
        facs = linear_factors(hp[n])
        if facs is None or any(mu % a for _, mu in facs):
            continue
        H = 1
        for L, mu in facs:
            H = sp.expand(H*L**(mu//a))
        if tdeg(H) != g:
            continue
        c = sp.simplify(hp[n]/H**a)
        cp = sp.simplify(hq[m]/H**b)
        H0, dd = primitive_part(H)
        g0 = tdeg(H0)
        if a >= b:
            Pn1, Qm1 = hp.get(n-1, 0), hq.get(m-1, 0)
            quo, rem = sp.div(sp.Poly(Pn1, x, y), sp.Poly(H**(a-b), x, y)) \
                if Pn1 != 0 else (sp.Poly(0, x, y), sp.Poly(0, x, y))
            div_ok = (sp.expand(rem.as_expr()) == 0)
            A1 = sp.expand(quo.as_expr())
            Theta = sp.expand(c*a*Qm1 - cp*b*A1)
            Theta = sp.expand(Theta)
            if g0 == 1 and (m-1) % g0 == 0:
                # Theta must be a scalar multiple of H0^(m-1)
                if Theta == 0:
                    lam_ok = True
                else:
                    q2, r2 = sp.div(sp.Poly(Theta, x, y),
                                    sp.Poly(sp.expand(H0**(m-1)), x, y))
                    lam_ok = (sp.expand(r2.as_expr()) == 0 and
                              sp.expand(q2.as_expr()).is_number)
            else:
                lam_ok = (Theta == 0)
            checks.append((div_ok, lam_ok, g0, n, m, a, b))
            samples.append((n, m, g, a, b, g0))
    ok = all(d and l for d, l, *_ in checks)
    print("   samples : %d   (n,m,g,a,b,g0) : %s" % (len(samples), sorted(set(samples))))
    print("   H^(a-b) | P_(n-1) everywhere        :", all(d for d, l, *_ in checks))
    print("   T2 residue is the allowed H0-power  :", all(l for d, l, *_ in checks))
    print("   PART3b PASS =", ok)
    RESULTS.append(("PART3b T2 determination relation", ok))
    return ok


# ----------------------------------------------------------------------------
# PART 4 -- hull-matching arithmetic
# ----------------------------------------------------------------------------
def exact_quot(F, D):
    """F/D as a polynomial; returns None if the division is not exact."""
    F = sp.expand(F)
    if F == 0:
        return sp.Integer(0)
    q, r = sp.div(sp.Poly(F, x, y), sp.Poly(sp.expand(D), x, y))
    if sp.expand(r.as_expr()) != 0:
        return None
    return sp.expand(q.as_expr())


def part3c():
    """the explicit reduced order-j identities (R_j), j = 1,2,3."""
    print()
    print("PART 3c  the REDUCED tower (R_j): order j divided by H^(a-(j-1)b-1)")
    print("   (R_1)  c a Jac(H,Q_(m-1)) + c' b Jac(P_(n-1),H)/H^(a-b-... ) = 0")
    print("   (R_2)  c a H^b Jac(H,Q_(m-2)) + Jac(P_(n-1),Q_(m-1))/H^(a-b-1)")
    print("            + c' b Jac(P_(n-2),H)/H^(a-2b) = 0")
    print("   every term is a form of degree j(m-1)+g-2; valid for")
    print("   1 <= j <= floor((a-1)/b)+1  (= the first partial quotient of a/b).")
    rng = random.Random(101)
    rows, ok = [], True
    tries = 0
    while len(rows) < 18 and tries < 300:
        tries += 1
        P, Q = build_auto(rng, nfac=rng.choice([2, 3]), maxdeg=rng.choice([2, 3]))
        hp, hq = hparts(P), hparts(Q)
        if not hp or not hq:
            continue
        n, m = max(hp), max(hq)
        if n+m < 8 or n+m > 30 or min(n, m) < 2:
            continue
        Jc = jac(P, Q)
        if not (Jc.is_number and Jc != 0):
            continue
        g = sp.igcd(n, m)
        a, b = n//g, m//g
        if a < b:
            P, Q, n, m, a, b = Q, P, m, n, b, a
            hp, hq = hq, hp
        facs = linear_factors(hp[n])
        if facs is None or any(mu % a for _, mu in facs):
            continue
        H = 1
        for L, mu in facs:
            H = sp.expand(H*L**(mu//a))
        if tdeg(H) != g:
            continue
        c = sp.simplify(hp[n]/H**a)
        cp = sp.simplify(hq[m]/H**b)
        jmaxR = (a-1)//b + 1
        good = True
        for j in range(1, min(jmaxR, n+m-3)+1):
            div = sp.expand(H**(a-(j-1)*b-1))
            terms = []
            for i in range(0, j+1):
                k = j-i
                A, B = hp.get(n-i, 0), hq.get(m-k, 0)
                if A == 0 or B == 0:
                    continue
                t = exact_quot(jac(A, B), div)
                if t is None:
                    good = False
                    print("   NON-EXACT DIVISION", n, m, a, b, j, i, k)
                    break
                terms.append(t)
            if not good:
                break
            tot = sp.expand(sum(terms))
            degs = set(tdeg(t) for t in terms if t != 0)
            if tot != 0:
                good = False
                print("   (R_j) FAIL", n, m, a, b, j, tot)
            if degs and degs != {j*(m-1)+g-2}:
                good = False
                print("   DEGREE FAIL", n, m, a, b, j, degs, j*(m-1)+g-2)
        ok &= good
        rows.append((n, m, g, a, b, jmaxR, good))
    print("   %-24s %s" % ("(n,m,g,a,b,jmax_R)", "reduced identities hold"))
    for r in rows:
        print("   (%3d,%3d,%2d,%2d,%2d,%2d)        %s" % r)
    print("   PART3c PASS =", ok)
    RESULTS.append(("PART3c reduced tower (R_j)", ok))
    return ok


def part4():
    print()
    print("PART 4  single-slope hull matching  =>  sigma = e/t with 1<=t<=g-1")
    print("   For each L|H with L not dividing P_1 nor Q_1 the two hulls share")
    print("   their descent rate sigma (PART 2 'sigma_scoped').  If both hulls")
    print("   are single-edge, sigma = a e / alpha = b e / beta with alpha,beta")
    print("   integers, alpha<=n-1, beta<=m-1.  gcd(a,b)=1 forces alpha=a t,")
    print("   beta = b t, sigma = e/t, and alpha<=ga-1 forces t<=g-1.")
    ok = True
    checked = 0
    for g in range(1, 13):
        for a in range(2, 13):
            for b in range(2, 13):
                if sp.igcd(a, b) != 1 or a == b:
                    continue
                for e in range(1, g+1):
                    n, m = g*a, g*b
                    sols = []
                    for alpha in range(1, n):
                        for beta in range(1, m):
                            # a e / alpha == b e / beta  <=>  a beta == b alpha
                            if a*beta == b*alpha:
                                sols.append((alpha, beta))
                    checked += 1
                    derived = [(a*t, b*t) for t in range(1, g)
                               if a*t <= n-1 and b*t <= m-1]
                    if sorted(sols) != sorted(derived):
                        ok = False
                        print("   MISMATCH", g, a, b, e, sols, derived)
                    if g == 1 and sols:
                        ok = False
                        print("   g=1 should be empty:", g, a, b, sols)
    print("   parameter triples checked : %d" % checked)
    print("   'single-slope solutions = {(at,bt): 1<=t<=g-1}' holds everywhere :", ok)
    print("   consequence: single-slope hull matching is IMPOSSIBLE when g=1.")
    print("   PART4 PASS =", ok)
    RESULTS.append(("PART4 hull arithmetic", ok))
    return ok


# ----------------------------------------------------------------------------
# PART 5 -- the search order
# ----------------------------------------------------------------------------
def cf(a, b):
    """continued fraction of a/b (a>b>0), as list of partial quotients."""
    out = []
    while b:
        out.append(a//b)
        a, b = b, a % b
    return out


def fib_list(N=40):
    F = [1, 1]
    while len(F) < N:
        F.append(F[-1]+F[-2])
    return F


def metallic_pair(a, b, q):
    """is (a,b) a pair of consecutive terms of a_k = q a_{k-1} + a_{k-2}?"""
    s, t = 1, 0          # a_1=1, a_0=0
    seen = []
    while s <= a:
        seen.append((s, t))
        s, t = q*s+t, s
    return (a, b) in seen


def part5(nrows=44):
    print()
    print("PART 5  search order over coprime exponent pairs (a,b), a,b>=2")
    F = set(fib_list(30))
    fibpairs = set()
    Fl = fib_list(30)
    for i in range(2, len(Fl)-1):
        fibpairs.add((Fl[i+1], Fl[i]))
    cands = []
    for s in range(5, 60):
        for b in range(2, s//2+1):
            a = s-b
            if a <= b or sp.igcd(a, b) != 1:
                continue
            q = cf(a, b)
            depth = len(q)                     # number of division steps
            sub = sum(q)                       # subtractive depth
            qmax = max(q)
            qmin = min(q)
            met = None
            for qq in range(1, 7):
                if metallic_pair(a, b, qq):
                    met = qq
                    break
            cands.append(dict(a=a, b=b, s=s, cf=q, depth=depth, sub=sub,
                              qmax=qmax, qmin=qmin, met=met,
                              fib=((a, b) in fibpairs)))
    # search order: hardest corner first == smallest Zaremba quantity
    # z = max partial quotient, then smallest a+b.
    cands.sort(key=lambda c: (c['qmax'], c['s'], c['a']))
    print("   ordering key = (z = max partial quotient ASC, then a+b ASC).")
    print("   Rationale: the PROVED forced chain  H^(a-ib) | P_(n-i),")
    print("   0<=i<=floor(a/b), is longest when a partial quotient is large, so")
    print("   large-quotient pairs are the MOST constrained and the small-z")
    print("   (golden / metallic) pairs are the hard corner.  z=2 with all")
    print("   leading quotients 1 is exactly the Fibonacci ray.")
    print("   'chain' = floor(a/b); 'sub' = sum of partial quotients = total")
    print("   number of subtractive Euclid steps the tower must run.")
    print("   g_min: JC(2) is classical for gcd(n,m)=1 (and PART 4 kills g=1 for")
    print("   single-slope hulls); CITED classical filters gcd<=8 and Moh's")
    print("   max(deg)<=100 give the last two columns.")
    print()
    hdr = ("  #    a    b  a/b CF          z  dep  sub chain metal fib  "
           "minimal (n,m):  g=2        g=9         Moh")
    print(hdr)
    print("  " + "-"*(len(hdr)-2))
    for i, c in enumerate(cands[:nrows]):
        a, b = c['a'], c['b']
        gmoh = max(9, -(-101//a))
        met = ("%d" % c['met']) if c['met'] else "-"
        print("  %-3d %3d %4d  %-13s %2d %3d %4d %4d  %-4s %-4s "
              "  (%3d,%3d) (%4d,%4d) g=%-3d(%4d,%4d)"
              % (i+1, a, b, "[" + ";".join(str(t) for t in c['cf']) + "]",
                 c['qmax'], c['depth'], c['sub'], a//b, met,
                 "Y" if c['fib'] else "-",
                 2*a, 2*b, 9*a, 9*b, gmoh, gmoh*a, gmoh*b))
    print()
    print("   the golden (Fibonacci) ray -- maximal CF depth per size, z=2:")
    Fl2 = fib_list(14)
    for i in range(3, len(Fl2)-1):
        a, b = Fl2[i+1], Fl2[i]
        if sp.igcd(a, b) != 1 or b < 2:
            continue
        q = cf(a, b)
        gmoh = max(9, -(-101//a))
        print("      (a,b)=(%3d,%3d)  CF=%-22s depth=%d sub=%d  Moh-minimal "
              "(n,m)=(%d,%d) at g=%d"
              % (a, b, "[" + ";".join(str(t) for t in q) + "]", len(q), sum(q),
                 gmoh*a, gmoh*b, gmoh))
    print()
    print("   metallic rays a_k = q a_(k-1) + a_(k-2)  (THM-3010 alternation):")
    for qq in range(1, 5):
        s, t = 1, 0
        seq = []
        while s <= 200:
            seq.append(s)
            s, t = qq*s+t, s
        pairs = [(seq[i+1], seq[i]) for i in range(1, len(seq)-1)
                 if seq[i] >= 2 and sp.igcd(seq[i+1], seq[i]) == 1]
        print("      q=%d (%s): %s" %
              (qq, {1: 'golden', 2: 'silver', 3: 'bronze'}.get(qq, 'q-metallic'),
               pairs[:7]))
    RESULTS.append(("PART5 search order table", True))
    return cands


if __name__ == '__main__':
    part1()
    part2()
    part3()
    part3b()
    part3c()
    part4()
    part5()
    print()
    print("=== SUMMARY ===")
    allok = True
    for name, ok in RESULTS:
        print("   %-42s %s" % (name, ok))
        allok &= bool(ok)
    print("   ALL =", allok)
