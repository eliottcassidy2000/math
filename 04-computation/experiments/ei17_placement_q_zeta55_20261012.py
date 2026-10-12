"""Can the 17-vertex Exoo-Ismailescu graph EI17 (HoG 51375) be placed in Q(zeta55)?

Run from the repository root:  python3 -B <script>            (fast part, ~1 minute)
                               python3 -B <script> --census 96 (adds the real-realization census, ~10 min)
Dependencies: numpy, scipy, sympy (with python-flint for fast integer LLL), mpmath.

Contents.
1. Graph data: the House of Graphs adjacency list (API, 2026-10-12) and the edge list of
   arXiv:2607.19995 (Costa Silva Filho, July 2026), with the isomorphism between the two
   labellings (|Aut(EI17)| = 2).
2. The certified realization of that paper (Appendix A, 24 digits, gauge v17 = 0, v11 = -1):
   all 31 edges unit; Newton refinement to 3000 digits in the full 30-equation system.
3. Structure: the nine unit 4-cycles are parallelograms in every faithful realization, so the
   31 edges carry 13 unit directions with 6 independent integer relations (6 free angles, no
   sign branches).
4. The paper's degree-20 polynomial p (Appendix B, coefficients parsed from the PDF text):
   irreducible over Q, Frobenius cycle types (20) at 3 and (3,17) at 71, hence Gal(p) = S20
   (exact facts about p) -- but p does NOT vanish at the certified x_v1: its nearest real root
   differs at the eighth significant digit.  Integer LLL at 3000 digits finds no polynomial of
   degree <= 60 vanishing at x_v1 below the random height; in particular no degree-20
   relation of height <= 1e128.  So p is not the minimal polynomial of the coordinate.
5. Placement test: in the gauge, EI17 lies in L = Q(zeta55) iff all 13 directions lie in L.
   Integer LLL on the 40-dimensional power basis at 2900 digits finds, for each of the twelve
   non-gauge directions, only vectors of the random height ~1e69: no placement with
   coordinate height below about 1e63.
6. (--census g1,g2,...) Real faithful realizations by the two-angle apex construction over all
   2^13 sign branches, seeded from grids over the two angles and Newton-polished in the full
   system, merged over the grids by distance multiset (labelling-independent).  With the
   seeding threshold 0.12, grids 96 and 160 find ten each and eleven in the union; the count
   grows with the seeding and is a lower bound: at least eleven non-congruent real faithful
   realizations against the paper's "more than one".  Run time about 35 minutes for 96,160.
"""
import itertools, json, sys, time
from fractions import Fraction
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix
from scipy.optimize import least_squares

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


# ---------------------------------------------------------------- 1. graph data
HOG = [[11, 12, 15], [10, 12, 16], [11, 13, 15], [10, 14, 16], [9, 13, 15], [8, 14, 16], [7, 9, 14, 15],
       [6, 8, 13, 16], [5, 7, 10, 12], [4, 6, 11, 12], [1, 3, 8, 11], [0, 2, 9, 10], [0, 1, 8, 9],
       [2, 4, 7, 14], [3, 5, 6, 13], [0, 2, 4, 6], [1, 3, 5, 7]]
PE = [(1, 4), (1, 7), (1, 8), (2, 3), (2, 6), (2, 10), (3, 5), (3, 7), (3, 12), (4, 6), (4, 13), (4, 14), (5, 8),
      (5, 9), (5, 14), (6, 11), (6, 16), (7, 9), (7, 15), (8, 13), (8, 15), (9, 10), (9, 11), (10, 12), (10, 16),
      (11, 12), (11, 17), (13, 17), (14, 15), (14, 17), (16, 17)]
# Appendix A of arXiv:2607.19995 (24 significant digits; gauge v17 = (0,0), v11 = (-1,0))
TABLE = {1: ('-2.498806168621063073781172', '-0.549687325764897172971352'),
         2: ('-2.385349777286040586649497', '-0.087042289116627046343162'),
         3: ('-1.976688116468964503756166', '-0.999728241116672384265386'),
         4: ('-1.855664803158276574864903', '0.216060141766444058530457'),
         5: ('-1.760489090853103061928075', '-0.023378927949119585938183'),
         6: ('-1.647529953519117279782757', '-0.762039998487959826837054'),
         7: ('-1.624860686432937524721422', '-0.063663361167507460404979'),
         8: ('-1.612262559956085560912846', '-1.012332362413167299599543'),
         9: ('-1.408661660817076082893331', '0.912685952000045337922224'),
         10: ('-1.385349777286040586649497', '-0.087042289116627046343162'),
         11: ('-1', '0'),
         12: ('-0.976688116468964503756166', '-0.999728241116672384265386'),
         13: ('-0.969121194493299061996577', '-0.246584894881826068097732'),
         14: ('-0.886543608664977512868325', '0.462645036648270126628190'),
         15: ('-0.738317077767960011853096', '-0.526308397815777587033169'),
         16: ('-0.647529953519117279782757', '-0.762039998487959826837054'),
         17: ('0', '0')}
# Appendix B of arXiv:2607.19995, parsed from the PDF text (two transcriptions agree)
P20 = {20: 10846277150822505740639962570558213942763119983043461674845448343918212271056651745480014057747467,
       19: 19180475212539477649458780252280766370785104165204168446975900051667169039146647895248310235378771,
       18: -26070936888654617384259635465342351757876050482308610020563950635118095256174246419408065866234597,
       17: -3962009254586100531018427473429740363965658407149928691646435897841056035927607923182216019350940,
       16: 18166859496891545508566917411230114971902727740304713885900161305285905262350527321353616634759110,
       15: -7607431700536231966080908557596836755687226593135326666832499865701187554080177241394973984057350,
       14: 45351587309956624267644739061100632322504856432748648077436547828253855443955203136258124315668464,
       13: -9408123675708568099217732631464064702259192795720169395209699138778275488455547776691534292051591,
       12: 14682673355886588168945506517818998475845345062830936736552073931495949256607496993509908842889085,
       11: 15143426019910260147689787064208205709627636321450044137594530613944920637911243716728461040656088,
       10: 65968618518965545764431180875089192219762433065476009078893452297224262471375559826743563580939,
       9: 6933435324110256159086227291350057737956789050172029740276872060533462997921572243016990516167855,
       8: -9790517400576717511906222722622431562355763122979780641991852086843847078118813130413048942134112,
       7: 6244397981079102575268653077249642024975028707259943515274717143990817911425984104440579854365874,
       6: -23586454334648278938656158865770452705897027029874210958571764632770769116482786067681998817533516,
       5: -53392439028840081073638800831628946914310140652692035611147137997008514386412828676786323411211286,
       4: -26449839798917337714527054796537315591594669043182637188470470315247453343648348898405567316167209,
       3: -67496405170016508230749607293466456667836751966863744185909998225730057226496683094528284790742964,
       2: 29801792463858164957001377196524717774830567182899091859788291738590738838397086829441065048179864,
       1: 0,
       0: -1312711393788876800678365780809286777142912632161468883999549147019329458866534946330375157169324}
COEFFS = [P20[k] for k in range(20, -1, -1)]  # high -> low


def isomorphisms():
    padj = {v: set() for v in range(1, 18)}
    for i, j in PE:
        padj[i].add(j); padj[j].add(i)
    hadj = {i: set(n) for i, n in enumerate(HOG)}
    order = sorted(range(1, 18), key=lambda v: -len(padj[v]))
    found = []

    def extend(m):
        if len(m) == 17:
            found.append(dict(m)); return
        v = order[len(m)]
        for h in range(17):
            if h in m.values() or len(hadj[h]) != len(padj[v]):
                continue
            if all((m[w] in hadj[h]) == (w in padj[v]) for w in m):
                m[v] = h; extend(m); del m[v]
    extend({})
    return found


# ------------------------------------------------------- 2. certified realization
FREEV = [v for v in range(1, 18) if v not in (11, 17)]
E30 = [e for e in PE if e != (11, 17)]


def unpack(X):
    z = {11: (mp.mpf(-1), mp.mpf(0)), 17: (mp.mpf(0), mp.mpf(0))}
    for k, v in enumerate(FREEV):
        z[v] = (X[2*k], X[2*k+1])
    return z


def residual(X):
    z = unpack(X)
    return mp.matrix([(z[i][0]-z[j][0])**2+(z[i][1]-z[j][1])**2-1 for i, j in E30])


def jacobian(X):
    z = unpack(X); M = mp.matrix(30, 30); col = {v: k for k, v in enumerate(FREEV)}
    for r, (i, j) in enumerate(E30):
        dx = 2*(z[i][0]-z[j][0]); dy = 2*(z[i][1]-z[j][1])
        if i in col:
            M[r, 2*col[i]] += dx; M[r, 2*col[i]+1] += dy
        if j in col:
            M[r, 2*col[j]] -= dx; M[r, 2*col[j]+1] -= dy
    return M


def refine(digits):
    mp.mp.dps = digits+50
    X = mp.matrix([mp.mpf(TABLE[v][c]) for v in FREEV for c in (0, 1)])
    for _ in range(25):
        r = residual(X)
        if mp.norm(r) < mp.mpf(10)**(-digits-20):
            break
        X = X-mp.lu_solve(jacobian(X), r)
    check(mp.norm(residual(X)) < mp.mpf(10)**(-digits-20), 'Newton refinement converged')
    z = unpack(X)
    for v in TABLE:  # the 24 printed digits are reproduced
        for c in (0, 1):
            check(abs(z[v][c]-mp.mpf(TABLE[v][c])) < mp.mpf(10)**-23, 'Appendix A digits reproduced')
    return z


# ---------------------------------------------------------- 3. direction structure
def four_cycles():
    adj = {v: set() for v in range(1, 18)}
    for i, j in PE:
        adj[i].add(j); adj[j].add(i)
    cyc = set()
    for a in range(1, 18):
        for b in adj[a]:
            for c in adj[b]:
                if c == a:
                    continue
                for d in adj[c]:
                    if d in (a, b) or a not in adj[d]:
                        continue
                    t = (a, b, c, d)
                    rots = [t[i:]+t[:i] for i in range(4)]
                    rots += [tuple(reversed(r)) for r in rots]
                    cyc.add(min(rots))
    return sorted(cyc)


def direction_structure(z):
    """Edge direction classes (up to sign) of the certified realization and their integer relations."""
    P = {v: mp.mpc(z[v][0], z[v][1]) for v in z}
    reps, cls = [], {}
    for i, j in PE:
        u = P[j]-P[i]
        for k, r in enumerate(reps):
            if abs(u-r) < mp.mpf(10)**-40:
                cls[(i, j)] = (k, 1); break
            if abs(u+r) < mp.mpf(10)**-40:
                cls[(i, j)] = (k, -1); break
        else:
            reps.append(u); cls[(i, j)] = (len(reps)-1, 1)
    n = len(reps)
    U = sp.symbols('u1:%d' % (n+1))
    pos = {17: sp.Integer(0)}; adj = {v: [] for v in range(1, 18)}
    for i, j in PE:
        adj[i].append(j); adj[j].append(i)
    tree = set()
    while len(pos) < 17:
        for v in list(pos):
            for w in adj[v]:
                if w in pos:
                    continue
                e = (v, w) if (v, w) in cls else (w, v)
                k, s = cls[e]
                pos[w] = sp.expand(pos[v]+(s*U[k] if e == (v, w) else -s*U[k])); tree.add(e)
    rels = [sp.expand(pos[j]-pos[i]-cls[(i, j)][1]*U[cls[(i, j)][0]]) for (i, j) in PE if (i, j) not in tree]
    rels = [r for r in rels if r != 0]
    M = sp.Matrix([[r.coeff(u) for u in U] for r in rels])
    rank = M.rank()
    gauge = [k for k, r in enumerate(reps) if abs(r-1) < mp.mpf(10)**-40]
    return reps, cls, rank, gauge


# ------------------------------------------------------- 4. the paper's polynomial
def polymod(a, m, p):
    a = [x % p for x in a]; dm = len(m)-1
    for i in range(len(a)-1, dm-1, -1):
        c = a[i]
        if c:
            for j in range(dm+1):
                a[i-dm+j] = (a[i-dm+j]-c*m[j]) % p
    a = a[:dm]
    while a and a[-1] == 0:
        a.pop()
    return a


def polymul(a, b, p):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] = (r[i+j]+x*y) % p
    return r


def powmod(base, e, m, p):
    r = [1]; b = polymod(base, m, p)
    while e:
        if e & 1:
            r = polymod(polymul(r, b, p), m, p)
        b = polymod(polymul(b, b, p), m, p); e >>= 1
    return r


def polygcd(a, b, p):
    a = [x % p for x in a]; b = [x % p for x in b]
    while True:
        while b and b[-1] == 0:
            b.pop()
        if not b:
            break
        inv = pow(b[-1], -1, p); bm = [(x*inv) % p for x in b]
        a = polymod(a, bm, p) if len(a) >= len(bm) else a
        a, b = b, a
    while a and a[-1] == 0:
        a.pop()
    return a


def cycle_type(coeffs_high_to_low, p):
    """Degrees of the irreducible factors mod p (None if p is bad), by distinct-degree factorization."""
    f = [c % p for c in coeffs_high_to_low[::-1]]
    if f[-1] == 0:
        return None
    inv = pow(f[-1], -1, p); f = [(x*inv) % p for x in f]
    df = [(i*f[i]) % p for i in range(1, len(f))]
    if len(polygcd(f, df, p)) > 1:
        return None
    degs, g, h, k = [], f[:], [0, 1], 0
    while len(g) > 1 and k < len(f):
        k += 1
        h = powmod(h, p, f, p)
        hk = h[:]+[0]*(2-len(h)) if len(h) < 2 else h[:]
        hk[1] = (hk[1]-1) % p
        d = polygcd(g, hk, p)
        if len(d) > 1:
            degs += [k]*((len(d)-1)//k)
            dinv = pow(d[-1], -1, p); dm = [(x*dinv) % p for x in d]
            q = [0]*(len(g)-len(dm)+1); r = g[:]
            for i in range(len(g)-1, len(dm)-2, -1):
                c = r[i]
                if c:
                    q[i-len(dm)+1] = c
                    for j in range(len(dm)):
                        r[i-len(dm)+1+j] = (r[i-len(dm)+1+j]-c*dm[j]) % p
            g = q
            while g and g[-1] == 0:
                g.pop()
    return tuple(sorted(degs))


def roots_mod(coeffs_high_to_low, p):
    return sum(1 for x in range(p) if not sum(c % p*pow(x, k, p) for k, c in enumerate(coeffs_high_to_low[::-1])) % p)


def lll_relation(values, digits):
    """Shortest LLL vector of the lattice (e_i, round(values_i * 10^digits)); returns coefficients."""
    S = mp.mpf(10)**digits
    n = len(values)
    rows = []
    for i in range(n):
        r = [0]*n; r[i] = 1; r.append(int(mp.nint(values[i]*S))); rows.append(r)
    M = DomainMatrix([[ZZ(v) for v in r] for r in rows], (n, n+1), ZZ)
    return [int(c) for c in M.lll().to_list()[0]][:n]


def paper_polynomial(z):
    x = sp.symbols('x')
    P = sp.Poly(COEFFS, x)
    check(P.degree() == 20 and sp.gcd_list(COEFFS) == 1, 'degree 20 primitive')
    fl = sp.factor_list(P)[1]
    check(len(fl) == 1 and fl[0][0].degree() == 20, 'irreducible over Q')
    types = {}
    for p in (3, 5, 7, 29, 71, 109):
        t = cycle_type(COEFFS, p)
        types[p] = t
        check(t is not None and sum(t) == 20, 'good prime with squarefree reduction')
        check(t.count(1) == roots_mod(COEFFS, p), 'root count agrees with the distinct-degree factorization')
    check(types[3] == (20,) and types[71] == (3, 17), 'a 20-cycle at 3 and a 17-cycle at 71')
    # Jordan: transitive (irreducible) + 17-cycle with 17 <= 20-3 => contains A20; the 20-cycle is odd => S20.
    mp.mp.dps = 90
    roots = mp.polyroots(COEFFS, maxsteps=500, extraprec=3000)
    real = sorted(r.real for r in roots if abs(r.imag) < mp.mpf(10)**-60)
    x1 = z[1][0]
    nearest = min(real, key=lambda r: abs(r-x1))
    gap = abs(nearest-x1)
    check(mp.mpf('8e-8') < gap < mp.mpf('9e-8'), 'nearest real root of p is 8.8e-8 from the certified x_v1')
    rel = abs(sum(c*x1**(20-k) for k, c in enumerate(COEFFS)))/abs(COEFFS[0]*x1**20)
    check(mp.mpf('3e-8') < rel < mp.mpf('5e-8'), 'relative |p(x_v1)| about 3.8e-8, not zero')
    # all 34 table coordinates and simple transforms: none is a root of p
    rels = []
    for v in TABLE:
        for val in (z[v][0], z[v][1], -z[v][0], -z[v][1]):
            if val == 0:
                continue
            rels.append(float(abs(sum(c*val**(20-k) for k, c in enumerate(COEFFS)))/(abs(COEFFS[0]*val**20)+1)))
    check(min(rels) > 1e-9, 'no coordinate of the realization is a root of p')
    return dict(irreducible=True, cycle_types={p: list(t) for p, t in types.items()}, galois_group='S20 (exact)',
                real_roots=[mp.nstr(r, 30) for r in real], nearest_root_gap=mp.nstr(gap, 5),
                relative_p_at_x_v1=mp.nstr(rel, 5), min_relative_over_coordinates=min(rels))


def algdep_search(z):
    """Integer LLL for polynomials vanishing at x_v1 at 2900 digits, degrees 20..60."""
    mp.mp.dps = 3010
    x1 = z[1][0]
    out = {}
    for d in (20, 24, 30, 40, 48, 60):
        pw = [mp.mpf(1)]
        for _ in range(d):
            pw.append(pw[-1]*x1)
        v = lll_relation(pw, 2900)
        h = max(abs(c) for c in v)
        out[d] = dict(height_digits=len(str(h))-1, random_expectation_digits=2900//(d+2))
        check(out[d]['height_digits'] >= out[d]['random_expectation_digits']-2, 'no relation below the random height')
    return out


# ------------------------------------------------------------- 5. placement test
def placement_test(z, reps, cls):
    mp.mp.dps = 3010
    zeta = mp.exp(2*mp.pi*1j/55); ZP = [zeta**k for k in range(40)]; lam = mp.pi
    basis = [mp.re(w)+lam*mp.im(w) for w in ZP]
    out = {}
    for k, d in enumerate(reps):
        v = lll_relation(basis+[mp.re(d)+lam*mp.im(d)], 2900)
        h = max(abs(c) for c in v)
        out['u%d' % (k+1)] = len(str(h))-1
    # the gauge direction 1 is trivially in L; every other direction has only random-height vectors
    check(sum(1 for k in out if out[k] == 0) == 1, 'exactly the gauge direction has an exact relation')
    check(all(out[k] >= 2900//42-2 for k in out if out[k] > 0), 'no non-gauge direction has a relation below the random height')
    return out


# ------------------------------------------------------------ 6. real census
ORDER = [(16, 6, 17), (10, 2, 16), (12, 10, 11), (3, 2, 12), (9, 10, 11), (7, 3, 9), (5, 3, 9), (14, 5, 17),
         (4, 6, 14), (13, 4, 17), (8, 5, 13), (1, 4, 7), (15, 7, 8)]
CLOSURES = [(1, 8), (14, 15)]


def tower(t4, t5, signs):
    z = {17: np.zeros_like(t4, dtype=complex), 11: -np.ones_like(t4, dtype=complex)}
    z[6] = z[11]+np.exp(1j*t4); z[2] = z[6]+np.exp(1j*t5)
    for (v, a, b), s in zip(ORDER, signs):
        d = z[b]-z[a]; D = np.abs(d)**2
        h = np.sqrt(np.maximum(1-D/4, 0)); n = 1j*d/np.sqrt(np.where(D > 0, D, 1))
        z[v] = np.where(D > 4, np.nan+0j, z[a]+d/2+s*h*n)
    return z


def np_unpack(x):
    z = {11: np.array([-1.0, 0.0]), 17: np.array([0.0, 0.0])}
    for k, v in enumerate(FREEV):
        z[v] = x[2*k:2*k+2]
    return z


def np_resid(x):
    z = np_unpack(x)
    return np.array([((z[i]-z[j])**2).sum()-1 for i, j in E30])


def census(ngrid, found=None, threshold=0.12):
    found = {} if found is None else found
    t = np.linspace(0, 2*np.pi, ngrid, endpoint=False)
    T4, T5 = np.meshgrid(t, t, indexing='ij'); T4 = T4.ravel(); T5 = T5.ravel()
    for signs in itertools.product((1, -1), repeat=13):
        z = tower(T4, T5, signs)
        c1 = np.abs(z[8]-z[1])**2-1; c2 = np.abs(z[15]-z[14])**2-1
        ok = np.isfinite(c1) & np.isfinite(c2)
        for k in np.where(ok & (np.abs(c1) < threshold) & (np.abs(c2) < threshold))[0]:
            zz = tower(np.array([T4[k]]), np.array([T5[k]]), signs)
            x0 = np.array([[zz[v][0].real, zz[v][0].imag] for v in FREEV]).ravel()
            if not np.all(np.isfinite(x0)):
                continue
            r = least_squares(np_resid, x0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=200)
            if r.cost < 1e-22:
                zs = np_unpack(r.x)
                P = np.array([zs[v] for v in range(1, 18)])
                D = np.sqrt(((P[:, None, :]-P[None, :, :])**2).sum(2)); iu = np.triu_indices(17, 1)
                if D[iu].min() < 1e-6:
                    continue
                sig = tuple(np.round(np.sort(D[iu]), 6))
                if sig not in found:
                    found[sig] = dict(hits={}, v1=zs[1].tolist(), mindist=float(D[iu].min()), signs=signs,
                                      coordinates={v: zs[v].tolist() for v in zs})
                found[sig]['hits'][ngrid] = found[sig]['hits'].get(ngrid, 0)+1
    return found


def main(census_grid):
    out = {}
    isos = isomorphisms()
    check(len(isos) == 2, 'exactly two isomorphisms paper -> HoG (|Aut| = 2)')
    out['isomorphisms'] = isos
    z = refine(3000)
    out['refined_digits'] = 3000
    cycles = four_cycles()
    check(len(cycles) == 9, 'nine unit 4-cycles')
    P = {v: mp.mpc(z[v][0], z[v][1]) for v in z}
    for a, b, c, d in cycles:
        check(abs((P[b]-P[a])-(P[c]-P[d])) < mp.mpf(10)**-2900, 'each 4-cycle is a parallelogram')
    reps, cls, rank, gauge = direction_structure(z)
    check(len(reps) == 13 and rank == 6 and len(gauge) == 1, '13 directions, 6 relations, one gauge direction')
    out['four_cycles'] = cycles
    out['directions'] = dict(classes=13, independent_relations=rank, free_angles=13-1-rank)
    out['paper_polynomial'] = paper_polynomial(z)
    out['algdep_lll_3000'] = algdep_search(z)
    out['placement_lll_2900'] = placement_test(z, reps, cls)
    if census_grid:
        found = {}
        for g in census_grid:
            census(g, found)
        out['census'] = dict(grids=census_grid, non_congruent_real_faithful_at_least=len(found),
                             realizations=[dict(hits={str(k): v for k, v in f['hits'].items()}, v1=f['v1'], mindist=f['mindist'],
                                                signs=list(f['signs'])) for f in found.values()])
        # labelling-independent membership test: the certified point set's distance multiset
        Pc = np.array([[float(TABLE[v][0]), float(TABLE[v][1])] for v in range(1, 18)])
        Dc = np.sqrt(((Pc[:, None, :]-Pc[None, :, :])**2).sum(2))
        sig_c = tuple(np.round(np.sort(Dc[np.triu_indices(17, 1)]), 6))
        check(sig_c in found, 'the certified realization is among them (by distance multiset)')
        out['census']['certified_entry_hits'] = {str(k): v for k, v in found[sig_c]['hits'].items()}
    out['checks'] = CHECKS
    return out


if __name__ == '__main__':
    grids = [int(g) for g in sys.argv[sys.argv.index('--census')+1].split(',')] if '--census' in sys.argv else []
    print(json.dumps(main(grids), indent=1, sort_keys=True, default=str))
