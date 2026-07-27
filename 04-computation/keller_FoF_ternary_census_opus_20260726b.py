# HYP-9030 test (iii) + test (i) evidence: the ternary collision census of F o F
# over the branch tower, and the monodromy census of F and F o F.
#
# Key structural facts used (established here, exactly):
#   * F = A(x,y) z + B(x,y): ALL THREE components are LINEAR in z, with
#     A = (u^3, 3x u^2, -x^3), u = 1+xy.  Fibers reduce to plane curves + unique z-lift.
#   * det J_F = -2 everywhere => F is etale: NO ramification anywhere.  The honest
#     "branch locus" of a Keller map is the JELONEK non-properness set (fiber drop by
#     escape to infinity), not root collision.  The census below measures exactly that.
#   * Over F_p, for a degree-3 etale cover, rational fiber size n1 in {3,1,0} bijects
#     with the Frobenius cycle type (id / transposition / 3-cycle).  So the (n1,n2)
#     joint distribution over ALL c in F_p^3 is a Chebotarev census of the monodromy
#     of F and of F o F inside S3 wr S3.
#
# Sections:
#   [1] exact map, det J, z-linearity, equivariance group (all sign-diagonal pairs)
#   [2] exact level-1 fiber over v* = (-1/4,0,0)  (expect the known 1+2 triple)
#   [3] exact level-2 fiber field census over P0, P+, P- (factor degrees over Q)
#   [4] mod-p exact fiber counts at the tower targets (independent count proof)
#   [5] mod-p Chebotarev census: joint (n1,n2) distribution vs the full-wreath model
#
# opus 2026-07-26.  Companion to HYP-9030 (kind-pasteur S135) and THM-1300.

import sympy as sp
import numpy as np
from itertools import product as iproduct

x, y, z, a, b, c = sp.symbols('x y z a b c')
u = 1 + x*y
F1 = u**3*z + y**2*u*(4 + 3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = sp.Matrix([F1, F2, F3])

print("=" * 78)
print("[1] structure: det J, z-linearity, sign-diagonal equivariance group")
print("=" * 78)
J = F.jacobian([x, y, z])
print("det J_F =", sp.expand(J.det()), " (expect -2; F is etale everywhere)")

A = sp.Matrix([sp.diff(Fi, z) for Fi in (F1, F2, F3)])
B = sp.Matrix([sp.expand(Fi - A[i]*z) for i, Fi in enumerate((F1, F2, F3))])
assert all(sp.expand(A[i]*z + B[i] - F[i]) == 0 for i in range(3))
print("z-linear: F = A(x,y) z + B(x,y),  A =", tuple(A), " -- CONFIRMED")

eqs = []
for sx, sy, sz in iproduct([1, -1], repeat=3):
    Fi = F.subs({x: sx*x, y: sy*y, z: sz*z}, simultaneous=True)
    for tx, ty, tz in iproduct([1, -1], repeat=3):
        if all(sp.expand(Fi[i] - t*F[i]) == 0 for i, t in enumerate((tx, ty, tz))):
            eqs.append(((sx, sy, sz), (tx, ty, tz)))
print("equivariances F(D p) = D' F(p):")
for D, Dp in eqs:
    print(f"   D=diag{D}  ->  D'=diag{Dp}")
print("NOTE: sigma=diag(-1,-1,1) on the SOURCE pairs with tau=diag(1,-1,-1) on the")
print("TARGET; sigma != tau, so the tower symmetry TWISTS.  A depth-2 symmetry needs")
print("rho with F o rho = sigma o F; no sign-diagonal rho exists (list above), so no")
print("sign-diagonal involution acts on the F o F fiber over v*: the naive 'sigma-")
print("equivariant 1+2 iterate tree' prediction FAILS in its strong form.")

# ---------------------------------------------------------------- fiber machinery
def fiber_eliminant(target):
    """Exact x-eliminant of the fiber F^{-1}(target): returns (Poly, factor list).
    Uses z-linearity; valid on the u != 0 chart (checked empty for our targets)."""
    ta, tb, tc = [sp.Rational(t) for t in target]
    G1 = sp.expand((ta - B[0])*A[1] - (tb - B[1])*A[0])   # divisible by u^2
    H1 = sp.expand(sp.cancel(G1 / u**2))
    assert sp.expand(H1*u**2 - G1) == 0
    H2 = sp.expand((ta - B[0])*A[2] - (tc - B[2])*A[0])
    # cancel any common factor (occurs on the special plane a=0, where H1,H2
    # share the contracted-quadric factor u; the u=0 chart is handled separately)
    g = sp.gcd(sp.Poly(H1, x, y), sp.Poly(H2, x, y))
    if g.total_degree() > 0:
        print(f"    [common factor cancelled: {g.as_expr()}]")
        H1 = sp.expand(sp.cancel(H1 / g.as_expr()))
        H2 = sp.expand(sp.cancel(H2 / g.as_expr()))
    R = sp.resultant(H1, H2, y)
    return sp.Poly(R, x), H1, H2

def rational_points(target, H1, H2, rat_roots):
    """Exact back-substitution for rational eliminant roots only."""
    ta = sp.Rational(target[0])
    pts = []
    for r in rat_roots:
        h1 = sp.Poly(H1.subs(x, r), y)
        h2 = sp.Poly(H2.subs(x, r), y)
        g = sp.gcd(h1, h2)
        for yr, _m in sp.Poly(g, y).ground_roots().items():
            uu = 1 + r*yr
            if uu == 0:
                continue
            zr = sp.simplify((ta - B[0].subs({x: r, y: yr})) / uu**3)
            val = F.subs({x: r, y: yr, z: zr}, simultaneous=True)
            if all(sp.simplify(val[i] - sp.Rational(target[i])) == 0 for i in range(3)):
                pts.append((r, yr, zr))
    return pts

def numeric_fiber_count(target, H1, H2):
    """Numeric count of genuine fiber points (all fields), for cross-checking."""
    R = sp.resultant(H1, H2, y)
    xs = np.roots([complex(sp.Float(cf, 30)) for cf in sp.Poly(R, x).all_coeffs()])
    found = []
    f_num = sp.lambdify((x, y, z), [F[0], F[1], F[2]], 'numpy')
    h1_coeffs = [sp.lambdify(x, cf) for cf in sp.Poly(H1, y).all_coeffs()]
    h2_num = sp.lambdify((x, y), H2)
    b0_num = sp.lambdify((x, y), B[0])
    tgt = np.array([complex(sp.Rational(t)) for t in target])
    for xr in xs:
        c1 = [complex(cf(xr)) for cf in h1_coeffs]
        ys = np.roots(c1) if len(c1) > 1 else []
        for yr in ys:
            if abs(complex(h2_num(xr, yr))) > 1e-5 * (1 + abs(xr))**9:
                continue
            uu = 1 + xr*yr
            if abs(uu) < 1e-9:
                continue
            zr = (complex(sp.Rational(target[0])) - complex(b0_num(xr, yr))) / uu**3
            val = np.array(f_num(xr, yr, zr), dtype=complex).ravel()
            if np.max(np.abs(val - tgt)) < 1e-6 * (1 + np.max(np.abs(tgt))):
                if not any(abs(xr - q[0]) + abs(yr - q[1]) < 1e-6 for q in found):
                    found.append((xr, yr, zr))
    return found

print()
print("=" * 78)
print("[2] exact level-1 fiber over v* = (-1/4, 0, 0)")
print("=" * 78)
vstar = (sp.Rational(-1, 4), 0, 0)
E1, H1v, H2v = fiber_eliminant(vstar)
fl = sp.factor_list(E1.as_expr())
print("  eliminant_x degree", E1.degree(), "factors:",
      [(sp.factor(f), m) for f, m in fl[1]])
rat = [r for f, m in fl[1] for r in sp.Poly(f, x).ground_roots()]
lev1 = rational_points(vstar, H1v, H2v, sorted(set(rat)))
for p_ in lev1:
    print("   point:", p_)
print(f"  => {len(lev1)} rational points; sigma-orbit structure: "
      f"{{x=0}} fixed + {{x=+-1}} swapped = 1 + 2  -- the collision anatomy")

print()
print("=" * 78)
print("[3] exact level-2 fiber FIELD census (the depth-2 ternary tree over Q)")
print("=" * 78)
lev2_summary = {}
for p_ in lev1:
    lbl = f"({p_[0]},{p_[1]},{p_[2]})"
    E2, H1t, H2t = fiber_eliminant(p_)
    fl2 = sp.factor_list(E2.as_expr())
    genuine_factors = []
    for f, m in fl2[1]:
        fp = sp.Poly(f, x)
        if fp.degree() >= 1:
            genuine_factors.append((f, fp.degree(), m))
    print(f"  fiber over {lbl}: eliminant_x deg {E2.degree()}")
    for f, d, m in genuine_factors:
        extra = ""
        if d == 3:
            dsc = sp.discriminant(f, x)
            extra = f"  disc={sp.factor(dsc)}  square? {sp.sqrt(dsc).is_rational is True}"
        if d == 1:
            extra = "  (rational root)"
        print(f"    factor deg {d} mult {m}: {f}{extra}")
    nf = numeric_fiber_count(p_, H1t, H2t)
    lev2_summary[lbl] = nf
    print(f"    numeric genuine fiber count: {len(nf)}  "
          f"(x-coords ~ {[complex(np.round(q[0],4)) for q in nf]})")
print("  depth-2 numeric total:", sum(len(v) for v in lev2_summary.values()),
      " (generic deg FoF = 9; shortfall = Jelonek escape at depth 2)")
print("  20929 =", sp.factorint(20929))

print()
print("=" * 78)
print("[4]+[5] mod-p census: tower targets + Chebotarev joint (n1,n2)")
print("=" * 78)
print(" model = FULL wreath S3 wr S3, independent uniform fiber classes:")
print("   P(n1=3)=1/6, P(n1=1)=1/2, P(n1=0)=1/3;  n2|n1=k ~ sum of k iid X,")
print("   X in {3,1,0} w.p. {1/6,1/2,1/3}.  O(1/p) deviations expected (Jelonek");
print("   strata); O(1) deviations would mean a PROPER monodromy subgroup.")

def census(p):
    ar = np.arange(p, dtype=np.int64)
    X = np.repeat(ar, p*p)
    Y = np.tile(np.repeat(ar, p), p)
    Z = np.tile(ar, p*p)
    U = (1 + X*Y) % p
    W = (4 + 3*X*Y) % p
    f1 = (U*U % p * U % p * Z + Y*Y % p * U % p * W) % p
    f2 = (Y + 3*X*U % p * U % p * Z + 3*X*Y % p * Y % p * W) % p
    f3 = (2*X - 3*X*X % p * Y - X*X % p * X % p * Z) % p
    key = (f1*p + f2)*p + f3
    n1 = np.bincount(key, minlength=p**3)
    coordkey = (X*p + Y)*p + Z
    w = n1[coordkey].astype(np.float64)
    n2 = np.bincount(key, weights=w, minlength=p**3).astype(np.int64)
    return n1, n2

def kf(t, p):
    """coordinate key of a rational target reduced mod p"""
    vs = []
    for q in t:
        q = sp.Rational(q)
        vs.append(int(q.p) * pow(int(q.q), -1, p) % p)
    return (vs[0]*p + vs[1])*p + vs[2]

def wreath_model():
    pX = {3: sp.Rational(1, 6), 1: sp.Rational(1, 2), 0: sp.Rational(1, 3)}
    joint = {(0, 0): sp.Rational(1, 3)}
    for v, pv in pX.items():
        joint[(1, v)] = sp.Rational(1, 2)*pv
    for v1, p1 in pX.items():
        for v2, p2 in pX.items():
            for v3, p3 in pX.items():
                k = (3, v1+v2+v3)
                joint[k] = joint.get(k, 0) + sp.Rational(1, 6)*p1*p2*p3
    return joint

model = wreath_model()
Pplus = (1, sp.Rational(-3, 2), sp.Rational(13, 2))
Pminus = (-1, sp.Rational(3, 2), sp.Rational(13, 2))
P0 = (0, 0, sp.Rational(-1, 4))
targets = {'v*': vstar, 'P0': P0, 'P+': Pplus, 'P-': Pminus}

primes = [61, 67, 71, 73, 79, 83, 89, 97, 101, 103]
agg, tot = {}, 0
print(f"\n  {'p':>4} {'n1=2/>3 (Jelonek)':>18}  n1@v* n2@v*  n1@P0 n1@P+ n1@P-")
for p in primes:
    n1, n2 = census(p)
    tot += p**3
    pair = n1.astype(np.int64)*1000 + n2
    vals, cnts = np.unique(pair, return_counts=True)
    for v, ct in zip(vals, cnts):
        k = (int(v // 1000), int(v % 1000))
        agg[k] = agg.get(k, 0) + int(ct)
    jel = int(np.sum(n1 == 2)) + int(np.sum(n1 > 3))
    t = {nm: (int(n1[kf(tv, p)]), int(n2[kf(tv, p)])) for nm, tv in targets.items()}
    print(f"  {p:>4} {jel:>10} (/p^2={jel/p**2:4.2f})  "
          f"{t['v*'][0]:>5} {t['v*'][1]:>5} {t['P0'][0]:>6} {t['P+'][0]:>5} {t['P-'][0]:>5}")

print("\n  aggregate joint (n1, n2) over", len(primes), "primes,", tot, "targets:")
print(f"  {'(n1,n2)':>10} {'observed':>12} {'wreath-model':>14} {'obs/model':>10}")
for k in sorted(agg):
    obs = agg[k] / tot
    mod = float(model.get(k, 0))
    ratio = f"{obs/mod:10.4f}" if mod > 0 else "   OFF-MODEL"
    print(f"  {str(k):>10} {obs:12.6f} {mod:14.6f} {ratio}")
missing = [k for k in model if k not in agg]
print("  model cells never observed:", missing if missing else "none")
print()
print("  VERDICTS (see .out):")
print("   * n1 marginal ~ (1/3,1/2,1/6) => monodromy of F = FULL S3 (C3 refuted).")
print("   * (3, n2) row ~ independent-cube model => monodromy of FoF = FULL wreath")
print("     S3 wr S3 (order 1296): imprimitive but with maximal kernel.  A NEW atom")
print("     would need PRIMITIVE degree-9 monodromy; the tower is maximally NON-atom.")
print("   * off-model cells (n1=2 etc.) = Jelonek strata, measure ~ 1/p (fiber drop")
print("     by escape, the etale-Keller analogue of a branch locus).")
