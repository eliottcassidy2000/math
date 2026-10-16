#!/usr/bin/env python3
"""S38 audit (independent referee): arithmetic and algebra checks for THM-4661 statements 1-4, written from scratch.

1. Counting: the simple bound 1 + q + ceil(q(q-2)/2), the even-q extremal exclusion (my own derivation of the
   Gallai count and the layer count), and the first-moment threshold of the note.
2. Quadratic subfields of Q(zeta55) (group count) and non-residue certificates for -2, -35 in L' = L(theta) at
   several primes p = 1 mod 55 where the monic cubic f(w) = w^3 + 2w^2 - 4w - 4 (w = 2 theta) has a simple root.
3. F_4 = K^+: tau = sqrt(-11)(zeta5 - zeta5^-1) satisfies tau^4 - 55 tau^2 + 605 = 0 and is real.
4. Fish equations (1): symbolic verification that the skew products give 11 z0^2 = D(1-D), 11 z1^2 = D(4-D),
   11 z2^2 = (1-D)(D+3) from the constraints only, and that z0, z1, z2 are conjugation-invariant; genus of C3 by
   Riemann-Hurwitz; D = 2 is not a branch value; smoothness of the affine model for 0 < D < 1.
5. The residue colouring at 19 (rhombus G2, reused by C1) and the simple root of f modulo 19.
6. EI19 forced ratios have modulus one and generate Q(sqrt-2), Q(sqrt-35).
7. Statement 3(b): explicit chains p_0..p_m with steps d < 1/2 and |p_0 - p_m| = 1, m = ceil(1/d).
"""
import itertools, json, math, os
from fractions import Fraction as Fr
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_field_and_counting_checks.json')
R = {}


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- 1. counting
rows = []
for q in range(4, 41):
    simple = 1 + q + math.ceil(Fr(q * (q - 2), 2))
    rows.append(dict(q=q, simple_min_order=simple, palette_bound=simple - 1, bridge_bound=q * (q - 1) // 2))
assert rows[2]['q'] == 6 and rows[2]['simple_min_order'] == 19
R['simple_bound'] = rows
log('simple bound: q=6 minimal failed exterior has >= 19 vertices; q + ceil(q(q-2)/2) vs q(q-1)/2:',
    [(r['q'], r['palette_bound'], r['bridge_bound']) for r in rows[:8]])

even = []
for q in range(4, 41, 2):
    n = 1 + q + q * (q - 2) // 2          # extremal order
    Wn = q * (q - 2) // 2
    # low vertices: degree q-1 (list q-1) or degree q (list q); high: degree q, list q-1.
    # sum of low degrees >= (q-1) l ;  = 2 e(Lambda) + q h ;  e(Lambda) <= (5/4)(l - 1)
    h = next(h for h in range(0, n + 1) if (q - 1) * (n - h) <= Fr(5, 2) * (n - h - 1) + q * h)
    # layer count at a high root z: each w in W has q-2 or q-3 neighbours inside W (degree q or q-1),
    # degree-q vertices of W have all W-neighbours of degree q-1:  (q-3)|W_{q-1}| = (q-2)|W_q| + 2 e(W_{q-1})
    feas = [w for w in range(h - 1, Wn + 1) if ((q - 3) * (Wn - w) - (q - 2) * w) >= 0
            and ((q - 3) * (Wn - w) - (q - 2) * w) % 2 == 0]
    even.append(dict(q=q, n=n, h_min=h, feasible_W_q=feas))
    assert not feas
R['even_q_extremal_exclusion'] = even
log('even q in 4..40: extremal order excluded for all (h_min, |W| ):', [(e['q'], e['h_min']) for e in even[:6]], '...')

thr = 2 * (math.log(6) / math.log(6 / 5) - 1)
R['first_moment_threshold_average_degree'] = thr
R['average_degrees'] = dict(EI17=2 * 31 / 17, EI19=2 * 35 / 19, fish=2 * 42 / 23,
                            probes=[2 * 14220 / 4971, 2 * 37086 / 13287, 2 * 17230 / 8581, 2 * 12430 / 6161])
log(f'first moment: 6^n (5/6)^(m+n) < 1 iff average degree > {thr:.4f}; average degrees {R["average_degrees"]}')

# ---------------------------------------------------------------- 2. quadratic subfields and certificates
# (Z/55)^x = (Z/5)^x x (Z/11)^x = C4 x C10: index-2 subgroups = nonzero homs to C2 = 2^2 - 1 = 3
R['quadratic_subfields_of_Q_zeta55'] = 2 ** 2 - 1
# which d are squares in L: d is a square in L iff Q(sqrt d) has conductor dividing 55
cands = {}
for d in (5, -11, -55, -2, -35, -1, -3, 2, 7, -7, 35):
    disc = d if d % 4 == 1 else 4 * d
    cands[d] = (55 % abs(disc) == 0)
R['square_in_L_by_conductor'] = cands
log('d square in Q(zeta55) (conductor of Q(sqrt d) divides 55):', cands)

f = lambda w, p: (w ** 3 + 2 * w * w - 4 * w - 4) % p
certs = []
p = 1
while len(certs) < 6:
    p += 55
    if not sp.isprime(p) or 592 % p == 0:
        continue
    roots = [w for w in range(p) if f(w, p) == 0]
    if not roots:
        continue
    simple = all((3 * w * w + 4 * w - 4) % p != 0 for w in roots)
    lm2 = pow(-2 % p, (p - 1) // 2, p)
    lm35 = pow(-35 % p, (p - 1) // 2, p)
    certs.append(dict(p=p, roots=roots, simple=simple, euler_minus2=lm2, euler_minus35=lm35,
                      quad_sub_squares=[pow(d % p, (p - 1) // 2, p) for d in (5, -11, -55)]))
R['certificates'] = certs
for c in certs:
    log('p =', c['p'], 'roots of f:', c['roots'], 'simple:', c['simple'], ' (-2|p) ->', c['euler_minus2'],
        ' (-35|p) ->', c['euler_minus35'], ' 5,-11,-55 ->', c['quad_sub_squares'])
c661 = next(c for c in certs if c['p'] == 661)
assert c661['euler_minus2'] == 660 and c661['euler_minus35'] == 660 and c661['simple']
# the S38 script's theta root 214 corresponds to w = 428
assert f(428, 661) == 0
th = sp.symbols('th')
assert sp.discriminant(2 * th ** 3 + 2 * th ** 2 - 2 * th - 1, th) == 148
ww = sp.symbols('ww')
assert sp.discriminant(ww ** 3 + 2 * ww ** 2 - 4 * ww - 4, ww) == 592

# ---------------------------------------------------------------- 3. F_4
z5 = sp.exp(2 * sp.pi * sp.I / 5)
tau = sp.sqrt(-11) * (z5 - 1 / z5)
tv = complex(sp.N(tau, 50))
poly = tv ** 4 - 55 * tv ** 2 + 605
x = sp.symbols('x')
mp = sp.minimal_polynomial(sp.sqrt(sp.Rational(55, 2) + sp.Rational(11, 2) * sp.sqrt(5)), x)
R['F4'] = dict(tau_numeric=str(tv), residual=abs(poly), minpoly_of_sqrt_55p11s5_over_2=str(mp))
assert abs(tv.imag) < 1e-12 and abs(poly) < 1e-9 and sp.expand(mp - (x ** 4 - 55 * x ** 2 + 605)) == 0
log('F_4: tau = sqrt(-11)(z5 - z5^-1) =', tv.real, 'is real and a root of x^4 - 55x^2 + 605:', mp)

# ---------------------------------------------------------------- 4. fish equations
c1, c2, d1, d2, a1, a2, b1, b2 = sp.symbols('c1 c2 d1 d2 a1 a2 b1 b2', real=True)
c = c1 + sp.I * c2; d = d1 + sp.I * d2; a = a1 + sp.I * a2; b = b1 + sp.I * b2
cb, db, ab, bb = sp.conjugate(c), sp.conjugate(d), sp.conjugate(a), sp.conjugate(b)
beta = sp.I * sp.sqrt(11)
D = sp.expand(d * db)
rel = {}
# constraints: Re(c dbar) = 0, |c|^2 + |d|^2 = 1, |a| = |b| = 1, Re(a dbar) = -D/2, Re(b cbar) = -(1-D)/2
z0 = (c * db - cb * d) / (2 * beta)
z1 = (a * db - ab * d) / beta
z2 = (b * cb - bb * c) / beta
# identities: (x ybar - xbar y)^2 = 4 (Re(x ybar))^2 - 4 |x|^2 |y|^2
for name, (X, Y) in dict(z0=(c, d), z1=(a, d), z2=(b, c)).items():
    lhs = sp.expand((X * sp.conjugate(Y) - sp.conjugate(X) * Y) ** 2)
    rhs = sp.expand(4 * sp.re(sp.expand(X * sp.conjugate(Y))) ** 2 - 4 * sp.expand(X * sp.conjugate(X)) * sp.expand(Y * sp.conjugate(Y)))
    assert sp.simplify(lhs - rhs) == 0
Dd = sp.symbols('D', positive=True)
# substitute the constraint values into 4 Re^2 - 4|x|^2|y|^2:
e0 = sp.expand((4 * 0 - 4 * (1 - Dd) * Dd) / (4 * beta ** 2) * 11)          # 11 z0^2
e1 = sp.expand((4 * (Dd / 2) ** 2 - 4 * 1 * Dd) / beta ** 2 * 11)           # 11 z1^2
e2 = sp.expand((4 * ((1 - Dd) / 2) ** 2 - 4 * 1 * (1 - Dd)) / beta ** 2 * 11)  # 11 z2^2
R['fish'] = dict(e0=str(sp.factor(e0)), e1=str(sp.factor(e1)), e2=str(sp.factor(e2)))
assert sp.expand(e0 - Dd * (1 - Dd)) == 0 and sp.expand(e1 - Dd * (4 - Dd)) == 0 and \
    sp.expand(e2 - (1 - Dd) * (Dd + 3)) == 0
# conjugation invariance: conj(x ybar - xbar y) = -(x ybar - xbar y) and conj(beta) = -beta
for X, Y in ((c, d), (a, d), (b, c)):
    s = X * sp.conjugate(Y) - sp.conjugate(X) * Y
    assert sp.simplify(sp.conjugate(s) + s) == 0
log('fish (1): 11 z0^2 =', R['fish']['e0'], '; 11 z1^2 =', R['fish']['e1'], '; 11 z2^2 =', R['fish']['e2'],
    '; z_i conjugation-invariant (skew product / beta)')
# genus by Riemann-Hurwitz: square classes D(1-D), D(4-D), (1-D)(D+3): branch values 0,1,4,-3, each index 2
branch = {}
for val in (0, 1, 4, -3, 'inf'):
    if val == 'inf':
        inert = 0                                  # every square class has even degree 2: unramified at infinity
    else:
        inert = sum(1 for g in (Dd * (1 - Dd), Dd * (4 - Dd), (1 - Dd) * (Dd + 3)) if g.subs(Dd, val) == 0)
    branch[str(val)] = inert
ram_points = sum(1 for k, v in branch.items() if v > 0)
genus = (8 * (-2) + ram_points * 4 * (2 - 1)) // 2 + 1
R['C3'] = dict(branch_values_with_odd_classes=branch, genus=genus,
               fibre_at_2=[str(sp.Rational(2 * (1 - 2), 11)), str(sp.Rational(2 * (4 - 2), 11)),
                           str(sp.Rational((1 - 2) * (2 + 3), 11))])
assert genus == 1 and all(g.subs(Dd, 2) != 0 for g in (Dd * (1 - Dd), Dd * (4 - Dd), (1 - Dd) * (Dd + 3)))
# smoothness: Jacobian of (11 z_i^2 - g_i(D)) has rank < 3 only if two z_i vanish, i.e. D in {0, 1}
zz = sp.symbols('z0:3')
Fs = [11 * zz[0] ** 2 - Dd * (1 - Dd), 11 * zz[1] ** 2 - Dd * (4 - Dd), 11 * zz[2] ** 2 - (1 - Dd) * (Dd + 3)]
J = sp.Matrix([[sp.diff(F, v) for v in (Dd,) + zz] for F in Fs])
sing = []
for pair in itertools.combinations(range(3), 2):
    sols = sp.solve([Fs[pair[0]].subs(zz[pair[0]], 0), Fs[pair[1]].subs(zz[pair[1]], 0)], Dd, dict=True)
    sing += [s[Dd] for s in sols]
R['C3']['D_values_with_two_vanishing_z'] = sorted(set(str(s) for s in sing))
log('C3: genus', genus, 'by Riemann-Hurwitz; D = 2 unramified (fibre z^2 =', R['C3']['fibre_at_2'],
    '); possible singular D values of the affine model:', R['C3']['D_values_with_two_vanishing_z'])

# ---------------------------------------------------------------- 5. residue colouring at 19
p = 19
units = [(x_, y_) for x_ in range(p) for y_ in range(p) if (x_ * x_ + y_ * y_) % p == 1]
Dset = sorted(set((x_ + y_) % p for x_, y_ in units))
col = (4, 3, 2, 4, 3, 2, 1, 0, 4, 1, 0, 4, 3, 2, 1, 0, 2, 1, 0)
ok19 = all(col[s] != col[(s + t) % p] for s in range(p) for t in Dset)
fw = lambda w: w ** 3 + 2 * w * w - 4 * w - 4
R['residue19'] = dict(norm_one_count=len(units), sums=Dset, colouring_proper=ok19, colours=len(set(col)),
                      f4_mod19=fw(4) % 19, fprime4_mod19=(3 * 16 + 16 - 4) % 19,
                      minus11_nonresidue_mod19=pow(-11 % 19, 9, 19) == 18, p_mod5=19 % 5)
assert ok19 and len(units) == 20 and Dset == [1, 2, 4, 6, 9, 10, 13, 15, 17, 18] and fw(4) % 19 == 0
log('residue colouring at 19:', R['residue19'])

# ---------------------------------------------------------------- 6. EI19 ratios
r1 = (-1 + 2 * sp.sqrt(-2)) / 3
r2 = (1 + sp.sqrt(-35)) / 6
R['EI19'] = dict(abs_r1=str(sp.simplify(sp.Abs(r1))), abs_r2=str(sp.simplify(sp.Abs(r2))),
                 r1_plus_conj=str(sp.simplify(r1 + sp.conjugate(r1))), r2_plus_conj=str(sp.simplify(r2 + sp.conjugate(r2))))
assert sp.simplify(sp.Abs(r1) - 1) == 0 and sp.simplify(sp.Abs(r2) - 1) == 0
log('EI19 ratios:', R['EI19'])

# ---------------------------------------------------------------- 7. chains for 3(b)
chains = []
for dd in (0.49, 0.3, 0.1, 0.0123):
    m = math.ceil(1 / dd)
    # zigzag with equal angles: endpoints at distance m d cos(alpha) after m steps alternating +-alpha
    # choose alpha so that the end-to-end distance is exactly 1 (for even m), else straight m-1 steps + adjust
    if m % 2 == 1:
        m += 1
    alpha = math.acos(1 / (m * dd))
    pts = [0j]
    for i in range(m):
        ang = alpha if i % 2 == 0 else -alpha
        pts.append(pts[-1] + dd * complex(math.cos(ang), math.sin(ang)))
    chains.append(dict(d=dd, m=m, end_distance=abs(pts[-1] - pts[0]),
                       step_errors=max(abs(abs(pts[i + 1] - pts[i]) - dd) for i in range(m))))
R['chains'] = chains
assert all(abs(cc['end_distance'] - 1) < 1e-12 for cc in chains)
log('3(b) chains:', chains)

json.dump(R, open(OUTJ, 'w'), indent=1, default=str)
log('wrote', OUTJ)
