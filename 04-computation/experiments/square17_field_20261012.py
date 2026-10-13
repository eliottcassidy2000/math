"""The field of the seventeen-square packing constant (sq17-lean, 2026) against the eleven-square octic.

The Lean formalization github.com/ofitzharding/sq17-lean (read 2026-10-12; its optimality theorem is
CITED, not re-verified here) defines the optimal side as  side = (6 + 4 t) / (1 + 2 t - t^2)  where
(t, beta) is the unique root in a 1e-12 box around
  t0    = 0.362043893314323383453621978739944164218293243036741520638379,
  beta0 = 0.330948674722739910429808533282915634174672003999685385633441
of the two integer polynomials below (Sq17Lean/EndpointRootData.lean; transcribed from the source and
validated by the root).  This script derives the minimal polynomials of t and of the side, proves the
Galois group is S_18 (Dedekind cycle types plus Jordan), factors the discriminant, and compares with the
eleven-square octic of square11_octic_langlands_collatz_20261006.md.
Run from the repository root: python3 -B <script>   (about two minutes)
"""
import collections, json
import sympy as sp
import mpmath as mp

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


x, y, s = sp.symbols('x y s')
P0 = x**5*y - x**4*y**2 - x**4*y - x**3*y**2 - x**4 + x**2*y**2 - x**3 - 3*x*y**2 - 3*x**2 - x*y + x + y
P1 = (x**5*y**3 + x**4*y**3 + 2*x**2*y**4 - 2*x**4*y**2 + x**5*y + 2*x*y**4 + 4*x**2*y**3 + 2*x**3*y**2 + x**4*y
      - 5*x*y**3 - 2*x**2*y**2 - 2*x**4 - y**3 - 2*x*y**2 - 4*x**2*y - 2*x**3 + 4*y**2 + 3*x*y - y)
T0 = sp.Rational(362043893314323383453621978739944164218293243036741520638379, 10**60)
B0 = sp.Rational(330948674722739910429808533282915634174672003999685385633441, 10**60)
check(abs(sp.N(P0.subs({x: T0, y: B0}), 70)) < 1e-58 and abs(sp.N(P1.subs({x: T0, y: B0}), 70)) < 1e-58,
      'the Lean midpoint is a root of both polynomials to 1e-58 (validates the transcription)')
mp.mp.dps = 80
F0 = sp.lambdify((x, y), P0, 'mpmath'); F1 = sp.lambdify((x, y), P1, 'mpmath')
root = mp.findroot(lambda u, v: [F0(u, v), F1(u, v)], (mp.mpf(str(T0)), mp.mpf(str(B0))))
t, beta = root[0], root[1]
side = (6+4*t)/(1+2*t-t*t)
check(abs(side-mp.mpf('4.675530093604551')) < 1e-15, 'side 4.675530093604551 (README)')

# minimal polynomial of t: the degree-18 factor of the resultant
R = sp.Poly(sp.expand(sp.resultant(P0, P1, y)), x)
check(R.degree() == 28, 'resultant degree 28')
factors = sp.factor_list(R.as_expr())[1]
check(sorted((sp.Poly(f, x).degree(), m) for f, m in factors) == [(1, 2), (1, 2), (1, 2), (2, 2), (18, 1)], 'factor pattern 1,1,1,2 (squared) and 18')
mt = sp.Poly(min((f for f, _ in factors), key=lambda f: abs(sp.N(f.subs(x, sp.Float(str(t), 60)), 60))), x)
check(mt.degree() == 18 and mt.is_irreducible and abs(sp.N(mt.as_expr().subs(x, sp.Float(str(t), 60)), 60)) < 1e-50, 'degree-18 minimal polynomial of t')
mt_coeffs = [int(c) for c in mt.all_coeffs()]
# minimal polynomial of the side
RT = sp.Poly(sp.expand(sp.resultant(mt.as_expr(), (1+2*x-x**2)*s-(6+4*x), x)), s)
fT = sp.factor_list(RT.as_expr())[1]
check(len(fT) == 1 and sp.Poly(fT[0][0], s).degree() == 18, 'side polynomial irreducible of degree 18')
mT = sp.Poly(fT[0][0], s)
check(abs(sp.N(mT.as_expr().subs(s, sp.Float(str(side), 60)), 60)) < 1e-30, 'side is a root')
mT_coeffs = [int(c) for c in mT.all_coeffs()]
check(mT_coeffs[0] == 4775 and mT_coeffs[-1] == -2631254953, 'leading 4775 = 5^2 191, constant -2631254953')


# Galois group by Dedekind cycle types (distinct-degree factorization mod p) and Jordan's theorem
def polymod(a_, m_, p):
    a_ = [c % p for c in a_]; dm = len(m_)-1
    for i in range(len(a_)-1, dm-1, -1):
        cc = a_[i]
        if cc:
            for j in range(dm+1):
                a_[i-dm+j] = (a_[i-dm+j]-cc*m_[j]) % p
    a_ = a_[:dm]
    while a_ and a_[-1] == 0:
        a_.pop()
    return a_


def polymul(a_, b_, p):
    r = [0]*(len(a_)+len(b_)-1)
    for i, u in enumerate(a_):
        if u:
            for j, w in enumerate(b_):
                r[i+j] = (r[i+j]+u*w) % p
    return r


def powmod(base, e, m_, p):
    r = [1]; bb = polymod(base, m_, p)
    while e:
        if e & 1:
            r = polymod(polymul(r, bb, p), m_, p)
        bb = polymod(polymul(bb, bb, p), m_, p); e >>= 1
    return r


def polygcd(a_, b_, p):
    a_ = [c % p for c in a_]; b_ = [c % p for c in b_]
    while True:
        while b_ and b_[-1] == 0:
            b_.pop()
        if not b_:
            break
        inv = pow(b_[-1], -1, p); bm = [(c*inv) % p for c in b_]
        a_ = polymod(a_, bm, p) if len(a_) >= len(bm) else a_
        a_, b_ = b_, a_
    while a_ and a_[-1] == 0:
        a_.pop()
    return a_


def cycle_type(coeffs_high_to_low, p):
    f = [c % p for c in coeffs_high_to_low[::-1]]
    if f[-1] == 0:
        return None
    inv = pow(f[-1], -1, p); f = [(c*inv) % p for c in f]
    df = [(i*f[i]) % p for i in range(1, len(f))]
    if len(polygcd(f, df, p)) > 1:
        return None
    degs, g, h, k = [], f[:], [0, 1], 0
    while len(g) > 1 and k < len(f):
        k += 1; h = powmod(h, p, f, p)
        hk = h[:]+[0]*(2-len(h)) if len(h) < 2 else h[:]
        hk[1] = (hk[1]-1) % p
        d = polygcd(g, hk, p)
        if len(d) > 1:
            degs += [k]*((len(d)-1)//k)
            dinv = pow(d[-1], -1, p); dm = [(c*dinv) % p for c in d]
            q = [0]*(len(g)-len(dm)+1); r = g[:]
            for i in range(len(g)-1, len(dm)-2, -1):
                cc = r[i]
                if cc:
                    q[i-len(dm)+1] = cc
                    for j in range(len(dm)):
                        r[i-len(dm)+1+j] = (r[i-len(dm)+1+j]-cc*dm[j]) % p
            g = q
            while g and g[-1] == 0:
                g.pop()
    return tuple(sorted(degs))


def roots_mod(coeffs_high_to_low, p):
    return sum(1 for v in range(p) if not sum(c % p*pow(v, k, p) for k, c in enumerate(coeffs_high_to_low[::-1])) % p)


types = collections.Counter(); first = {}
p = 2
while p < 400:
    ct = cycle_type(mt_coeffs, p)
    if ct is not None:
        types[ct] += 1; first.setdefault(ct, p)
        check(ct.count(1) == roots_mod(mt_coeffs, p), 'root count agrees with the distinct-degree factorization')
    p = sp.nextprime(p)
check(types[(5, 13)] >= 1 and first[(5, 13)] == 3, 'a 13-cycle at p = 3 (prime cycle of length > 9, Jordan: A_18)')
check(types[(18,)] >= 1 and first[(18,)] == 7, 'an 18-cycle at p = 7 (odd permutation): S_18')
disc_t = sp.factorint(mt.discriminant())
check(disc_t == {-1: 1, 2: 181, 5: 1, 154943: 2, 388183: 1, 193319963: 1}, 'discriminant of the t-polynomial')
nroots = mt.nroots(n=40)
real_roots = sorted(float(r) for r in nroots if abs(sp.im(r)) < 1e-30)
check(len(real_roots) == 4, 'signature (4, 7)')
# the eleven-square octic (square11_octic_langlands_collatz_20261006.md): u = tan(theta/2), T = (6u+4)/(1+2u-u^2)
P11 = [1, -20, 178, -842, 1923, -496, -6754, 12420, -6865]
disc11 = sp.factorint(sp.Poly(P11, x).discriminant())
shared = sorted(set(disc_t)-{-1} & set(disc11)-{-1})
u = mp.mpf('0.365769307604677293388545')
angles = dict(trump_deg=float(2*mp.atan(u)*180/mp.pi), bidwell_deg=float(2*mp.atan(t)*180/mp.pi))
check(abs(angles['trump_deg']-40.18193729) < 1e-6 and abs(angles['bidwell_deg']-39.80495898) < 1e-6, 'tilt angles')
# the mirrored side formulas: s11 = (2 + 2 cos + 3 sin)/(cos + sin), s17 = (3 + 3 cos + 2 sin)/(cos + sin)
th = 2*mp.atan(t)
check(abs((3+3*mp.cos(th)+2*mp.sin(th))/(mp.cos(th)+mp.sin(th))-side) < 1e-60, 's17 = (3 + 3 cos + 2 sin)/(cos + sin)')
thu = 2*mp.atan(u)
check(abs((2+2*mp.cos(thu)+3*mp.sin(thu))/(mp.cos(thu)+mp.sin(thu))-(6*u+4)/(1+2*u-u*u)) < 1e-20, 's11 = (2 + 2 cos + 3 sin)/(cos + sin)')
out = dict(t=mp.nstr(t, 50), beta=mp.nstr(beta, 50), side=mp.nstr(side, 50), min_poly_t=mt_coeffs, min_poly_side=mT_coeffs,
           cycle_types={str(k): v for k, v in sorted(types.items(), key=lambda kv: -kv[1])[:12]}, first_primes={str(k): v for k, v in first.items() if k in ((5, 13), (18,))},
           galois_group='S18 (exact: irreducible, 13-cycle at 3, 18-cycle at 7)', discriminant_t={str(k): v for k, v in disc_t.items()},
           real_roots=real_roots, shared_discriminant_primes_with_octic=shared, angles=angles, checks=CHECKS)
print(json.dumps(out, indent=1))
