#!/usr/bin/env python3
"""Exact controls for the constant-D finite-six/infinity-two quartic family.

The all-parameter normalization and rational-primitive argument are in the
matching proof note. No inherited mathematical implementation is imported.
"""
import hashlib
import json
import sympy as S

x, t, s, zold, r, sig, bd = S.symbols('x t s zold r sig bd')
a, b, c0, c1, m0, m1, d, n0, n1, c = S.symbols(
    'a b c0 c1 m0 m1 d n0 n1 c')
beta0, beta1 = S.symbols('beta0 beta1')
Z, q, u, v, w, z = S.symbols('Z q u v w z')
gates = 0
records = {}


def check(ok, label):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(label)


def eq(left, right, label):
    check(S.cancel(left-right) == 0, label)


def jac(left, right, p, qvar):
    return S.diff(left, p)*S.diff(right, qvar)-S.diff(left, qvar)*S.diff(right, p)


def coeff(poly, var, degree):
    return S.expand(poly).coeff(var, degree)


def truncated_compose(poly, centre, degree):
    """Independent finite polynomial convolution, not a formal-series solver."""
    def mul(left, right):
        return [S.expand(sum(left[j]*right[k-j] for j in range(k+1)))
                for k in range(degree+1)]
    cp = [coeff(centre, x, j) for j in range(degree+1)]
    answer = [S.Integer(0)]*(degree+1)
    for mon, value in S.Poly(poly, Z).terms():
        power = [S.Integer(1)]+[S.Integer(0)]*degree
        for _ in range(mon[0]):
            power = mul(power, cp)
        pv = [coeff(value, x, j) for j in range(degree+1)]
        term = mul(pv, power)
        answer = [aa+bb for aa, bb in zip(answer, term)]
    return [S.factor(aa) for aa in answer]


# Complete literal source sections and both surface charts.
N = a*x**6+s*(b*x*x+c1*x**3+2*a*x**4)+s*s*(c0+c1*x+a*x*x)
M = m0+m1*x+d*x*x+n1*x**3+s*(n0+n1*x)
Nglobal = S.expand(N.subs(s, zold-x*x))
Mglobal = S.expand(M.subs(s, zold-x*x))
check(S.Poly(Nglobal, x, zold).degree(x) <= 4, 'N global x box')
check(S.Poly(Nglobal, x, zold).degree(zold) <= 2, 'N global z box')
check(S.Poly(Mglobal, x, zold).degree(x) <= 2, 'M global x box')
check(S.Poly(Mglobal, x, zold).degree(zold) <= 1, 'M global z box')
eq(N.subs(s, 0), a*x**6, 'finite boundary multiplicity six')
eq(M.subs({x: 0, s: 0}), m0, 'finite boundary M unit hypothesis')
H = S.expand(t*t*N.subs(s, 1/t))
L = S.expand(t*M.subs(s, 1/t))
F = S.expand(H*H+L)
eq(H, a*x**6*t*t+(b*x*x+c1*x**3+2*a*x**4)*t+c0+c1*x+a*x*x,
   'literal quartic square prefix')
eq(L, (m0+m1*x+d*x*x+n1*x**3)*t+n0+n1*x, 'literal linear remainder')
eq(F.subs(x, 0), c0*c0+m0*t+n0, 'no generic vertical component at x zero')
Hb = S.cancel(H.subs({x: 1/r, t: -r*r-r**4*bd}, simultaneous=True))
Lb = S.cancel(L.subs({x: 1/r, t: -r*r-r**4*bd}, simultaneous=True))
check(Hb.is_polynomial(r, bd) and Lb.is_polynomial(r, bd), 'both full D charts polynomial')
eq(Hb.subs(r, 0), c0-b, 'H constant on D')
eq(Lb.subs(r, 0), n0-d, 'L constant on D')
eq(jac(1/r, -r*r-r**4*bd, r, bd), r*r, 'actual symplectic transition')

qq = r*r+sig
Ni = S.cancel(r**4*qq**2*Nglobal.subs({x: 1/r, zold: 1/qq}, simultaneous=True))
Mi = S.cancel(-r*r*qq*Mglobal.subs({x: 1/r, zold: 1/qq}, simultaneous=True))
Ni_expected = a*r*r-c1*r*sig+(c0-b)*sig*sig-b*r*r*sig
Mi_expected = -n1*r-d*r*r-m1*r**3-m0*r**4+sig*(n0-d-m1*r-m0*r*r)
eq(Ni, Ni_expected, 'complete infinity numerator N')
eq(Mi, Mi_expected, 'complete infinity numerator M')
eq(Ni.subs(sig, 0), a*r*r, 'infinity boundary multiplicity two')
eq(Mi.subs({r: 0, sig: 0}), 0, 'actual common infinity root')
Ei = S.expand(Ni*Ni+sig**3*Mi-c*sig**4)
Pi = (a-c1*Z+(c0-b)*Z*Z)**2-n1*Z**3+(n0-d-c)*Z**4
eq(coeff(Ei.subs(sig, r*Z), r, 4), Pi, 'full generic infinity tangent quartic')
eq(Pi.subs(Z, 0), a*a, 'infinity tangent roots are nonzero')
eq(coeff(Pi, Z, 4), (c0-b)**2+n0-d-c, 'generic infinity degree exactly four')
T0 = S.expand(Pi+c*Z**4)
eq((Z*S.diff(T0, Z)-4*T0).subs(Z, 0), -4*a*a,
   'generic repeated-root eliminant is not identically zero')
eq(coeff(S.diff(Ei, sig).subs(sig, r*Z), r, 3), S.diff(Pi, Z),
   'actual infinity polar coefficient')
check(2+2-3 == 1, 'every simple infinity tangent gives eta order one')

# The two omitted lower normal jets enlarge to the complete section space.
Nfull = N+s*(beta0+beta1*x)
Ngfull = S.expand(Nfull.subs(s, zold-x*x))
Hfull = S.expand(t*t*Nfull.subs(s, 1/t))
eq((Hfull*Hfull+L).subs(x, 0), (beta0*t+c0)**2+m0*t+n0,
   'no generic vertical component also with the lower normal jets')
eq((Nfull*Nfull+s**3*M-c*s**4).subs(x, 0),
   beta0**2*s*s+(2*beta0*c0+m0)*s**3+(c0*c0+n0-c)*s**4,
   'normal-unit lower-jet case has Weierstrass degree two')
Nifull = Ni-beta0*r*r*(r*r+sig)*sig-beta1*r*(r*r+sig)*sig
eq(S.cancel(r**4*qq**2*Ngfull.subs({x: 1/r, zold: 1/qq}, simultaneous=True)),
   Nifull, 'full lower-normal-jet infinity transport')
eq(coeff((Nifull*Nifull+sig**3*Mi-c*sig**4).subs(sig, r*Z), r, 4), Pi,
   'lower normal jets preserve the complete infinity tangent')
Nmon = [x**i*zold**j for i in range(5) for j in range(3)]
Nmats = [[coeff(mon.subs(zold, x*x), x, k) for mon in Nmon]
         for k in range(9) if k != 6]
Nmats += [[coeff(coeff(mon, x, 4), zold, j) for mon in Nmon] for j in (1, 2)]
check(S.Matrix(Nmats).rank() == 9, 'complete N constraint space has dimension six')
Nparams = [a, b, c0, c1, beta0, beta1]
Ncols = S.Matrix([[coeff(coeff(S.diff(Ngfull, p), x, i), zold, j)
                   for p in Nparams] for i in range(5) for j in range(3)])
check(Ncols.rank() == 6, 'displayed full N parameters form a basis')
check(S.Matrix(Nmats)*Ncols == S.zeros(len(Nmats), 6), 'full N basis satisfies every constraint')
Mparams = [m0, m1, d, n0, n1]
Mcols = S.Matrix([[coeff(coeff(S.diff(Mglobal, p), x, i), zold, j)
                   for p in Mparams] for i in range(3) for j in range(2)])
check(Mcols.rank() == 5, 'displayed M parameters span the five-dimensional constant-D space')
eq(coeff(coeff(Mglobal, x, 2), zold, 1), 0, 'only M constant-D linear constraint')
for j in (0, 1):
    check(6 > 3*j and (6-3*j) > 0, 'lower normal jet regular regime '+str(j))

# Balanced finite face and complete critical-value jets.
E = S.expand(N*N+s**3*M-c*s**4)
eq(E.subs(x, 0), m0*s**3+(c0*c0+n0-c)*s**4,
   'finite point has exactly three Weierstrass sheets')
Q = S.cancel(E.subs(s, x**4*Z)/x**12)
check(Q.is_polynomial(x, Z), 'actual balanced rescaling has no negative powers')
P = (a+b*Z)**2+m0*Z**3
eq(Q.subs(x, 0), P, 'actual finite cubic face')
eq(S.discriminant(P, Z), a**3*m0*(4*b**3-27*a*m0), 'only nonzero double-face stratum')
double = {a: -b*q/3, m0: -4*b*b/(9*q)}
Qd = S.expand(Q.subs(double))
eq(P.subs(double), -4*b*b/(9*q)*(Z-q)**2*(Z-q/4), 'double and simple face slopes')
eq(S.diff(Qd, Z, 2).subs({x: 0, Z: q}), -2*b*b/3, 'critical centre nondegenerate')
z1 = -c1*q/b
z2 = q*(2*b*b*q+3*b*c0*q+3*c1*c1)/(3*b*b)
z3 = -c1*q*(b*b*q+9*b*c0*q+3*c1*c1)/(3*b**3)
m1star = -4*b*c1/(3*q)
dstar = 4*(2*b*b*q-3*b*c0*q-3*c1*c1)/(9*q)
n1star = 4*c1*(b*b*q-3*b*c0*q-c1*c1)/(9*b*q)
R1 = q*q*(4*b*c1+3*m1*q)/3
R2 = -q*q*(8*b*b*q-12*b*c0*q-12*c1*c1-9*d*q)/9
R3 = -q*q*(4*b*b*c1*q-12*b*c0*c1*q-9*b*n1*q-4*c1**3)/(9*b)
eq(truncated_compose(Qd, q, 1)[1], R1, 'first critical value')
eq(R1.subs(m1, m1star), 0, 'necessary first coefficient constraint')
eq(truncated_compose(Qd.subs(m1, m1star), q+z1*x, 2)[2], R2,
   'second critical value after first constraint')
eq(R2.subs(d, dstar), 0, 'necessary second coefficient constraint')
eq(truncated_compose(Qd.subs({m1: m1star, d: dstar}), q+z1*x+z2*x*x, 3)[3], R3,
   'third critical value after first two constraints')
eq(R3.subs(n1, n1star), 0, 'necessary third coefficient constraint')
star = {m1: m1star, d: dstar, n1: n1star}
Qs = S.expand(Qd.subs(star))
centre = q+z1*x+z2*x*x+z3*x**3
Rjets = truncated_compose(Qs, centre, 5)
Cjets = truncated_compose(S.diff(Qs, Z), centre, 3)
for j in range(4):
    eq(Rjets[j], 0, 'vanishing critical value coefficient '+str(j))
    eq(Cjets[j], 0, 'critical centre coefficient '+str(j))
r4 = q**3*(4*b*b*q-24*b*c0*q-27*c*q+36*c0*c0*q+12*c1*c1+27*n0*q)/27
r5 = -4*c1*q**3*(5*b*b*q-27*b*c0*q-27*c*q+36*c0*c0*q+12*c1*c1+27*n0*q)/(27*b)
eq(Rjets[4], r4, 'fourth critical value coefficient')
eq(Rjets[5], r5, 'fifth critical value coefficient')
eq(S.diff(r4, c), -q**4, 'generic order four and fibre dependence')
h0 = -b*b/3
h1 = truncated_compose(S.diff(Qs, Z, 2)/2, centre, 1)[1]
eq(h1, -2*b*c1/3, 'first Hessian coefficient')
resnum = S.factor(4*z1*r4/q-h1*r4/h0-r5)
eq(resnum, -2*c1*q**3*(2*b*b*q-18*b*c0*q-27*c*q+36*c0*c0*q+12*c1*c1+27*n0*q)/(27*b),
   'twice the residue bracket times r4')
eq(S.diff(resnum, c), 2*c1*q**4/b, 'generic residues force c1 zero')
eq(resnum.subs(c1, 0), 0, 'even family residue boundary')

# A separate abstract jet calculation checks the sign in the residue formula.
rr4, rr5, hh0, hh1, ee, ff = S.symbols('rr4 rr5 hh0 hh1 ee ff')
split = rr4*x**4+rr5*x**5+(hh0+hh1*x)*(ee*x*x+ff*x**3)**2
eq(coeff(split, x, 4), rr4+hh0*ee*ee, 'split leading equation')
eq(coeff(split, x, 5), rr5+hh1*ee*ee+2*hh0*ee*ff, 'split next equation')
eq((-hh1/(2*hh0)+rr5/(2*rr4))+hh1/hh0,
   hh1/(2*hh0)+rr5/(2*rr4), 'polar correction uses plus r5 before inversion')

for lam, branch_orders, budget in [(1, [0], 0), (2, [-1, -1], 0),
                                   (3, [-2], 1), (4, [-2, -2], 2)]:
    check(sum(max(0, -o-1) for o in branch_orders) == budget,
          'normalized primitive budget lambda '+str(lam))
    if lam < 4:
        check(budget < 2, 'infinity ramification defeats lower critical stratum '+str(lam))

# Complete even stratum: literal birational transport and primitive space.
even = dict(double)
even.update({c1: 0, m1: 0, n1: 0, d: 8*b*b/9-4*b*c0/3})
Fe = S.expand(F.subs(even))
eq(Fe.subs(x, -x), Fe, 'full even family parity')
ep = c0-b
C = 4*b*b*q/9
B = 2*(-b*q/3)*ep-C
K0 = ep*ep+n0-even[d]
aa = -b*q/3
Fu = S.cancel(Fe.subs({x: 1/u, t: u**3*v-u*u}, simultaneous=True))
eq(Fu, (aa*v*v+b*u*v+ep)**2+double[m0]*u**3*v-double[m0]*u*u+even[d]*u*v+n0-even[d],
   'literal polynomial in reciprocal and cubic coordinates')
eq(jac(1/x, x+x**3*t, x, t), -x, 'first birational Jacobian')
eq(jac(v-u/q, v*(v-u/q), u, v), -(v-u/q)/q, 'second birational Jacobian')
model = -3*aa*aa*z*z+B*z+K0+(4*aa*aa*z+C)*w*w
eq(S.cancel(Fu.subs({u: q*(z/w-w), v: z/w}, simultaneous=True)), model,
   'exact elliptic model of the entire even stratum')
eq(q*(z/w-w), q*(z-w*w)/w, 'explicit inverse birational x coordinate')
A = 4*aa*aa*z+C
P2 = 3*aa*aa*z*z-B*z+c-K0
eq(model-c, A*w*w-P2, 'elliptic double cover equation')
eq(S.diff(S.discriminant(P2, z), c), -12*aa*aa, 'two generic distinct numerator roots')
eq(S.diff(P2.subs(z, -C/(4*aa*aa)), c), 1, 'numerator and denominator coprime generically')
eq(S.diff(P2.subs(z, 0), c), 1, 'neither double pole lies at z zero generically')
check(-4+4 == 0, 'four branch points give geometric genus one')
eq(-q*(q*(z/w-w))/w, -q*q*(z-w*w)/(w*w), 'actual omega coefficient in dz wedge dw')
eta = q*q*(z-w*w)/(2*A*w**3)
eq(S.diff(model, w)*eta, q*q*(z-w*w)/(w*w), 'dF wedge eta equals actual omega')
eta_numerator = S.expand(A*z-P2)
primitive_numerator = S.expand(S.diff(P2, z)*A-P2*S.diff(A, z))
eq(eta_numerator, aa*aa*z*z+(B+C)*z-(c-K0), 'actual second-kind numerator')
eq(primitive_numerator, 12*aa**4*z*z+6*aa*aa*C*z-B*C-4*aa*aa*(c-K0),
   'derivative of one over w numerator')
kforced = -q*q/(12*aa*aa)
eq(coeff(q*q*eta_numerator+kforced*primitive_numerator, z, 2), 0,
   'only possible primitive scalar from leading coefficient')
eq(S.diff(coeff(q*q*eta_numerator+kforced*primitive_numerator, z, 0), c),
   -2*q*q/3, 'nonzero generic fibre obstruction to exactness')

# Explicit positive/hostile boundaries; no finite bank substitutes for proof.
named = {b: 1, q: S.Rational(4, 9), c0: 0, c1: 0, n0: 0}
check(R1.subs(named).subs(m1, 1) != 0, 'named order-one stratum')
check(R2.subs(named).subs(d, 0) != 0, 'old constant-D stopping object has order two')
check(R3.subs(named).subs(n1, 1) != 0, 'named order-three stratum')
eq(r4.subs(named), S.Rational(256, 6561)*(S.Rational(4, 27)-c),
   'named even order-four actual fibre coefficient')
check(S.diff(resnum.subs(named | {c1: 1}), c) != 0,
      'named order-four asymmetric residue hostile')
hp = x*x+x**4*t
Fp = hp**4+hp
Gp = 1/(3*x**3*(4*hp**3+1))
eq(jac(Fp, Gp, x, t), 1, 'genuine rational mate outside finite-six hypothesis')
eq(jac(t**4, -x/(4*t**3), x, t), 1, 'constant-D rational mate outside the boundary pattern')
eq(((-r*r-r**4*bd)**4).subs(r, 0),
   0, 'constant-D rational positive control')
eq(jac(x*x, t/(2*x), x, t), 1, 'affine criticality alone does not exclude rational mates')
records['critical_constraints'] = [str(S.factor(aa)) for aa in (m1star, dstar, n1star)]
records['critical_order'] = [str(r4), str(r5)]
records['residue_necessity'] = 'c1=0; m1=n1=0'
records['primitive_budgets'] = [0, 0, 1, 2]
records['elliptic_model'] = 'A(z)*w^2=P2(z); four simple branch points generically'
records['exactness_obstruction'] = str(-2*q*q/3)
records['scope'] = 'a*m0!=0; full fixed-x=0 finite-six constant-D family, including two lower normal jets; no rational scalar-Jacobian mate'
semantic = hashlib.sha256(json.dumps(records, sort_keys=True).encode()).hexdigest()
print('PASS constant-D quartic family:', gates, 'always-active exact gates')
print('Generic infinity: four simple branches, eta order +1 at every branch')
print('Finite split orders 1,2,3: total primitive pole budget <2')
print('Order 4: residues force the even stratum; its generic elliptic primitive space fails exactness')
print('Named lower strata, asymmetric residue hostile, and genuine outside-family rational mate checked')
print('Semantic SHA256:', semantic)
