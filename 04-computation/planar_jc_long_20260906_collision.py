#!/usr/bin/env python3
"""Exact controls for the all-order moving-collision period theorem.

No repository modules or output files are imported.  The proof is analytic;
the finite universe is four characteristic-zero coefficient strata, each
solved through s^6, all three target two-form slots, and five formal clocks.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sympy as sp

CHECKS = 0
ORDER = 6


def check(ok, name):
    global CHECKS
    if not ok:
        raise AssertionError(name)
    CHECKS += 1


# Coefficients in Q[alpha]/F4.  Rational controls use its rational subfield.
# General field division is unnecessary: every formal denominator below has
# nonzero rational constant term.
def k(x=0):
    return (F(x), F(0), F(0), F(0))


ZERO, ONE = k(), k(1)
ALPHA = (F(0), F(1), F(0), F(0))
REDUCE4 = tuple(map(F, (1276420, -7849770, 28419741, 77822208)))
REDUCE4 = tuple(a / 72783360 for a in REDUCE4)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(a, b):
    return tuple(x * b for x in a)


def neg(a):
    return scale(a, -1)


def mul(a, b):
    if not any(a) or not any(b):
        return ZERO
    if not any(a[1:]):
        return scale(b, a[0])
    if not any(b[1:]):
        return scale(a, b[0])
    z = [F(0)] * 7
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            z[i + j] += x * y
    for d in range(6, 3, -1):
        for j in range(4):
            z[d - 4 + j] += z[d] * REDUCE4[j]
    return tuple(z[:4])


def total(values):
    out = ZERO
    for a in values:
        out = add(out, a)
    return out


def const(a):
    return [a] + [ZERO] * ORDER


def plus(a, b):
    return [add(x, y) for x, y in zip(a, b)]


def times(a, b):
    return [total(mul(a[i], b[n-i]) for i in range(n+1))
            for n in range(ORDER+1)]


def scal(a, c):
    return [scale(x, c) for x in a]


def minus(a, b):
    return plus(a, scal(b, -1))


def power(a, n):
    out = const(ONE)
    for _ in range(n):
        out = times(out, a)
    return out


def evaluate(coefficients, a):
    out = const(ZERO)
    for c in reversed(coefficients):
        out = plus(times(out, a), const(c))
    return out


def derivative(a):
    return [scale(a[i+1], i+1) for i in range(ORDER)] + [ZERO]


def inverse(a):
    check(a[0][0] != 0 and not any(a[0][1:]), "rational unit constant")
    out = [k(1/a[0][0])] + [ZERO] * ORDER
    for n in range(1, ORDER+1):
        out[n] = scale(total(mul(a[i], out[n-i]) for i in range(1,n+1)),
                       -1/a[0][0])
    check(times(a, out) == const(ONE), "formal series inverse")
    return out


def det(a, b):
    return minus(times(a[0], b[1]), times(a[1], b[0]))


def fmt(a):
    return str(a[0]) if not any(a[1:]) else [str(x) for x in a]


def emit(tag, *values):
    print(tag, json.dumps(values, separators=(",", ":")))


x, q, r = sp.symbols("x q r")
D = 1+x*x*q
B = (D-1)*(D+2)**2
C = x*D*(D+2)
E = q*(D+3)
jac = sp.diff(C,x)*sp.diff(E,q)-sp.diff(C,q)*sp.diff(E,x)
check(sp.expand(jac-6*(B+2)) == 0, "exact compiler Jacobian")
check(sp.expand(C*C*E-B*(B+4)) == 0, "surface equation")
f4 = 72783360*r**4-77822208*r**3-28419741*r**2+7849770*r-1276420
g3 = 104569906176*r**3+87403912326*r**2-19653044202*r+4049817787
check(sp.gcd(f4,g3) == 1, "no common root of fourth and fifth coefficients")
check(sp.Poly(f4,r).is_irreducible, "exceptional quotient is a field")
bez_a, bez_b, bez_gcd = sp.gcdex(f4,g3,r)
check(sp.expand(bez_a*f4+bez_b*g3) == 1 and bez_gcd == 1,
      "Bezout identity certifies every algebraic conjugate")
emit("SYMBOLIC", "det(C,E)/det(x,q)=6(B+2)", "gcd(F4,G3)=1",
     "F4 irreducible over Q")

J = sp.Matrix([[3,0,0,-1,0,-2],[-9,0,0,0,-1,2],
               [0,3,0,-1,0,0],[0,4,0,0,-1,0],
               [0,0,3,-1,0,-2],[0,0,9,0,-1,-2]])
check(J.det() == -288, "inherited implicit Jacobian determinant")
JI = [[F(a) for a in row] for row in J.inv().tolist()]


def q_abc(a,b,c):
    out = list(map(k, [F(-3,4),1,F(-27,4),-2,F(9,2),1,0,0,0]))
    for d, z in enumerate([a,b,c]):
        for degree, coefficient in [(2,1),(4,-2),(6,1)]:
            out[d+degree] = add(out[d+degree], scale(z,coefficient))
    return out


def compiler(Q, z, chi):
    ss = const(ZERO)
    ss[1] = ONE
    qz = plus(plus(evaluate(Q,z),ss),times(chi,z))
    dz = plus(const(ONE), times(power(z,2),qz))
    bz = times(minus(dz,const(ONE)),power(plus(dz,const(k(2))),2))
    cz = times(times(z,dz),plus(dz,const(k(2))))
    ez = times(qz,plus(dz,const(k(3))))
    qx = plus(evaluate([scale(Q[i],i) for i in range(1,len(Q))],z),chi)
    dx = plus(scal(times(z,qz),2),times(power(z,2),qx))
    cx = plus(times(dz,plus(dz,const(k(2)))),
              scal(times(times(z,plus(dz,const(ONE))),dx),2))
    ex = plus(times(qx,plus(dz,const(k(3)))),times(qz,dx))
    cq = scal(times(power(z,3),plus(dz,const(ONE))),2)
    eq = scal(plus(dz,const(ONE)),2)
    return bz,cz,ez,(cx,ex),(cq,eq)


def residual(Q, state):
    vals = [compiler(Q,z,state[5]) for z in state[:3]]
    return [minus(v[j],state[j+2]) for v in vals for j in (1,2)]


def solve_case(name,a,b,c,first,expected):
    Q = q_abc(a,b,c)
    state = [const(k(z)) for z in (-1,0,1,0,-3,0)]
    check(all(v[0] == ZERO for v in residual(Q,state)),name+" initial collision")
    for n in range(1,ORDER+1):
        residual_n = [v[n] for v in residual(Q,state)]
        delta = [neg(total(scale(z,coefficient) for z,coefficient
                           in zip(residual_n,row))) for row in JI]
        for j in range(6):
            state[j][n] = add(state[j][n],delta[j])
        check(all(not any(v[n]) for v in residual(Q,state)),
              name+" implicit coefficient "+str(n))
    check(all(v == const(ZERO) for v in residual(Q,state)),name+" full collision")
    chi = state[5]
    check(all(z == ZERO for z in chi[:first]) and chi[first] == expected,
          name+" first compensator coefficient")
    vals = [compiler(Q,z,chi) for z in state[:3]]
    check(vals[0][0] == vals[1][0] == vals[2][0],name+" B also collides")
    tangents = [v[3] for v in vals]
    raw = [det(tangents[1],tangents[2]),det(tangents[2],tangents[0]),
           det(tangents[0],tangents[1])]
    check([z[0] for z in raw] == list(map(k,[15,-54,39])),name+" raw orientation")
    lambdas = [times(z,inverse(scal(raw[1],-1))) for z in raw]
    pi = [total(z[n] for z in lambdas) for n in range(ORDER+1)]
    xi = [total(z[n] for z in [times(l,z) for l,z in zip(lambdas,state[:3])])
          for n in range(ORDER+1)]
    check(xi[0] == k(F(4,9)),name+" Xi is a unit")
    equation = plus(pi,times(derivative(chi),xi))
    check(all(z == ZERO for z in equation[:ORDER]),name+" all-order identity through cutoff")
    check(all(z == ZERO for z in pi[:first-1]) and
          pi[first-1] == scale(expected,F(-4*first,9)),name+" sharp period valuation")
    for slot in range(3):
        pullbacks = []
        for z,v in zip(state[:3],vals):
            rho = plus(const(ONE),times(derivative(chi),z))
            fixed_s = [times(vq,rho) for vq in v[4]]
            jj = det(v[3],fixed_s) if slot == 0 else v[3][slot-1]
            pullbacks.append(jj)
        period = [total(z[n] for z in [times(l,jj) for l,jj in zip(lambdas,pullbacks)])
                  for n in range(ORDER+1)]
        check(all(z == ZERO for z in period[:ORDER]),name+" descended slot "+str(slot))
    denominator = scal(plus(vals[0][0],const(k(2))),6)
    for z,v in zip(state[:3],vals):
        rho = plus(const(ONE),times(derivative(chi),z))
        fixed_s = [times(vq,rho) for vq in v[4]]
        density = times(det(v[3],fixed_s),inverse(denominator))
        check(density[:ORDER] == rho[:ORDER],name+" localized compatible unit density")
    # Polynomial-in-t clocks are controls for the analytic arbitrary-clock result.
    clocks = [(0,1),(0,0,1),(0,0,0,1),(0,1,1),(0,0,2,-1)]
    for coefficients in clocks:
        tt = sp.Symbol("t")
        phi = sum(sp.Integer(c)*tt**i for i,c in enumerate(coefficients))
        phi_order = next(i for i,c in enumerate(coefficients) if c)
        first_coeff = coefficients[phi_order]
        # Only coefficients up to the first live term are needed; later Pi
        # coefficients cannot cancel it in characteristic zero.
        expansion = sp.Poly(sp.expand(phi**(first-1)),tt)
        target_order = phi_order*(first-1)
        check(min(monom[0] for monom,_ in expansion.terms()) == target_order and
              expansion.nth(target_order) == first_coeff**(first-1),
              name+" clock valuation "+str(coefficients))
    emit("STRATUM",name,"chi-order",first,"chi-leading",fmt(expected),
         "Pi-order",first-1,"Pi-leading",fmt(pi[first-1]),"cutoff",ORDER)
    return chi,pi


solve_case("generic",ZERO,ZERO,ZERO,2,k(F(-259,9)))
solve_case("zero-chi2",k(F(-259,36)),ZERO,ZERO,3,k(F(16*5717,2187)))
p0 = F(-5717,729)
solve_case("zero-chi2-chi3",k(F(-259,36)+p0),ZERO,k(-p0),4,
           k(F(-16*1276420,177147)))
palpha = add(add(k(F(-5717,729)),scale(ALPHA,F(-1688,81))),
             scale(mul(ALPHA,ALPHA),F(520,9)))
a = add(add(k(F(-259,36)),palpha),scale(ALPHA,4))
b,c = scale(ALPHA,-9),neg(palpha)
g_alpha = total([k(4049817787),scale(ALPHA,-19653044202),
                 scale(mul(ALPHA,ALPHA),87403912326),
                 scale(mul(mul(ALPHA,ALPHA),ALPHA),104569906176)])
chi5 = scale(g_alpha,F(64,3**17))
exception_chi,exception_pi = solve_case("exceptional-field",a,b,c,5,chi5)
kappa = scale(chi5,F(-20,9))
check(exception_pi[4] == kappa and kappa != ZERO,"exceptional period equals inherited kappa")

# The density/source-coordinate boundary is an exact symbolic identity,
# independent of every finite coefficient control.
h = sp.Symbol("h")
X = x+h*x*x/2
check(sp.diff(X,x) == 1+h*x,"formal source Jacobian equals compatible density")
# An earlier source-coordinate change alters the sum of tangent weights.
l1,l2,l3,z1,z2,z3 = sp.symbols("l1 l2 l3 z1 z2 z3")
P = l1+l2+l3
Xi = l1*z1+l2*z2+l3*z3
transformed_sum = l1*(1+h*z1)+l2*(1+h*z2)+l3*(1+h*z3)
check(sp.expand(transformed_sum-(P+h*Xi)) == 0,"coordinate-covariant period")
emit("SCOPE", "all target 2-form slots tested", "formal unit density is allowed",
     "source shear changes fixed dx wedge dt benchmark", "no polynomial termination claim")
emit("PASS",CHECKS)
emit("SOURCE_SHA256",hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
