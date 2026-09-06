"""Exact identities and controls for simplicity in the full anchored model.

The all-parameter proof is in the companion report. No producer imports,
sampling, floating roots, or optimization are used. Ordinary D/B residue
moments are not Newton moments of D. Raw stdout is explicitly LF.
"""
from math import comb
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
GATES = 0


def require(ok, label):
    global GATES
    if not ok:
        raise RuntimeError(label)
    GATES += 1


def same(a, b, label):
    require(S.expand(a-b) == 0, label)


v, u, r, x, y, z, s = S.symbols("v u r x y z s")
B = v**5-13*v**4+55*v**3-x*v**2+y*v-z
C = v**4-12*v**3+45*v**2-S.Rational(2, 3)*x*v+S.Rational(3, 7)*y
D = v**4-11*v**3+36*v**2-S.Rational(5, 12)*x*v+y/7
den = [1, -13, 55, -x, y, -z]


def residue_moments(A, order):
    num_poly = S.Poly(S.expand(u**4*A.subs(v, 1/u)), u)
    num = [num_poly.nth(k) for k in range(order+1)]
    out = []
    for k in range(order+1):
        out.append(S.expand(num[k]-sum(den[j]*out[k-j]
                                     for j in range(1, min(5, k)+1))))
    for k in range(order+1):
        same(sum(den[j]*out[k-j] for j in range(min(5, k)+1)),
             num[k], "complete division product coefficient")
    return out


mu, nu = residue_moments(C, 8), residue_moments(D, 8)
expected_C = [1, 1, 3, x/3-16, 16*x/3-373-4*y/7,
              54*x-59*y/7+z-3969]
expected_D = [1, 2, 7, 7*x/12-19, 115*x/12-632-6*y/7,
              199*x/2-92*y/7+z-7171]
for k in range(6):
    same(mu[k], expected_C[k], "C/B formal residue moment")
    same(nu[k], expected_D[k], "D/B formal residue moment")
HC2 = S.Matrix(2, 2, lambda i, j:mu[i+j+1]).det()
HC3 = S.factor(S.Matrix(3, 3, lambda i, j:mu[i+j+1]).det())
HD3 = S.factor(S.Matrix(3, 3, lambda i, j:nu[i+j]).det())
same(HC2, (x-75)/3, "C shifted order-two determinant")
same(HC3, -x*(x-75)**2/27+S.Rational(15, 7)*(x-75)*y
     -S.Rational(16, 49)*y*y+(x-75)*z/3, "C shifted order-three determinant")
same(HD3, -S.Rational(49, 144)*x*x+S.Rational(269, 4)*x
     -S.Rational(18, 7)*y-3132, "D/B ordinary order-three determinant")
same(S.diff(B, v).subs(v, 0), y, "double zero forces y=0")
same(B.subs(v, 0), -z, "zero forces z=0")
same(HC3.subs({y:0, z:0}), -x*(x-75)**2/27, "double-zero C obstruction")
same(HD3.subs({x:75, y:0}), -S.Rational(37, 16), "double-zero D obstruction")
print("Residue division through moment 8 and all three Gram identities pass.")
print("Double zero: C forces x>=75 and then x=75; D/B ordinary H3 determinant=-37/16.")

# A repeated positive B root must be a root of both quartic interlacers.
coeff = S.Matrix([[S.diff(A, x).subs(v, r), S.diff(A, y).subs(v, r)] for A in (C, D)])
same(coeff.det(), r/12, "invertible common-root linear equations for r>0")
solution = S.solve([C.subs(v, r), D.subs(v, r)], [x, y])
X = S.Rational(24, 7)*r**3-36*r*r+108*r
Y = 3*r**4-28*r**3+63*r*r
same(solution[x], X, "forced repeated-root third coefficient")
same(solution[y], Y, "forced repeated-root fourth coefficient")
Z = S.expand((B+z).subs({v:r, x:X, y:Y}))
same(Z, S.Rational(4, 7)*r**5-5*r**4+10*r**3, "forced constant coefficient")
condition = 2*r*r-14*r+21
same(S.diff(B, v).subs({v:r, x:X, y:Y}),
     S.Rational(4, 7)*r*r*condition, "only two positive repeated-root candidates")


def reduced(poly):
    return S.rem(S.expand(poly), condition, r)


same(reduced(X), 126-12*r, "candidate x reduction")
same(reduced(Y), S.Rational(735, 4)-49*r, "candidate y reduction")
same(reduced(Z), S.Rational(441, 4)-42*r, "candidate z reduction")
same(reduced(HD3.subs({x:X, y:Y})), 5*r-S.Rational(75, 4), "candidate D/B Gram obstruction")
for eps in (-1, 1):
    rr = (7+eps*S.sqrt(7))/2
    require(S.simplify(condition.subs(r, rr)) == 0, "complete quadratic roots")
    yy = S.simplify((S.Rational(735, 4)-49*r).subs(r, rr))
    dd = S.simplify((5*r-S.Rational(75, 4)).subs(r, rr))
    require(rr.is_positive is True, "positive candidate root")
    if eps == 1:
        require(yy.is_negative is True, "upper candidate has negative elementary coefficient")
    else:
        require(dd.is_negative is True, "lower candidate has negative D/B H3")
    print(f"r=(7{'+sqrt(7)' if eps == 1 else '-sqrt(7)'})/2: y={yy}; D/B H3={dd}.")

# Exact hostile: C alone permits repeated roots; its D predicate fails.
hostile = {x:75, y:0, z:0}
same(B.subs(hostile), v*v*(v-3)*(v-5)**2, "C-only repeated beta hostile")
require(S.cancel(C.subs(hostile)/B.subs(hostile)-S.Rational(2, 3)/v
                 -S.Rational(1, 3)/(v-3)) == 0, "positive C-only atomic measure")
same(D.subs(hostile).subs(v, 5), -S.Rational(25, 4), "D fails the repeated-root contact")
same(HD3.subs(hostile), -S.Rational(37, 16), "hostile D ordinary moment failure")
print("C-only hostile B=v^2(v-3)(v-5)^2 retained: C/B=2/(3v)+1/(3(v-3)), D(5)=-25/4.")

# The entire simple-zero boundary has two fixed atoms and degree-four tails.
G = v**4-13*v**3+55*v*v-x*v+y
Ct = v**3-S.Rational(45, 4)*v*v+S.Rational(75, 2)*v-S.Rational(5, 12)*x
Dt = v**3-S.Rational(32, 3)*v*v+S.Rational(197, 6)*v-S.Rational(23, 72)*x
same(B.subs(z, 0), v*G, "zero-root factorization")
same(C, S.Rational(3, 7)*G+S.Rational(4, 7)*v*Ct, "fixed C zero-atom peeling")
same(D, G/7+S.Rational(6, 7)*v*Dt, "fixed D zero-atom peeling")
same(G.subs(v, 0), y, "nonzero remaining constant coefficient")
for vals in ({x:84, y:35, z:0}, {x:84, y:35, z:1}):
    BB = S.Poly(B.subs(vals), v)
    require(S.gcd(BB, BB.diff()).degree() == 0, "positive-model simple beta control")
    if vals[z] == 0:
        require(S.Poly(G.subs(vals), v).count_roots(0, S.oo) == 4, "four strictly positive residual beta roots")
    else:
        require(BB.count_roots(0, S.oo) == 5, "five strictly positive beta roots")
    for label, moments in (("C", mu), ("D", nu)):
        H = S.Matrix(5, 5, lambda i, j:moments[i+j].subs(vals) if hasattr(moments[i+j], "subs") else moments[i+j])
        minors = [H[:k, :k].det() for k in range(1, 6)]
        for minor in minors:
            require(minor > 0, "literal full positive residue model")
        print(f"(x,y,z)=(84,35,{vals[z]}), {label}/B leading Gram minors: {minors}.")

# Reconstruct the original carriers and retain them on the zero-root boundary.
def product(a, b):
    out = {}
    for i, ai in a.items():
        for j, bj in b.items():
            out[i+j] = out.get(i+j, 0)+ai*bj
    return {i:S.expand(a) for i, a in out.items()}


def add(a, b):
    return {i:S.expand(a.get(i, 0)+b.get(i, 0)) for i in set(a)|set(b)}


def shift(a, n, scale=1):
    return {i+n:scale*ai for i, ai in a.items()}


O = {j:comb(14, 2*j+1) for j in range(7)}
E = {j:comb(14, 2*j) for j in range(8)}
beta = dict(zip(range(-1, 5), [1, 13, 55, x, y, z]))
cr = dict(zip(range(-1, 4), [1, 12, 45, 2*x/3, 3*y/7]))
dr = dict(zip(range(-1, 4), [1, 11, 36, 5*x/12, y/7]))
PP = {i:S.expand(O[i]*beta[i]) for i in O.keys() & beta.keys()}
QQa = add(product(O, O), shift(product(E, E), -1))
QQb = add(product(beta, beta), shift(product(cr, dr), 1, 2))
QQ = {i:S.expand(QQa[i]*QQb[i]) for i in QQa.keys() & QQb.keys()}
same(sum(pi*(-s)**i for i, pi in PP.items()),
     182-20020*s+2002*x*s*s-3432*y*s**3+2002*z*s**4, "unchanged original phase equation")
same(QQ[-1], 28, "unchanged inverse response carry")
same(QQ[8].subs(z, 0), 0, "zero-root top response coefficient disappears")
require(S.expand(QQ[7].subs(z, 0)) != 0, "mixed response remains at exponent seven")
require(S.Poly(sum(pi.subs(z, 0)*(-s)**i for i, pi in PP.items()), s).degree() == 3,
        "original zero-boundary phase degree")
print("Zero-root atoms: C weights 3/7 and 4/7; D weights 1/7 and 6/7. Original O14/E14 and inverse carry28 retained.")
print("Compactness supplies a positive mutual beta-root gap, with no explicit constant and no positive distance from zero.")
print("No B-only simplicity, universal strict reference, boundary sign, or general Laurent noncancellation claimed.")
print(f"PASS: {GATES} always-active exact gates.")
