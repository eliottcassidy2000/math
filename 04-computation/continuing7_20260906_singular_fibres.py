"""Exact original-response sign on every shared-root singular fibre.

Standalone source: native binomial convolution, formal residue moments,
quartic quotient arithmetic, rational intervals, and always-active gates.
No producer imports. All certificates and stdout are actual UTF-8 LF.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb, floor, ceil
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
GATES = 0


def require(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(label)


def rat(a):
    return F(int(S.numer(a)), int(S.denom(a)))


def box(a, b=None):
    return F(a), F(a if b is None else b)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def mul(a, b):
    z = [p * q for p in a for q in b]
    return min(z), max(z)


def inv(a):
    require(a[0] * a[1] > 0, "nonzero interval denominator")
    return 1 / a[1], 1 / a[0]


def ev(poly, var, interval):
    ans = box(0)
    for c in S.Poly(poly, var).all_coeffs():
        ans = add(mul(ans, interval), box(rat(c)))
    return ans


def ev2(poly, rbox, xbox):
    ans = box(0)
    for c in S.Poly(poly, x).all_coeffs():
        ans = add(mul(ans, xbox), ev(c, r, rbox))
    return ans


def rational_ev(expr, interval):
    n, d = S.fraction(S.cancel(expr))
    return mul(ev(n, r, interval), inv(ev(d, r, interval)))


def reduction(expr, P):
    n, d = S.fraction(S.cancel(expr))
    require(S.degree(S.gcd(d, P), r) == 0, "quartic denominator is a unit")
    out = S.rem(S.rem(n, P, r) * S.invert(d, P, r), P, r)
    require(S.rem(n - d * out, P, r) == 0, "exact quotient coefficient identity")
    return S.expand(out)


def wire(poly, var=None):
    if var is None:
        var = r
    return [str(c) for c in S.Poly(poly, var).all_coeffs()]


def strbox(b):
    return list(map(str, b))


def signed(expr, interval, direction, label):
    b = ev(expr, r, interval)
    require(b[1] < 0 if direction < 0 else b[0] > 0, label)
    return strbox(b)


r, u, x, y, z, s, v = S.symbols("r u x y z s v")
B = v**5 - 13*v**4 + 55*v**3 - x*v**2 + y*v - z
C = v**4 - 12*v**3 + 45*v**2 - S.Rational(2, 3)*x*v + S.Rational(3, 7)*y
D = v**4 - 11*v**3 + 36*v**2 - S.Rational(5, 12)*x*v + S.Rational(1, 7)*y
phase = z*s**4 - S.Rational(12, 7)*y*s**3 + x*s**2 - 10*s + S.Rational(1, 11)
cap = S.Rational(7, 72)*(x-75)*(135-x)
YC = S.Rational(14, 9)*x*r - S.Rational(7, 3)*r**2*(r**2-12*r+45)
ZC = r**2*(S.Rational(5, 9)*x-r*(4*r**2-45*r+150)/3)
YD = S.Rational(35, 12)*x*r-7*r**2*(r**2-11*r+36)
ZD = r**2*(S.Rational(23, 12)*x-r*(6*r**2-64*r+197))
AC, AD = 5*u**2-24*u+9, 23*u**2-60*u+12
HC = 44*r**2*u**4-132*r**2*u**3-495*r*u**4+1584*r*u**3-3*r+1650*u**4-5940*u**3+330*u
HD = 66*r**2*u**4-132*r**2*u**3-704*r*u**4+1452*r*u**3-r+2167*u**4-4752*u**3+110*u
PC = 705672*r**4-15079284*r**3+115450055*r**2-375696750*r+441317250
PD = 120434688*r**4-2310536448*r**3+15142183583*r**2-39453584808*r+35022385200
UC = (11286*r**2-109085*r+245025)/(27126*r**2-261360*r+586025)
UD = (446688*r**2-4023865*r+7744176)/(1978416*r**2-17739216*r+34339426)
charts = {
    "C": (C, YC, ZC, AC, HC, 9, 33, PC, UC, [3528453826, 3550824804, 6081020695, 8208387541]),
    "D": (D, YD, ZD, AD, HD, 12, 11, PD, UD, [2049574992, 2831681376, 6130326664, 8173391713]),
}

# Reconstruct the ORIGINAL two rows, including both carries. A separate
# binomial addition identity verifies every convolution weight.
beta = {-1: 1, 0: 13, 1: 55, 2: x, 3: y, 4: z}
cc = {-1: 1, 0: 12, 1: 45, 2: S.Rational(2, 3)*x, 3: S.Rational(3, 7)*y}
dd = {-1: 1, 0: 11, 1: 36, 2: S.Rational(5, 12)*x, 3: S.Rational(1, 7)*y}
raw, weights = {}, {}
for a, aa in beta.items():
    for b, bb in beta.items():
        raw[a+b] = raw.get(a+b, 0) + aa*bb
for a, aa in cc.items():
    for b, bb in dd.items():
        raw[a+b+1] = raw.get(a+b+1, 0) + 2*aa*bb
for a in range(7):
    for b in range(7):
        weights[a+b] = weights.get(a+b, 0) + comb(14, 2*a+1)*comb(14, 2*b+1)
for a in range(8):
    for b in range(8):
        weights[a+b-1] = weights.get(a+b-1, 0) + comb(14, 2*a)*comb(14, 2*b)
for j, weight in weights.items():
    require(weight == comb(28, 2*j+2), "independent closed binomial convolution weight")
q = {j: S.expand(c*weights[j]) for j, c in raw.items() if j in weights}
Q = S.expand(sum(c*(-s)**j for j, c in q.items()))
require(sorted(q) == list(range(-1, 9)) and q[-1] == 28, "full original support and inverse carry")
require(S.expand(sum(beta[j]*comb(14, 2*j+1)*(-s)**j for j in range(5))-2002*phase) == 0,
        "literal original first-row normalization")


def moments(top):
    den, ans = [1, -13, 55, -x, y, -z], []
    for k in range(7):
        ans.append(S.expand((top[k] if k < len(top) else 0)
                            - sum(den[j]*ans[k-j] for j in range(1, min(k, 5)+1))))
    return ans


cm = moments([1, -12, 45, -S.Rational(2, 3)*x, S.Rational(3, 7)*y])
dm = moments([1, -11, 36, -S.Rational(5, 12)*x, S.Rational(1, 7)*y])


def gram(m, n, shift=0):
    return S.expand(S.det(S.Matrix(n, n, lambda i, j: m[i+j+shift])))


H4 = gram(cm, 4)
require(gram(cm, 2, 1) == (x-75)/3, "C shifted two determinant")
require(S.expand(gram(cm, 3) - (S.Rational(1, 9)*(x-75)*(135-x)-S.Rational(8, 7)*y)) == 0,
        "C ordinary three determinant and curved packet")

certificate = {"scope": "all singular shared-root original phases; general off-wall sign OPEN",
               "original_q": {str(j): str(c) for j, c in sorted(q.items())},
               "C_moments": list(map(str, cm)), "D_moments": list(map(str, dm)),
               "C_ordinary_four": str(H4), "charts": {}, "fibres": {}, "controls": {}}
root_boxes, reduced_q = {}, {}
for name, (A, Y, Z, a, h, k, den, P, U, initials) in charts.items():
    require(S.expand(B.subs({v: r, y: Y, z: Z})) == 0 and S.expand(A.subs({v: r, y: Y})) == 0,
            "literal shared-root chart " + name)
    require(S.factor(phase.subs({y: Y, z: Z, s: u/r}) - x*u*u*a/(k*r*r)+h/(den*r)) == 0,
            "singular original-phase affine identity " + name)
    require(S.degree(S.cancel(S.resultant(a, h, u)/P), r) == 0, "complete singular resultant " + name)
    require(reduction(a.subs(u, U), P) == reduction(h.subs(u, U), P) == 0,
            "rational reconstruction of every singular pair " + name)
    require(S.Poly(P, r).count_roots(0, S.oo) == 4, "all four positive singular roots " + name)
    qchart = S.Poly(Q.subs({y: Y, z: Z, s: u/r}), x)
    require(qchart.degree() == 2, "singular response remains quadratic " + name)
    qa, qb, qc = [reduction(c.subs(u, U), P) for c in qchart.all_coeffs()]
    reduced_q[name] = (qa, qb, qc)
    certificate["charts"][name] = {"P": str(P), "U": str(U), "Y": str(Y), "Z": str(Z),
                                   "response_coefficients_descending_x": [wire(c) for c in (qa, qb, qc)]}
    for j, initial in enumerate(initials, 1):
        aa, bb = S.Rational(initial, 10**9), S.Rational(initial+1, 10**9)
        require(S.Poly(P, r).count_roots(aa, bb) == 1, "inherited interval contains one root")
        refined = S.Poly(P, r).refine_root(aa, bb, eps=S.Rational(1, 10**34))
        low = floor(refined[0]*10**28)
        lo, hi = S.Rational(low, 10**28), S.Rational(low+1, 10**28)
        require(lo < refined[0] <= refined[1] < hi and S.Poly(P, r).count_roots(lo, hi) == 1,
                "exact refined decimal isolator")
        R = box(rat(lo), rat(hi))
        root_boxes[name+str(j)] = R
        ub = rational_ev(U, R)
        require(ub[0] > 0, "positive original u coordinate")
        record = {"r_interval": strbox(R), "u_interval": strbox(ub),
                  "s_interval": strbox(mul(ub, inv(R)))}
        certificate["fibres"][name+str(j)] = record
        if j == 4:
            require(R[0] > F(71, 10), "largest pair violates actual beta-root maximum")
            record["result"] = "excluded by r<71/10"

require(F(71, 10)**2 + (13-F(71, 10))**2/4 > 59,
        "exact Cauchy beta-root ceiling")

# Four degree-four packet caps. Negative derivative and concavity make each
# endpoint sign an exclusion of the ENTIRE x-ray beyond that endpoint.
limits = {"C1": 94, "C3": 89, "D1": 102, "D2": 98, "D3": 95}
for key, limit in limits.items():
    R = root_boxes[key]
    Y, Z = charts[key[0]][1:3]
    gap = S.expand(cap-Y)
    record = certificate["fibres"][key]
    record["x_upper_strict"] = limit
    record["cap_at_upper"] = signed(gap.subs(x, limit), R, -1, key+" curved cap is negative")
    record["cap_derivative_at_upper"] = signed(S.diff(gap, x).subs(x, limit), R, -1,
                                               key+" cap decreases on excluded ray")
    require(S.diff(gap, x, 2) == -S.Rational(7, 36), "strict global concavity of curved cap")
    if key == "D3":
        record["z_at_upper"] = signed(Z.subs(x, limit), R, -1, "D3 z negative below x=95")
        record["z_slope"] = signed(S.diff(Z, x), R, 1, "D3 z strictly increasing")
        record["result"] = "empty already under C degree-four packet and z>=0"

# The only failed degree-four response projection is repaired by the next
# ordinary C Hankel matrix, using moments only through degree six.
R = root_boxes["C2"]
shifted = S.Poly(S.expand(H4.subs({y: YC, z: ZC}).subs(x, x+88)), x)
shifted_coeffs = [S.rem(c, PC, r) for c in shifted.all_coeffs()]
require(shifted.degree() == 4, "degree-four shifted Hankel determinant")
record = certificate["fibres"]["C2"]
record["x_upper_strict"] = 88
record["H4_shifted_coefficients_descending_x"] = [wire(c) for c in shifted_coeffs]
record["H4_shifted_coefficient_intervals"] = [signed(c, R, -1, "negative coefficient of H4(88+X)")
                                                for c in shifted_coeffs]

# All five surviving necessary fibres: negative curvature plus positive
# derivative at the upper cap puts the whole feasible interval on the
# increasing, negative part of the same original response parabola.
for key, limit in {"C1": 94, "C2": 88, "C3": 89, "D1": 102, "D2": 98}.items():
    a, b, c = reduced_q[key[0]]
    R = root_boxes[key]
    rec = certificate["fibres"][key]
    rec["Q_x2"] = signed(a, R, -1, key+" negative response curvature")
    rec["Q_derivative_at_upper"] = signed(2*a*limit+b, R, 1, key+" response increases up to cap")
    rec["Q_at_upper"] = signed(a*limit**2+b*limit+c, R, -1, key+" negative original response at cap")
    rec["result"] = "Q(-s)<0 for every admissible x on the complete singular fibre"
    print("WHOLE_FIBRE", key, "x<", limit, "Q_upper<", ceil(F(rec["Q_at_upper"][1])))

# Cheap hostile: retain original phase, shared root, both order-five
# Stieltjes packets and the curved cap, but delete the C degree-six gate.
hostile = {}
for label, expression in [("x_minus75", x-75), ("y", YC), ("z", ZC), ("C_curved_cap", cap-YC)]:
    hostile[label] = signed(expression.subs(x, 92), R := root_boxes["C2"], 1, "hostile passes " + label)
for label, mm in [("C", cm), ("D", dm)]:
    for shift in (0, 1):
        for size in range(1, 4):
            expr = gram(mm, size, shift).subs({x: 92, y: YC.subs(x, 92), z: ZC.subs(x, 92)})
            hostile[label+str(shift)+str(size)] = signed(S.rem(S.expand(expr), PC, r), R, 1,
                                                        "hostile positive definite order-five packet")
a, b, c = reduced_q["C"]
hostile["Q"] = signed(a*92**2+b*92+c, R, 1, "hostile has positive ORIGINAL response")
hostile["H4"] = signed(S.rem(S.expand(H4.subs({x: 92, y: YC.subs(x, 92), z: ZC.subs(x, 92)})), PC, r),
                         R, -1, "degree-six C determinant rejects hostile")
certificate["controls"]["C2_x92_order_five_hostile"] = hostile

# Nonvacuity on an explicit interval, not only an isolated model point.
# Remove the known shared factor symbolically; then preserve all remaining
# real roots and their ordering uniformly over r and x intervals.
XBOX = box(F(86999, 1000), F(87001, 1000))
R = root_boxes["C2"]
BB = S.expand(B.subs({y: YC, z: ZC}))
CC = S.expand(C.subs(y, YC))
DD = S.expand(D.subs(y, YC))
Bquot, br = S.div(BB, v-r, v)
Cquot, cr = S.div(CC, v-r, v)
require(br == cr == 0, "exact shared factor without any x division")
brackets = {"B": [(80,100), (560,610), (2440,2470), (6320,6340)],
            "C": [(400,410), (1910,1930), (6120,6140)],
            "D": [(180,200), (1470,1500), (3380,3400), (5930,5950)]}
nonvacuous = {"x_interval": strbox(XBOX), "r_interval": strbox(R), "brackets": []}
for label, poly in [("B", Bquot), ("C", Cquot), ("D", DD)]:
    for aa, bb in brackets[label]:
        lval = ev2(S.expand(poly.subs(v, S.Rational(aa, 1000))), R, XBOX)
        rval = ev2(S.expand(poly.subs(v, S.Rational(bb, 1000))), R, XBOX)
        require(lval[1] < 0 < rval[0] or rval[1] < 0 < lval[0], "uniform actual root bracket " + label)
        nonvacuous["brackets"].append([label, aa, bb, strbox(lval), strbox(rval)])
require(F(100,1000) < F(180,1000) < F(200,1000) < F(400,1000) < F(410,1000) < F(560,1000),
        "first beta gap contains both first interlacer roots")
require(F(610,1000) < F(1470,1000) < F(1500,1000) < F(1910,1000) < F(1930,1000) < F(2440,1000),
        "second beta gap contains both second interlacer roots")
require(F(2470,1000) < F(3380,1000) < F(3400,1000) < R[0], "third gap ends at known shared C root")
require(R[1] < F(5930,1000) < F(5950,1000) < F(6120,1000) < F(6140,1000) < F(6320,1000),
        "last beta gap contains both final interlacer roots")
require(ev2(YC, R, XBOX)[0] > 0 and ev2(ZC, R, XBOX)[0] > 0, "positive y,z throughout actual interval")
certificate["controls"]["C2_actual_full_model_interval"] = nonvacuous

here = Path(__file__).resolve().parent
dest = here.parent / "05-knowledge/results" if here.name == "04-computation" else here
path = dest / (Path(__file__).stem + "_certificate.json")
path.write_bytes((json.dumps(certificate, indent=2, sort_keys=True)+"\n").encode("utf-8"))
print("EXCLUDED C4,D4 beta maximum; D3 curved C cap and z>=0")
print("HOSTILE C2,x92: original positive response, both degree-five Stieltjes packets PD; C degree-six H4 negative")
print("NONVACUOUS C2 full-model interval x in [86999/1000,87001/1000]; original singular phase fixed")
print("SCOPE every singular shared-positive-root fibre is negative; general off-wall sign OPEN")
print("CERTIFICATE_SHA256", hashlib.sha256(path.read_bytes()).hexdigest())
print("PASS", GATES, "always-active exact gates; actual LF")
