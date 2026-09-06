"""Independent exact referee: symbolic quotient algebra and literal source products.

No producer import. Requires SymPy; all gates remain active under python -O.
Run beside the producer, or from 04-computation with producer outputs in results.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
BASE = Path(__file__).resolve().parent
STEM = "continuing4_20260906_moments_packet"
PINS = {
    ".py": "1b7d02c20c631f95b8a539af8948bf92a9b6dc1504ddc05d5a00b6b13892a9bc",
    ".out": "54e693133895b249e31aeeed325d786ee4ed6c46d8343fa082e802adcbf6932e",
    "_certificate.json": "e396e221c4cd8e4b668f0f07a572ed5314b25cb7c66a8ac68eb439d84979170f",
}
GATES = 0


def gate(value, label):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(label)


def locate(name):
    for p in (BASE / name, BASE.parent / "05-knowledge" / "results" / name):
        if p.exists():
            return p
    raise FileNotFoundError(name)


for ext, expected in PINS.items():
    gate(hashlib.sha256(locate(STEM + ext).read_bytes()).hexdigest() == expected,
         "frozen producer " + ext)
cert = json.loads(locate(STEM + "_certificate.json").read_text(encoding="utf-8"))
v, q, t, s = S.symbols("v q t s")
x, y, z = S.symbols("x y z", real=True)
Q = S.Rational
B = v**5 - 13*v**4 + 55*v**3 - x*v**2 + y*v - z
C = v**4 - 12*v**3 + 45*v**2 - Q(2, 3)*x*v + Q(3, 7)*y
D = v**4 - 11*v**3 + 36*v**2 - Q(5, 12)*x*v + Q(1, 7)*y


def eq(a, b, label):
    if isinstance(a, S.MatrixBase) or isinstance(b, S.MatrixBase):
        gate(all(S.cancel(entry) == 0 for entry in a-b), label)
    else:
        gate(S.cancel(a - b) == 0, label)


def coeffs(poly):
    p = S.Poly(poly, v)
    return [p.nth(i) for i in range(p.degree() + 1)]


def from_coeffs(cs):
    return sum(Q(c) * v**i for i, c in enumerate(cs))


def expansion(A, denominator, n=10):
    # Literal expansion of A(1/q)/[q B(1/q)], not a recurrence implementation.
    transformed = S.cancel(A.subs(v, 1/q) / (q*denominator.subs(v, 1/q)))
    series = S.series(transformed, q, 0, n).removeO().expand()
    return [series.coeff(q, j) for j in range(n)]


def hankel(ms, size, shift=0):
    return S.Matrix(size, size, lambda i, j: ms[i+j+shift])


def all_principal_signs(M, strict, label):
    for size in range(1, M.rows + 1):
        for inds in combinations(range(M.rows), size):
            d = M.extract(inds, inds).det(method="domain-ge")
            gate(d > 0 if strict else d >= 0, label + str(inds))


def real_count(poly, lo=-S.oo, hi=S.oo):
    return S.Poly(poly, v).count_roots(lo, hi)


def simple(poly):
    return S.Poly(S.gcd(poly, S.diff(poly, v)), v).degree() == 0


# UNIVERSAL identities, separate from all finite controls.
symbolic = {}
for label, A, starts in [
    ("C", C, [1, 1, 3, x/3-16, 16*x/3-373-4*y/7]),
    ("D", D, [1, 2, 7, 7*x/12-19, 115*x/12-632-6*y/7]),
]:
    ms = expansion(A, B)
    symbolic[label] = ms
    for j in range(5):
        eq(ms[j], starts[j], label + " literal initial coefficient")
    for j in range(5, 10):
        eq(ms[j], 13*ms[j-1]-55*ms[j-2]+x*ms[j-3]-y*ms[j-4]+z*ms[j-5],
           label + " denominator recurrence")
    H = hankel(ms, 5)
    # Columns are v times the five monomials, reduced modulo B.
    T = S.zeros(5)
    for j in range(4):
        T[j+1, j] = 1
    T[:, 4] = S.Matrix([z, -y, x, -55, 13])
    for entry in H*T-T.T*H:
        eq(entry, 0, label + " symbolic selfadjoint companion")
    for j in range(10):
        eq((H*T**j)[0, 0], ms[j], label + " quotient reproduces formal moments")

# Parent's additional fibre corollary: fixed x,y gives two affine PSD pencils.
for label, slopes in [
    ("C", [1, 14, 130, 904+4*x/3]),
    ("D", [1, 15, 147, 1067+19*x/12]),
]:
    for j in range(9):
        eq(S.diff(symbolic[label][j], z), 0 if j < 5 else slopes[j-5],
           label + " exact affine z coefficient")
        eq(S.diff(symbolic[label][j], z, 2), 0, label + " no hidden z curvature through degree8")

a, b = S.symbols("a b", real=True)
r = a + S.I*b
xr = Q(24, 7)*r**3 - 36*r**2 + 108*r
yr = 3*r**4 - 28*r**3 + 63*r**2
eq(C.subs({v: r, x: xr, y: yr}), 0, "common zero parameterization C")
eq(D.subs({v: r, x: xr, y: yr}), 0, "common zero parameterization D")
eq(S.im(xr).expand()/b, Q(24, 7)*(3*a*a-b*b-21*a+Q(63, 2)), "imaginary x")
b2 = 3*a*a - 21*a + Q(63, 2)
eq(S.im(yr).expand()/b,
   12*a**3-12*a*b*b-84*a*a+28*b*b+126*a, "imaginary y before elimination")
eq(S.cancel(S.im(yr).expand()/b).subs(b*b, b2),
   -6*(2*a-7)*(2*a*a-14*a+21), "common nonreal resultant condition")
eq(b2.subs(a, Q(7, 2)), -Q(21, 4), "linear branch impossible")
eq(b2, Q(3, 2)*(2*a*a-14*a+21), "quadratic branch forces b zero")
u = S.symbols("u", positive=True)
eq(B.subs(v, -u), -(u**5+13*u**4+55*u**3+x*u*u+y*u+z), "support nonnegativity")

# A PSD quotient alone does not recover an arbitrary denominator.
generic_B = (v-1)*(v*v+1)**2
generic_A = (v*v+1)**2
generic_ms = expansion(generic_A, generic_B)
gate(generic_ms == [1]*10, "cancellation hostile moments")
gate(hankel(generic_ms, 5).rank() == 1 and real_count(generic_B) == 1,
     "one positive Hankel can hide four nonreal roots")

# DEGREE-SEVEN HOSTILE, recovered from exact original rational polynomials.
vals = {x: Q(78071, 1000), y: Q(601, 50), z: Q(127473806477, 203203019250)}
s0 = Q(57, 2)
cap = Q(707, 100)
for name, val in [("x", vals[x]), ("y", vals[y]), ("z", vals[z]), ("s", s0), ("support_cap", cap)]:
    eq(val, Q(cert["witness"][name]), "witness pin " + name)
Bb = B.subs(vals)
gate(simple(Bb), "hostile beta is simple")
gate(real_count(Bb) == 3 and real_count(Bb, 0, S.oo) == 3, "hostile beta three positive real roots")
eq(vals[z], 12*vals[y]/(7*s0)-vals[x]/s0**2+10/s0**3-1/(11*s0**4), "exact original phase equation")

O = sum(S.binomial(14, k)*t**((k-1)//2) for k in range(1, 15, 2))
E = sum(S.binomial(14, k)*t**(k//2) for k in range(0, 15, 2))
beta = (1+13*t+55*t*t+x*t**3+y*t**4+z*t**5)/t
raw_C = (1+12*t+45*t*t+Q(2, 3)*x*t**3+Q(3, 7)*y*t**4)/t
raw_D = (1+11*t+36*t*t+Q(5, 12)*x*t**3+Q(1, 7)*y*t**4)/t
left = S.expand(O*O + E*E/t)
right = S.expand(beta*beta + 2*t*raw_C*raw_D)
# Independent ordinary SymPy multiplication and literal coefficient extraction.
P = sum(S.expand(O).coeff(t, j)*S.expand(beta).coeff(t, j)*t**j for j in range(-1, 7))
response = sum(left.coeff(t, j)*right.coeff(t, j)*t**j for j in range(-2, 15))
Pphase = S.expand(P.subs(vals).subs(t, -v))
eq(Pphase.subs(v, s0), 0, "literal first source zero")
eq(Pphase, from_coeffs(cert["P_minus_s"]), "first polynomial entire coefficient pin")
gate(simple(Pphase) and real_count(Pphase) == 4 and real_count(Pphase, 0, S.oo) == 4,
     "four positive original phases")
resp_w = S.expand(response.subs(vals))
for j in range(-2, 15):
    eq(resp_w.coeff(t, j), Q(cert["Q_coefficients"].get(str(j), "0")), "every literal Q coefficient")
eq(resp_w.coeff(t, -1), 28, "lowest nonzero carrier")
for j in range(-1, 9):
    gate(resp_w.coeff(t, j) != 0, "all response carriers retained")
Qvalue = S.cancel(resp_w.subs(t, -s0))
eq(Qvalue, Q(cert["Q_at_minus_s"]), "positive full Q value")
gate(Qvalue > 0, "actual carried response is positive")
without_cross = sum(left.coeff(t, j)*S.expand(beta*beta).coeff(t, j)*t**j for j in range(-2, 15))
gate(S.cancel((response-without_cross).subs(vals).subs(t, -s0)) != 0, "cross carry is material")

denominators = []
for label, A in [("C", C), ("D", D)]:
    ms = expansion(A.subs(vals), Bb)
    frozen = cert[label]
    for actual, stored in zip(ms, frozen["moments"]):
        eq(actual, Q(stored), label + " hostile literal moment")
    H, K = hankel(ms, 4), hankel(ms, 4, 1)
    matrices = {"ordinary4": H, "shifted4": K, "upper4": cap*H-K,
                "quadratic_localizer3": cap*hankel(ms, 3, 1)-hankel(ms, 3, 2)}
    for name, M in matrices.items():
        all_principal_signs(M, True, label + name)
        for j in range(1, M.rows+1):
            eq(M[:j, :j].det(), Q(frozen["leading_minors"][name][j-1]), label + name + " frozen leading minor")
    det5 = hankel(ms, 5).det(method="domain-ge")
    gate(det5 < 0, label + " first omitted square fails")
    eq(det5, Q(frozen["ordinary5_determinant"]), label + " negative determinant exact pin")
    T = H.inv()*K
    gate(H*T == T.T*H, label + " quadrature selfadjoint")
    basis = S.eye(4)
    for j in range(3):
        gate(T*basis[:, j] == basis[:, j+1], label + " cyclic monomial shift")
    # Different reconstruction from producer moment linear system:
    # characteristic polynomial, spectral resolvent via adjugate, powers of T.
    quadrature_B = T.charpoly(v).as_expr()
    quadrature_A = S.expand((H*(v*S.eye(4)-T).adjugate())[0, 0])
    eq(quadrature_B, from_coeffs(frozen["quadrature_denominator"]), label + " characteristic polynomial")
    eq(quadrature_A, from_coeffs(frozen["quadrature_numerator"]), label + " resolvent numerator")
    gate(simple(quadrature_B) and real_count(quadrature_B, 0, cap) == 4,
         label + " four distinct interior quadrature nodes")
    gate(S.degree(S.gcd(quadrature_A, quadrature_B), v) == 0, label + " no quadrature zero weights")
    for j in range(8):
        eq((H*T**j)[0, 0], ms[j], label + " quadrature exact degree-seven representation")
    m8quad = (H*T**8)[0, 0]
    eq(m8quad, Q(frozen["quadrature_eighth_moment"]), label + " quadrature moment8")
    gate(m8quad > ms[8], label + " native eighth moment too small")
    eq(ms[8]-m8quad, det5/H.det(), label + " missing orthogonal square Schur complement")
    denominators.append(quadrature_B)
    print(label + " det(H5) = " + str(det5))
    print(label + " four-atom measure has moment8 strictly above native moment8.")
gate(S.degree(S.gcd(*denominators), v) == 0, "two canonical quadrature supports disjoint")

elementary = [1, 13, 55, vals[x], vals[y], vals[z]]
normalized = [S.cancel(elementary[k]/S.binomial(5, k)) for k in range(6)]
for k in range(1, 5):
    margin = normalized[k]**2-normalized[k-1]*normalized[k+1]
    gate(margin > 0, "strict Newton inequality")
    eq(margin, Q(cert["newton_margins"][k-1]), "Newton margin pin")
gate(cap > Q(13, 5) and 5*cap**2-26*cap-67 < 0, "support cap below exact endpoint")
gate(5*Q(71, 10)**2-26*Q(71, 10)-67 > 0, "exact endpoint below 7.1")
gate(vals[y] > Q(161875, 888583) and vals[x] > 75, "inherited residue floors")
gate(vals[y] > Q(8750, 8241)*(vals[x]-75), "inherited residue slope")
gate(5*vals[z] < cap*vals[y], "five z endpoint consequence")

# FULL-MODEL POSITIVE AND SINGULAR CONTROLS.
for control in cert["controls"]:
    name = control["name"]
    cv = {symbol: Q(control[str(symbol)]) for symbol in (x, y, z)}
    cb = B.subs(cv)
    ca = [C.subs(cv), D.subs(cv)]
    common = S.gcd(cb, S.gcd(*ca))
    eq(common, from_coeffs(control["gcd"]), name + " exact common factor")
    gate(simple(cb) and real_count(cb, 0, S.oo) == 5, name + " all five beta roots positive distinct")
    for label, numerator in zip(("C", "D"), ca):
        ms = expansion(numerator, cb)
        H5 = hankel(ms, 5)
        all_principal_signs(H5, name == "strict", name + label)
        for j in range(1, 6):
            eq(H5[:j, :j].det(), Q(control["leading_minors"][label][j-1]), name + label + " leading minors pin")
        gate(H5.rank() == (5 if name == "strict" else 4), name + label + " exact rank")
        if name == "strict":
            nodes = [Q(k, 5) for k in (1, 3, 9, 22, 30)]
            eq(cb, S.prod(v-r for r in nodes), "strict literal roots")
            residues = [S.cancel(numerator.subs(v, r)/S.diff(cb, v).subs(v, r)) for r in nodes]
            gate(all(w > 0 for w in residues), label + " literal positive residues")
            eq(sum(residues), 1, label + " probability normalization")
        else:
            eq(H5*S.Matrix([-Q(9, 7), Q(123, 7), -25, 10, -1]), S.zeros(5, 1),
               label + " reduced beta spans Hankel radical")
    if name == "common_node":
        root_boxes = [(Q(8, 100), Q(9, 100)), (Q(104, 100), Q(105, 100)),
                      (Q(225, 100), Q(226, 100)), (Q(662, 100), Q(663, 100))]
        reduced = S.cancel(cb/(v-3))
        for lo, hi in root_boxes:
            gate(real_count(reduced, lo, hi) == 1, "singular reduced beta root box")
        for label, numerator, boxes in [
            ("C", ca[0], [(Q(5, 10), Q(7, 10)), (Q(20, 10), Q(21, 10)), (Q(63, 10), Q(64, 10))]),
            ("D", ca[1], [(Q(2, 10), Q(3, 10)), (Q(16, 10), Q(17, 10)), (Q(61, 10), Q(62, 10))]),
        ]:
            for j, (lo, hi) in enumerate(boxes):
                gate(root_boxes[j][1] < lo < hi < root_boxes[j+1][0], label + " strict reduced interlacing order")
                gate(real_count(S.cancel(numerator/(v-3)), lo, hi) == 1, label + " reduced numerator box")

cv = {x: 75, y: 0, z: 0}
eq(B.subs(cv), v*v*(v-3)*(v-5)**2, "C-only repeated boundary")
msC, msD = expansion(C.subs(cv), B.subs(cv)), expansion(D.subs(cv), B.subs(cv))
eq(S.cancel(C.subs(cv)/B.subs(cv)), Q(2, 3)/v+Q(1, 3)/(v-3), "C-only positive measure")
all_principal_signs(hankel(msC, 5), False, "C-only singular PSD")
eq(hankel(msD, 3).det(), -Q(37, 16), "second channel rejects repeated hostile")

print("Literal original Q(-57/2) = " + str(Qvalue))
print("Universal recurrence/selfadjoint identities and no-common-nonreal-zero obstruction: PASS.")
print("Fixed x,y: both degree-eight Hankels affine in z; decoded admissible z fibre is a compact interval (or empty).")
print("All principal minors, exact Sturm counts, two spectral quadratures, boundaries and carries: PASS.")
print("Full degree-eight geometry decoder accepted; finite-phase response sign remains OPEN.")
print("PASS: " + str(GATES) + " always-active exact gates.")
