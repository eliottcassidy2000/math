"""Exact joint-D certificate excluding the closed nonpositive circuit octant.

Native B,C,D moments and coefficients are retained. The full proof is in the
companion report. The finite sign certificate contains EVERY tensor Bernstein
coefficient of two polynomials on the unit square (133+40 coefficients).
No floating-point computation or prior producer import is used.
"""
from pathlib import Path
from hashlib import sha256
from itertools import product
from math import comb
import json
import sys
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", newline="\n")
HERE = Path(__file__).resolve().parent
DEST = HERE.parent/"05-knowledge/results" if HERE.name == "04-computation" else HERE
STEM = Path(__file__).stem
x, y, z, u, v = sp.symbols("x y z u v")
gates = 0


def need(test, label):
    global gates
    gates += 1
    if not bool(test):
        raise RuntimeError(label)


def native_moments(which, count=9):
    numerator = ([1, -12, 45, -sp.Rational(2, 3)*x, sp.Rational(3, 7)*y]
                 if which == "C" else
                 [1, -11, 36, -sp.Rational(5, 12)*x, sp.Rational(1, 7)*y])
    denominator = [1, -13, 55, -x, y, -z]
    moments = []
    for j in range(count):
        moments.append(sp.expand((numerator[j] if j < len(numerator) else 0)
                       - sum(denominator[i]*moments[j-i]
                             for i in range(1, min(5, j)+1))))
    return moments


def hankel(mm, size):
    return sp.Matrix([[mm[i+j] for j in range(size)] for i in range(size)])


md = native_moments("D")
mc = native_moments("C")
H3 = sp.factor(hankel(md, 3).det())
H4 = sp.factor(hankel(md, 4).det())
need(sp.expand(H3 + (343*x**2-67788*x+2592*y+3157056)/sp.Integer(1008)) == 0,
     "native D ordinary3 identity")
need(sp.diff(H4, z, 2) == -6, "strict quadratic concavity in z")
need(sp.Poly(H4, x, y, z).length() == 15, "all15 native determinant terms")

xc = sp.Rational(831875, 8788)
yc = sp.Rational(3460080078125, 52206766144)
zc = sp.Rational(791547381622314453125, 52414354446428935168)
q = sp.Rational(275, 338)
ylo = sp.Rational(13, 166375)*x**3
yhi = (-343*x**2+67788*x-3157056)/sp.Integer(2592)
zlo = 44*y**3/x**3
need(ylo.subs(x, xc) == yc and zlo.subs({x: xc, y: yc}) == zc,
     "circuit tie normalization")
need(0 < xc < 99, "x chart endpoints")
g = sp.factor(343*x**2-67788*x+3157056+2592*ylo)
need(g.subs(x, 99) > 0, "x99 excluded by joint y bounds")
shifted_derivative = sp.Poly(sp.expand(sp.diff(g, x).subs(x, 99+u)), u)
need(all(co > 0 for co in shifted_derivative.all_coeffs()),
     "g strictly increasing for x at least99")
need(sp.expand(H3 + sp.Rational(18, 7)*(y-yhi)) == 0, "H3 is the upper y constraint")


def primitive_integer(expression, variables):
    pp = sp.Poly(sp.cancel(expression), *variables, domain=sp.QQ)
    den, integral = pp.clear_denoms()
    content, primitive = integral.primitive()
    need(den > 0 and content > 0, "positive normalization")
    need(sp.expand(expression - content*primitive.as_expr()/den) == 0,
         "primitive numerator identity")
    return primitive, den/content


P0, scale0 = primitive_integer(-x**6*H4.subs(z, zlo), (x, y))
P1, scale1 = primitive_integer(-x**3*sp.diff(H4, z).subs(z, zlo), (x, y))
need(P0.degree_list() == (10, 6) and P1.degree_list() == (5, 3),
     "native sign polynomial bidegrees")


def chart(pp):
    return sp.Poly(sp.expand(pp.as_expr().subs(y, ylo+(yhi-ylo)*v)
                           .subs(x, xc+(99-xc)*u)), u, v, domain=sp.QQ)


def bernstein(pp):
    degrees = pp.degree_list()
    power = dict(pp.terms())
    records = []
    for idx in product(*(range(d+1) for d in degrees)):
        coefficient = sum(co*sp.prod(sp.Rational(comb(idx[j], exp[j]), comb(degrees[j], exp[j]))
                              for j in range(2))
                          for exp, co in power.items() if all(exp[j] <= idx[j] for j in range(2)))
        need(coefficient > 0, "every Bernstein coefficient strictly positive")
        records.append([list(idx), coefficient])
    # Reconstruct by evaluation at a complete degree-bounded tensor grid.
    # This supplements the closed basis-change identity with a distinct
    # finite-polynomial-identity path and keeps every coefficient in the JSON.
    for i in range(degrees[0]+1):
        U = sp.Rational(i, max(1, degrees[0]))
        for j in range(degrees[1]+1):
            V = sp.Rational(j, max(1, degrees[1]))
            rebuilt = sum(co*comb(degrees[0], ii)*U**ii*(1-U)**(degrees[0]-ii)
                          *comb(degrees[1], jj)*V**jj*(1-V)**(degrees[1]-jj)
                          for (ii, jj), co in records)
            need(rebuilt == pp.eval({u: U, v: V}), "complete tensor-grid identity")
    return {"degrees": list(degrees), "power_terms": pp.terms(),
            "bernstein": records, "minimum": min(co for idx, co in records)}


Q0, Q1 = chart(P0), chart(P1)
need(Q0.degree_list() == (18, 6) and Q1.degree_list() == (9, 3),
     "complete transformed bidegrees")
C0, C1 = bernstein(Q0), bernstein(Q1)
need(len(C0["bernstein"]) == 133 and len(C1["bernstein"]) == 40,
     "complete173-coefficient universe")


def point(xx, yy, zz):
    subst = {x: xx, y: yy, z: zz}
    circuits = [sp.Rational(831875, 8788)/xx, 13*xx**3/(166375*yy),
                44*yy**3/(xx**3*zz)]
    row = {"xyz": [xx, yy, zz], "circuits": circuits,
           "R": [1/q, circuits[0]/q, circuits[0]*circuits[1]/q,
                 sp.prod(circuits)/q]}
    for A, mm in [("C", mc), ("D", md)]:
        numeric = [m.subs(subst) if isinstance(m, sp.Basic) else m for m in mm]
        row[A+"_minors"] = [hankel(numeric, j).det() for j in range(1, 6)]
    return row


# The canonical correction: H4 alone can be positive in the strict
# all-negative circuit/Newton region, but H3 rejects that surrogate.
bad = point(sp.Rational(4159375, 38207),
            sp.Rational(4325100097656250, 42626885365649),
            sp.Rational(494717113513946533203125000000, 13159366797968048415097695043))
need(bad["circuits"] == [sp.Rational(2939, 3380), sp.Rational(29201, 29390),
                        sp.Rational(276701, 292010)], "exact inherited determinant hostile")
need(all(0 < c < 1 for c in bad["circuits"]) and all(r > 1 for r in bad["R"]),
     "hostile is in strict negative-circuit Newton region")
need(bad["D_minors"][2] < 0 < bad["D_minors"][3],
     "H3 rejects while H4 determinant alone survives")

# Conversely, H3 alone survives an actual positive-root strict C-only point
# in the negative octant; H4 supplies the missing rejection.
r = sp.Rational(2047, 2048)
converse = point(xc/r, yc/r**4, zc/r**10)
need(converse["circuits"] == [r, r, r], "strict C-only negative octant")
need(all(a > 0 for a in converse["C_minors"]), "strict C-only positive residue geometry")
need(converse["D_minors"][2] > 0 > converse["D_minors"][3],
     "H3 alone survives; H4 rejects")
center = point(xc, yc, zc)
need(center["circuits"] == [1, 1, 1] and center["D_minors"][2] > 0 > center["D_minors"][3],
     "closed tie boundary also rejected")

# Two-negative pair extensions are not automatic. These exact surrogates
# retain all Newton inequalities and the full ordinary D packet through6.
pair_survivors = []
for word, coords in [("--+", (95, 68, 1)), ("+--", (86, 50, 9))]:
    rec = point(*map(sp.Integer, coords))
    need([int(sp.sign(c-1)) for c in rec["circuits"]] == [1 if ch=="+" else -1 for ch in word],
         "exact two-negative survivor word")
    need(all(r > 1 for r in rec["R"]) and all(m > 0 for m in rec["D_minors"][:4]),
         "two-negative survivor passes full D ordinary degree6 and Newton")
    need(rec["D_minors"][4] < 0, "native D degree8 rejects pair survivor")
    if word == "+--":
        need(all(m > 0 for m in rec["C_minors"]), "C-only pair survivor retains full C geometry")
        Btest = sp.Poly(x**5-13*x**4+55*x**3-86*x**2+50*x-9, x)
        need(Btest.count_roots(0, sp.oo) == 5 and sp.degree(sp.gcd(Btest, Btest.diff())) == 0,
             "literal positive simple beta roots for C-only pair survivor")
    rec["word"] = word
    pair_survivors.append(rec)

# Nonvacuity and scope: exact strict full-model points in four other words.
positives = []
for word, coords in [("+++", ("77.454", "8.902", ".02558")),
                     ("++-", ("77.3613", "8.6001", ".17694")),
                     ("+-+", ("86.2333", "51.3919", "8.6469")),
                     ("-++", ("97.7028", "70.6021", "14.5020"))]:
    rec = point(*map(sp.Rational, coords))
    need(all(a > 0 for a in rec["C_minors"]+rec["D_minors"]), "strict full model positive control")
    need([int(sp.sign(c-1)) for c in rec["circuits"]] == [1 if ch=="+" else -1 for ch in word],
         "exact positive control circuit word")
    rec["word"] = word
    positives.append(rec)


def encode(value):
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, sp.Basic):
        return str(value)
    return value


certificate = {"status": "Joint-D nonpositive-octant proof candidate pending independent audit",
               "D_moments": md[:7], "D_H3": H3, "D_H4": H4,
               "x0": xc, "y_lower": ylo, "y_upper": yhi, "z_lower": zlo,
               "g": g, "g_at99": g.subs(x, 99), "g_derivative_shift": shifted_derivative.as_expr(),
               "P0": P0.as_expr(), "P1": P1.as_expr(), "scale0": scale0, "scale1": scale1,
               "chart0": C0, "chart1": C1,
               "H4_only_hostile": bad, "H3_only_hostile": converse,
               "tie_center": center, "two_negative_survivors": pair_survivors,
               "positive_controls": positives, "gates": gates}
payload = (json.dumps(encode(certificate), sort_keys=True, indent=2)+"\n").encode("utf-8")
(DEST/(STEM+"_certificate.json")).write_bytes(payload)
print("NATIVE D-H3 and15-term D-H4 reconstructed by formal division; dzzH4=-6")
print("DOMAIN x0=831875/8788 <= x <99; y_lower <= y <= y_upper; z >=44y^3/x^3")
print("BERNSTEIN P0 degree(18,6):133 strictly positive coefficients; full tensor identity PASS")
print("BERNSTEIN P1 degree(9,3):40 strictly positive coefficients; full tensor identity PASS")
print("CONCLUSION joint D-H3,H4 positivity excludes every C2<=1,C3<=1,C4<=1 point")
print("HOSTILES H4-only positive/H3 negative; H3-only positive/H4 negative; ties excluded")
print("PAIR_SURVIVORS --+ at95,68,1 and +-- at86,50,9 pass D through6; D through8 rejects both; latter retains strict C")
print("POSITIVE full-model controls +++,++-,+-+,-++; two-negative classification remains OPEN")
print("CERTIFICATE_SHA256", sha256(payload).hexdigest())
print("PASS", gates, "always-active exact gates; raw LF")
