"""Independent exact controls for the seven-point tetrahedron classification.

All checks remain active under python -O. No producer code is imported.
"""
import hashlib
import json
from pathlib import Path

import sympy as S

V = S.Matrix
ZERO = V([0, 0, 0])
HALF = S.Rational(1, 2)
QUARTER = S.Rational(1, 4)


def require(truth, message):
    if not bool(truth):
        raise RuntimeError(message)


def same(a, b):
    d = a - b
    if isinstance(d, S.MatrixBase):
        return all(S.simplify(x) == 0 for x in d)
    return S.simplify(d) == 0


def geometry(vertices):
    """Derive the construction independently from face and line equations."""
    require(all(same(x.dot(x), 1) for x in vertices), "Not unit vertices")
    v, a, b, c = vertices
    tetra_det = S.simplify(S.det(V.hstack(a-v, b-v, c-v)))
    n = (b-a).cross(c-a)
    require(n.dot(n) != 0, "Collinear face")
    w = S.simplify(n.dot(a) * n / n.dot(n))
    s = a+b+c
    g = (v+s)/4
    d = g-v
    A = S.simplify(v-2*v.dot(d)*d/d.dot(d))
    points = {"A": A, "C": v, "D": 4*g,
              "F": 4*w, "H": s-2*w}
    if g.dot(g) != 0:
        points["E"] = S.simplify(g/g.dot(g))
    if w.dot(w) != 0:
        # Solve the three tangent plane equations, rather than use w/||w||^2.
        tangent_matrix = V.vstack(a.T, b.T, c.T)
        points["B"] = S.simplify(tangent_matrix.inv()*V([1, 1, 1]))
        require(same(points["B"], w/w.dot(w)), "Tangent intersection formula")
    require(same(A.dot(A), 1), "Second intersection is not on the original sphere")
    require(not same(A, v), "Second intersection unexpectedly tangent")
    for x, y, z in [(a,b,c), (b,c,a), (c,a,b)]:
        require(same((points["H"]-x).dot(y-z), 0), "Orthocenter altitude")
    return {"s": s, "w": w, "g": g, "points": points,
            "tetra_det": tetra_det}


def sphere_ranks(points):
    rows = [[2*x[0], 2*x[1], 2*x[2], -1] for x in points]
    M = V(rows)
    b = V([x.dot(x) for x in points])
    return M.rank(), M.row_join(b).rank()


def verify_sphere(points, center, radius2):
    return all(same((x-center).dot(x-center), radius2) for x in points)


report = {}

# Generic rational counterexample: all definitions finite, no common sphere.
generic = geometry([V(x) for x in [(0,0,-1),(1,0,0),(0,1,0),(0,0,1)]])
P = generic["points"]
require(generic["tetra_det"] != 0 and len(P) == 7, "Counterexample assumptions")
require(sphere_ranks([P[x] for x in ["A","C","D","E"]]) == (3,4),
        "Counterexample four-point sphere system must be inconsistent")
circle_det = S.factor(V([[x[0], x[2], x.dot(x), 1]
                        for x in [P[k] for k in ["A","C","D","E"]]]).det())
require(circle_det == -S.Rational(56,9), "Counterexample circle determinant")
require(same((P["B"]-P["H"]).cross(P["F"]-P["H"]), ZERO), "Collinearity")
require(len({tuple(P[k]) for k in ["B","H","F"]}) == 3, "Distinct collinear points")
require(verify_sphere([P[k] for k in ["A","C","E"]]+[ZERO],
                      V([1,1,-HALF]), S.Rational(9,4)),
        "Nearby A,C,O,E circle theorem control")
report["generic_counterexample"] = {"four_point_ranks": [3,4],
                                    "circle_determinant": str(circle_det),
                                    "three_distinct_collinear_labels": ["H","B","F"]}

# Valid exceptional example, and exact absence of any additional coincidence.
r = S.sqrt(3)/2
face = [V([r,0,HALF]), V([0,r,HALF]), V([0,-r,HALF])]
valid = geometry([V([0,1,0])]+face)
P = valid["points"]
center = V([0,HALF,1])
require(valid["tetra_det"] != 0, "Valid tetrahedron degenerates")
require(verify_sphere(list(P.values()), center, S.Rational(5,4)), "Valid sphere residuals")
require(sphere_ranks(list(P.values())) == (4,4), "Unique sphere rank")
groups = {}
for label, x in P.items():
    groups.setdefault(tuple(x), []).append(label)
coincidences = sorted([sorted(g) for g in groups.values() if len(g)>1])
require(coincidences == [["B","F"],["D","E"]] and len(groups) == 5,
        "Unexpected coincidence in valid example")
affine_det = S.simplify(V.hstack(P["A"]-P["C"], P["H"]-P["C"], P["D"]-P["C"]).det())
require(affine_det != 0, "Valid exceptional points must not be concyclic")
report["valid_exception"] = {"tetrahedron_determinant": str(valid["tetra_det"]),
                            "center": list(map(str,center)), "radius_squared": "5/4",
                            "sphere_ranks": [4,4], "distinct_points": 5,
                            "coincidences": coincidences,
                            "four_point_affine_determinant": str(affine_det)}

# Rank-collapse endpoint: the three scalar conditions alone do NOT suffice.
endpoint = geometry([V([r,0,-HALF])]+face)
s,w,g,P = endpoint["s"],endpoint["w"],endpoint["g"],endpoint["points"]
v = P["C"]
require(endpoint["tetra_det"] != 0, "Negative endpoint should be a genuine tetrahedron")
require(same(16*g.dot(g),4) and same(w.dot(w),QUARTER) and same(s.dot(s),3),
        "Negative endpoint scalar conditions")
require(same(v,s-4*w) and same(v.dot(w),-QUARTER), "Negative endpoint dependence")
require(sphere_ranks(list(P.values())) == (3,4), "Negative endpoint inconsistency")
report["negative_endpoint"] = {"tetrahedron_determinant": str(endpoint["tetra_det"]),
                              "all_three_scalar_conditions_hold": True,
                              "sphere_ranks": [3,4], "no_finite_sphere": True}

# The positive endpoint is not a tetrahedron; if admitted, uniqueness fails.
coplanar = geometry([V([-r,0,HALF])]+face)
P = coplanar["points"]
require(coplanar["tetra_det"] == 0, "Positive endpoint must be coplanar")
require(sphere_ranks(list(P.values())) == (3,3), "Positive endpoint underdetermination")
require(verify_sphere(list(P.values()),V([0,0,1]),1), "Coplanar sphere 1")
require(verify_sphere(list(P.values()),V([0,1,1]),2), "Coplanar sphere 2")
report["positive_endpoint"] = {"tetrahedron_determinant": "0", "sphere_ranks": [3,3],
                              "nonunique_spheres_if_degenerate_inputs_allowed": True}

# Definition degeneracies: regular tetrahedron centroid, and great-circle face.
q = S.sqrt(3)/3
centered = geometry([V(x)*q for x in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]])
require(centered["tetra_det"] != 0 and same(centered["g"],ZERO), "Centered tetrahedron")
require("E" not in centered["points"], "G=O inversion must be undefined")
greatface = geometry([V(x) for x in [(0,0,1),(1,0,0),(-1,0,0),(0,1,0)]])
require(greatface["tetra_det"] != 0 and same(greatface["w"],ZERO), "Great-circle face")
require("B" not in greatface["points"], "w=O tangent intersection must be undefined")
faceM = V([[1,0,0],[-1,0,0],[0,1,0]])
require(faceM.rank() == 2 and faceM.row_join(V([1,1,1])).rank() == 3,
        "Great-circle tangent planes must be inconsistent, not merely nonunique")
report["undefined_controls"] = {"G_equals_O": "E undefined",
                                "w_equals_O": "B undefined; tangent system ranks 2,3"}

# Coordinate-free symbolic controls using the exceptional Gram matrix.
p = S.symbols("p", real=True)
Gram = V([[1,0,p],[0,3,S.Rational(3,4)],[p,S.Rational(3,4),QUARTER]])
require(same(Gram.det(),3*(S.Rational(1,16)-p*p)), "Gram determinant formula")
coeff = V([1,4*p,4])/(2*(1+4*p))
require(same(Gram*coeff,V([HALF,S.Rational(3,2),HALF])), "Center scalar equations")
require(same((coeff.T*Gram*coeff)[0],(5+12*p)/(4*(1+4*p))), "Radius formula")
affine_coeff = V([[-S.Rational(3,2),-1,0],[HALF,1,1],[0,-2,0]])
require(affine_coeff.det() == -3, "Affine determinant factor")
T = S.symbols("T", positive=True)
c = (T-4)/3
require(same(16/T-16*(1+c)/T+c,(T-4)*(T-16)/(3*T)), "Four-point residual")
report["symbolic_controls"] = {"gram_determinant": str(S.factor(Gram.det())),
                              "center_and_radius_formulas": "verified",
                              "affine_determinant_factor": -3,
                              "four_point_residual": "(T-4)(T-16)/(3T)"}

# Independently audit the producer's subsequently added generic six-point core.
core = geometry([V([0,0,-1]), V([1,0,0]), V([0,1,0]),
                 V([S.Rational(3,5),0,S.Rational(4,5)])])
s,w,g,P = core["s"],core["w"],core["g"],core["points"]
v = P["C"]
T0 = 16*g.dot(g)
k0 = w.dot(w)
delta = v.dot(s.cross(w))
require(delta != 0, "Generic core independence")
# Obtain the center by solving the actual sphere equations, not the cross formula.
M = V.vstack((2*v).T,(2*s).T,(2*w).T)
center0 = S.simplify(M.inv()*V([5,T0-1,1+4*k0]))
cross_center = S.simplify((5*s.cross(w)+(T0-1)*w.cross(v)+(1+4*k0)*v.cross(s))/(2*delta))
require(same(center0,cross_center), "Reciprocal-basis center")
require(same(center0,V([-S.Rational(43,24),S.Rational(37,6),-S.Rational(5,2)])),
        "Generic exact center")
rho0 = S.simplify(center0.dot(center0)-4)
require(rho0 == S.Rational(25049,576), "Generic exact radius")
core_points = [P[k] for k in ["B","C","D","E","F"]]+[4*v]
require(verify_sphere(core_points,center0,rho0), "Generic six-point residuals")
require(sphere_ranks(core_points) == (4,4), "Generic unique sphere ranks")
outside = {label: S.simplify(P[label].dot(P[label])-2*center0.dot(P[label])+4)
           for label in ["A","H"]}
require(outside == {"A":S.Rational(1178,225),"H":S.Rational(38,5)}, "Strict outside values")
for label1,label2 in [("D","E"),("B","F")]:
    require(same(4*P[label1]/P[label1].dot(P[label1]),P[label2]), "Radius-2 inverse pair")
# Symbolic residual and inversion checks.
T,k,vs,mu,rr,dot = S.symbols("T k vs mu rr dot", nonzero=True)
residuals = [1-5+4,
             T-(T+4)+4,
             16/T-4*(T+4)/T+4,
             1/k-(1+4*k)/k+4,
             16*k-4*(1+4*k)+4,
             16-4*5+4]
require(all(S.simplify(x)==0 for x in residuals), "Symbolic six-point residuals")
require(same((T-1-2*vs)-8*k-((T-1)-2*(1+4*k))+4,6-2*vs),
        "Symbolic H residual")
require(same(-mu*((T+4)-4*5),mu*(16-T)), "Symbolic A residual")
require(same(16/rr-8*dot/rr+4,4*(rr-2*dot+4)/rr), "Radius-2 inversion invariance")
report["generic_six_point_core"] = {"delta":str(delta),
                                    "center":list(map(str,center0)),
                                    "radius_squared":str(rho0),
                                    "sphere_ranks":[4,4],
                                    "outside_residuals":{k:str(v) for k,v in outside.items()},
                                    "symbolic_residuals_and_inversion":"verified"}

# Review precisely identifies which producer artifacts were inspected.
base = Path(__file__).resolve().parent.parent / "session23b-geometry"
report["reviewed_sources"] = {}
for name in ["geometry_report.md","verify_geometry.py","verification.txt"]:
    path = base / name
    if path.exists():
        report["reviewed_sources"][name] = hashlib.sha256(path.read_bytes()).hexdigest()
report["status"] = "PASS"
print(json.dumps(report,indent=2))
