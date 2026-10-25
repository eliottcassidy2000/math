"""Exact controls for the seven-support Cremona repair and actual hourglass.

Run with Python 3 + SymPy, normally and under -O. No producer imports.
The finite universes and geometric limits are specified in RESULT.md.
"""
from itertools import combinations, permutations, product
from collections import Counter
import hashlib
import json
from pathlib import Path

import sympy as s

checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def norm(v):
    for x in v:
        if x:
            return tuple(s.simplify(y / x) for y in v)
    raise ValueError("zero projective vector")


def pg2(q):
    return [(1, a, b) for a in range(q) for b in range(q)] + [
        (0, 1, a) for a in range(q)
    ] + [(0, 0, 1)]


supports = list(range(1, 8))
graph_supports = [(a, b) for a in supports for b in supports
                  if a & b == 0 or a == b == 7]
fixed_strata = [(a, b) for a, b in graph_supports
                if a.bit_count() == b.bit_count() == 1]
seven = [p for p in graph_supports if p not in fixed_strata]
require(len(graph_supports) == 13, "all Cremona graph strata")
require(len(fixed_strata) == 6 and len(seven) == 7, "rank refinement counts")
require(sorted(a for a, _ in seven) == supports, "blowdown A1 bijection")
for a, b in seven:
    require(b == (7 if a == 7 else 7 ^ a), "modified complement action")
require((1, 6) in seven and (3, 4) in seven, "specialization hostile states")
require(1 & 3 == 1 and not (6 & 4 == 6), "blowdown does not reflect order")

finite_counts = {}
for q in (2, 3, 5, 7):
    points = pg2(q)
    graph = [(a, b) for a in points for b in points
             if len({(a[i] * b[i]) % q for i in range(3)}) == 1]
    removed = [(a, b) for a, b in graph
               if sum(x != 0 for x in a) == sum(x != 0 for x in b) == 1]
    require(len(graph) == q*q + 4*q + 1, "projective graph count")
    require(len(removed) == 6, "six rational fixed strata")
    require(len(graph)-len(removed) == q*q+4*q-5, "open graph count")
    finite_counts[q] = {"P2": len(points), "Gamma": len(graph),
                        "V": len(graph)-len(removed)}
require(finite_counts[2]["P2"] == finite_counts[2]["V"] == 7, "F2 coincidence")
require(finite_counts[3]["V"] != finite_counts[3]["P2"], "F3 hostile")

# Coordinate flags are the six torus-fixed points. They carry the native
# transitive-tournament complement, independently of the six boundary labels.
orders = list(permutations(range(3)))
edges = list(combinations(range(3), 2))
def mask(order):
    return tuple(order.index(i) < order.index(j) for i, j in edges)
def flag(order):
    return (1 << order[0], (1 << order[0]) | (1 << order[1]))
adjacencies = []
for a, b in combinations(orders, 2):
    one_arc = sum(x != y for x, y in zip(mask(a), mask(b))) == 1
    one_flag = sum(x != y for x, y in zip(flag(a), flag(b))) == 1
    require(one_arc == one_flag, "flag/tournament adjacency")
    if one_arc:
        adjacencies.append((a, b))
for a in orders:
    reverse = a[::-1]
    require(mask(reverse) == tuple(not bit for bit in mask(a)), "converse")
    require(flag(reverse) == (7 ^ flag(a)[1], 7 ^ flag(a)[0]), "fan inversion")
require(len(adjacencies) == 6, "hexagon edges")

fan_rays = [(1,0),(1,1),(0,1),(-1,0),(-1,-1),(0,-1)]
for i, ray in enumerate(fan_rays):
    nxt = fan_rays[(i+1) % 6]
    require(ray[0]*nxt[1]-ray[1]*nxt[0] == 1, "smooth hourglass slope fan")
    require(tuple(-v for v in ray) == fan_rays[(i+3) % 6], "opposition")
require(set(fan_rays) == {(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)},
        "actual peripheral slope labels")

# Anticanonical monomials are cubics vanishing at all coordinate points.
monomials = [e for e in product(range(4), repeat=3) if sum(e) == 3
             and max(e) < 3]
require(len(monomials) == 7, "seven anticanonical sections")
quadratic_weights = {tuple(a[i]+b[i] for i in range(3))
                     for a, b in combinations(monomials, 2)}
quadratic_weights |= {tuple(2*x for x in a) for a in monomials}
require(len(quadratic_weights) == 19, "antibicanonical dimension")
require(28-len(quadratic_weights) == 9, "nine toric quadratic relations")

E, x, z, X, Y, Z, lam = s.symbols("E x z X Y Z lam")
old = [E+x, E-x, E+z, E-z, E+x+z, E-x-z]
B = s.Matrix([[1,-1,-1],[-1,-1,0],[-1,2,-1]])
require(B.det() == 5, "noncollinear Cremona centers")
base_sub = dict(zip((E, x, z), B * s.Matrix([X,Y,Z])))
new = [s.expand(f.subs(base_sub, simultaneous=True)) for f in old]
require(new == [-2*Y-Z, 2*X-Z, Y-2*Z, 2*X-3*Y,
                -X-2*Z, 3*X-2*Y], "actual transformed hourglass labels")
pairing = [(0,2), (1,4), (3,5)]
c = {X:Y*Z, Y:X*Z, Z:-X*Y}
F = s.prod(new)
require(s.expand(F.subs(c, simultaneous=True)-(X*Y*Z)**2*F) == 0,
        "actual sextic Cremona identity")
cc = [s.expand(c[t].subs(c, simultaneous=True)) for t in (X,Y,Z)]
require(cc == [-X*X*Y*Z, -X*Y*Y*Z, -X*Y*Z*Z], "projective involution")
for (i,j), factor, scalar in zip(pairing, (X,Y,Z), (-1,-1,1)):
    require(s.expand((new[i]*new[j]).subs(c, simultaneous=True)
                     -scalar*factor**2*new[i]*new[j]) == 0, "paired branch pullback")

rows = [s.Matrix([[s.diff(f,t) for t in (X,Y,Z)]]) for f in new]
nodes = {}
for i,j in combinations(range(6),2):
    point = norm(rows[i].T.cross(rows[j].T))
    incidence = tuple(k for k,row in enumerate(rows) if (row*s.Matrix(point))[0] == 0)
    require(len(incidence) == 2, "six branch lines have no triples")
    nodes[point] = incidence
require(len(nodes) == 15, "all fifteen double-plane nodes")
centers = {norm(s.eye(3).col(i)) for i in range(3)}
remaining = [p for p in nodes if p not in centers]
require(len(remaining) == 12, "twelve nodes after blowing up centers")
node_orbits = set()
for point in remaining:
    a,b,d = point
    image = norm((b*d,a*d,-a*b))
    require(image in nodes and image != point, "Cremona pairs all remaining nodes")
    node_orbits.add(tuple(sorted((point,image), key=str)))
require(len(node_orbits) == 6, "six two-node orbits")

fixed_values = []
for u,v in product((1,-1), (s.I,-s.I)):
    value = s.simplify(F.subs({X:1,Y:u,Z:v})/(u*v)**2)
    require(value != 0, "fixed base point avoids branch")
    fixed_values.append(str(value))

# A new coordinate axis is not an original branch line. The pullback of
# L0/L1 has odd valuation at X=0: this blocks a lift to the full root field.
require(all(row[0,1] != 0 or row[0,2] != 0 for row in rows), "X not branch")
ratio = s.factor((new[0]/new[1]).subs(c, simultaneous=True))
require(s.simplify(ratio - X*(Y-2*Z)/(Y*(X+2*Z))) == 0,
        "odd new-axis valuation witness")

# Auxiliary symmetric family. Its exact critical parameter values are
# derived in RESULT.md; this independently checks the elimination identities.
a,b = s.symbols("a b", nonzero=True)
g = a+1/a+b+1/b+a/b+b/a
da, db = s.expand(a*s.diff(g,a)), s.expand(b*s.diff(g,b))
require(s.factor(da+db) == (a+b)*(a*b-1)/(a*b), "critical split")
require(s.factor(da.subs(b,1/a)) == (a-1)*(a+1)*(a*a+a+1)/a**2,
        "critical branch ab=1")
require(s.simplify(da.subs(b,-a)-(a-1/a)) == 0, "critical branch a+b=0")
critical = {(1,1):-6, (-1,-1):2, (1,-1):2, (-1,1):2}
for (aa,bb), parameter in critical.items():
    require(da.subs({a:aa,b:bb}) == db.subs({a:aa,b:bb}) == 0,
            "rational critical point")
    require(g.subs({a:aa,b:bb})+parameter == 0, "rational critical value")
cyclo = a*a+a+1
require(s.rem(s.together((g.subs(b,1/a)+3)).as_numer_denom()[0],cyclo,a)==0,
        "primitive cube-root critical value")
sym = sum(t**4*u**2 for t,u in permutations((X,Y,Z),2)) + lam*X*X*Y*Y*Z*Z
require(s.expand(sym.subs(lam,2)-(X*X+Y*Y)*(X*X+Z*Z)*(Y*Y+Z*Z))==0,
        "auxiliary six-line degeneration")

receipt = {
    "status":"FINITE-EXACT; geometric proofs are separate",
    "checks":checks,
    "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "support_strata":{"Gamma":13,"V":7,"removed_rank_zero":6},
    "finite_field_counts":finite_counts,
    "transitive_tournament_vertices":6,
    "transitive_tournament_edges":len(adjacencies),
    "hourglass_fan_rays":fan_rays,
    "anticanonical_dimensions":[7,19],
    "actual_hourglass":{"basis_matrix":str(B),"determinant":5,
                        "line_pairing_zero_based":pairing,"nodes":15,
                        "remaining_node_pairs":len(node_orbits),
                        "fixed_base_branch_values":fixed_values,
                        "nonlift_ratio":str(ratio)},
    "auxiliary_bad_parameters":[-6,2,3]
}
print(json.dumps(receipt, indent=2, sort_keys=True))
