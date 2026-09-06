"""Independent exact referee. No producer imports; complete slope pencils.

Native triangle extraction uses every least-row vertex, with no modular-conic
partition. Parameter identities use independent symbolic polynomial algebra.
"""
from pathlib import Path
from itertools import combinations
from math import gcd
from collections import defaultdict
import hashlib
import json
import sys
import argparse
import sympy as S

sys.stdout.reconfigure(newline="\n")
N = 0


def need(ok, why):
    global N
    N += 1
    if not ok:
        raise ArithmeticError(why)


here = Path(__file__).resolve().parent
stem = "continuing7_20260906_reciprocal_second_move"
parser = argparse.ArgumentParser()
parser.add_argument("--producer-dir", type=Path, default=None)
args = parser.parse_args()
producer = args.producer_dir or (here if (here / (stem + ".py")).exists()
                                 else here.parent / "continuing7_20260906_wildcard")
src = producer / (stem + ".py")
dest = producer.parent / "05-knowledge/results" if producer.name == "04-computation" else producer
certpath = dest / (stem + "_certificate.json")
outpath = dest / (stem + ".out")
for path, expected in [
    (src, "1e1571440d0ea569b7e3f65ffdc1ebe9ae61b2ae3d1f11cc39a7a6667d4b708c"),
    (certpath, "f5117590d335c3df7579b8bb483d2a0a9d0906ab67c06bbb1cfd01ff7a579c3f"),
    (outpath, "060493fb23ff746151420891fd9cce0df5eb3e0a7ebd7b4e0a56d3a84f2d328c"),
]:
    need(hashlib.sha256(path.read_bytes()).hexdigest() == expected, "frozen producer pin")
certificate = json.loads(certpath.read_bytes())


def board(p):
    return [1, 0] + [next(y for y in range(1, p) if x*y % p == 1) for x in range(2, p)]


def triples(g):
    """A triple is found exactly once, at its vertex with the least row."""
    ans = set()
    for a in range(len(g)):
        lines = defaultdict(list)
        for b in range(a+1, len(g)):
            dx, dy = b-a, g[b]-g[a]
            divisor = gcd(dx, dy)
            lines[dx//divisor, dy//divisor].append(b)
        for line in lines.values():
            for b, c in combinations(line, 2):
                ans.add((a, b, c))
    return ans


def move(g, u, v):
    h = g.copy()
    a, b = h.index(u), h.index(v)
    h[a], h[b] = h[b], h[a]
    return h


def collinear(g, tri):
    a, b, c = tri
    return (b-a)*(g[c]-g[a]) == (c-a)*(g[b]-g[a])


literal_moves = 0
for item in certificate["finite"]:
    p = item["p"]
    g = board(p)
    need([g[g[i]] for i in range(p)] == list(range(p)), "involutive initial board")
    old = triples(g)
    need(old == {tuple(t) for t in item["old"]}, "independent complete old triangle set")
    pack = []
    for tri in sorted(t for t in old if 0 in t):
        xx = set(tri) - {0}
        yy = {g[i] for i in xx}
        need(len(xx) == len(yy) == 2 and not xx & yy, "native packet disjoint shores")
        pack.append((xx, yy))
    need([[sorted(xx), sorted(yy)] for xx, yy in pack] == item["packets"], "complete packet certificate")
    need(len(old) == 2*len(pack), "ordered anchor pairing normalization")
    support = {t: {g[i] for i in t} for t in old}
    count_rows = {(u, v): count for u, v, count in item["move_counts"]}
    declared_covers = {(u, v): {tuple(t) for t in ts} for u, v, ts in item["covers"]}
    safe, covers = [], []
    for u, v in combinations(range(p), 2):
        missed = next((t for t, labels in support.items() if not {u, v} & labels), None)
        h = move(g, u, v)
        if missed is not None:
            need(all(h[i] == g[i] for i in missed) and collinear(h, missed), "literal unchanged old witness")
        else:
            covers.append((u, v))
        if (u, v) in count_rows:
            now = triples(h)
            literal_moves += 1
            need(len(now) == count_rows[u, v], "every declared move count by complete slope pencils")
            if missed is None:
                need(now == declared_covers[u, v], "complete created and surviving triangles at each cover")
            if not now:
                safe.append([u, v])
    need(set(covers) == set(declared_covers), "all-cover universe")
    need(safe == item["safe"], "exact successful second-move set")
    need([move(g, *uv) for uv in safe] == item["safe_permutations"], "complete safe permutations")
    if len(pack) >= 2:
        need(not safe and covers == [(0, 1)], "multiple-packet supplied controls")
    for u, v in covers:
        h = move(g, u, v)
        invlabels = (g[u], g[v])
        ht = move(g, *invlabels)
        need(all(ht[h[i]] == i for i in range(p)), "transpose is conjugated output move")
    print("FINITE", p, "packets", len(pack), "covers", len(covers), "safe", safe)

# Pure incidence referee for the only subtle universal case N=2. This
# universe deliberately includes a non-anchor cover, unlike the named prime
# controls, and verifies that it must be an inverse-pair diagonal trap.
g = {0: 1, 1: 0, 10: 10, **{i: i+1 if i % 2 == 0 else i-1 for i in range(2, 10)}}
possible = []
for xx in combinations(range(2, 10), 2):
    xx = frozenset(xx)
    yy = frozenset(g[i] for i in xx)
    if not xx & yy:
        possible.append((xx, yy))
abstract = cross_covers = 0
for (xx, yy), (aa, bb) in combinations(possible, 2):
    if xx & aa or yy & bb or len(xx & bb) > 1 or len(yy & aa) > 1:
        continue
    abstract += 1
    supports = [{1} | set(yy), {0} | set(xx), {1} | set(bb), {0} | set(aa)]
    for uv in combinations(range(11), 2):
        if not all(set(uv) & ss for ss in supports):
            continue
        need(uv == (0, 1) or (0 not in uv and 1 not in uv and g[uv[0]] == uv[1]),
             "all abstract two-packet covers are undo or inverse pairs")
        if uv != (0, 1):
            cross_covers += 1
            h = g.copy()
            a, b = [next(i for i in g if g[i] == value) for value in uv]
            h[a], h[b] = h[b], h[a]
            need(h[a] == a and h[b] == b and h[10] == 10, "universal newly created diagonal")
need(cross_covers > 0, "abstract controls exercise non-anchor covers")
print("INCIDENCE", abstract, "two-packet patterns;", cross_covers, "non-anchor inverse-pair covers")

# Rebuild all affine parameter witnesses with a symbolic determinant and
# Euclidean polynomial division, rather than the producer's coefficient engine.
r = S.symbols("r", integer=True, nonnegative=True)
p = 37+360*r
X1, Y1, X2, Y2 = 24+216*r, 17+160*r, 30+270*r, 21+200*r


def det(points):
    return S.expand(S.det(S.Matrix([[1, a, b] for a, b in points])))


def positive(poly):
    coeff = S.Poly(poly, r).all_coeffs()
    need(all(c >= 0 for c in coeff) and S.Poly(poly, r).eval(0) > 0,
         "uniform positive polynomial on nonnegative parameter ray")


def decode(poly):
    return sum(S.Rational(v)*r**i for i, v in enumerate(poly))


def reciprocal(point):
    a, b = point
    quotient, remainder = S.div(S.expand(a*b-1), p, r)
    need(remainder == 0 and all(c.q == 1 for c in S.Poly(quotient, r).all_coeffs()),
         "uniform native reciprocal quotient")
    positive(a-1); positive(b-1); positive(p-a); positive(p-b)
    return quotient


order = [0, 1, Y1, Y2, X1, X2, p-1]
for a, b in zip(order, order[1:]):
    positive(b-a)
need(det([(0, 1), (X1, Y1), (X2, Y2)]) == 0, "uniform old packet")
reciprocal((X1, Y1)); reciprocal((X2, Y2))
families = [((1, X1), ((Y1, 1), ((p+1)/2, 2), ((8*p+1)/9, 9))),
            ((1, X2), ((Y2, 1), ((2*p+1)/3, 3), ((5*p+1)/6, 6))),
            ((Y1, X2), ((X1, X2), ((3*p+4)/5, (3*p+5)/4), ((p-5)/4, (2*p-4)/5)))]
inverse = {S.sympify(0): S.Integer(1), S.sympify(1): S.Integer(0),
           X1: Y1, Y1: X1, X2: Y2, Y2: X2}
covered = set()
for (u, v), points in families:
    record = certificate["uniform_family"][len(covered)//2]
    need([S.expand(decode(a)) for a in record["labels"]] == list(map(S.sympify, (u, v))), "family labels")
    saved_points = [tuple(S.expand(decode(a)) for a in point) for point in record["points"]]
    need(saved_points == [tuple(map(S.expand, point)) for point in points], "family native points")
    need(det(points) == 0, "identically zero created integer determinant")
    for a, b in combinations([S.sympify(q[0]) for q in points], 2):
        roots = S.solve(S.Eq(a, b), r)
        need(not any(root.is_integer is True and root >= 0 for root in roots), "distinct native rows for every integer parameter")
    old_y = u if points[0][1] == v else v
    reciprocal((points[0][0], old_y))
    quotients = []
    for point in points[1:]:
        quotients.append(reciprocal(point))
        for label in (u, v):
            roots = S.solve(S.Eq(point[1], label), r)
            need(not any(root.is_integer is True and root >= 0 for root in roots), "other triangle points remain unchanged")
    need(quotients == [S.expand(decode(q)) for q in record["product_quotients"]], "independent exact product quotients")
    trans = (inverse[S.sympify(u)], inverse[S.sympify(v)])
    need([decode(a) for a in record["transpose_labels"]] == list(trans), "family conjugate labels")
    need([tuple(decode(a) for a in point) for point in record["transpose_points"]] == [(b, a) for a, b in points], "family transposed witnesses")
    covered.add(frozenset(map(S.sympify, (u, v)))); covered.add(frozenset(trans))
allcovers = {frozenset((a, b)) for a in [S.Integer(1), Y1, Y2] for b in [S.Integer(0), X1, X2]}
diagonal = {frozenset((S.Integer(0), S.Integer(1))), frozenset((X1, Y1)), frozenset((X2, Y2))}
need(len(allcovers) == 9 and covered | diagonal == allcovers and not covered & diagonal,
     "entire nine-cover family classified")
need(gcd(37, 360) == 1, "Dirichlet progression is coprime")
print("FAMILY all integer r>=0; all9 covers blocked, allowing extra native packets")
print("GEOMETRY", literal_moves, "moved boards reconstructed by complete slope pencils")
print("PASS", N, "always-active independent exact gates; actual LF")
