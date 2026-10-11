"""Exact omitted-unit-edge tests in the Moser field.

The all-size field obstruction is inherited mathematics, implemented with two
F4 residue frames.  No seven-chromatic graph or six-equality gadget is claimed.
"""
from collections import deque
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json

import finite_algebraic_witness_20261007 as proof


CHECKS = 0
ZERO = (F(0),) * 4
ONE = (F(1), F(0), F(0), F(0))
A = (F(0), F(1), F(0), F(0))
B = (F(0), F(0), F(1), F(0))
AB = (F(0), F(0), F(0), F(1))


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def element(z):
    need(type(z) is tuple and len(z) == 4
         and all(type(c) in (int, F) for c in z), "four exact rational coefficients")
    return tuple(F(c) for c in z)


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def scale(c, z):
    return tuple(c * x for x in z)


def sub(x, y):
    return add(x, scale(-1, y))


def times(x, y):
    """Basis 1,a,b,ab, with a²=a−1 and b²=(5/3)b−1."""
    out = [F(0)] * 4
    for i, u in enumerate(x):
        if not u:
            continue
        for j, v in enumerate(y):
            if not v:
                continue
            ae, be = (i & 1) + (j & 1), (i >> 1) + (j >> 1)
            aa = ((0, -1), (1, 1)) if ae == 2 else ((ae, 1),)
            bb = ((0, -1), (1, F(5, 3))) if be == 2 else ((be, 1),)
            for ap, ac in aa:
                for bp, bc in bb:
                    out[ap + 2 * bp] += u * v * ac * bc
    return tuple(out)


def conjugate(x):
    u, v, w, z = x
    return (u + v + F(5, 3) * (w + z), -v - F(5, 3) * z, -w-z, z)


def norm(x):
    return times(x, conjugate(x))


def power(x, n):
    need(type(n) is int and n >= 0, "natural exponent")
    out = ONE
    for _ in range(n):
        out = times(out, x)
    return out


def binary_order(q):
    need(type(q) is int and q > 0, "positive exact integer")
    return (q & -q).bit_length() - 1


def denominator_depth(z):
    z = element(z)
    return max(binary_order(c.denominator) for c in z)


def f4times(x, y):
    need(type(x) is int and type(y) is int and 0 <= x < 4 and 0 <= y < 4,
         "F4 labels are exact integers0..3")
    raw = 0
    for i in range(2):
        for j in range(2):
            if (x >> i & 1) and (y >> j & 1):
                raw ^= 1 << (i+j)
    if raw & 4:
        raw ^= 7
    return raw


def residue(z, frame=0):
    z = element(z)
    need(type(frame) is int and frame in (0, 1), "exact residue-frame bit")
    need(denominator_depth(z) == 0, "reduce a source offset with odd denominators")
    beta = 2 if frame == 0 else 3
    basis = (1, 2, beta, f4times(2, beta))
    out = 0
    for c, v in zip(z, basis):
        if c.numerator & 1:
            out ^= v
    return out


def unit_denominator_certificate(z):
    """Authenticate the elementary leading-denominator obstruction at this z."""
    z = element(z)
    need(norm(z) == ONE, "actual unit vector required")
    d = denominator_depth(z)
    w = scale(1 << d, z)
    if d:
        # This branch is mathematically impossible; it states the exact gates.
        need(any(residue(w, f) for f in (0, 1)), "minimal denominator keeps parity")
        need(all(f4times(residue(w, f), residue(conjugate(w), f)) == 0
                 for f in (0, 1)), "cleared unit norm has zero residue")
        raise ValueError("inconsistent unit vector")
    need(all(residue(z, f) != 0 for f in (0, 1)), "unit vectors have nonzero residues")
    return d


def validate_points(points):
    need(type(points) is tuple, "source points must be a tuple")
    points = tuple(element(p) for p in points)
    need(len(set(points)) == len(points), "coincident vertices rejected")
    return points


def unit_edges(points):
    points = validate_points(points)
    return _unit_edges(points)


@lru_cache(None)
def _unit_edges(points):
    return frozenset((i, j) for i, j in combinations(range(len(points)), 2)
                     if norm(sub(points[i], points[j])) == ONE)


def graph_edges(n, edges):
    need(type(n) is int and n >= 0, "exact vertex count")
    need(type(edges) is frozenset, "canonical edge set")
    for e in edges:
        need(type(e) is tuple and len(e) == 2
             and all(type(i) is int for i in e)
             and 0 <= e[0] < e[1] < n, "exact ordered edge endpoints")
    return edges


def component_origins(n, edges):
    graph_edges(n, edges)
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    origins = [-1] * n
    for root in range(n):
        if origins[root] != -1:
            continue
        origins[root] = root
        todo = [root]
        while todo:
            u = todo.pop()
            for v in adj[u]:
                if origins[v] == -1:
                    origins[v] = root
                    todo.append(v)
    return tuple(origins)


@dataclass(frozen=True)
class FieldColour:
    points: tuple
    frame: int
    origins: tuple
    colours: tuple


def field_colour(points, frame=0):
    points = validate_points(points)
    need(type(frame) is int and frame in (0, 1), "exact residue-frame bit")
    edges = unit_edges(points)
    origins = component_origins(len(points), edges)
    colours = tuple(residue(sub(p, points[o]), frame) for p, o in zip(points, origins))
    result = FieldColour(points, frame, origins, colours)
    verify_field_colour(result)
    return result


def verify_field_colour(packet):
    need(type(packet) is FieldColour, "FieldColour packet")
    points = validate_points(packet.points)
    need(type(packet.frame) is int and packet.frame in (0, 1), "exact frame")
    need(type(packet.origins) is tuple and len(packet.origins) == len(points)
         and all(type(i) is int for i in packet.origins), "exact component markers")
    edges = unit_edges(points)
    expected = component_origins(len(points), edges)
    need(packet.origins == expected, "retain the full unit-graph components")
    need(type(packet.colours) is tuple and len(packet.colours) == len(points)
         and all(type(c) is int and 0 <= c < 4 for c in packet.colours), "four-colour tuple")
    need(packet.colours == tuple(residue(sub(p, points[o]), packet.frame)
                                 for p, o in zip(points, expected)), "source-owned colour evaluation")
    proof.verify_colouring(len(points), edges, 4, packet.colours)
    return edges


@dataclass(frozen=True)
class Fragment:
    points: tuple
    edges: frozenset
    ports: tuple
    palette: int


def verify_fragment(fragment):
    need(type(fragment) is Fragment, "Fragment packet")
    points = validate_points(fragment.points)
    graph_edges(len(points), fragment.edges)
    need(type(fragment.ports) is tuple and len(fragment.ports) == 2
         and all(type(v) is int for v in fragment.ports)
         and 0 <= fragment.ports[0] < fragment.ports[1] < len(points), "exact distinct ports")
    need(type(fragment.palette) is int and fragment.palette >= 1, "exact palette")
    actual = unit_edges(points)
    need(fragment.edges <= actual, "retained edge must be an actual unit edge")
    need(fragment.ports in actual and fragment.ports not in fragment.edges,
         "omitted pair must still be exactly one unit apart")
    return actual


def separate(fragment, packet):
    """Consume the same source's full-unit-graph colouring, padded to k>=4."""
    verify_fragment(fragment)
    verify_field_colour(packet)
    need(fragment.palette >= 4, "this direct separator uses four colours")
    need(fragment.points == packet.points, "exact source coordinates must match")
    proof.verify_colouring(len(fragment.points), fragment.edges, fragment.palette, packet.colours)
    u, v = fragment.ports
    need(packet.colours[u] != packet.colours[v], "ports must separate")
    return packet.colours


@dataclass(frozen=True)
class Equality:
    fragment: Fragment
    colouring: tuple
    refutation: tuple


def verify_equality(packet):
    need(type(packet) is Equality, "Equality packet")
    f = packet.fragment
    verify_fragment(f)
    proof.verify_colouring(len(f.points), f.edges, f.palette, packet.colouring)
    nodes = proof.verify_refutation(len(f.points), f.edges | {f.ports},
                                    f.palette, packet.refutation)
    return nodes


def unique(points):
    # Candidate generation may identify coincident algebraic construction paths;
    # the resulting graph vertices are the distinct exact points.
    return tuple(sorted(set(points)))


def spindle():
    return (ZERO, ONE, A, add(ONE, A), B, AB, add(B, AB))


def triangular(radius):
    need(type(radius) is int and radius >= 0, "natural shell radius")
    return unique(add(scale(i, ONE), scale(j, A))
                  for i in range(-radius, radius+1) for j in range(-radius, radius+1)
                  if max(abs(i), abs(j), abs(i+j)) <= radius)


def minkowski(left, right):
    return unique(add(x, y) for x in left for y in right)


def candidates():
    wheel = triangular(1)
    rotated = tuple(times(B, p) for p in wheel)
    triangle = (ZERO, B, AB)
    s = spindle()
    return (
        ("triangular_radius1", wheel),
        ("triangular_radius2", triangular(2)),
        ("moser_spindle", s),
        ("two_triangular_shells", unique(wheel + rotated)),
        ("wheel_plus_rotated_triangle", minkowski(wheel, triangle)),
        ("wheel_plus_rotated_wheel", minkowski(wheel, rotated)),
        ("six_rotated_spindles", unique(times(power(B, j), p) for j in range(6) for p in s)),
        ("coefficient_box_minus1_to1", tuple(element(c) for c in product((-1,0,1), repeat=4))),
        ("half_translated_spindle", tuple(add(scale(F(1,2), ONE), p) for p in s)),
        ("halved_spindle", tuple(scale(F(1,2), p) for p in s)),
    )


def rejects(f):
    try:
        f()
    except (ValueError, TypeError):
        return True
    return False


def main():
    check(times(A,A) == sub(A,ONE), "first exact defining relation")
    check(times(B,B) == sub(scale(F(5,3),B),ONE), "second exact defining relation")
    check(norm(A) == norm(B) == ONE, "both generators have actual unit norm")
    parities = tuple(element(c) for c in product(range(2),repeat=4))
    paired = {(residue(z,0),residue(z,1)) for z in parities}
    check(paired == set(product(range(4),repeat=2)), "R/2R is F4 times F4")
    for x in parities:
        for f in (0,1):
            check(residue(conjugate(x),f) == f4times(residue(x,f),residue(x,f)),
                  "conjugation is Frobenius in each frame")
        for y in parities:
            for f in (0,1):
                check(residue(add(x,y),f) == residue(x,f)^residue(y,f), "additive reduction")
                check(residue(times(x,y),f) == f4times(residue(x,f),residue(y,f)),
                      "multiplicative reduction")
    # The leading 2-denominator argument needs both frames, not just one.
    kernel = sub(A,B)
    check(residue(kernel,0)==0 and residue(kernel,1)!=0 and denominator_depth(kernel)==0,
          "one residue frame loses a nonzero leading parity")
    check(rejects(lambda: residue(scale(F(1,2),kernel),0)),
          "vanishing in one frame does not authorize division by two")
    # Exact rational-parameter norm-one rotations, including even denominators
    # in the raw formula whose factors must cancel before canonical reduction.
    rotations=0
    for m in range(-6,7):
        for n in range(-6,7):
            if m == n == 0:
                continue
            den = m*m + 11*n*n
            # (m+n sqrt(-11))/(m-n sqrt(-11)); sqrt(-11)=6b-5.
            z = add(scale(F(m*m-11*n*n-10*m*n,den),ONE),
                    scale(F(12*m*n,den),B))
            check(norm(z)==ONE and unit_denominator_certificate(z)==0, "unit denominator lemma")
            rotations += 1

    table=[]
    total_edges=0
    for name, pts in candidates():
        edges=unit_edges(pts)
        colours=[field_colour(pts,f) for f in (0,1)]
        for c in colours:
            check(verify_field_colour(c)==edges, "full exact unit incidence authenticated")
        for e in sorted(edges):
            f=Fragment(pts, edges-{e}, e, 6)
            for c in colours:
                result=separate(f,c)
                check(result[e[0]]!=result[e[1]], "actual six-colour separation for every omission")
        table.append({"name":name,"vertices":len(pts),"unit_edges":len(edges),
                      "components":len(set(colours[0].origins)),
                      "palette_used":[len(set(c.colours)) for c in colours],
                      "six_forced_equal_omissions":0})
        total_edges += len(edges)

    # Nonvacuous lower-palette positives: all eleven spindle edge deletions.
    pts=spindle(); edges=unit_edges(pts)
    check(len(pts)==7 and len(edges)==11, "exact spindle control")
    no_colour, refutation=proof.colour_or_refute(7,edges,3)
    check(no_colour is None, "spindle is not three-colourable")
    ref_nodes=proof.verify_refutation(7,edges,3,refutation)
    counts=[]
    for e in sorted(edges):
        frag=Fragment(pts,edges-{e},e,3)
        col, _=proof.colour_or_refute(7,frag.edges,3)
        check(col is not None, "fragment nonvacuity")
        check(verify_equality(Equality(frag,col,refutation))==ref_nodes, "forced-equality certificate")
        solutions=equal=0
        for c in product(range(3),repeat=7):
            if all(c[u]!=c[v] for u,v in frag.edges):
                solutions += 1
                equal += c[e[0]]==c[e[1]]
        check(solutions>0 and solutions==equal, "independent literal all-colourings check")
        counts.append((e,solutions))
        for k in (4,6):
            check(separate(replace(frag,palette=k),field_colour(pts))[e[0]]
                  != separate(replace(frag,palette=k),field_colour(pts))[e[1]],
                  "forcing does not survive spare palette")

    # Coordinate translations and fragment replacement cannot be silently erased.
    translated=tuple(add(scale(F(1,2),ONE),p) for p in pts)
    check(rejects(lambda: residue(translated[0])), "absolute reduction hostile")
    tc=field_colour(translated)
    check(tc.colours==field_colour(pts).colours, "component-origin translation repair")
    e=min(edges); frag=Fragment(pts,edges-{e},e,6)
    check(rejects(lambda: separate(frag,tc)), "isomorphic translated source is not the same receipt")
    # A nonunit omitted pair is not a geometric equality target.
    nonunit=next(e for e in combinations(range(7),2) if e not in edges)
    bads=(
        lambda: element((True,0,0,0)),
        lambda: element((0.0,0,0,0)),
        lambda: unit_edges((ZERO,ZERO)),
        lambda: field_colour(pts,True),
        lambda: verify_field_colour(replace(field_colour(pts),colours=(False,)+field_colour(pts).colours[1:])),
        lambda: verify_field_colour(replace(field_colour(pts),origins=(False,)*7)),
        lambda: verify_fragment(replace(frag,ports=(False,1))),
        lambda: verify_fragment(replace(frag,ports=nonunit)),
        lambda: verify_fragment(replace(frag,edges=edges)),
        lambda: verify_fragment(replace(frag,palette=True)),
        lambda: verify_fragment(replace(frag,edges=frozenset({(False,1)}))),
    )
    for f in bads:
        check(rejects(f), "malformed receipt rejected")
    print(json.dumps({"scope":"exact candidates inside the common rational Moser field",
                      "checks":CHECKS,"leading_parity_states":len(parities),
                      "unit_rotation_controls":rotations,
                      "candidate_sets":table,"all_unit_edge_omissions":total_edges,
                      "spindle_three_colour_equalities":counts,
                      "shared_refutation_nodes":ref_nodes,
                      "six_equality_found":False,
                      "stopping_reason":"explicit full-unit-graph four-colour certificate",
                      "plane_six_or_seven":"not decided"},indent=2))


if __name__=="__main__":
    main()
