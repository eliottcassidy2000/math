#!/usr/bin/env python3
"""Independent exact audits of overlap legality, roots, lattices and resultants.

No repository mathematical implementation is imported. All decisions are exact.
Run normally and with python3 -O; the JSON receipts must agree.
"""
from fractions import Fraction
from itertools import product
from math import comb, gcd
import json
import sympy as sp

X = sp.Symbol("X")
checks = 0


def need(ok, message):
    global checks
    checks += 1
    if not ok:
        raise ArithmeticError(message)


def carry(word):
    out = 0
    for i, bit in enumerate(word):
        out = (3 if bit else 1) * out + (bit << i)
    return out


def source(p, r):
    q = p + r
    return [((i + 1) * p) // q - (i * p) // q for i in range(q)]


def labels(p, r):
    q = p + r
    return {(-j * r) % p: (j + (j * r + p - 1) // p - 1) % q
            for j in range(p)}


def edited(p, r, u, delta):
    w = source(p, r)
    q = len(w)
    pos = labels(p, r)
    residues = (u, u + r) if delta == 1 else (u + r, u)
    direction = 1 if delta == 1 else -1
    edges = []
    for residue in residues:
        left = pos[residue]
        right = (left + direction) % q
        need(w[left] == 1 and w[right] == 0, "illegal prescribed overlap")
        w[left], w[right] = w[right], w[left]
        pos[residue] = right
        edges.append((left, right))
    need(bool(set(edges[0]) & set(edges[1])), "prescribed overlap is disjoint")
    return w, edges


def egcd(a, b):
    if not b:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def root(p, r, g):
    _, a, b = egcd(p, r)
    return pow(2, a, g) * pow(3 * pow(2, -1, g) % g, b, g) % g


def independent(v, w):
    return any(v[i] * w[j] != v[j] * w[i]
               for i, j in ((0, 1), (0, 2), (1, 2)))


def normalized(v):
    b, c, m = v
    if c < 0 or (c == 0 and (b < 0 or (b == 0 and m < 0))):
        return tuple(-x for x in v)
    return v


def short_pair(p, r, u):
    first = None
    for h in range(1, p + 1):
        for b, c, m in product(range(-h, h + 1), repeat=3):
            v = (b, c, m)
            if max(map(abs, v)) != h or normalized(v) != v:
                continue
            if (m - b * r - c * u) % p:
                continue
            if first is None:
                first = v
            elif independent(first, v):
                return first, v
    raise ArithmeticError("no short independent pair")


def polynomial(p, r, u, delta, v):
    b, c, m = v
    need((m - b * r - c * u) % p == 0, "row not in lattice")
    a = (m - b * r - c * u) // p
    e = a - b + delta * c
    sign = -1 if delta == 1 and c % 2 else 1
    numerator = sign * 2 ** max(e, 0) * 3 ** max(b, 0) * 5 ** max(-c, 0)
    denominator = 2 ** max(-e, 0) * 3 ** max(-b, 0) * 5 ** max(c, 0)
    f = sp.Poly(denominator * X ** max(m, 0) * (X - 1) ** max(c, 0)
                - numerator * X ** max(-m, 0) * (X - 1) ** max(-c, 0), X)
    h = max(map(abs, v))
    need(abs(a) <= h, "A height bound")
    need(not f.is_zero, "zero polynomial")
    need(f.eval(0) != 0 and f.eval(1) != 0, "spurious root at zero or one")
    need((f.degree() == 0) == (m == c == 0), "constant degeneracy")
    need(f.degree() <= 2 * h, "degree bound")
    need(sum(abs(int(x)) for x in f.all_coeffs()) <= (38 if delta == 1 else 46) ** h,
         "improved coefficient bound")
    return f


# Independent actual-edge enumeration establishes the exact legality ranges.
legality_shapes = 0
legality_paths = 0
all_two_paths = 0
for p in range(2, 19):
    for r in range(1, p):
        if gcd(p, r) != 1:
            continue
        q = p + r
        bitword = source(p, r)
        pos = labels(p, r)
        tagged = [None] * q
        for a, loc in pos.items():
            tagged[loc] = a
        found = set()
        for edge1 in range(q):
            j1 = (edge1 + 1) % q
            if (tagged[edge1] is None) == (tagged[j1] is None):
                continue
            t = tagged[:]
            first_label = t[edge1] if t[edge1] is not None else t[j1]
            direction1 = 1 if t[edge1] is not None else -1
            t[edge1], t[j1] = t[j1], t[edge1]
            for edge2 in range(q):
                j2 = (edge2 + 1) % q
                if (t[edge2] is None) == (t[j2] is None):
                    continue
                all_two_paths += 1
                second_label = t[edge2] if t[edge2] is not None else t[j2]
                direction2 = 1 if t[edge2] is not None else -1
                if first_label == second_label:
                    continue
                overlap = bool({edge1, j1} & {edge2, j2})
                if not overlap:
                    continue
                need(direction1 == direction2, "mixed distinct-label overlap")
                u = min(first_label, second_label)
                need(max(first_label, second_label) == u + r, "wrong residue pair")
                delta = 1 if direction1 == 1 else 2
                need(first_label == (u if delta == 1 else u + r), "wrong leader")
                found.add((delta, u))
                legality_paths += 1
        expected = {(1, u) for u in range(min(r, p - r))}
        expected |= {(2, u) for u in range(max(0, p - 2 * r), p - r)}
        need(found == expected, "all-r legality mismatch")
        legality_shapes += 1

# Large exact actual-carry/root sweep. It deliberately includes both clock signs.
root_shapes = root_paths = nontrivial = cut_crossings = 0
negative_shapes = 0
nonendpoint_square_controls = []
for p in range(2, 121):
    for r in range(1, p):
        if gcd(p, r) != 1:
            continue
        root_shapes += 1
        D = 2 ** (p + r) - 3 ** p
        negative_shapes += int(D < 0)
        for delta, low, high in ((1, 0, min(r, p - r)),
                                 (2, max(0, p - 2 * r), p - r)):
            for u in range(low, high):
                word, edges = edited(p, r, u, delta)
                g = gcd(abs(D), carry(word))
                root_paths += 1
                cut_crossings += int(any(abs(a - b) > 1 for a, b in edges))
                need(gcd(g, 30) == 1, "bad prime in content")
                if g > 1:
                    nontrivial += 1
                    z = root(p, r, g)
                    need(pow(z, p, g) == 2 % g and 2 * pow(z, r, g) % g == 3 % g,
                         "clock root")
                    need(gcd(z * (z - 1), g) == 1, "nonunit root coordinate")
                    expected = -2 if delta == 1 else 4
                    need((5 * pow(z, u, g) * (z - 1) - expected) % g == 0,
                         "overlap root constant")
                    endpoint = u == (p - 1 - r if delta == 1 else 0)
                    if not endpoint and len(nonendpoint_square_controls) < 8:
                        square_primes = [ell for ell in (7, 11, 13, 17, 19, 23, 29, 31)
                                         if g % (ell * ell) == 0]
                        if square_primes:
                            nonendpoint_square_controls.append(
                                dict(p=p, r=r, u=u, delta=delta, content=g,
                                     squared_primes=square_primes))

# Every u in the whole algebraic range, not only legal physical overlaps.
lattice_instances = resultant_tests = dependent_projections = constant_cases = 0
largest_minima_product = {"product": 0}
for p in range(2, 37):
    for r in range(1, p):
        if gcd(p, r) != 1:
            continue
        for u in range(p - r):
            pair = short_pair(p, r, u)
            heights = [max(map(abs, v)) for v in pair]
            prod = heights[0] * heights[1]
            need(prod ** 3 <= p ** 2, "Minkowski product bound")
            need(independent(*pair), "dependent lattice rows")
            lattice_instances += 1
            dependent_projections += int(pair[0][2] * pair[1][1] ==
                                         pair[1][2] * pair[0][1])
            if prod > largest_minima_product["product"]:
                largest_minima_product = dict(p=p, r=r, u=u, pair=pair,
                                               heights=heights, product=prod)
            for delta in (1, 2):
                f, h = [polynomial(p, r, u, delta, v) for v in pair]
                need(sp.gcd(f, h).degree() == 0, "common Q-root")
                R = int(sp.resultant(f, h))
                need(R != 0, "vanishing resultant")
                need(abs(R) <= (38 if delta == 1 else 46) ** (4 * prod),
                     "resultant size bound")
                legal = (u < r if delta == 1 else u >= p - 2 * r)
                if legal:
                    word, _ = edited(p, r, u, delta)
                    g = gcd(abs(2 ** (p + r) - 3 ** p), carry(word))
                    need(R % g == 0, "content does not divide resultant")
                    if g > 1:
                        z = root(p, r, g)
                        need(f.eval(z) % g == 0 and h.eval(z) % g == 0,
                             "cleared relations fail modulo full content")
                constant_cases += int(f.degree() == 0 or h.degree() == 0)
                resultant_tests += 1

# Force the dependent-projection and constant cases, including negative rows.
hostile_controls = []
for p, r, u in ((2, 1, 0), (5, 1, 3), (11, 7, 2), (13, 5, 4)):
    for delta in (1, 2):
        for tag, pair in (
                ("dependent_projection", ((1, 0, r), (0, 0, p))),
                ("constant", ((p, 0, 0), (0, 1, u))),
                ("negative_constant", ((-p, 0, 0), (0, -1, -u)))):
            need(independent(*pair), "hostile rows dependent")
            f, h = [polynomial(p, r, u, delta, v) for v in pair]
            need(sp.gcd(f, h).degree() == 0, "hostile common root")
            need(sp.resultant(f, h) != 0, "hostile zero resultant")
            hostile_controls.append(dict(p=p, r=r, u=u, delta=delta,
                                         tag=tag, degrees=[f.degree(), h.degree()]))

# Outside-main-universe prime-power controls from the exact endpoint packet.
primepower_controls = []
for p, r, u, delta, expected in ((418, 237, 180, 1, 11849),
                                (203, 165, 0, 2, 169)):
    word, _ = edited(p, r, u, delta)
    g = gcd(abs(2 ** (p + r) - 3 ** p), carry(word))
    need(g == expected, "known prime-power content")
    pair = short_pair(p, r, u)
    f, h = [polynomial(p, r, u, delta, v) for v in pair]
    R = int(sp.resultant(f, h))
    need(R != 0 and R % g == 0, "prime-power resultant depth")
    z = root(p, r, g)
    need(f.eval(z) % g == 0 and h.eval(z) % g == 0, "prime-power root depth")
    primepower_controls.append(dict(p=p, r=r, u=u, delta=delta, content=g,
                                   pair=pair, resultant=R))

# A rational series certificate for the final uniform exponent.
need(sum(Fraction(4 ** j, sp.factorial(j)) for j in range(7)) > 46,
     "exp(4)>46 certificate")

print(json.dumps(dict(
    status="PASS_EXACT",
    checks=checks,
    legality_universe="2<=p<=18, 1<=r<p, gcd(p,r)=1, every ordered actual two-edge path",
    legality_shapes=legality_shapes,
    all_two_edge_paths=all_two_paths,
    distinct_one_overlapping_paths=legality_paths,
    root_universe="2<=p<=120, 1<=r<p, gcd(p,r)=1, every legal same-sign overlap",
    root_shapes=root_shapes,
    root_paths=root_paths,
    negative_clock_shapes=negative_shapes,
    cut_crossing_paths=cut_crossings,
    nontrivial_content_paths=nontrivial,
    nonendpoint_square_controls=nonendpoint_square_controls,
    lattice_universe="2<=p<=36, 1<=r<p, gcd(p,r)=1, every 0<=u<p-r, both constants",
    lattice_instances=lattice_instances,
    resultant_tests=resultant_tests,
    selected_dependent_projection_pairs=dependent_projections,
    selected_constant_cases=constant_cases,
    largest_selected_minima_product=largest_minima_product,
    hostile_controls=hostile_controls,
    primepower_controls=primepower_controls,
), indent=2))
