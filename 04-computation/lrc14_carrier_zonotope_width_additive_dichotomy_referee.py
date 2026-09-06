#!/usr/bin/env python3
"""Exact checks for the carrier zonotope area, width, and layer dictionary.

The analytic identities are proved in REPORT.md.  This standalone program
audits their integer algebra, the sharp layer controls, and the half-body
sumset map without importing repository mathematics.
"""

from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0


class AuditFailure(RuntimeError):
    pass


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AuditFailure(payload)


def cross(x, y):
    return (
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    )


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def neg(x):
    return tuple(-a for a in x)


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def strict_floor(num, den):
    return (num - 1) // den


def carriers(w):
    box = tuple(strict_floor(3 * (sum(w) - w[i]), 14) for i in range(3))
    out = set()
    for x in range(-box[0], box[0] + 1):
        for y in range(-box[1], box[1] + 1):
            nz = -w[0] * x - w[1] * y
            if nz % w[2]:
                continue
            z = nz // w[2]
            C = (x, y, z)
            if abs(z) <= box[2] and all(c % 3 for c in C):
                out.add(C)
    return out


def primitive(v):
    d = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
    q = tuple(x // d for x in v)
    if next(x for x in q if x) < 0:
        q = neg(q)
    return q


def relation_minimum(w):
    seeds = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        d = gcd(w[i], w[j])
        v = [0, 0, 0]
        v[i], v[j] = w[j] // d, -w[i] // d
        seeds.append((sum(abs(x) for x in v), tuple(v)))
    bound, best = min(seeds)
    initial = bound
    for x in range(-initial, initial + 1):
        for y in range(-initial, initial + 1):
            nz = -w[0] * x - w[1] * y
            if nz % w[2]:
                continue
            z = nz // w[2]
            norm = abs(x) + abs(y) + abs(z)
            if norm and norm < bound:
                bound, best = norm, (x, y, z)
    return bound, primitive(best)


def egcd(a, b):
    if b == 0:
        return (abs(a), 1 if a > 0 else -1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)


def cross_preimage(w, n):
    """Find integral z with w cross z=n."""
    d, p, q = egcd(w[1], w[2])
    need(n[0] % d == 0, ("cross-divisibility", w, n, d))
    z3 = p * (n[0] // d)
    z2 = -q * (n[0] // d)
    for t in range(d):
        zz3 = z3 + (w[2] // d) * t
        numerator = n[1] + w[0] * zz3
        if numerator % w[2]:
            continue
        z = (numerator // w[2], z2 + (w[1] // d) * t, zz3)
        need(cross(w, z) == n, ("cross-preimage", w, n, z, cross(w, z)))
        return z
    raise AuditFailure(("no-cross-preimage", w, n))


def affine_rank(points):
    points = tuple(points)
    if len(points) <= 1:
        return 0
    base = points[0]
    ds = [tuple(x[i] - base[i] for i in range(3)) for x in points[1:]]
    first = next((x for x in ds if x != (0, 0, 0)), None)
    if first is None:
        return 0
    return 2 if any(cross(first, x) != (0, 0, 0) for x in ds) else 1


def area_and_residue_audit():
    residues = 0
    for wr in product((1, 2), repeat=3):
        live = {
            C for C in product(range(3), repeat=3)
            if all(C) and dot(wr, C) % 3 == 0
        }
        expected = {tuple(x % 3 for x in wr), tuple((-x) % 3 for x in wr)}
        need(live == expected, ("live-cosets", wr, live, expected))
        residues += len(live)

    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    triples = 0
    values = tuple(x for x in range(1, 32, 2) if x % 3)
    for w in combinations(values, 3):
        triples += 1
        generators = tuple(neg(cross(w, e)) for e in basis)
        coefficients = []
        for i, j in combinations(range(3), 2):
            k = 3 - i - j
            determinant = cross(generators[i], generators[j])
            expected = tuple(w[k] * x for x in w)
            need(determinant == expected or determinant == neg(expected),
                 ("zonotope-minor", w, i, j, determinant, expected))
            coefficients.append(w[k])
        need(sum(coefficients) == sum(w), ("area-coefficients", w, coefficients))
    return triples, residues


def layer_control(w, chosen, expected_mu, expected_layers, expected_directions):
    mu, minimizer = relation_minimum(w)
    chosen = primitive(chosen)
    need(mu == expected_mu and dot(w, chosen) == 0,
         ("minimum", w, mu, minimizer, chosen))
    need(sum(abs(x) for x in chosen) == mu, ("chosen-norm", w, chosen, mu))
    z = cross_preimage(w, chosen)
    rows = carriers(w)
    layers = {dot(z, C) for C in rows}
    need(layers == set(expected_layers), ("layers", w, layers, expected_layers))
    need(all(Q(abs(dot(z, C))) < Q(3 * mu, 14) for C in rows),
         ("width-support", w, mu, layers))
    unit = all(x % 3 for x in chosen)
    expected_residues = {0} if unit else {1, 2}
    need({t % 3 for t in layers} == expected_residues,
         ("layer-residue", w, chosen, layers, expected_residues))
    ds = {primitive(C) for C in rows}
    need(len(ds) == expected_directions,
         ("direction-count", w, len(ds), expected_directions, ds))
    return mu, unit, len(rows), tuple(sorted(layers)), len(ds), chosen, z


def half_body_audit(height=79):
    values = tuple(x for x in range(1, height + 1, 2) if x % 3)
    triples = nonempty = rank2_core = 0
    max_core = (Q(0), None)
    for w in combinations(values, 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        triples += 1
        rows = carriers(w)
        if not rows:
            continue
        nonempty += 1
        positive_class = {C for C in rows if tuple(x % 3 for x in C) == tuple(x % 3 for x in w)}
        need(len(positive_class) * 2 == len(rows), ("class-balance", w, len(rows), positive_class))
        half = {
            C for C in positive_class
            if all(28 * abs(C[i]) < 3 * (sum(w) - w[i]) for i in range(3))
        }
        sums = {add(x, y) for x in half for y in half}
        need(sums <= rows, ("half-sum-map", w, sums - rows))
        h, m = len(half), len(positive_class)
        if h:
            need(len(sums) >= 2 * h - 1, ("torsionfree-sum-bound", w, h, len(sums)))
            need(2 * h - 1 <= m, ("half-core-general", w, h, m))
            if affine_rank(half) == 2:
                rank2_core += 1
                need(len(sums) >= 3 * h - 3, ("freiman-d2-control", w, h, len(sums)))
                need(3 * h - 3 <= m, ("half-core-rank2", w, h, m))
            ratio = Q(h, m)
            if ratio > max_core[0]:
                max_core = (ratio, w, h, m, affine_rank(half))
    return triples, nonempty, rank2_core, max_core


def dense_bulk_audit(height=79):
    values = tuple(x for x in range(1, height + 1, 2) if x % 3)
    dense = 0
    smallest_excess = None
    for w in combinations(values, 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        N = len(carriers(w))
        if 11 * N <= 2 * w[2]:
            continue
        dense += 1
        bulk = Q(2 * sum(w), 49)
        excess = Q(N) - bulk
        forced = Q(32 * w[2] + 132, 539)
        need(excess > forced, ("linear-excess", w, N, excess, forced))
        row = (excess, w, N, bulk, forced)
        if smallest_excess is None or row < smallest_excess:
            smallest_excess = row
    need(dense == 114, ("dense-H79", dense))
    return dense, smallest_excess


def main():
    area_triples, residue_points = area_and_residue_audit()
    controls = (
        ((1, 5, 11), (1, 2, -1), 4, (0,), 1),
        ((1, 19, 23), (4, 1, -1), 6, (0,), 1),
        ((1, 5, 61), (5, -1, 0), 6, (-1, 1), 2),
        ((17, 23, 25), (1, -4, 3), 8, (-1, 1), 2),
        ((19, 23, 29), (3, -5, 2), 10, (-2, -1, 1, 2), 3),
        ((1, 17, 131), (5, -8, 1), 14, (0,), 1),
        ((71, 95, 97), (8, -7, 1), 16, (-3, 0, 3), 5),
    )
    layer_rows = tuple(layer_control(*row) for row in controls)

    a2 = carriers((19, 23, 29))
    u, v = (1, 8, -7), (10, -7, -1)
    need(a2 == {u, neg(u), v, neg(v), add(u, v), neg(add(u, v))},
         ("A2-circuit", a2))

    # The first rank-two half-body core is a two-shell A2 saturation.
    shell_w = (85, 97, 107)
    shell_rows = carriers(shell_w)
    shell_A = {C for C in shell_rows if tuple(x % 3 for x in C) == tuple(x % 3 for x in shell_w)}
    shell_P = {
        C for C in shell_A
        if all(28 * abs(C[i]) < 3 * (sum(shell_w) - shell_w[i]) for i in range(3))
    }
    need(len(shell_P) == 3 and affine_rank(shell_P) == 2,
         ("two-shell-core", shell_P))
    need(tuple(sum(C[i] for C in shell_P) for i in range(3)) == (0, 0, 0),
         ("two-shell-zero-sum", shell_P))
    shell_sumset = {add(x, y) for x in shell_P for y in shell_P}
    need(shell_sumset == {neg(x) for x in shell_A},
         ("two-shell-saturation", shell_sumset, shell_A))
    need(shell_A == shell_P | {tuple(-2 * x for x in C) for C in shell_P},
         ("two-shell-dilation", shell_A, shell_P))

    half = half_body_audit()
    dense, smallest = dense_bulk_audit()
    print("LRC CARRIER ZONOTOPE WIDTH + ADDITIVE DICHOTOMY REFEREE")
    print("status=PASS proved-identities+finite-controls; universal_projection=OPEN; LRC14=OPEN")
    print("support=Z_w=-w_cross[-3/14,3/14]^3")
    print("area=9*norm(w)*(a+b+c)/49 covol(Lambda)=norm(w) covol(3Lambda)=9*norm(w)")
    print("two_live_cosets_bulk=2*(a+b+c)/49 area_triples=%d residue_points=%d" %
          (area_triples, residue_points))
    print("width_Lambda=3*mu1/7 width_3Lambda=mu1/7 layers=abs(t)<3*mu1/14")
    for row, data in zip(controls, layer_rows):
        print("layer_control w=%s mu=%d type=%s N=%d layers=%s directions=%d n=%s z=%s" %
              (row[0], data[0], "unit" if data[1] else "one-zero", data[2],
               data[3], data[4], data[5], data[6]))
    print("A2_exact=w:(19,23,29),u:%s,v:%s,u+v:%s" % (u, v, add(u, v)))
    print("two_shell_A2=w:%s,N:%d,inner_triangle:%s,A=P_union_minus2P,P_plus_P=minusA" %
          (shell_w, len(shell_rows), tuple(sorted(shell_P))))
    print("half_body_H79 triples=%d nonempty=%d rank2_inner_cores=%d max_core_ratio=%s row=%s h=%d m=%d rank=%d" %
          (half[0], half[1], half[2], half[3][0], half[3][1], half[3][2], half[3][3], half[3][4]))
    print("dense_bulk_H79 rows=%d smallest_excess=%s" % (dense, smallest))
    print("checks=%d verdict=PASS" % CHECKS)


if __name__ == "__main__":
    main()
