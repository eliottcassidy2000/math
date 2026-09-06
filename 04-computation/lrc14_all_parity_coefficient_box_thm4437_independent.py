#!/usr/bin/env python3
"""Clean-room exact audit of the parity-free short-relation slope box.

No project module is imported.  Cube sections are rebuilt from intersections
of the cutting plane with all twelve cube edges, while speed-domain vertices
are rebuilt from the plane v.w=0 and the unit cube.
"""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement, permutations, product
from math import gcd


R = Q(3, 14)
TARGET = Q(15, 98)
LOW = {(1, 1, 1), (1, 1, 2), (1, 2, 2)}
CHECKS = 0


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(payload)


def gcd3(p):
    return gcd(gcd(p[0], p[1]), p[2])


def eligible_pattern(p):
    return (
        sum(x != 0 for x in p) >= 2
        and gcd3(p) == 1
        and sum(x % 3 == 0 for x in p) <= 1
        and p != (0, 1, 1)
    )


def signed_vectors(p):
    """Distinct signed permutations modulo simultaneous sign reversal."""
    answer = set()
    for u in set(permutations(p)):
        nz = next(i for i, x in enumerate(u) if x)
        for signs in product((-1, 1), repeat=3):
            v = tuple(u[i] * signs[i] for i in range(3))
            if v[nz] <= 0:
                continue
            if min(v) >= 0 or max(v) <= 0:
                continue
            answer.add(v)
    return sorted(answer)


def speed_vertices(v):
    """Vertices of {w in [0,1]^3:v.w=0}, excluding the origin."""
    points = set()
    for free in range(3):
        fixed = [j for j in range(3) if j != free]
        for bits in product((Q(0), Q(1)), repeat=2):
            rhs = -(v[fixed[0]] * bits[0] + v[fixed[1]] * bits[1])
            if v[free]:
                x = rhs / v[free]
                if Q(0) <= x <= Q(1):
                    w = [Q(0), Q(0), Q(0)]
                    w[fixed[0]], w[fixed[1]], w[free] = bits[0], bits[1], x
                    if any(w):
                        points.add(tuple(w))
            elif rhs == 0:
                for x in (Q(0), Q(1)):
                    w = [Q(0), Q(0), Q(0)]
                    w[fixed[0]], w[fixed[1]], w[free] = bits[0], bits[1], x
                    if any(w):
                        points.add(tuple(w))
    # Every nonzero vertex has a unit coordinate.  This independently checks
    # that division by max(w) has already been normalized.
    need(all(max(w) == 1 and sum(v[i] * w[i] for i in range(3)) == 0 for w in points),
         ("bad speed vertex", v, points))
    return sorted(points)


def cube_section_width(v, w, delta):
    """Width of e -> (w cross e)_i/v_i on v.e=delta, e in [-R,R]^3."""
    axis = max(range(3), key=lambda i: abs(v[i]))
    need(v[axis] != 0, ("zero chart", v))
    vertices = set()
    for free in range(3):
        fixed = [j for j in range(3) if j != free]
        for signs in product((-1, 1), repeat=2):
            e = [Q(0), Q(0), Q(0)]
            e[fixed[0]] = signs[0] * R
            e[fixed[1]] = signs[1] * R
            residual = delta - v[fixed[0]] * e[fixed[0]] - v[fixed[1]] * e[fixed[1]]
            if v[free]:
                e[free] = residual / v[free]
                if -R <= e[free] <= R:
                    vertices.add(tuple(e))
            elif residual == 0:
                for endpoint in (-R, R):
                    e[free] = endpoint
                    vertices.add(tuple(e))
    if not vertices:
        return Q(0)
    # Coordinate of w cross e.  Use cyclic indices to avoid importing any
    # geometry implementation from the repository.
    j, k = (axis + 1) % 3, (axis + 2) % 3
    values = [(w[j] * e[k] - w[k] * e[j]) / v[axis] for e in vertices]
    return max(values) - min(values)


def slope(v, w):
    s = sum(abs(x) for x in v)
    defects = range(-(3 * s // 14 + 1), 3 * s // 14 + 2)
    widths = {
        d: cube_section_width(v, w, Q(d))
        for d in defects if 14 * abs(d) < 3 * s
    }
    zero_residues = sum(x % 3 == 0 for x in v)
    need(zero_residues <= 1, ("bad residue case", v))
    if zero_residues == 0:
        return Q(2, 3) * sum(width for d, width in widths.items() if d % 3 == 0)
    return Q(1, 3) * sum(width for d, width in widths.items() if d % 3 != 0)


def compile_pattern(p):
    best = Q(-1)
    witnesses = []
    sectors = 0
    vertices = 0
    for v in signed_vectors(p):
        ws = speed_vertices(v)
        if not ws:
            continue
        sectors += 1
        for w in ws:
            vertices += 1
            value = slope(v, w)
            if value > best:
                best = value
                witnesses = [(v, w)]
            elif value == best:
                witnesses.append((v, w))
    need(best >= 0, ("no positive speed sector", p))
    return best, tuple(witnesses), sectors, vertices


def main():
    patterns = [
        p for p in combinations_with_replacement(range(19), 3)
        if eligible_pattern(p)
    ]
    need(len(patterns) == len(set(patterns)), "duplicate patterns")
    support = {2: 0, 3: 0}
    for p in patterns:
        support[sum(x != 0 for x in p)] += 1

    digest = sha256()
    rows = []
    low = {}
    total_sectors = total_vertices = 0
    for p in patterns:
        value, witnesses, sectors, vertices = compile_pattern(p)
        total_sectors += sectors
        total_vertices += vertices
        digest.update(repr((p, value, witnesses, sectors, vertices)).encode("ascii"))
        if p in LOW:
            low[p] = (value, witnesses)
        else:
            need(value <= TARGET, ("slope hostile", p, value, witnesses))
            rows.append((value, p, witnesses))

    rows.sort(key=lambda row: (row[0], row[1]), reverse=True)
    need(rows[0][0] == TARGET and rows[0][1] == (1, 7, 8),
         ("wrong sharp leader", rows[0]))
    unique_sharp = [row[1] for row in rows if row[0] == TARGET]
    need(unique_sharp == [(1, 7, 8)], ("nonunique sharp pattern", unique_sharp))
    need(set(low) == LOW, ("missing low circuit", low))

    print("CLEANROOM ALL-PARITY COEFFICIENT BOX")
    print("UNIVERSE", len(patterns), "SUPPORT2", support[2], "SUPPORT3", support[3])
    print("SIGNED_SECTORS", total_sectors, "SPEED_VERTICES", total_vertices)
    for p in sorted(low):
        print("LOW", p, "MAX", low[p][0], "WITNESSES", low[p][1][:4])
    print("TOP_NONLOW")
    for value, p, witnesses in rows[:12]:
        print(p, value, witnesses[:2])
    print("MAX_NONLOW", rows[0][0], "UNIQUE_PATTERN", rows[0][1])
    print("SEMANTIC_SHA256", digest.hexdigest())
    print("CHECKS", CHECKS)
    print("PASS")


if __name__ == "__main__":
    main()
