"""Exact flat-star census and boundary-aware hierarchical tiling controls.

No third-party dependencies; assertions deliberately remain active under -O.
The unbounded statements are proved in the companion note, not inferred here.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import permutations
from math import gcd, lcm
import argparse
from pathlib import Path


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def reciprocal_partitions(k):
    """Complete Egyptian-fraction search; denominator bound is derived."""
    def rec(left, remaining, lower, word):
        if left == 1:
            if remaining.numerator == 1 and remaining.denominator >= lower:
                yield word + (remaining.denominator,)
            return
        # All future denominators >= p imply remaining <= left/p.
        for p in range(lower, int(left / remaining) + 1):
            rest = remaining - Q(1, p)
            if rest > 0:
                yield from rec(left - 1, rest, p, word + (p,))
    return tuple(rec(k, Q(k - 2, 2), 3, ()))


def bracelet(word):
    w = tuple(word)
    return min(v[i:] + v[:i] for v in (w, w[::-1])
               for i in range(len(w)))


def angle(p):
    return Q(p - 2, 2 * p)  # turns


def act(v, rotations, reflection, modulus):
    a, b = v
    if reflection:
        a, b = a + b, -b
    for _ in range(rotations):
        a, b = -b, a + b
    return a % modulus, b % modulus


def orbit_count(n):
    unseen = {(a, b) for a in range(n) for b in range(n)} - {(0, 0)}
    count = 0
    while unseen:
        v = min(unseen)
        orbit = {act(v, j, r, n) for j in range(6) for r in (0, 1)}
        require(orbit <= unseen, "D6 orbits form a partition")
        unseen -= orbit
        count += 1
    return count


def draw_figure(path):
    """Unit triangular faces, with whole six-triangle clusters merged.

    Rendering alone uses square roots; none of the exact checks do.
    The rhombus is a finite crop of each full-plane tiling.
    """
    from math import sqrt
    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="720" viewBox="0 0 1400 720">',
        '<rect width="1400" height="720" fill="#f8fafc"/>',
        '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#152a3a}.title{font-size:27px;font-weight:600}.sub{font-size:19px}.caption{font-size:17px}</style>',
        '<text x="40" y="43" class="title">The same two local stars can carry very different global addresses</text>',
        '<text x="40" y="76" class="sub">Six unit triangles merge into each hexagon. All shared edges fit exactly.</text>',
    ]
    directions = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))
    for panel, kind in enumerate(("periodic", "hierarchical")):
        ox, oy, scale = 40 + panel * 700, 126, 17
        def point(v):
            a, b = v
            return ox + (a + b / 2 + 2) * scale, oy + 400 - (sqrt(3) * b / 2 + 2) * scale
        def points(vertices):
            return " ".join(f"{x:.4f},{y:.4f}" for x, y in map(point, vertices))
        def chosen(a, b):
            if kind == "periodic":
                return a % 4 == 0 and b % 4 == 0
            return a >= 0 and b >= 0 and a % 3 == b % 3 == 0 and not ((a // 3) & (b // 3))
        boundary = ((-1, -1), (22, -1), (22, 22), (-1, 22))
        elements += [f'<text x="{ox}" y="118" class="sub">'
                     + ("Periodic centers: 4Λ · exactly 3 vertex orbits" if panel == 0
                        else "Hierarchical centers: 3(i,j), i AND j = 0") + '</text>',
                     f'<defs><clipPath id="clip{panel}"><polygon points="{points(boundary)}"/></clipPath></defs>',
                     f'<g clip-path="url(#clip{panel})">']
        for a in range(-3, 26):
            for b in range(-3, 26):
                triangles = (((a, b), (a + 1, b), (a, b + 1)),
                             ((a + 1, b), (a, b + 1), (a + 1, b + 1)))
                for tri in triangles:
                    if not any(chosen(*v) for v in tri):
                        elements.append(f'<polygon points="{points(tri)}" fill="#edf3f6" stroke="#778994" stroke-width="0.55"/>')
        for a in range(-3, 26):
            for b in range(-3, 26):
                if chosen(a, b):
                    hexagon = tuple((a + da, b + db) for da, db in directions)
                    elements.append(f'<polygon points="{points(hexagon)}" fill="#58b6a6" stroke="#135d57" stroke-width="1.3"/>')
        elements += ['</g>', f'<polygon points="{points(boundary)}" fill="none" stroke="#243d4d" stroke-width="1.5"/>']
    elements += [
        '<text x="40" y="572" class="caption">Local stars in both tilings: six triangles; or four triangles and one hexagon.</text>',
        '<text x="40" y="605" class="caption">Left: a finite crop of a periodic full-plane tiling. More widely spaced centers give unbounded k.</text>',
        '<text x="40" y="638" class="caption">Right: 27 selected centers in the 8 × 8 address window; Q = 2Q + {(0,0),(1,0),(0,1)}.</text>',
        '<text x="40" y="671" class="caption">The right full-plane tiling has infinitely many vertex orbits. Its center-selection hierarchy is fractal; all tiles stay regular.</text>',
        '</svg>',
    ]
    Path(path).write_text("\n".join(elements) + "\n", encoding="utf-8")


def main():
    multisets = {k: reciprocal_partitions(k) for k in range(3, 7)}
    stars = {k: sorted({bracelet(w) for m in multisets[k]
                       for w in permutations(m)}) for k in multisets}
    require([len(multisets[k]) for k in range(6, 2, -1)] == [1, 2, 4, 10],
            "17 multisets")
    require([len(stars[k]) for k in range(6, 2, -1)] == [1, 3, 7, 10],
            "21 dihedral stars")
    # Independent integer-denominator census, after the complete search has
    # established the maximum denominator 42 (not an inherited assumption).
    require(max(p for row in multisets.values() for w in row for p in w) == 42,
            "derived maximum denominator")
    scale = lcm(*range(3, 43))
    masses = {p: scale * (p - 2) // p for p in range(3, 43)}
    def integer_angles(left, start, total, word):
        if left == 0:
            if total == 2 * scale:
                yield word
            return
        for p in range(start, 43):
            if total + left * masses[p] > 2 * scale:
                break
            if total + masses[p] + (left - 1) * masses[42] < 2 * scale:
                continue
            yield from integer_angles(left - 1, p, total + masses[p], word + (p,))
    for k in multisets:
        brute = set(integer_angles(k, 3, 0, ()))
        require(brute == set(multisets[k]), "independent complete census")
    print("Complete angle census: 17 multisets, 21 dihedral cyclic stars")
    for k in range(6, 2, -1):
        print(f"valency {k}: multisets={len(multisets[k])}, stars={len(stars[k])}")
        print("  " + "; ".join(".".join(map(str, w)) for w in stars[k]))

    all_multisets = sum((list(v) for v in multisets.values()), [])
    # For each hostile star, choose an odd central polygon.  Once one
    # neighbor is present, the census forces alternating distinct neighbors.
    hostile = [((3, 7, 42), 7, 3, 42), ((3, 8, 24), 3, 8, 24),
               ((3, 9, 18), 9, 3, 18), ((3, 10, 15), 15, 3, 10),
               ((4, 5, 20), 5, 4, 20), ((5, 5, 10), 5, 5, 10)]
    for star, p, q, r in hostile:
        require(p % 2 == 1 and q != r, "odd boundary alternation")
        for neighbor, forced in ((q, r), (r, q)):
            matching = []
            for w in all_multisets:
                c = Counter(w)
                c.subtract([p, neighbor])
                if min(c.values()) >= 0:
                    matching.append(tuple(sorted(c.elements())))
            require(matching == [(forced,)], "all completions forced by census")
        require((q, r)[p % 2] != q, "odd-length word fails closure")
        print(f"Nonextension certificate {star}: central {p}, neighbors {q}<->{r}")

    all_stars = {w for row in stars.values() for w in row}
    edges = set()
    for w in all_stars:
        for j, p in enumerate(w):
            splits = ((3, 3),) if p == 6 else ((3, 4), (4, 3)) if p == 12 else ()
            for split in splits:
                target = bracelet(w[:j] + split + w[j + 1:])
                require(angle(p) == sum(map(angle, split)), "angle-preserving split")
                require(target in all_stars, "flat-star refinement remains in census")
                require(len(target) == len(w) + 1, "strict DAG grading")
                edges.add((w, target))
    print(f"Flat-star refinement DAG: 21 vertices, {len(edges)} directed edges")
    for a, b in sorted(edges):
        print(f"  {a} -> {b}")

    spherical = [(p, q) for p in range(3, 7) for q in range(3, 7)
                 if (p - 2) * (q - 2) < 4]
    euclidean = [(p, q) for p in range(3, 7) for q in range(3, 7)
                 if (p - 2) * (q - 2) == 4]
    require(spherical == [(3, 3), (3, 4), (3, 5), (4, 3), (5, 3)], "five solids")
    require(euclidean == [(3, 6), (4, 4), (6, 3)], "three regular plane tilings")
    print(f"Positive/zero regular curvature: {spherical} / {euclidean}")

    # One merged hexagon per N^2-cell triangular-lattice fundamental domain.
    for n in range(3, 61):
        numerator = n*n + 6*n - 10 + 2*gcd(n, 3) + gcd(n, 2)**2
        require(numerator % 12 == 0, "Burnside integrality")
        require(orbit_count(n) == numerator // 12, "independent D6 orbit census")
    print("Periodic two-star tilings: k_N exact for N=3..60")
    print("  k_3..k_20 = " + ",".join(str(orbit_count(n)) for n in range(3, 21)))
    print("  k_N=(N^2+6N-10+2*gcd(N,3)+gcd(N,2)^2)/12")
    # Hierarchical centre selection.  This is a decoration hierarchy, not a
    # claim that the full unit-sided tiling is invariant under scaling.
    current = {(0, 0)}
    for r in range(0, 9):
        direct = {(a, b) for a in range(2**r) for b in range(2**r) if not (a & b)}
        require(current == direct, "independent bit and substitution construction")
        require(len(current) == 3**r, "three-branch hierarchy")
        current = {(2*a + i, 2*b + j) for a, b in current
                   for i, j in ((0, 0), (1, 0), (0, 1))}
    print("Sierpinski centre recursion: bit predicate and 3-branch substitution agree, r=0..8")
    print("All exact checks passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure", metavar="PATH", help="write the optional standalone SVG")
    args = parser.parse_args()
    if args.figure:
        draw_figure(args.figure)
    else:
        main()
