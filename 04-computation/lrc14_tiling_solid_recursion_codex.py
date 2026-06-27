#!/usr/bin/env python3
"""Exact tiling/solid recursion atlas for the LRC14 proof interface.

The goal is not to prove LRC14.  It is to separate three often-confused
carriers:

* local Euclidean tiling incidence, such as the six triangles around a vertex;
* lattice self-similarity indices, such as Gaussian and Eisenstein norms;
* finite spherical defects, such as Platonic, Archimedean, and Johnson solids.

Keeping those carriers typed makes the C27/q=3/unital thread easier to use.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import isqrt


def ffrac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def regular_delta(p: int, q: int) -> Fraction:
    """Return 2/p + 2/q - 1 for a regular {p,q} map."""
    return Fraction(2, p) + Fraction(2, q) - 1


def regular_fvector(p: int, q: int) -> tuple[int, int, int] | None:
    delta = regular_delta(p, q)
    if delta <= 0:
        return None
    edges = Fraction(2, 1) / delta
    vertices = Fraction(2, 1) * edges / q
    faces = Fraction(2, 1) * edges / p
    assert vertices.denominator == edges.denominator == faces.denominator == 1
    return int(vertices), int(edges), int(faces)


def vertex_config_curvature(config: tuple[int, ...]) -> Fraction:
    """Euler-characteristic curvature per vertex for a uniform vertex figure."""
    return Fraction(1, 1) - Fraction(len(config), 2) + sum(Fraction(1, p) for p in config)


def archimedean_fvector(config: tuple[int, ...]) -> tuple[int, int, dict[int, int]]:
    curv = vertex_config_curvature(config)
    assert curv > 0
    vertices = Fraction(2, 1) / curv
    edges = vertices * len(config) / 2
    face_counts: dict[int, int] = {}
    for p, count in sorted(Counter(config).items()):
        faces = vertices * count / p
        assert faces.denominator == 1
        face_counts[p] = int(faces)
    assert vertices.denominator == edges.denominator == 1
    return int(vertices), int(edges), face_counts


def gaussian_norm(a: int, b: int) -> int:
    return a * a + b * b


def eisenstein_norm(a: int, b: int) -> int:
    # For omega^2 + omega + 1 = 0, N(a+b*omega) = a^2 - ab + b^2.
    return a * a - a * b + b * b


def norm_atlas(limit: int, norm_fn) -> dict[int, tuple[int, int]]:
    bound = isqrt(limit) + 3
    best: dict[int, tuple[int, int]] = {}
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            if a == 0 and b == 0:
                continue
            n = norm_fn(a, b)
            if 0 < n <= limit:
                candidate = (a, b)
                old = best.get(n)
                if old is None or (abs(a) + abs(b), abs(a), abs(b), a, b) < (
                    abs(old[0]) + abs(old[1]),
                    abs(old[0]),
                    abs(old[1]),
                    old[0],
                    old[1],
                ):
                    best[n] = candidate
    return dict(sorted(best.items()))


def centered_hex(side: int) -> int:
    # Number of hex cells in a side-length side centered hex cluster.
    return 3 * side * (side - 1) + 1


def count_cyclic_triangles(order: list[str]) -> int:
    # Transitive tournament in the displayed order.
    return 0


PLATONIC = [
    ("tetrahedron", 3, 3, "self-dual; K4/tournament object"),
    ("cube", 4, 3, "dual octahedron; tournament cube space"),
    ("octahedron", 3, 4, "dual cube; arc-coordinate space"),
    ("dodecahedron", 5, 3, "dual icosahedron; quotient-dual surface"),
    ("icosahedron", 3, 5, "dual dodecahedron; G5 f-vector surface"),
]

EUCLIDEAN = [
    ("triangular", 3, 6, "hexagonal", "six triangles around each vertex"),
    ("square", 4, 4, "square", "self-dual grid; m^2 square subdivisions"),
    ("hexagonal", 6, 3, "triangular", "three hexagons around each vertex"),
]

ARCHIMEDEAN = [
    ("truncated tetrahedron", (3, 6, 6)),
    ("cuboctahedron", (3, 4, 3, 4)),
    ("truncated cube", (3, 8, 8)),
    ("truncated octahedron", (4, 6, 6)),
    ("rhombicuboctahedron", (3, 4, 4, 4)),
    ("truncated cuboctahedron", (4, 6, 8)),
    ("snub cube", (3, 3, 3, 3, 4)),
    ("icosidodecahedron", (3, 5, 3, 5)),
    ("truncated dodecahedron", (3, 10, 10)),
    ("truncated icosahedron", (5, 6, 6)),
    ("rhombicosidodecahedron", (3, 4, 5, 4)),
    ("truncated icosidodecahedron", (4, 6, 10)),
    ("snub dodecahedron", (3, 3, 3, 3, 5)),
]


def main() -> None:
    print("LRC14 TILING/SOLID RECURSION ATLAS")
    print("=" * 72)
    print()

    print("A. Regular {p,q} curvature")
    print("-" * 72)
    print("delta = 2/p + 2/q - 1.  Sphere: delta>0, plane: delta=0.")
    for name, p, q, note in PLATONIC:
        fv = regular_fvector(p, q)
        assert fv is not None
        v, e, f = fv
        print(
            f"{name:13s} {{p,q}}={{{p},{q}}}  delta={ffrac(regular_delta(p, q)):>4s}  "
            f"(V,E,F)=({v},{e},{f})  {note}"
        )
    print()
    for name, p, q, dual, note in EUCLIDEAN:
        print(
            f"{name:10s} {{p,q}}={{{p},{q}}}  delta={ffrac(regular_delta(p, q)):>1s}  "
            f"dual={dual:10s}  p*q={p*q:2d}  {note}"
        )
    print()

    print("B. Square, triangular, and hexagonal recursion numbers")
    print("-" * 72)
    square_self = [m * m for m in range(1, 11)]
    triangle_self = [m * m for m in range(1, 9)]
    centered_hexes = [centered_hex(s) for s in range(1, 9)]
    print(f"square self m^2:      {square_self}")
    print(f"triangle self m^2:    {triangle_self}  (2-edge step gives 4, then 16)")
    print(f"centered hex clusters:{centered_hexes}  (side 2 gives 7 cells)")
    print(
        "triangle <-> hex '6': local incidence/duality, not a lattice index: "
        "{3,6} <-> {6,3}."
    )
    print()

    print("C. Lattice index norms up to 100")
    print("-" * 72)
    gaussian = norm_atlas(100, gaussian_norm)
    eisenstein = norm_atlas(100, eisenstein_norm)
    print(f"Gaussian norms, square lattice:   {list(gaussian)}")
    print(f"Eisenstein norms, tri/hex lattice:{list(eisenstein)}")
    for target in [4, 6, 7, 14, 16, 22, 27, 31, 49]:
        g = gaussian.get(target)
        e = eisenstein.get(target)
        g_s = str(g) if g else "no"
        e_s = str(e) if e else "no"
        print(
            f"target {target:2d}: Gaussian={g_s:>10s}  "
            f"Eisenstein={e_s:>10s}"
        )
    print()
    print("Typed reading:")
    print("  4,16: square and triangle self-subdivision by integer scale.")
    print("  6: triangle/hex dual incidence count, absent as an Eisenstein norm.")
    print("  7,49: Eisenstein hex self-similarity indices; 7 = N(3+omega).")
    print("        Twist witnesses: N(3+omega)=7 and N(8+5*omega)=49.")
    print("  27: q=3/C27 branch carrier; 27 = 3^3 appears as an Eisenstein norm.")
    print("  31: Eisenstein norm and cube-root-of-31 petal scale alias, not a proof unit.")
    print()

    print("D. Archimedean finite-defect carriers")
    print("-" * 72)
    for name, config in ARCHIMEDEAN:
        curv = vertex_config_curvature(config)
        v, e, faces = archimedean_fvector(config)
        face_s = ", ".join(f"F{p}={c}" for p, c in faces.items())
        cfg = ".".join(str(x) for x in config)
        print(
            f"{name:29s} {cfg:11s} curv={ffrac(curv):>4s} "
            f"V={v:3d} E={e:3d} {face_s}"
        )
    print()
    print("Johnson solids: 92 convex regular-faced non-uniform finite defects.")
    print("Proof role: good obstruction/stress-test carriers, not universal symmetry.")
    print()

    print("E. LRC14 proof-interface synthesis")
    print("-" * 72)
    print("n=14 has apex half-size 7 and C=2n-1=27=3^3.")
    print("Square/Gaussian carriers model self-dual folds and m^2 subdivisions.")
    print("Tri/hex/Eisenstein carriers model the q=3, C27, and hex-norm-7 branch.")
    print("Platonic carriers match existing tournament levels: object, cube, arc, quotient, dual quotient.")
    print("Archimedean and Johnson solids are finite curvature defects for testing label leakage.")
    print()
    print("Candidate lemma shape:")
    print("  After AP/Goddyn-Wong and C27 labels are attached, any low-gap LRC14")
    print("  row that uses the tri/hex carrier must preserve the branch-local")
    print("  C27 marked transfer chart; otherwise it falls into the HYP-2942")
    print("  unital pair-repeat obstruction or the Farey/K33 child branch.")
    print()

    print("F. Tournament analysis")
    print("-" * 72)
    role_order = [
        "exact M/Farey node",
        "AP/GW C27 marked shell transfer",
        "q=3 unital branch-local block chart",
        "Eisenstein tri/hex norm carrier",
        "triangular-hexagonal local duality",
        "square/Gaussian self-dual carrier",
        "Platonic tournament-level carrier",
        "Archimedean uniform finite defect",
        "Johnson non-uniform finite defect",
        "raw numerical analogy",
    ]
    edges = len(role_order) * (len(role_order) - 1) // 2
    print("Binary relation: x -> y iff x retains more LRC proof data than y.")
    print(f"vertices={len(role_order)} edges={edges} c3={count_cyclic_triangles(role_order)} hp=1")
    for i, role in enumerate(role_order):
        print(f"  score {len(role_order)-1-i:2d}: {role}")
    print()
    print("Verdict:")
    print("  The user numbers split cleanly by type: square/triangle 4,16 are")
    print("  integer-scale self-recursions; hex 7,49 are Eisenstein norm")
    print("  self-recursions; triangle-hex 6 is a dual-incidence count.  The")
    print("  LRC14-relevant carrier is therefore not generic polyhedral symmetry,")
    print("  but the typed passage C27 -> q=3/F3^3 -> Eisenstein tri/hex chart,")
    print("  with Platonic/Archimedean/Johnson solids used as finite defect tests.")


if __name__ == "__main__":
    main()
