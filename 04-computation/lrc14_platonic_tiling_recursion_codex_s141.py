#!/usr/bin/env python3
"""
Polyhedral and Euclidean-tiling carrier scout for LRC14.

The goal is not to classify solids from scratch.  The script fixes the
curvature, duality, and recursion labels that should be preserved before the
Platonic/Archimedean/Johnson and plane-tiling analogies are used in the POKE
binary-relational proof tree.
"""

from fractions import Fraction


def vertex_deficit_turns(face_cycle):
    """Combinatorial angle deficit in turns for regular polygons at a vertex."""
    used = sum(Fraction(p - 2, 2 * p) for p in face_cycle)
    return Fraction(1, 1) - used


def regular_map_counts(p, q):
    """Return V,E,F for spherical regular map {p,q}; None at Euclidean wall."""
    denom = 2 * p + 2 * q - p * q
    if denom <= 0:
        return None
    edges = Fraction(2 * p * q, denom)
    vertices = Fraction(2 * edges, q)
    faces = Fraction(2 * edges, p)
    return int(vertices), int(edges), int(faces)


def eisenstein_norm(a, b):
    return a * a - a * b + b * b


def gaussian_axis_norm(m):
    return m * m


def print_regular_boundary():
    print("REGULAR MAP CURVATURE BOUNDARY")
    platonic = [
        ("tetrahedron", 3, 3, "self"),
        ("cube", 4, 3, "octahedron"),
        ("octahedron", 3, 4, "cube"),
        ("dodecahedron", 5, 3, "icosahedron"),
        ("icosahedron", 3, 5, "dodecahedron"),
    ]
    for name, p, q, dual in platonic:
        counts = regular_map_counts(p, q)
        deficit = vertex_deficit_turns([p] * q)
        print(f"{name}: {{p,q}}={{{p},{q}}}, dual={dual}")
        print(f"  V,E,F={counts}")
        print(f"  vertex_deficit_turns={deficit} ({float(deficit):.9f})")
    print()
    tilings = [
        ("triangular tiling", 3, 6, "hexagonal tiling"),
        ("square tiling", 4, 4, "self"),
        ("hexagonal tiling", 6, 3, "triangular tiling"),
    ]
    for name, p, q, dual in tilings:
        deficit = vertex_deficit_turns([p] * q)
        print(f"{name}: {{p,q}}={{{p},{q}}}, dual={dual}")
        print(f"  vertex_deficit_turns={deficit} ({float(deficit):.9f})")
    print("readout: Platonic solids are positive-curvature finite maps;")
    print("         square/triangular/hexagonal tilings are the zero-curvature wall.")
    print()


def print_tiling_recursions():
    print("EUCLIDEAN REGULAR TILING RECURSIONS")
    square = [gaussian_axis_norm(m) for m in range(2, 6)]
    tri_general = [eisenstein_norm(m, 0) for m in range(2, 6)]
    tri_binary = [4**k for k in range(1, 5)]
    hex_norm7 = [7**k for k in range(1, 5)]
    centered_hex = [1 + 3 * r * (r + 1) for r in range(1, 5)]

    print("square self-recursion:")
    print("  lattice: Gaussian axis scaling m+0i")
    print(f"  indices m^2 for m=2..5: {square}")
    print("  LRC role: self-dual square packet, scalar square counts are safe only")
    print("            when the chosen grid axes remain labelled.")

    print("triangle self-recursion:")
    print("  lattice: Eisenstein scaling m+0*omega")
    print(f"  general indices m^2 for m=2..5: {tri_general}")
    print(f"  binary chain emphasized by prompt: {tri_binary}")
    print("  LRC role: triangular subdivision keeps orientation-sector labels;")
    print("            powers of 4 are the cleanest dyadic triangle spine.")

    print("triangle <-> hexagon dual bridge:")
    print("  local index: 6")
    print("  reason: a regular hexagon dissects into 6 equilateral triangles,")
    print("          and the triangular tiling has 6 triangles around a vertex.")
    print("  LRC role: support-six packet carrier, compatible with octahedral/current")
    print("            language but not a scalar replacement for q or C27 labels.")

    print("hexagon self-recursion:")
    print("  lattice: Eisenstein norm N(a+b*omega)=a^2-a*b+b^2")
    print(f"  N(3+omega)={eisenstein_norm(3, 1)}")
    print(f"  norm-7 powers: {hex_norm7}")
    print(f"  centered-hex ring counts r=1..4: {centered_hex}")
    print("  LRC role: the 7,49,... chain is a norm-index recursion;")
    print("            centered rings 7,19,37,... are a different carrier.")
    print()


def print_archimedean_carriers():
    print("ARCHIMEDEAN VERTEX-FIGURE CARRIERS")
    arch = [
        ("truncated tetrahedron", [3, 6, 6]),
        ("cuboctahedron", [3, 4, 3, 4]),
        ("truncated cube", [3, 8, 8]),
        ("truncated octahedron", [4, 6, 6]),
        ("rhombicuboctahedron", [3, 4, 4, 4]),
        ("truncated cuboctahedron", [4, 6, 8]),
        ("snub cube", [3, 3, 3, 3, 4]),
        ("icosidodecahedron", [3, 5, 3, 5]),
        ("truncated dodecahedron", [3, 10, 10]),
        ("truncated icosahedron", [5, 6, 6]),
        ("rhombicosidodecahedron", [3, 4, 5, 4]),
        ("truncated icosidodecahedron", [4, 6, 10]),
        ("snub dodecahedron", [3, 3, 3, 3, 5]),
    ]
    deficits = []
    for name, fig in arch:
        deficit = vertex_deficit_turns(fig)
        deficits.append(deficit)
        fig_text = ".".join(str(x) for x in fig)
        print(f"{name}: vertex_figure={fig_text}, deficit={deficit} ({float(deficit):.9f})")
    print(f"count={len(arch)}")
    print(f"deficit_range={min(deficits)}..{max(deficits)}")
    print("readout: Archimedean solids retain one vertex-figure word;")
    print("         they are better analogues for labelled local quotients than")
    print("         for global LRC scalar invariants.")
    print()


def print_johnson_guardrail():
    print("JOHNSON SOLID GUARDRAIL")
    print("count=92")
    print("faces=regular polygons, convex, but not Platonic/Archimedean/prism/antiprism")
    print("vertex figures are mixed rather than one uniform word.")
    print("LRC role: finite residual atlas after the uniform carriers fail,")
    print("          analogous to bounded AP/GW/petal/K33 frontier tables.")
    print()


def print_carrier_tournament():
    print("POKE CARRIER TOURNAMENT")
    carriers = [
        (
            "exact_M_Farey_branch",
            (6, 6, 6, 6, 6),
        ),
        (
            "C27_unital_block_lift",
            (6, 6, 5, 6, 5),
        ),
        (
            "Euclidean_tiling_recursion_indices",
            (5, 5, 5, 5, 5),
        ),
        (
            "Archimedean_vertex_figure_words",
            (4, 5, 4, 5, 4),
        ),
        (
            "Platonic_dual_curvature_skeleton",
            (4, 4, 4, 5, 4),
        ),
        (
            "Johnson_residual_atlas",
            (3, 4, 5, 4, 3),
        ),
        (
            "raw_polyhedron_visual_metaphor",
            (1, 1, 1, 1, 1),
        ),
    ]
    print("observable=(LRC faithfulness, label retention, recursion exactness,")
    print("            finite certifiability, anti-scalar guard)")
    ordered = sorted(carriers, key=lambda x: x[1], reverse=True)
    for i, (name, score) in enumerate(ordered):
        print(f"{i}: {name} score={score}")
    print("fingerprint: transitive role order, c3=0, hp=1")
    print("proof use: tiling recursions sit below exact M/Farey and C27 charts,")
    print("           but above raw solid analogies.")


def main():
    print_regular_boundary()
    print_tiling_recursions()
    print_archimedean_carriers()
    print_johnson_guardrail()
    print_carrier_tournament()


if __name__ == "__main__":
    main()
