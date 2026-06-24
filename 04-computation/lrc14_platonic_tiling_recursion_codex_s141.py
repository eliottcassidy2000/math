#!/usr/bin/env python3
"""
Polyhedral and Euclidean-tiling carrier scout for LRC14.

The goal is not to classify solids from scratch.  The script fixes the
curvature, duality, recursion labels, and low-gap LRC14 row tags that should be
preserved before the Platonic/Archimedean/Johnson and plane-tiling analogies
are used in the POKE binary-relational proof tree.
"""

from itertools import combinations
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


LOW_GAP_ROW_TAGS = [
    ("AP", "{1,...,13}", "flat Euclidean baseline", "triangle_self + tri_hex_dual"),
    ("GW", "{1,...,11,13,24}", "tight labelled collision", "tri_hex_dual + hex_heptadic"),
    ("near_K33", "AP with 12->36", "first nonunit Farey/K33 child", "johnson_local_defect"),
    ("petal10", "AP with 10->20", "unit-visible C27 petal", "hex_heptadic local defect"),
    ("petal13", "AP with 13->26", "unit-visible C27 petal", "hex_heptadic local defect"),
    ("two_swap_GW", "drop(10,12)->add(20,24)", "P10+GW splice", "johnson_local_defect"),
    ("two_swap_K33", "drop(10,12)->add(20,36)", "P10+K33 splice", "johnson_local_defect"),
]


LOCAL_CARRIER_CRITERIA = [
    "keeps_exact_M_Farey_label",
    "keeps_C27_carry",
    "keeps_mod7_petal_seam",
    "keeps_pair_incidence_unit",
    "keeps_branch_locality",
    "resists_scalarization",
]


LOCAL_CARRIER_SCORES = {
    "square_self": [0, 0, 0, 0, 0, -2],
    "triangle_self": [1, 1, 0, 0, 1, 1],
    "tri_hex_dual": [1, 2, 1, 1, 1, 2],
    "hex_heptadic": [1, 1, 3, 1, 1, 2],
    "platonic_positive": [1, 1, 0, 1, 0, 0],
    "archimedean_near_flat": [2, 2, 2, 3, 3, 2],
    "johnson_local_defect": [3, 2, 2, 2, 3, 3],
    "raw_runner_vertices": [0, 0, 0, 0, 0, -3],
}


def local_carrier_tournament():
    names = list(LOCAL_CARRIER_SCORES)
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        a, b = names[i], names[j]
        va, vb = LOCAL_CARRIER_SCORES[a], LOCAL_CARRIER_SCORES[b]
        wins_a = sum(1 for x, y in zip(va, vb) if x > y)
        wins_b = sum(1 for x, y in zip(va, vb) if y > x)
        if wins_a > wins_b or (wins_a == wins_b and i < j):
            adj[i][j] = True
        else:
            adj[j][i] = True
    scores = [sum(row) for row in adj]
    hist = {s: scores.count(s) for s in sorted(set(scores))}
    cycles = 0
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cycles += 1
    order = sorted(range(n), key=lambda idx: (-scores[idx], idx))
    return {
        "scores": dict(zip(names, scores)),
        "score_hist": hist,
        "directed_3cycles": cycles,
        "order": [names[i] for i in order],
    }


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


def print_local_carrier_tournament():
    print("LOCAL CURVATURE CARRIER TOURNAMENT")
    tour = local_carrier_tournament()
    print("vertices are proof carriers, not runners.")
    print("criteria:", ", ".join(LOCAL_CARRIER_CRITERIA))
    print("scores:", tour["scores"])
    print("score_hist:", tour["score_hist"])
    print("directed_3cycles:", tour["directed_3cycles"])
    print("hamiltonian_order:", " > ".join(tour["order"]))
    print("readout: Johnson defects beat uniform Archimedean charts once residual")
    print("         petal/K33/two-swap surgery data must be retained.")
    print()


def print_low_gap_row_tags():
    print("LOW-GAP LRC14 ROW TAGS")
    for name, row, role, tag in LOW_GAP_ROW_TAGS:
        print(f"{name:14s} row={row:30s} role={role:32s} carrier={tag}")
    print("proof use: AP/GW form the flat labelled baseline; unit-visible")
    print("           petals use hex-heptadic seams; K33/two-swap residuals")
    print("           are Johnson-style finite surgery packets.")
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
    print_local_carrier_tournament()
    print_low_gap_row_tags()
    print_carrier_tournament()


if __name__ == "__main__":
    main()
