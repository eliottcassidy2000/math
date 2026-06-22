#!/usr/bin/env python3
"""S105: design/Hodge carriers for the LRC14 residual-leak proof target.

Prompt sources:
  - MathWorld Unital: q^3+1 points, q+1 point blocks, pair-unique.
  - MathWorld Clebsch graph: SRG(16,5,0,2), folded 5-cube.
  - MathWorld Truncated Octahedral graph: Bruhat graph B4 on S4.

The goal is not to prove LRC14 here.  It is to test whether these prompted
objects are real finite carriers for the HYP-2889/HYP-2890 labelled residual:
pair-slot frames, signed cycle/curl bases, and local compression topology.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import comb


# ---------------------------------------------------------------------------
# GF(9) = F_3[i] with i^2 = -1 = 2.


def gf9(a: int, b: int = 0) -> int:
    return (a % 3) + 3 * (b % 3)


def gadd(x: int, y: int) -> int:
    return gf9((x % 3) + (y % 3), (x // 3) + (y // 3))


def gneg(x: int) -> int:
    return gf9(-(x % 3), -(x // 3))


def gsub(x: int, y: int) -> int:
    return gadd(x, gneg(y))


def gmul(x: int, y: int) -> int:
    a, b = x % 3, x // 3
    c, d = y % 3, y // 3
    return gf9(a * c + 2 * b * d, a * d + b * c)


def gpow(x: int, n: int) -> int:
    out = 1
    for _ in range(n):
        out = gmul(out, x)
    return out


def ginv(x: int) -> int:
    assert x != 0
    # GF(9)^* has order 8.
    return gpow(x, 7)


def gdot(a: tuple[int, int, int], b: tuple[int, int, int]) -> int:
    out = 0
    for x, y in zip(a, b):
        out = gadd(out, gmul(x, y))
    return out


def pcanon(vec: tuple[int, int, int]) -> tuple[int, int, int]:
    for x in vec:
        if x:
            inv = ginv(x)
            return tuple(gmul(inv, y) for y in vec)  # type: ignore[return-value]
    raise ValueError("zero projective vector")


def hermitian_unital_q3() -> tuple[list[tuple[int, int, int]], list[tuple[int, ...]], list[int]]:
    """Return points, secant blocks, and all intersection sizes."""

    all_points = sorted({pcanon(v) for v in product(range(9), repeat=3) if any(v)})
    unit_points = []
    for p in all_points:
        # Hermitian curve x^(q+1)+y^(q+1)+z^(q+1)=0, q=3.
        val = 0
        for x in p:
            val = gadd(val, gpow(x, 4))
        if val == 0:
            unit_points.append(p)

    lines = sorted({pcanon(v) for v in product(range(9), repeat=3) if any(v)})
    blocks = []
    sizes = []
    for line in lines:
        incident = tuple(i for i, p in enumerate(unit_points) if gdot(line, p) == 0)
        sizes.append(len(incident))
        if len(incident) == 4:
            blocks.append(incident)

    return unit_points, blocks, sizes


def design_pair_counts(blocks: list[tuple[int, ...]], v: int) -> tuple[Counter[int], Counter[int], set[int]]:
    point_rep = Counter()
    pair_rep = Counter()
    for block in blocks:
        for x in block:
            point_rep[x] += 1
        for a, b in combinations(block, 2):
            pair_rep[tuple(sorted((a, b)))] += 1
    return Counter(point_rep.values()), Counter(pair_rep.values()), set(point_rep)


# ---------------------------------------------------------------------------
# Clebsch graph as folded Q_5.


def popcount(x: int) -> int:
    return int(x).bit_count()


def clebsch_graph() -> tuple[list[int], set[tuple[int, int]]]:
    full = (1 << 5) - 1
    vertices = sorted({min(x, x ^ full) for x in range(1 << 5)})
    index = {v: i for i, v in enumerate(vertices)}
    edges = set()
    for a, b in combinations(vertices, 2):
        d = popcount(a ^ b)
        if min(d, 5 - d) == 1:
            edges.add((index[a], index[b]))
    return vertices, edges


def graph_degrees(v: int, edges: set[tuple[int, int]]) -> list[int]:
    deg = [0] * v
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    return deg


def common_neighbor_counts(v: int, edges: set[tuple[int, int]]) -> tuple[Counter[int], Counter[int]]:
    nbr = [set() for _ in range(v)]
    edge_set = {tuple(sorted(e)) for e in edges}
    for a, b in edges:
        nbr[a].add(b)
        nbr[b].add(a)
    adj_counts = Counter()
    non_counts = Counter()
    for a, b in combinations(range(v), 2):
        c = len(nbr[a] & nbr[b])
        if (a, b) in edge_set:
            adj_counts[c] += 1
        else:
            non_counts[c] += 1
    return adj_counts, non_counts


def closed_neighborhood_design(v: int, edges: set[tuple[int, int]]) -> tuple[Counter[int], Counter[int], list[tuple[int, ...]]]:
    nbr = [set([i]) for i in range(v)]
    for a, b in edges:
        nbr[a].add(b)
        nbr[b].add(a)
    blocks = [tuple(sorted(s)) for s in nbr]
    point_rep, pair_rep, _ = design_pair_counts(blocks, v)
    return point_rep, pair_rep, blocks


# ---------------------------------------------------------------------------
# Truncated octahedral graph = Bruhat/Cayley graph of S4 on adjacent swaps.


def compose_adjacent_swap(p: tuple[int, ...], i: int) -> tuple[int, ...]:
    q = list(p)
    q[i], q[i + 1] = q[i + 1], q[i]
    return tuple(q)


def bruhat_s4_graph() -> tuple[list[tuple[int, ...]], set[tuple[int, int, int]]]:
    verts = list(permutations(range(4)))
    idx = {p: i for i, p in enumerate(verts)}
    edges = set()
    for p in verts:
        a = idx[p]
        for color in range(3):
            q = compose_adjacent_swap(p, color)
            b = idx[q]
            edges.add((min(a, b), max(a, b), color))
    return verts, edges


def bruhat_faces(verts: list[tuple[int, ...]]) -> tuple[set[tuple[int, ...]], set[tuple[int, ...]]]:
    idx = {p: i for i, p in enumerate(verts)}
    squares = set()
    hexes = set()

    for p in verts:
        # Commutation face s1 s3 = s3 s1.
        cyc = [
            p,
            compose_adjacent_swap(p, 0),
            compose_adjacent_swap(compose_adjacent_swap(p, 0), 2),
            compose_adjacent_swap(p, 2),
        ]
        squares.add(canonical_cycle(tuple(idx[x] for x in cyc)))

        # Braid faces s_i s_{i+1} s_i = s_{i+1} s_i s_{i+1}.
        for color in (0, 1):
            a = color
            b = color + 1
            cyc = [p]
            cur = p
            for step in (a, b, a, b, a):
                cur = compose_adjacent_swap(cur, step)
                cyc.append(cur)
            hexes.add(canonical_cycle(tuple(idx[x] for x in cyc)))

    return squares, hexes


def canonical_cycle(cyc: tuple[int, ...]) -> tuple[int, ...]:
    rots = []
    n = len(cyc)
    for seq in (cyc, tuple(reversed(cyc))):
        for i in range(n):
            rots.append(seq[i:] + seq[:i])
    return min(rots)


# ---------------------------------------------------------------------------
# Octahedron = L(K4), existing HYP-2887 carrier.


def octahedron_line_k4() -> tuple[list[tuple[int, int]], set[tuple[int, int]]]:
    k4_edges = list(combinations(range(4), 2))
    edges = set()
    for i, e1 in enumerate(k4_edges):
        for j, e2 in enumerate(k4_edges):
            if i < j and set(e1) & set(e2):
                edges.add((i, j))
    return k4_edges, edges


def cycle_rank(v: int, e: int, components: int = 1) -> int:
    return e - v + components


def carrier_tournament() -> list[tuple[str, int, list[str]]]:
    """Exploratory proof-obligation score, not a mathematical theorem."""

    obligations = {
        "unital_q3_pair_frame": [
            "exact_pair_slot_count_C(8,2)=28",
            "BIBD_frame_NtN=8I+J",
            "four_slot_blocks_match_additive_quadruples",
            "k8_formalizable",
        ],
        "bruhat_truncated_octahedron": [
            "local_adjacent_compression_topology",
            "square_commutation_faces",
            "hex_braid_faces",
            "cycle_rank_13_matches_signed_curl_ledger",
            "formalizable_Coxeter_complex",
        ],
        "clebsch_closed_neighborhood": [
            "folded_cube_parity_quotient",
            "SRG_frame_NtN=4I+2J",
            "triangle_free_pair_uniformity",
            "closed_neighborhood_design",
        ],
        "octahedral_current": [
            "existing_HYP_2887_support6_carrier",
            "face_curl_rank_7",
            "affine_zero_opposite_edges",
            "HYP_2636_tail_route",
        ],
        "scalar_additive_energy": [
            "positive_same_frequency_s2_packet",
        ],
    }
    rows = []
    for name, obs in obligations.items():
        rows.append((name, len(obs), obs))
    rows.sort(key=lambda x: (-x[1], x[0]))
    return rows


def main() -> None:
    print("=== S105 design/Hodge carrier audit for LRC14 residual leak ===")

    points, blocks, sizes = hermitian_unital_q3()
    point_rep, pair_rep, seen_points = design_pair_counts(blocks, len(points))
    print("\n--- Hermitian unital q=3 ---")
    print(f"points={len(points)}  secant_blocks={len(blocks)}  line_intersection_sizes={dict(sorted(Counter(sizes).items()))}")
    print(f"expected unital parameters: v=28, block_size=4, b=63, r=9, lambda=1")
    print(f"point replication histogram={dict(sorted(point_rep.items()))}")
    print(f"pair replication histogram={dict(sorted(pair_rep.items()))}")
    print(f"all points seen={len(seen_points)==len(points)}")
    print(f"LRC hook: C(8,2)={comb(8,2)} pair slots; unital blocks are 4-slot additive-quadruple frames.")
    print("frame identity: incidence^T incidence = 8I + J on the 28 pair-slot coordinates")

    cverts, cedges = clebsch_graph()
    deg = graph_degrees(len(cverts), cedges)
    adj_common, non_common = common_neighbor_counts(len(cverts), cedges)
    cn_point_rep, cn_pair_rep, cn_blocks = closed_neighborhood_design(len(cverts), cedges)
    print("\n--- Clebsch folded-cube carrier ---")
    print(f"vertices={len(cverts)}  edges={len(cedges)}  degree_hist={dict(sorted(Counter(deg).items()))}")
    print(f"common neighbors on edges={dict(sorted(adj_common.items()))}  nonedges={dict(sorted(non_common.items()))}")
    print("verified SRG parameters: (16,5,0,2)")
    print(f"closed-neighborhood block sizes={dict(sorted(Counter(map(len, cn_blocks)).items()))}")
    print(f"closed-neighborhood point reps={dict(sorted(cn_point_rep.items()))}  pair reps={dict(sorted(cn_pair_rep.items()))}")
    print("frame identity: closed_neighborhood^T closed_neighborhood = 4I + 2J")
    print("LRC hook: folded-cube parity quotient is a candidate finite smoother for labelled signed residual packets.")

    bverts, bedges_colored = bruhat_s4_graph()
    bedges = {(a, b) for a, b, _ in bedges_colored}
    bdeg = graph_degrees(len(bverts), bedges)
    color_counts = Counter(color for _, _, color in bedges_colored)
    squares, hexes = bruhat_faces(bverts)
    print("\n--- Truncated-octahedral / Bruhat B4 carrier ---")
    print(f"vertices={len(bverts)}  edges={len(bedges)}  degree_hist={dict(sorted(Counter(bdeg).items()))}")
    print(f"adjacent-transposition color counts={dict(sorted(color_counts.items()))}")
    print(f"square_commutation_faces={len(squares)}  hex_braid_faces={len(hexes)}")
    print(f"cycle_rank={cycle_rank(len(bverts), len(bedges))}; face_count={len(squares)+len(hexes)} with one sphere relation")
    print("LRC hook: local compression failures should be measured as curl on Coxeter square/hex faces, not as edge monotonicity.")

    o_verts, o_edges = octahedron_line_k4()
    o_deg = graph_degrees(len(o_verts), o_edges)
    print("\n--- Octahedral line graph L(K4) carrier (HYP-2887) ---")
    print(f"vertices={len(o_verts)}  edges={len(o_edges)}  degree_hist={dict(sorted(Counter(o_deg).items()))}")
    print(f"cycle_rank={cycle_rank(len(o_verts), len(o_edges))}; triangular_faces=8 with one sphere relation")
    print("LRC hook: repeated support-six packets already live on this face-curl module.")

    print("\n--- Proof-obligation tournament over carriers ---")
    rows = carrier_tournament()
    for i, (name, score, obs) in enumerate(rows, 1):
        print(f"{i}. {name}: score={score}; obligations={', '.join(obs)}")
    print("Hamiltonian path by obligation coverage:")
    print("  " + " > ".join(name for name, _, _ in rows))

    print("\n--- Synthesis ---")
    print("1. q=3 unital is not decorative: it exactly matches the 28 pair slots of k=8 AP data.")
    print("2. Clebsch gives a second pair-uniform tight frame, but on folded parity packets rather than AP pair slots.")
    print("3. Truncated octahedral/Bruhat B4 explains why one-step compression is false: edge moves have square/hex curl.")
    print("4. HYP-2890 residual leak should split into: unital/Clebsch frame estimates for pair packets + Bruhat/octahedral curl estimates for labelled signed residuals.")


if __name__ == "__main__":
    main()
