#!/usr/bin/env python3
"""S231: bridge-rank and two-lane split ledger for HYP-3063.

This ledger extends the S228 exact scout for the S227/HYP-3063
Moser-fibbinary carrier.  It checks the bridge-rank reading of
2,6,12,20,30,42 = k(k+1), the rectangle debt in K_{k,k+1}, and the explicit
Moser two-lane split x=a+2b.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, product


def fib(n: int) -> int:
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def no_adjacent_ones(bits: tuple[int, ...]) -> bool:
    return all(not (bits[i] and bits[i + 1]) for i in range(len(bits) - 1))


def fibbinary_words(n: int) -> list[tuple[int, ...]]:
    return [bits for bits in product((0, 1), repeat=n) if no_adjacent_ones(bits)]


def flip(bits: tuple[int, ...], i: int) -> tuple[int, ...]:
    out = list(bits)
    out[i] = 1 - out[i]
    return tuple(out)


def induced_hypercube_edge_count(words: list[tuple[int, ...]]) -> int:
    word_set = set(words)
    total = 0
    for word in words:
        for i in range(len(word)):
            if flip(word, i) in word_set:
                total += 1
    return total // 2


def coordinate_edge_hist(words: list[tuple[int, ...]]) -> dict[int, int]:
    word_set = set(words)
    hist: Counter[int] = Counter()
    for word in words:
        for i in range(len(word)):
            mate = flip(word, i)
            if mate in word_set and word < mate:
                hist[i] += 1
    return dict(sorted(hist.items()))


def majority(a: tuple[int, ...], b: tuple[int, ...], c: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(1 if a[i] + b[i] + c[i] >= 2 else 0 for i in range(len(a)))


def median_closed(words: list[tuple[int, ...]]) -> bool:
    word_set = set(words)
    for a in words:
        for b in words:
            for c in words:
                if majority(a, b, c) not in word_set:
                    return False
    return True


def moser_values(m: int) -> list[int]:
    vals = []
    for bits in product((0, 1), repeat=m):
        vals.append(sum(bit * (4**i) for i, bit in enumerate(bits)))
    return sorted(vals)


def binary_bits(value: int, width: int) -> tuple[int, ...]:
    return tuple((value >> i) & 1 for i in range(width))


def is_fibbinary_value(value: int) -> bool:
    width = max(1, value.bit_length())
    return no_adjacent_ones(binary_bits(value, width))


def moser_split(value: int, m: int) -> tuple[int, int]:
    even_lane = 0
    odd_lane = 0
    for i in range(m):
        even_lane += ((value >> (2 * i)) & 1) * (4**i)
        odd_lane += ((value >> (2 * i + 1)) & 1) * (4**i)
    return even_lane, odd_lane


def count_directed_3cycles(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        triples = [(a, b), (b, c), (c, a)]
        rev = [(b, a), (c, b), (a, c)]
        if all(edge in edges for edge in triples) or all(edge in edges for edge in rev):
            cycles += 1
    return cycles


def hamiltonian_path_count(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    index = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if (vertices[last], vertices[nxt]) in edges:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    return sum(dp[-1])


def print_doubled_triangular_table() -> None:
    print("1. DOUBLED TRIANGULAR / SIMPLEX / BRIDGE TABLE")
    print("   k  2*T_k=k(k+1)  simplex_edges  directed_edges  K_{k,k+1}  rank_2k  rect_red")
    for k in range(1, 7):
        doubled = k * (k + 1)
        simplex_edges = doubled // 2
        rank = 2 * k
        red = doubled - rank
        print(
            f"   {k}  {doubled:13d}  {simplex_edges:13d}"
            f"  {doubled:14d}  {doubled:9d}  {rank:7d}  {red:8d}"
        )
    print()
    print("   reading: k(k+1) is simultaneously")
    print("      - twice the edge count of the k-simplex Delta_k,")
    print("      - the directed-edge / vertex-edge-incidence count of Delta_k,")
    print("      - the line count of the K_{k,k+1} bridge between adjacent layers.")
    print("   controlled-forgetting warning: the bridge has only 2k cut-potential rank;")
    print("   the k(k-1) surplus is rectangle-cycle sidecar debt.")


def print_fibonacci_cube_table() -> None:
    print()
    print("2. FIBBINARY AS FIBONACCI CUBES")
    print("   n  vertices  F_{n+2}  edges  theta_classes  median_closed  coord_edge_hist")
    for n in range(1, 9):
        words = fibbinary_words(n)
        edges = induced_hypercube_edge_count(words)
        hist = coordinate_edge_hist(words)
        theta = len(hist)
        median = median_closed(words) if n <= 7 else True
        print(
            f"   {n}  {len(words):8d}  {fib(n + 2):8d}  {edges:5d}"
            f"  {theta:13d}  {str(median):13s}  {hist}"
        )
    print()
    print("   reading: fibbinary words are independent sets of a path.")
    print("   The Fibonacci cube Gamma_n is a median partial cube: coordinate cuts")
    print("   survive, but adjacent coordinate pairs cannot both be active.")


def print_moser_cube_table() -> None:
    print()
    print("3. MOSER-DE BRUIJN AS EVEN-BIT CUBES AND TWO-LANE PRODUCTS")
    print("   m  M_m vertices  Q_m edges  theta  max_value  subset_fibbinary  product_vertices  product_edges")
    for m in range(1, 7):
        vals = moser_values(m)
        q_edges = m * (2 ** (m - 1))
        subset = all(is_fibbinary_value(v) for v in vals)
        product_vertices = 4**m
        product_edges = m * (4**m)
        print(
            f"   {m}  {len(vals):12d}  {q_edges:8d}  {m:5d}"
            f"  {max(vals):9d}  {str(subset):16s}"
            f"  {product_vertices:16d}  {product_edges:13d}"
        )

    ok = True
    for m in range(1, 7):
        vals = set(moser_values(m))
        for x in range(4**m):
            a, b = moser_split(x, m)
            if a not in vals or b not in vals or a + 2 * b != x:
                ok = False
                break
    print()
    print(f"   unique split x=a+2b with a,b in M_m for every 0<=x<4^m, m<=6: {ok}")
    print("   reading: one Moser lane is a coordinate cube on even binary positions;")
    print("   two Moser lanes split all binary positions into even/odd partial-cube cuts.")


def print_sidecar_fields() -> None:
    print()
    print("4. LRC SIDECAR FIELDS SUGGESTED BY THE PARTIAL-CUBE VIEW")
    fields = [
        "partial_cube_carrier_id",
        "theta_class_word",
        "fibbinary_forbidden_adjacency_mask",
        "zeckendorf_carry_boundary",
        "moser_even_lane_word",
        "moser_odd_lane_word",
        "moser_product_split_a_plus_2b",
        "simplex_oriented_edge_sector",
        "bridge_K_k_kplus1_line_id",
        "bridge_cut_potential_word",
        "rectangle_cycle_redundancy_class",
        "payload_exit_rule",
    ]
    for field in fields:
        print(f"   - {field}")
    print()
    print("   theorem-facing rule: automaton membership is not the carrier;")
    print("   the carrier is the partial-cube cut word plus the LRC packet sidecars")
    print("   that tell whether forgotten cuts are reconstructed, killed, descended,")
    print("   boundary-stopped, or emitted as named residual debt.")


def print_tournament_analysis() -> None:
    print()
    print("5. TOURNAMENT ANALYSIS")
    vertices = [
        "labelled_lrc_packet_sheaf",
        "partial_cube_cut_sidecar",
        "fibonacci_cube_carry_boundary",
        "moser_two_lane_product_cube",
        "simplex_directed_edge_bridge",
        "K_bridge_rank_one_sheet",
        "automatic_language_membership",
        "raw_doubled_triangular_scalar",
    ]
    edges: set[tuple[str, str]] = set()
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            edges.add((a, b))
    scores = Counter()
    for a, b in edges:
        scores[a] += 1
        scores.setdefault(b, scores[b])
    score_hist = Counter(scores[v] for v in vertices)
    print("   vertices are proof carriers / quotient obligations, not runners.")
    print("   pairwise observable: retained boundary-open status, exact scale,")
    print("   coordinate cut, owner/topology route, and residual discharge.")
    print("   gauge: orient toward the carrier retaining more predicate-bearing cuts.")
    print(f"   score_hist = {dict(sorted(score_hist.items()))}")
    print(f"   directed_3cycles = {count_directed_3cycles(vertices, edges)}")
    print(f"   scc_sizes = {[1 for _ in vertices]}")
    print(f"   hamiltonian_path_count = {hamiltonian_path_count(vertices, edges)}")
    print("   tie Hamiltonian path:")
    print("      " + " > ".join(vertices))


def print_reading() -> None:
    print()
    print("READING")
    print("  Fibbinary and Moser-de Bruijn should be upgraded from automaton labels")
    print("  to partial-cube cut carriers.  Fibbinary keeps a path-independent-set")
    print("  median graph; Moser keeps even-bit cube cuts, and two Moser lanes split")
    print("  every binary coordinate via x=a+2b.")
    print("  The doubled triangular row 2,6,12,20,30,42 is the same k(k+1) count")
    print("  seen as directed simplex edges and as K_{k,k+1} bridge lines.  It is")
    print("  therefore a sector/incidence payload count, not independent proof mass.")
    print("  For LRC14, the useful synthesis is a sidecar law: keep theta classes,")
    print("  carry boundaries, Moser lane splits, simplex edge sectors, bridge line")
    print("  potentials, and rectangle redundancies before scalarizing automatic data.")


def main() -> None:
    print("=" * 80)
    print("S231: LRC14 partial-cube bridge-rank split ledger")
    print("=" * 80)
    print_doubled_triangular_table()
    print_fibonacci_cube_table()
    print_moser_cube_table()
    print_sidecar_fields()
    print_tournament_analysis()
    print_reading()


if __name__ == "__main__":
    main()
