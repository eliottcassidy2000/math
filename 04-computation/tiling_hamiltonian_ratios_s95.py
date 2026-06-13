#!/usr/bin/env python3
"""
tiling_hamiltonian_ratios_s95.py

Exact small-n investigation of the tournament-tiling-explorer fiber:

  fixed base Hamiltonian path P0 = 0 -> 1 -> ... -> n-1
  free tiles = all nonconsecutive arcs, m = C(n-1,2)
  tiling = one tournament containing P0 as a Hamiltonian path

For an isomorphism class C, let
  H(C) = Hamiltonian path count of any representative
  Aut(C) = automorphism group size
  F(C) = number of fixed-P0 tilings landing in C

The central theorem is:

  F(C) = H(C) / Aut(C)

This script verifies it exactly, then explores ratio spectra and pair sums.
"""

from collections import Counter, defaultdict
from itertools import permutations
from math import comb, factorial


def tournament_from_tiling(n, bits):
    A = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        A[i][i + 1] = 1
    idx = 0
    for i in range(n):
        for j in range(i + 2, n):
            if (bits >> idx) & 1:
                A[j][i] = 1
            else:
                A[i][j] = 1
            idx += 1
    return A


def score_partition_perms(A):
    n = len(A)
    scores = [sum(row) for row in A]
    blocks = defaultdict(list)
    for v, s in enumerate(scores):
        blocks[s].append(v)
    groups = [blocks[s] for s in sorted(blocks)]

    def rec(k):
        if k == len(groups):
            yield []
            return
        for p in permutations(groups[k]):
            for rest in rec(k + 1):
                yield list(p) + rest

    yield from rec(0)


def canonical(A):
    n = len(A)
    best = None
    for p in score_partition_perms(A):
        form = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best


def aut_size(A):
    n = len(A)
    scores = [sum(row) for row in A]
    blocks = defaultdict(list)
    for v, s in enumerate(scores):
        blocks[s].append(v)
    groups = list(blocks.values())
    total = 0
    p = [None] * n

    def rec(k):
        if k == len(groups):
            yield p[:]
            return
        group = groups[k]
        for image_group in permutations(group):
            for src, image in zip(group, image_group):
                p[src] = image
            yield from rec(k + 1)

    for p in rec(0):
        ok = True
        for i in range(n):
            for j in range(n):
                if A[p[i]][p[j]] != A[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += 1
    return total


def hamiltonian_paths(A):
    n = len(A)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if A[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[full])


def c3_count(A):
    n = len(A)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c += 1
    return c


def complement_canon(A):
    n = len(A)
    B = [[0 if i == j else 1 - A[i][j] for j in range(n)] for i in range(n)]
    return canonical(B)


def analyze(n):
    m = comb(n - 1, 2)
    total_tilings = 1 << m
    classes = {}
    tiling_class = []

    for bits in range(total_tilings):
        A = tournament_from_tiling(n, bits)
        can = canonical(A)
        if can not in classes:
            classes[can] = {
                "rep": A,
                "F": 0,
                "H": hamiltonian_paths(A),
                "aut": aut_size(A),
                "c3": c3_count(A),
            }
        classes[can]["F"] += 1
        tiling_class.append(can)

    cans = sorted(classes, key=lambda c: (classes[c]["H"], classes[c]["aut"], c))
    cid = {can: i for i, can in enumerate(cans)}
    for can in cans:
        info = classes[can]
        info["cid"] = cid[can]
        info["comp_cid"] = cid[complement_canon(info["rep"])]
        info["self_comp"] = info["comp_cid"] == info["cid"]

    mismatches = []
    for can in cans:
        info = classes[can]
        if info["F"] * info["aut"] != info["H"]:
            mismatches.append((info["cid"], info["F"], info["aut"], info["H"]))

    # Cube-edge quotient statistics on fixed-path tiling space.
    directed_edge_counts = Counter()
    unordered_cross_edges = Counter()
    silent_edges = 0
    for bits, can in enumerate(tiling_class):
        i = cid[can]
        for k in range(m):
            nb = bits ^ (1 << k)
            j = cid[tiling_class[nb]]
            directed_edge_counts[(i, j)] += 1
            if bits < nb:
                if i == j:
                    silent_edges += 1
                else:
                    unordered_cross_edges[tuple(sorted((i, j)))] += 1

    ratio_spectrum = Counter()
    h_fibers = defaultdict(list)
    for can in cans:
        info = classes[can]
        ratio_spectrum[f"1/{info['aut']}"] += 1
        h_fibers[info["H"]].append(info)

    sum_F = sum(classes[c]["F"] for c in cans)
    sum_H = sum(classes[c]["H"] for c in cans)
    sum_FH = sum(classes[c]["F"] * classes[c]["H"] for c in cans)
    sum_F2 = sum(classes[c]["F"] ** 2 for c in cans)
    same_unordered_pairs = sum(info["F"] * (info["F"] - 1) // 2 for info in classes.values())
    automorphic_defect = sum(info["H"] - info["F"] for info in classes.values())

    edge_tiling_energy = 0
    edge_h_energy = 0
    for (i, j), w in unordered_cross_edges.items():
        Fi, Fj = classes[cans[i]]["F"], classes[cans[j]]["F"]
        Hi, Hj = classes[cans[i]]["H"], classes[cans[j]]["H"]
        edge_tiling_energy += w * Fi * Fj
        edge_h_energy += w * Hi * Hj

    # H-collision table: same H can split into several automorphism/fiber types.
    collisions = []
    for H, infos in sorted(h_fibers.items()):
        if len(infos) > 1:
            collisions.append((
                H,
                len(infos),
                sum(x["F"] for x in infos),
                sorted((x["F"], x["aut"], x["c3"], x["self_comp"]) for x in infos),
            ))

    return {
        "n": n,
        "m": m,
        "classes": classes,
        "cans": cans,
        "total_tilings": total_tilings,
        "mismatches": mismatches,
        "ratio_spectrum": ratio_spectrum,
        "sum_F": sum_F,
        "sum_H": sum_H,
        "sum_FH": sum_FH,
        "sum_F2": sum_F2,
        "same_unordered_pairs": same_unordered_pairs,
        "automorphic_defect": automorphic_defect,
        "silent_edges": silent_edges,
        "cross_edge_weight": sum(unordered_cross_edges.values()),
        "cube_edges": m * (1 << (m - 1)) if m else 0,
        "edge_tiling_energy": edge_tiling_energy,
        "edge_h_energy": edge_h_energy,
        "collisions": collisions,
    }


def main():
    print("=" * 78)
    print("TILING-HAMILTONIAN RATIOS S95")
    print("=" * 78)
    print()
    print("Theorem under test: F(C) = H(C)/|Aut(C)| for the explorer fixed path.")
    print("Terminology: F(C) is the path-fiber; F/H = 1/|Aut| is the tiling ratio.")

    summaries = []
    for n in range(3, 8):
        data = analyze(n)
        summaries.append(data)
        print("\n" + "-" * 78)
        print(f"n={n}, m=C(n-1,2)={data['m']}, tilings=2^m={data['total_tilings']}, classes={len(data['cans'])}")
        print(f"  theorem mismatches: {len(data['mismatches'])}")
        print(f"  ratio spectrum F/H: {dict(sorted(data['ratio_spectrum'].items()))}")
        print(
            "  sums: "
            f"ΣF={data['sum_F']}, ΣH={data['sum_H']}, "
            f"ΣF·H={data['sum_FH']}, ΣF²={data['sum_F2']}, "
            f"ΣC(F,2)={data['same_unordered_pairs']}"
        )
        print(
            "  cube edges: "
            f"total={data['cube_edges']}, silent={data['silent_edges']}, "
            f"cross={data['cross_edge_weight']}, "
            f"edge_F_energy={data['edge_tiling_energy']}, "
            f"edge_H_energy={data['edge_h_energy']}"
        )
        print(f"  automorphic defect Σ(H-F)=ΣH(1-1/Aut)={data['automorphic_defect']}")
        if data["collisions"]:
            print("  H-collisions (H: class_count, total_F, [(F,Aut,c3,self_comp), ...]):")
            for H, count, total_F, rows in data["collisions"][:12]:
                print(f"    H={H}: {count} classes, total_F={total_F}, {rows}")

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    fields = [
        ("classes", lambda d: len(d["cans"])),
        ("ΣF = total tilings", lambda d: d["sum_F"]),
        ("ΣH over classes", lambda d: d["sum_H"]),
        ("ΣF·H = path-fiber H moment", lambda d: d["sum_FH"]),
        ("ΣF² = ordered same-class tiling pairs", lambda d: d["sum_F2"]),
        ("ΣC(F,2) = unordered same-class tiling pairs", lambda d: d["same_unordered_pairs"]),
        ("silent cube edges", lambda d: d["silent_edges"]),
        ("cross cube edges", lambda d: d["cross_edge_weight"]),
        ("edge tiling energy Σ_edges F_i F_j", lambda d: d["edge_tiling_energy"]),
        ("edge H energy Σ_edges H_i H_j", lambda d: d["edge_h_energy"]),
        ("automorphic defect Σ(H-F)", lambda d: d["automorphic_defect"]),
    ]
    for name, fn in fields:
        print(f"  {name}: {[fn(d) for d in summaries]}")

    print("\n" + "=" * 78)
    print("PROVED STATEMENTS")
    print("=" * 78)
    print("""
1. Fixed-path fiber theorem.
   For every isomorphism class C,
       F(C) = H(C)/|Aut(C)|.
   Proof: count pairs (labeled copy of C, Hamiltonian path in that copy).
   There are (n!/|Aut|)H such pairs. The symmetric group is transitive on
   the n! possible ordered vertex paths, so exactly H/|Aut| copies contain
   the explorer's fixed path P0.

2. Tiling ratio theorem.
       F(C)/H(C) = 1/|Aut(C)|.
   The only possible ratios are reciprocals of odd automorphism orders.
   Generic rigid classes have ratio 1; symmetric classes are compressed.

3. Total path-fiber theorem.
       Σ_C F(C) = Σ_C H(C)/|Aut(C)| = 2^C(n-1,2).
   This is just partitioning the explorer's staircase hypercube by class.

4. Same-class pair theorem.
       ordered same-class tiling pairs = Σ_C F(C)^2 = Σ_C H(C)^2/|Aut(C)|^2.
       unordered same-class pairs = Σ_C C(F(C),2).
   These are natural "second moment" sequences for the explorer quotient.

5. Cube-edge balance theorem.
   If W_ij counts directed single-tile flips from class i to class j, then
       Σ_j W_ij = m F_i
   and W is symmetric. Hence the quotient flip chain is reversible with
   stationary distribution proportional to F_i = H_i/|Aut_i|.
""")


if __name__ == "__main__":
    main()
