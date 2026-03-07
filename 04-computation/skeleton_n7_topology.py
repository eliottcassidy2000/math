#!/usr/bin/env python3
"""
Blue line skeleton topology at n=5 and n=7.

Analyze:
1. The GS flip graph (blue line skeleton) connecting SC classes
2. How NSC classes connect to the skeleton (one-sided vs two-sided)
3. Graph-theoretic properties: is it a tree? a cycle? Mobius strip?

Uses TILING model (backbone = standard path, flip = complement non-backbone bits).

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations
from collections import defaultdict
import sys

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def converse(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def is_SC(A, n):
    return canonical(A, n) == canonical(converse(A, n), n)

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs = []
    fixed = []
    seen = set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen:
            continue
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx)
            seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx))
            seen.add(idx)
            seen.add(tidx)
    return pairs, fixed

def is_gs_tiling(tiling_bits, n, pairs, fixed):
    for idx1, idx2 in pairs:
        if ((tiling_bits >> idx1) & 1) != ((tiling_bits >> idx2) & 1):
            return False
    return True

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

def analyze(n):
    m = num_tiling_bits(n)
    print(f"\n{'='*70}")
    print(f"BLUE LINE SKELETON TOPOLOGY at n={n}")
    print(f"{'='*70}")
    print(f"  Non-backbone edges (tiling bits): {m}")
    print(f"  Total tilings: {2**m}")

    pairs, fixed = tiling_transpose_pairs(n)
    gs_dof = len(pairs) + len(fixed)
    print(f"  GS degrees of freedom: {gs_dof}")
    print(f"  GS tilings: {2**gs_dof}")

    # Find all GS tilings by iterating over free bits
    def enumerate_gs_tilings():
        """Generate all GS tilings by setting free bits and mirroring."""
        free_indices = [idx for idx, _ in pairs] + fixed
        free_indices.sort()
        result = []
        for free_val in range(2**gs_dof):
            bits = 0
            # Set paired bits
            for k, (idx1, idx2) in enumerate(pairs):
                if (free_val >> k) & 1:
                    bits |= (1 << idx1) | (1 << idx2)
            # Set fixed bits
            for k, fidx in enumerate(fixed):
                if (free_val >> (len(pairs) + k)) & 1:
                    bits |= (1 << fidx)
            result.append(bits)
        return result

    gs_tilings = enumerate_gs_tilings()
    print(f"  Generated {len(gs_tilings)} GS tilings")

    # Build iso class database from ALL tilings
    canon_db = {}
    class_list = []

    print(f"  Building iso class database...", end="", flush=True)
    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        c = canonical(A, n)
        if c not in canon_db:
            canon_db[c] = len(class_list)
            class_list.append({
                'canon': c, 'rep': A, 'tilings': set(), 'sc': is_SC(A, n),
                'scores': score_seq(A, n), 'gs_tilings': set()
            })
        class_list[canon_db[c]]['tilings'].add(bits)

    for bits in gs_tilings:
        A = tournament_from_tiling(n, bits)
        c_idx = canon_db[canonical(A, n)]
        class_list[c_idx]['gs_tilings'].add(bits)

    num_classes = len(class_list)
    num_sc = sum(1 for c in class_list if c['sc'])
    sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
    nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]
    print(f" {num_classes} classes, {num_sc} SC, {len(nsc_indices)} NSC")

    # Build GS flip graph (blue line skeleton)
    gs_edges = defaultdict(int)
    same_class = 0
    cross_class = 0

    for bits in gs_tilings:
        A_from = tournament_from_tiling(n, bits)
        c_from = canon_db[canonical(A_from, n)]
        flipped = flip_tiling(bits, m)
        A_to = tournament_from_tiling(n, flipped)
        c_to = canon_db[canonical(A_to, n)]
        if c_from == c_to:
            same_class += 1
        else:
            cross_class += 1
            edge = (min(c_from, c_to), max(c_from, c_to))
            gs_edges[edge] += 1

    print(f"\n  BLUE LINE SKELETON:")
    print(f"    Same-class flips: {same_class}")
    print(f"    Cross-class flips: {cross_class}")
    if same_class + cross_class > 0:
        print(f"    Cross fraction: {cross_class/(same_class+cross_class):.1%}")

    # Build adjacency for SC classes
    sc_adj = defaultdict(set)
    for (i, j), w in gs_edges.items():
        if class_list[i]['sc'] and class_list[j]['sc']:
            sc_adj[i].add(j)
            sc_adj[j].add(i)

    print(f"\n    Blue edges (SC-SC):")
    for (i, j), w in sorted(gs_edges.items()):
        tag_i = "SC" if class_list[i]['sc'] else "NSC"
        tag_j = "SC" if class_list[j]['sc'] else "NSC"
        print(f"      {i}({tag_i},{class_list[i]['scores']}) <-[{w}]-> {j}({tag_j},{class_list[j]['scores']})")

    # Skeleton graph properties
    sc_in_skeleton = set()
    for (i, j) in gs_edges:
        sc_in_skeleton.add(i)
        sc_in_skeleton.add(j)

    print(f"\n    SC classes in skeleton: {len(sc_in_skeleton)} / {num_sc}")
    print(f"    SC classes NOT in skeleton: {num_sc - len(sc_in_skeleton)}")

    # Degree sequence
    degree = defaultdict(int)
    for (i, j) in gs_edges:
        degree[i] += 1
        degree[j] += 1

    print(f"    Degree sequence: {sorted(degree[i] for i in sc_in_skeleton)}")

    # Check if tree (|E| = |V| - 1)
    num_edges = len(gs_edges)
    num_verts = len(sc_in_skeleton)
    print(f"    Edges: {num_edges}, Vertices: {num_verts}")
    if num_edges == num_verts - 1:
        print(f"    ** TREE (|E| = |V|-1) **")
    elif num_edges == num_verts:
        print(f"    ** UNICYCLIC (|E| = |V|) **")
    else:
        print(f"    ** Genus = {num_edges - num_verts + 1} **")

    # NSC class connections to skeleton
    print(f"\n  NSC CLASS CONNECTIONS:")

    one_sided = 0
    two_sided = 0

    for i in nsc_indices:
        A_op = converse(class_list[i]['rep'], n)
        partner = canon_db[canonical(A_op, n)]

        # Where do tilings of this NSC class go under flip?
        flip_targets = defaultdict(int)
        for bits in class_list[i]['tilings']:
            flipped = flip_tiling(bits, m)
            A_to = tournament_from_tiling(n, flipped)
            c_to = canon_db[canonical(A_to, n)]
            flip_targets[c_to] += 1

        sc_targets = {k: v for k, v in flip_targets.items() if class_list[k]['sc']}
        nsc_targets = {k: v for k, v in flip_targets.items() if not class_list[k]['sc']}

        # Check sidedness: does this NSC class connect to both endpoints of a blue edge?
        sc_target_set = set(sc_targets.keys())

        # Find which blue edges this NSC class "straddles"
        straddles = []
        one_side_of = []
        for (a, b) in gs_edges:
            if a in sc_target_set and b in sc_target_set:
                straddles.append((a, b))
            elif a in sc_target_set or b in sc_target_set:
                side = a if a in sc_target_set else b
                one_side_of.append(((a, b), side))

        if i <= partner:  # Only print each pair once
            print(f"\n    NSC pair ({i}, {partner}), scores={class_list[i]['scores']}:")
            print(f"      Class {i} flip targets: SC={dict(sc_targets)}, NSC={dict(nsc_targets)}")

            # Partner's targets
            partner_targets = defaultdict(int)
            for bits in class_list[partner]['tilings']:
                flipped = flip_tiling(bits, m)
                A_to = tournament_from_tiling(n, flipped)
                c_to = canon_db[canonical(A_to, n)]
                partner_targets[c_to] += 1
            partner_sc = {k: v for k, v in partner_targets.items() if class_list[k]['sc']}
            partner_nsc = {k: v for k, v in partner_targets.items() if not class_list[k]['sc']}
            print(f"      Class {partner} flip targets: SC={dict(partner_sc)}, NSC={dict(partner_nsc)}")

            # Sidedness analysis
            combined_sc = set(sc_targets.keys()) | set(partner_sc.keys())
            shared_sc = set(sc_targets.keys()) & set(partner_sc.keys())
            only_i = set(sc_targets.keys()) - set(partner_sc.keys())
            only_partner = set(partner_sc.keys()) - set(sc_targets.keys())

            print(f"      SC targets: combined={combined_sc}")
            print(f"        shared={shared_sc}, only_{i}={only_i}, only_{partner}={only_partner}")

            if only_i or only_partner:
                one_sided += 1
                print(f"        ** ONE-SIDED: partners prefer different SC classes **")
            else:
                two_sided += 1
                print(f"        (two-sided: partners share all SC targets)")

    print(f"\n  SIDEDNESS SUMMARY:")
    print(f"    One-sided NSC pairs: {one_sided}")
    print(f"    Two-sided NSC pairs: {two_sided}")

    return class_list, gs_edges, sc_indices, nsc_indices

# Run for n=5
analyze(5)

print("\n\nDONE")
