#!/usr/bin/env python3
"""
Fast blue line skeleton at n=7 using optimized canonical form.

Strategy: use score-sequence bucketing + certificate comparison to
avoid full n! permutation search. For n=7, standard canonical form
is too slow (5040 perms * 32768 tournaments).

Instead:
1. Build class database using score-seq + refined invariant (sorted row sums of A^2)
2. Only compute full canonical form for representatives within each refined bucket
3. Use incremental matching for the rest

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations
from collections import defaultdict
import time

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

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def refined_invariant(A, n):
    """Score sequence + sorted A^2 row sums + sorted A^3 row sums."""
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    # A^2
    A2 = [[sum(A[i][k]*A[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    a2_scores = tuple(sorted(sum(A2[i]) for i in range(n)))
    # A^3 diagonal (number of directed 3-cycles through each vertex)
    a3_diag = tuple(sorted(A2[i][j]*A[j][i] for i in range(n) for j in range(n) if i != j))
    return (scores, a2_scores, a3_diag)

def canonical(A, n):
    """Full canonical form via permutation search."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def canonical_fast(A, n, known_canons=None):
    """Canonical form with early exit if we find a match in known_canons."""
    # First compute refined invariant for quick bucketing
    ri = refined_invariant(A, n)

    if known_canons is not None and ri in known_canons:
        # Try matching against known canonical forms with this invariant
        for known_canon in known_canons[ri]:
            # Check if A is isomorphic to this known canon
            for perm in permutations(range(n)):
                form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
                if form == known_canon:
                    return known_canon, ri
            # If we exhausted all perms without match, try next known
        # No match found among known — compute fresh

    c = canonical(A, n)
    return c, ri

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

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

def analyze_n7():
    n = 7
    m = num_tiling_bits(n)
    print(f"BLUE LINE SKELETON at n={n}")
    print(f"  Tiling bits: {m}")
    print(f"  Total tilings: {2**m}")

    pairs, fixed = tiling_transpose_pairs(n)
    gs_dof = len(pairs) + len(fixed)
    print(f"  Transpose pairs: {len(pairs)}, Fixed: {len(fixed)}")
    print(f"  GS DOF: {gs_dof}, GS tilings: {2**gs_dof}")

    # Generate GS tilings
    def gen_gs():
        result = []
        for free_val in range(2**gs_dof):
            bits = 0
            for k, (idx1, idx2) in enumerate(pairs):
                if (free_val >> k) & 1:
                    bits |= (1 << idx1) | (1 << idx2)
            for k, fidx in enumerate(fixed):
                if (free_val >> (len(pairs) + k)) & 1:
                    bits |= (1 << fidx)
            result.append(bits)
        return result

    gs_tilings = gen_gs()
    print(f"  Generated {len(gs_tilings)} GS tilings")

    # Phase 1: Build class database from ALL tilings using refined invariant
    # Group by refined invariant, then compute canonical forms within groups
    t0 = time.time()
    print(f"\n  Phase 1: Building class database from {2**m} tilings...")

    ri_groups = defaultdict(list)  # refined_invariant -> list of (bits, A)
    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        ri = refined_invariant(A, n)
        ri_groups[ri].append(bits)

    t1 = time.time()
    print(f"    Refined invariant grouping: {t1-t0:.1f}s, {len(ri_groups)} groups")

    # Now compute canonical forms within each group
    canon_db = {}  # canonical -> class_idx
    class_list = []
    bits_to_class = {}  # bits -> class_idx

    for ri, bits_list in ri_groups.items():
        # For this group, compute canonical forms
        group_canons = {}  # canonical -> class_idx (local)

        for bits in bits_list:
            A = tournament_from_tiling(n, bits)

            # Try to match against known canons in this group
            matched = False
            for known_c, known_idx in group_canons.items():
                # Quick check: try to find a perm that maps A to known_c
                for perm in permutations(range(n)):
                    form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
                    if form == known_c:
                        bits_to_class[bits] = known_idx
                        class_list[known_idx]['tilings'].add(bits)
                        matched = True
                        break
                if matched:
                    break

            if not matched:
                # Compute fresh canonical form
                c = canonical(A, n)
                if c in canon_db:
                    idx = canon_db[c]
                    bits_to_class[bits] = idx
                    class_list[idx]['tilings'].add(bits)
                    group_canons[c] = idx
                else:
                    idx = len(class_list)
                    canon_db[c] = idx
                    sc = is_SC(A, n)
                    class_list.append({
                        'canon': c, 'rep': A, 'tilings': {bits}, 'sc': sc,
                        'scores': score_seq(A, n), 'gs_tilings': set()
                    })
                    bits_to_class[bits] = idx
                    group_canons[c] = idx

    t2 = time.time()
    print(f"    Canonical form computation: {t2-t1:.1f}s")

    num_classes = len(class_list)
    num_sc = sum(1 for c in class_list if c['sc'])
    sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
    nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]
    print(f"    {num_classes} classes, {num_sc} SC, {len(nsc_indices)} NSC")

    # Mark GS tilings
    for bits in gs_tilings:
        idx = bits_to_class[bits]
        class_list[idx]['gs_tilings'].add(bits)

    # Phase 2: GS flip graph
    print(f"\n  Phase 2: GS flip graph...")
    gs_edges = defaultdict(int)
    same_class = 0
    cross_class = 0

    for bits in gs_tilings:
        c_from = bits_to_class[bits]
        flipped = flip_tiling(bits, m)
        c_to = bits_to_class[flipped]
        if c_from == c_to:
            same_class += 1
        else:
            cross_class += 1
            edge = (min(c_from, c_to), max(c_from, c_to))
            gs_edges[edge] += 1

    print(f"    Same-class: {same_class}, Cross-class: {cross_class}")
    if same_class + cross_class > 0:
        print(f"    Cross fraction: {cross_class/(same_class+cross_class):.1%}")

    print(f"\n    Blue edges:")
    for (i, j), w in sorted(gs_edges.items()):
        print(f"      {i}(scores={class_list[i]['scores']}) <-[{w}]-> {j}(scores={class_list[j]['scores']})")

    # Skeleton topology
    degree = defaultdict(int)
    sc_in_skel = set()
    for (i, j) in gs_edges:
        degree[i] += 1
        degree[j] += 1
        sc_in_skel.add(i)
        sc_in_skel.add(j)

    nv = len(sc_in_skel)
    ne = len(gs_edges)
    print(f"\n    Skeleton: {nv} vertices, {ne} edges")
    print(f"    Degree seq: {sorted(degree[i] for i in sc_in_skel)}")
    if ne == nv - 1:
        print(f"    ** TREE **")
    elif ne == nv:
        print(f"    ** UNICYCLIC **")
    else:
        print(f"    ** Genus {ne - nv + 1} **")

    # Phase 3: NSC sidedness
    print(f"\n  Phase 3: NSC sidedness analysis...")

    seen_pairs = set()
    one_sided = 0
    two_sided = 0

    for i in nsc_indices:
        A_op = converse(class_list[i]['rep'], n)
        c_op = canonical(A_op, n)
        partner = canon_db[c_op]

        pair_key = (min(i, partner), max(i, partner))
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        # Flip targets for class i
        targets_i = defaultdict(int)
        for bits in class_list[i]['tilings']:
            flipped = flip_tiling(bits, m)
            c_to = bits_to_class[flipped]
            targets_i[c_to] += 1

        sc_i = {k: v for k, v in targets_i.items() if class_list[k]['sc']}
        nsc_i = {k: v for k, v in targets_i.items() if not class_list[k]['sc']}

        # Flip targets for partner
        targets_p = defaultdict(int)
        for bits in class_list[partner]['tilings']:
            flipped = flip_tiling(bits, m)
            c_to = bits_to_class[flipped]
            targets_p[c_to] += 1

        sc_p = {k: v for k, v in targets_p.items() if class_list[k]['sc']}

        only_i = set(sc_i.keys()) - set(sc_p.keys())
        only_p = set(sc_p.keys()) - set(sc_i.keys())
        shared = set(sc_i.keys()) & set(sc_p.keys())

        is_one_sided = bool(only_i or only_p)
        if is_one_sided:
            one_sided += 1
        else:
            two_sided += 1

        tag = "ONE-SIDED" if is_one_sided else "two-sided"
        print(f"    Pair ({i},{partner}) scores={class_list[i]['scores']}: "
              f"SC_i={set(sc_i.keys())}, SC_p={set(sc_p.keys())} -> {tag}")
        if nsc_i:
            print(f"      NSC targets from {i}: {dict(nsc_i)}")

    print(f"\n  SIDEDNESS SUMMARY:")
    print(f"    One-sided: {one_sided}, Two-sided: {two_sided}")

    t3 = time.time()
    print(f"\n  Total time: {t3-t0:.1f}s")

analyze_n7()
