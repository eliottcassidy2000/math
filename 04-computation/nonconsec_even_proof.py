#!/usr/bin/env python3
"""
WHY are non-consecutive 3-cycle contributions even for GS tilings?

Hypothesis: the GS constraint pairs non-consecutive triples so that
each pair has even total contribution.

The GS constraint says: for each transpose pair (idx1, idx2) of
non-backbone edges, A[i1][j1] = A[i2][j2] (where edge idx1 is (i1,j1)
and its transpose edge idx2 is (i2,j2)).

For each non-consecutive triple {a,b,c}, check:
1. Does the staircase transpose map it to another triple?
2. If so, do their contributions pair up?

kind-pasteur-2026-03-06-S25h
"""
from itertools import combinations

def tiling_edge_list(n):
    """Non-backbone edges in order."""
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    return edges

def transpose_edge(i, j, n):
    """(i,j) -> (n-1-j, n-1-i)."""
    ti, tj = n-1-j, n-1-i
    if ti > tj:
        ti, tj = tj, ti
    return (ti, tj)

def transpose_vertex(v, n):
    return n - 1 - v

for n in [5, 7]:
    print(f"\nn={n}")

    # Non-consecutive triples (not of form {i, i+1, i+2})
    nonconsec = []
    for a, b, c in combinations(range(n), 3):
        if not (b == a+1 and c == a+2):
            nonconsec.append((a, b, c))

    print(f"  Non-consecutive triples: {len(nonconsec)}")

    # Check: does staircase transpose map triples to triples?
    # Transpose: v -> n-1-v. This reverses the order.
    # Triple {a,b,c} -> {n-1-c, n-1-b, n-1-a}
    # After sorting: same triple type.

    pairs = []
    fixed = []
    seen = set()
    for triple in nonconsec:
        if triple in seen:
            continue
        a, b, c = triple
        ta, tb, tc = n-1-c, n-1-b, n-1-a
        t_triple = tuple(sorted([ta, tb, tc]))
        if t_triple == triple:
            fixed.append(triple)
            seen.add(triple)
        else:
            if t_triple in seen:
                continue
            pairs.append((triple, t_triple))
            seen.add(triple)
            seen.add(t_triple)

    print(f"  Triple pairs: {len(pairs)}, Fixed triples: {len(fixed)}")

    # Are the fixed triples also non-consecutive?
    for t in fixed:
        a, b, c = t
        is_consec = (b == a+1 and c == a+2)
        print(f"    Fixed triple {t}: consecutive={is_consec}")

    # Check: for paired triples, does GS constraint force their
    # t3 contributions to be equal?

    # Under GS: A[i][j] = A[n-1-j][n-1-i] for all non-backbone edges.
    # For a triple {a,b,c}:
    #   3-cycle a->b->c->a: A[a][b]*A[b][c]*A[c][a]
    #   Under transpose: A[n-1-b][n-1-a]*A[n-1-c][n-1-b]*A[n-1-a][n-1-c]
    #   = A[ta][tc]*A[tb][ta]*A[tc][tb]  (where ta=n-1-a, etc.)
    #   This is the 3-cycle tc->tb->ta->tc
    #   Which is a 3-cycle on {ta,tb,tc} in REVERSE direction
    #   = counter-clockwise 3-cycle on the transpose triple

    # So: CW 3-cycle on {a,b,c} under GS <=> CCW 3-cycle on transpose triple
    # And: CCW on {a,b,c} under GS <=> CW on transpose triple

    # For a PAIRED triple (triple != transpose):
    #   t3 contribution of triple = #{3-cycles on {a,b,c}}
    #   t3 contribution of transpose = #{3-cycles on {ta,tb,tc}}
    #   Under GS: CW on original = CCW on transpose
    #   So: total = CW(orig) + CCW(orig) + CW(trans) + CCW(trans)
    #            = CW(orig) + CCW(orig) + CCW(orig) + CW(orig)  [by GS]
    #            = 2*(CW(orig) + CCW(orig))
    #   But CW(orig) + CCW(orig) = 1 (exactly one orientation is a 3-cycle)
    #   So total = 2 for each paired pair!

    # For a FIXED triple ({a,b,c} = its own transpose):
    #   GS says CW(abc) = CCW(abc). But exactly one of CW or CCW holds.
    #   So CW(abc) = CCW(abc) means... both must be 0? No, exactly one holds.
    #   Contradiction! Unless... the fixed triple is consecutive.

    # Wait, let me re-check. The GS constraint says A[i][j] = A[n-1-j][n-1-i]
    # for NON-BACKBONE edges only. If some edge in the triple IS backbone,
    # then GS doesn't constrain it.

    print(f"\n  PROOF ANALYSIS:")
    print(f"  For each paired non-consec triple pair:")
    print(f"    GS => CW on original = CCW on transpose")
    print(f"    total 3-cycles (both triples) = 2 * (CW + CCW on original) = 2 * 1 = 2")
    print(f"    ... but wait, need to check edge types")

    # Check which edges are backbone vs non-backbone for each triple
    for triple in nonconsec[:5]:
        a, b, c = triple
        edges = [(a,b), (b,c), (a,c)]
        backbone = [(x,y) for x,y in edges if abs(x-y) == 1]
        nonbb = [(x,y) for x,y in edges if abs(x-y) >= 2]
        print(f"    Triple {triple}: backbone={backbone}, nonbb={nonbb}")

    # Count backbone edges per non-consecutive triple
    bb_counts = []
    for a, b, c in nonconsec:
        bb = sum(1 for x, y in [(a,b),(b,c),(a,c)] if abs(x-y) == 1)
        bb_counts.append(bb)
    from collections import Counter
    print(f"  Backbone edges per non-consec triple: {Counter(bb_counts)}")

    # Key question: for fixed non-consec triples, how many backbone edges?
    for t in fixed:
        a, b, c = t
        bb = sum(1 for x, y in [(a,b),(b,c),(a,c)] if abs(x-y) == 1)
        print(f"    Fixed non-consec triple {t}: {bb} backbone edges")

    # Now let's verify the pairing argument computationally
    # For GS tilings, check that each pair contributes exactly 2 to total
    edges_list = tiling_edge_list(n)
    edge_to_idx = {e: i for i, e in enumerate(edges_list)}

    gs_pairs, gs_fixed = [], []
    seen2 = set()
    for idx, (i, j) in enumerate(edges_list):
        if idx in seen2:
            continue
        ti, tj = transpose_edge(i, j, n)
        tidx = edge_to_idx.get((ti, tj))
        if tidx == idx:
            gs_fixed.append(idx)
            seen2.add(idx)
        elif tidx is not None:
            gs_pairs.append((idx, tidx))
            seen2.add(idx)
            seen2.add(tidx)

    gs_dof = len(gs_pairs) + len(gs_fixed)

    # Generate a few GS tilings and check pair contributions
    for trial in range(min(8, 2**gs_dof)):
        bits = 0
        for k, (idx1, idx2) in enumerate(gs_pairs):
            if (trial >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(gs_fixed):
            if (trial >> (len(gs_pairs) + k)) & 1:
                bits |= (1 << fidx)

        def make_A(bits_val):
            A = [[0]*n for _ in range(n)]
            for i in range(n-1):
                A[i][i+1] = 1
            for idx, (i, j) in enumerate(edges_list):
                if (bits_val >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
            return A

        A = make_A(bits)
        Af = make_A(bits ^ ((1 << len(edges_list)) - 1))

        total_pair_sum = 0
        for (t1, t2) in pairs:
            a1,b1,c1 = t1
            a2,b2,c2 = t2
            # Count 3-cycles on t1 and t2, both before and after flip
            def count_on(M, triple):
                a,b,c = triple
                cnt = 0
                if M[a][b] and M[b][c] and M[c][a]: cnt += 1
                if M[a][c] and M[c][b] and M[b][a]: cnt += 1
                return cnt
            before_sum = count_on(A, t1) + count_on(A, t2)
            after_sum = count_on(Af, t1) + count_on(Af, t2)
            pair_total = before_sum + after_sum
            total_pair_sum += pair_total

        fixed_sum = 0
        for t in fixed:
            def count_on(M, triple):
                a,b,c = triple
                cnt = 0
                if M[a][b] and M[b][c] and M[c][a]: cnt += 1
                if M[a][c] and M[c][b] and M[b][a]: cnt += 1
                return cnt
            fixed_sum += count_on(A, t) + count_on(Af, t)

        if trial < 4:
            print(f"\n  Trial {trial}: pair_sum={total_pair_sum} (expected: {2*len(pairs)}), fixed_sum={fixed_sum}")
