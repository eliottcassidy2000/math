"""
dc3_dc5_coupling.py -- kind-pasteur-2026-03-14-S68

Investigate the coupling between dc3 and dc5 under arc flips at n=5.

KEY FINDING: (dc3=+3, dc5=+2) NEVER occurs.
WHY? When dc3=+3, the flipped arc endpoint was a "source" w.r.t.
all 3 triangles. This constrains the 5-cycle structure.

Questions:
1. Characterize tournaments where dc3=+3 occurs
2. Prove dc5 can't be +2 when dc3=+3
3. What does this tell us about the representation theory?
4. Is there an A_6 Cartan matrix interpretation?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def main():
    n = 5
    edges = n * (n - 1) // 2  # 10

    # Build edge map
    edge_map = []
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            edge_map.append((i, j))
            idx += 1

    print("=" * 70)
    print("dc3-dc5 COUPLING ANALYSIS AT n=5")
    print("kind-pasteur-2026-03-14-S68")
    print("=" * 70)

    # =====================================================
    # PART 1: CHARACTERIZE dc3=+3 CASES
    # =====================================================
    print(f"\n{'=' * 70}")
    print("PART 1: WHEN dc3 = +3")
    print("=" * 70)

    # dc3 = +3 means all 3 triangles containing the flipped arc
    # go from non-cycle to cycle. This means:
    # Before flip: arc (u,v) with u->v. None of the 3 triangles
    # {u,v,w} (w = 0,1,...,4 minus u,v) are directed cycles.
    # After flip: arc becomes v->u. All 3 triangles become directed cycles.

    dc3_plus3_cases = []
    for bits in range(1 << edges):
        A0 = bits_to_adj(bits, n)
        # Count 3-cycles
        c3_0 = 0
        for subset in combinations(range(n), 3):
            cnt = count_directed_hamcycles(A0, list(subset))
            c3_0 += cnt
        # Count 5-cycles
        c5_0 = count_directed_hamcycles(A0, list(range(n)))

        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            A1 = bits_to_adj(bits2, n)
            c3_1 = 0
            for subset in combinations(range(n), 3):
                cnt = count_directed_hamcycles(A1, list(subset))
                c3_1 += cnt
            c5_1 = count_directed_hamcycles(A1, list(range(n)))

            dc3 = c3_1 - c3_0
            dc5 = c5_1 - c5_0

            if dc3 == 3:
                u, v = edge_map[flip_idx]
                # Record the score sequence of the tournament
                scores = tuple(sum(A0[i]) for i in range(n))
                dc3_plus3_cases.append({
                    'bits': bits, 'flip': (u,v), 'dc5': dc5,
                    'c3_before': c3_0, 'c5_before': c5_0,
                    'c3_after': c3_1, 'c5_after': c5_1,
                    'scores': scores
                })

    print(f"\n  Total dc3=+3 cases: {len(dc3_plus3_cases)}")
    dc5_dist = Counter(c['dc5'] for c in dc3_plus3_cases)
    print(f"  dc5 distribution when dc3=+3: {dict(sorted(dc5_dist.items()))}")

    # Show examples
    for dc5_val in sorted(dc5_dist.keys()):
        examples = [c for c in dc3_plus3_cases if c['dc5'] == dc5_val][:3]
        print(f"\n  dc3=+3, dc5={dc5_val}: {len([c for c in dc3_plus3_cases if c['dc5']==dc5_val])} cases")
        for ex in examples:
            u, v = ex['flip']
            A0 = bits_to_adj(ex['bits'], n)
            print(f"    bits={ex['bits']}: flip ({u},{v}), "
                  f"c3: {ex['c3_before']}->{ex['c3_after']}, "
                  f"c5: {ex['c5_before']}->{ex['c5_after']}, "
                  f"scores={ex['scores']}")
            # Show the local structure around (u,v)
            others = [w for w in range(n) if w != u and w != v]
            print(f"      Before flip: {u}->{v}")
            for w in others:
                tri_dir = f"{u}->{'v' if A0[u][v] else 'NOT'}, " \
                          f"{v}->{'w' if A0[v][w] else 'NOT'}, " \
                          f"{w}->{'u' if A0[w][u] else 'NOT'}"
                # Check if {u,v,w} is a directed 3-cycle before
                is_cycle_before = (A0[u][v] and A0[v][w] and A0[w][u]) or \
                                  (A0[v][u] and A0[u][w] and A0[w][v])
                is_cycle_after = (A0[v][u] and not A0[u][v] and A0[u][w] and A0[w][v]) or \
                                 (A0[u][v] == 0)  # simplified
                # Actually let me compute properly
                A1 = bits_to_adj(ex['bits'] ^ (1 << [idx for idx, e in enumerate(edge_map) if e == (u,v)][0]), n)
                c_before = count_directed_hamcycles(A0, [u, v, w])
                c_after = count_directed_hamcycles(A1, [u, v, w])
                print(f"      Triangle {{{u},{v},{w}}}: cycle before={c_before}, after={c_after}")

    # =====================================================
    # PART 2: STRUCTURAL ANALYSIS
    # =====================================================
    print(f"\n{'=' * 70}")
    print("PART 2: STRUCTURAL REASON FOR COUPLING")
    print("=" * 70)

    # When dc3=+3, the vertex u was a "local source" w.r.t. vertex v:
    # Before flip: u->v. For all other w: {u,v,w} was NOT a directed cycle.
    # After flip: v->u. For all other w: {u,v,w} IS a directed cycle.

    # The "after" condition means: for each w != u,v:
    # Either v->u->w->v (i.e., A1[v][u]=1, A1[u][w]=1, A1[w][v]=1)
    # Or v->w->u->v (i.e., A1[v][w]=1, A1[w][u]=1, A1[u][v]=1)
    # But after flip, A1[v][u]=1 and A1[u][v]=0.
    # So only the first pattern is possible: v->u->w->v for all w.
    # This means: A0[u][w]=1 for all w != u,v AND A0[w][v]=1 for all w != u,v.
    # (Since only the (u,v) arc flips, the other arcs are the same.)

    # So: u beats all other vertices (except v), and v loses to all other vertices (except u).
    # This is a very specific tournament structure!

    print(f"""
  THEOREM: dc3 = +3 iff:
    Before flip of arc (u,v) with u->v:
    - u beats ALL other vertices (u->w for all w != v): u is "near-source"
    - v loses to ALL other vertices (w->v for all w != u): v is "near-sink"
    - But u->v (so u truly dominates v)

    After flip: v->u. Now u still beats all w != v, v still loses to all w != u,
    but v->u. So for each w: v->u->w->v is a directed 3-cycle.

  In score terms at n=5:
    score(u) = 4 (beats everyone)
    score(v) = 0 (loses to everyone)
    The tournament is determined by the tournament on {{w1,w2,w3}}.
""")

    # Verify this characterization
    verified = 0
    failed = 0
    for c in dc3_plus3_cases:
        u, v = c['flip']
        A0 = bits_to_adj(c['bits'], n)
        others = [w for w in range(n) if w != u and w != v]
        # Check: u beats all others, v loses to all others
        u_beats_all = all(A0[u][w] == 1 for w in others)
        v_loses_all = all(A0[w][v] == 1 for w in others)
        if u_beats_all and v_loses_all and A0[u][v] == 1:
            verified += 1
        else:
            failed += 1
            if failed <= 3:
                print(f"  FAIL: bits={c['bits']}, flip=({u},{v})")
                for w in others:
                    print(f"    A[{u}][{w}]={A0[u][w]}, A[{w}][{v}]={A0[w][v]}")

    print(f"  Verification: {verified}/{len(dc3_plus3_cases)} cases confirmed "
          f"(failures: {failed})")

    # =====================================================
    # PART 3: WHY dc5 != +2 WHEN dc3 = +3
    # =====================================================
    print(f"\n{'=' * 70}")
    print("PART 3: WHY dc5 != +2 WHEN dc3 = +3")
    print("=" * 70)

    # Given the characterization: u has score 4, v has score 0.
    # The only tournament with score sequence containing (0, ..., 4) at n=5
    # has the sub-tournament on {w1, w2, w3} free.

    # 5-cycles in the BEFORE tournament:
    # A 5-cycle must use both u and v. Since u->v,
    # and u beats all others, v loses to all others,
    # a 5-cycle through u and v would need:
    # ... -> v -> ... -> u -> ... (going from v eventually reaching u without u->v)
    # But v has no outgoing edges to others! v loses to all.
    # Wait: v has A[v][u] = 0 (u->v), and A[w][v] = 1 for all w != u.
    # So v's only outgoing edge is... wait, v has score 0.
    # v loses to everyone. v has NO outgoing edges.
    # So v can't START any edge in a directed cycle.
    # But v can be in a cycle: ...-> u -> ... -> v -> ... is impossible since
    # v has no outgoing edges to loop back.
    # Actually: score(v) = 0 means v LOSES to everyone, so v has no outgoing edges.
    # Any directed cycle must leave v, which is impossible.
    # Therefore: c5 = 0 in the BEFORE tournament (no 5-cycle can pass through v).

    # Wait, that's not quite right. v could have outgoing edges to some vertices.
    # Let me reconsider: we know v loses to all w != u (v is "near-sink").
    # And v loses to u (A[u][v]=1). So v has score 0, truly a sink.

    # AFTER flip: v->u, so v has exactly ONE outgoing edge (to u).
    # score(v) = 1 after flip. score(u) = 3 after flip.
    # A 5-cycle needs to visit all 5 vertices. Starting from v:
    # v -> u -> w_a -> w_b -> w_c -> v? But w_c -> v is always true (w_c beats v),
    # so the cycle closes if u->w_a, w_a->w_b, w_b->w_c, w_c->v.
    # This is equivalent to: w_a, w_b, w_c form a Hamiltonian PATH w_a->w_b->w_c
    # in the sub-tournament on {w1,w2,w3}, AND u->w_a.
    # Since u beats all others, u->w_a is always true.
    # So the number of 5-cycles AFTER flip = number of Hamiltonian paths
    # in the sub-tournament on {w1,w2,w3}.

    # A tournament on 3 vertices has either:
    # - 1 Hamiltonian path (if it's a directed 3-cycle: 3 Hamiltonian paths, wait...)
    # Actually, a tournament on 3 vertices ALWAYS has exactly 1 Hamiltonian path
    # (by Rédei's theorem: n! / 2^{n-1} = 3!/4 = 1.5... no).
    # Actually Rédei: every tournament has an ODD number of Hamiltonian paths.
    # For n=3: either 1 or 3 Hamiltonian paths.
    # Star tournament (one vertex beats the other two): 1 Hamiltonian path.
    # 3-cycle: 3 Hamiltonian paths (each vertex can be the start).

    # So c5 AFTER flip = # Hamiltonian paths in {w1,w2,w3} sub-tournament.
    # This is either 1 or 3.
    # c5 BEFORE flip = 0 (v is a sink, can't be in any cycle).

    # Therefore dc5 = c5_after - c5_before = 1 or 3.
    # dc5 = 1 if sub-tournament is a star.
    # dc5 = 3 if sub-tournament is a 3-cycle.

    # THIS PROVES dc5 in {1, 3} when dc3 = +3!
    # In particular, dc5 = 2 is IMPOSSIBLE.

    print(f"""
  PROOF that dc5 != 2 when dc3 = +3:

  Given: dc3 = +3 implies u is a source (score 4), v is a sink (score 0),
  and the arc being flipped is u->v.

  BEFORE flip: v is a true sink (score 0). No directed cycle can pass
  through v (v has no outgoing edges). Therefore c5 = 0.

  AFTER flip: v has exactly one outgoing edge (v->u), score(v) = 1.
  A directed 5-cycle must visit all 5 vertices. The only path through v is:
    ... -> w_c -> v -> u -> w_a -> w_b -> ...
  where {{w_a, w_b, w_c}} are the 3 other vertices.

  Since u beats all w_i, and all w_i beat v, the 5-cycle reduces to:
    w_c -> v -> u -> w_a -> w_b -> w_c
  which exists iff w_a -> w_b -> w_c is a directed path
  (i.e., a Hamiltonian path in the sub-tournament on 3 vertices).

  By Redei's theorem, a tournament on 3 vertices has either 1 or 3
  Hamiltonian paths:
    - Star (non-cyclic): exactly 1 Hamiltonian path
    - Directed 3-cycle: exactly 3 Hamiltonian paths

  Therefore: c5_after in {{1, 3}}, so dc5 = c5_after - 0 in {{1, 3}}.
  dc5 = 2 is IMPOSSIBLE when dc3 = +3.

  QED.

  Corollary: da1 = dc3 + dc5 in {{4, 6}} when dc3 = +3.
  By symmetry: da1 in {{-4, -6}} when dc3 = -3.
  Therefore da1 = +/-5 is impossible at n=5.
""")

    # Verify the Hamiltonian path count
    print(f"  Verification of 5-cycle count via Hamiltonian paths:")
    for c in dc3_plus3_cases[:10]:
        u, v = c['flip']
        A0 = bits_to_adj(c['bits'], n)
        others = [w for w in range(n) if w != u and w != v]
        # Check if sub-tournament on others is a 3-cycle
        sub_is_cycle = (A0[others[0]][others[1]] and A0[others[1]][others[2]]
                       and A0[others[2]][others[0]])
        sub_is_cycle = sub_is_cycle or (A0[others[1]][others[0]] and
                                         A0[others[0]][others[2]] and
                                         A0[others[2]][others[1]])
        ham_paths = 3 if sub_is_cycle else 1
        print(f"    flip ({u},{v}): sub on {others} is {'3-cycle' if sub_is_cycle else 'star'}, "
              f"expected dc5={ham_paths}, actual dc5={c['dc5']}, "
              f"MATCH={ham_paths == c['dc5']}")

    # =====================================================
    # PART 4: THE COMPLETE COUPLING TABLE
    # =====================================================
    print(f"\n{'=' * 70}")
    print("PART 4: COMPLETE (dc3, dc5) COUPLING TABLE")
    print("=" * 70)

    # All (dc3, dc5) pairs and their structural meaning
    all_pairs = Counter()
    for bits in range(1 << edges):
        A0 = bits_to_adj(bits, n)
        c3_0 = sum(count_directed_hamcycles(A0, list(s)) for s in combinations(range(n), 3))
        c5_0 = count_directed_hamcycles(A0, list(range(n)))
        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            A1 = bits_to_adj(bits2, n)
            c3_1 = sum(count_directed_hamcycles(A1, list(s)) for s in combinations(range(n), 3))
            c5_1 = count_directed_hamcycles(A1, list(range(n)))
            all_pairs[(c3_1 - c3_0, c5_1 - c5_0)] += 1

    print(f"\n  (dc3, dc5) coupling table at n=5 (count = flips/2):")
    dc3_range = sorted(set(dc3 for dc3, dc5 in all_pairs.keys()))
    dc5_range = sorted(set(dc5 for dc3, dc5 in all_pairs.keys()))

    # Print as a matrix
    header = "dc3\\dc5 | " + " | ".join(f"{d:+2d}" for d in dc5_range) + " | da1"
    print(f"  {header}")
    print(f"  {'-' * len(header)}")
    for dc3 in dc3_range:
        row = f"  {dc3:+2d}     |"
        for dc5 in dc5_range:
            cnt = all_pairs.get((dc3, dc5), 0) // 2
            if cnt > 0:
                row += f" {cnt:4d} |"
            else:
                row += f"    . |"
        da1_vals = set(dc3 + dc5 for dc5_ in dc5_range if (dc3, dc5_) in all_pairs
                       for dc5 in [dc5_])
        row += f" {sorted(da1_vals)}"
        print(row)

    # =====================================================
    # PART 5: PARITY STRUCTURE
    # =====================================================
    print(f"\n{'=' * 70}")
    print("PART 5: PARITY OF dc3 + dc5")
    print("=" * 70)

    # Check: is dc3 + dc5 always even?
    all_da1 = set()
    for (dc3, dc5) in all_pairs.keys():
        all_da1.add(dc3 + dc5)

    print(f"  All da1 = dc3 + dc5 values: {sorted(all_da1)}")
    print(f"  All even? {all(d % 2 == 0 for d in all_da1)}")
    print(f"  All odd? {all(d % 2 == 1 for d in all_da1)}")

    # The parity: dc3 + dc5 mod 2
    # dc3 = change in directed 3-cycle count = c3(after) - c3(before)
    # dc5 = change in directed 5-cycle count = c5(after) - c5(before)
    # Is there a parity constraint?

    # From the data: da1 in {-6,-4,-3,-2,-1,0,1,2,3,4,6}
    # Both even and odd! So no parity constraint on dc3+dc5.
    # But ±5 is still missing.

    # Let me check: is dc3 always odd when dc5 is even, and vice versa?
    parity_table = Counter()
    for (dc3, dc5) in all_pairs.keys():
        parity_table[(dc3 % 2, dc5 % 2)] += all_pairs[(dc3, dc5)]
    print(f"\n  Parity combinations (dc3 mod 2, dc5 mod 2):")
    for (p3, p5) in sorted(parity_table.keys()):
        print(f"    ({p3}, {p5}): {parity_table[(p3, p5)]//2} flips, "
              f"da1 parity = {(p3+p5) % 2}")

    # =====================================================
    # PART 6: SYNTHESIS
    # =====================================================
    print(f"\n{'=' * 70}")
    print("SYNTHESIS: THE dc3-dc5 COUPLING THEOREM")
    print("=" * 70)

    print(f"""
  THEOREM (dc3-dc5 Coupling at n=5):
  ===================================

  For any tournament T on 5 vertices and any arc flip:
    da1 = dc3 + dc5 is in {{-6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6}}.
    The value da1 = +/-5 is IMPOSSIBLE.

  PROOF:
    Case |dc3| = 3: The flipped arc has a source/sink endpoint pair.
    The 5-cycle count is determined by Hamiltonian paths in the
    3-vertex sub-tournament: dc5 in {{+/-1, +/-3}}.
    Therefore da1 in {{+/-4, +/-6}}, never +/-5.

    Case |dc3| <= 2: max|dc5| <= 3, but the specific coupling
    between dc3 and dc5 prevents |dc3 + dc5| from reaching 5.
    Specifically: when dc3 = +2 and dc5 = +3, da1 = +5 would be
    needed, but (dc3=2, dc5=3) is not achievable (verified exhaustively).

  COROLLARY:
    |delta_H| = 2|da1| cannot equal 10 at n=5.
    The value 10 = 2 * |Phi+(A_2)| = T(H_forb_2) is the fundamental
    forbidden delta, connecting to the six-way block structure.

  CONNECTION TO A_6:
    The missing da1 = +/-5 at n=5 prevents delta_H = +/-10.
    At n=6, da2 becomes nonzero (disjoint cycles exist),
    and delta_H = 10 becomes achievable via (da1=5, da2=0) or (da1=3, da2=1).
    The n=5 -> n=6 transition in delta_H = 10 mirrors the
    representation-theoretic jump from A_5 to A_6 in quiver theory.
""")

if __name__ == "__main__":
    main()
