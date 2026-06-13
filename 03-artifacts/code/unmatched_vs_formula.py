#!/usr/bin/env python3
"""
Test: #U_{T'} - #U_T = delta_I formula from THM-013.

Key findings from unmatched_decomposition.py:
- s_x=0 vertices never block (confirmed)
- s_x=-1 blocks T-paths (3-cycle a->i->j->a), s_x=+1 blocks T'-paths
- Per-vertex: uT_x != 2*H(B_x), identity is aggregate only
- Double-blocking occurs 23% of time at n=6

This script tests the aggregate identity and explores the
STRUCTURE of unmatched paths to find a proof approach.

Instance: kind-pasteur-2026-03-05-S6
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
from collections import defaultdict, Counter
import random


def flip_arc(T, i, j):
    n = len(T)
    Tp = [row[:] for row in T]
    Tp[i][j] = 1 - T[i][j]
    Tp[j][i] = 1 - T[j][i]
    return Tp


def compute_I(T):
    cycles = find_odd_cycles(T)
    if not cycles:
        return 1
    cg = conflict_graph(cycles)
    return independence_poly_at(cg, 2)


def ham_paths_list(T):
    """Return all Ham paths as tuples."""
    n = len(T)
    paths = []
    def bt(path, used):
        if len(path) == n:
            paths.append(tuple(path))
            return
        last = path[-1]
        for u in range(n):
            if not (used & (1 << u)) and T[last][u]:
                path.append(u)
                bt(path, used | (1 << u))
                path.pop()
    for v in range(n):
        bt([v], 1 << v)
    return paths


def count_adj_and_unmatched(T, Tp, i, j):
    """Count adj(i,j), adj'(j,i), and unmatched paths."""
    n = len(T)
    paths_T = ham_paths_list(T)
    paths_Tp = ham_paths_list(Tp)

    # adj(i,j): T-paths using i immediately before j
    adj_T_paths = []
    for p in paths_T:
        for k in range(n-1):
            if p[k] == i and p[k+1] == j:
                adj_T_paths.append((p, k))
                break

    # adj'(j,i): T'-paths using j immediately before i
    adj_Tp_paths = []
    for p in paths_Tp:
        for k in range(n-1):
            if p[k] == j and p[k+1] == i:
                adj_Tp_paths.append((p, k))
                break

    # Unmatched T-paths: can't swap i,j
    uT = 0
    for p, k in adj_T_paths:
        blocked = False
        if k > 0 and T[p[k-1]][j] == 0:
            blocked = True
        if k+1 < n-1 and T[i][p[k+2]] == 0:
            blocked = True
        if blocked:
            uT += 1

    # Unmatched T'-paths
    uTp = 0
    for p, k in adj_Tp_paths:
        blocked = False
        if k > 0 and T[p[k-1]][i] == 0:  # T'[p[k-1]][i] = T[p[k-1]][i]
            blocked = True
        if k+1 < n-1 and T[j][p[k+2]] == 0:  # T'[j][p[k+2]] = T[j][p[k+2]]
            blocked = True
        if blocked:
            uTp += 1

    return len(adj_T_paths), len(adj_Tp_paths), uT, uTp


def nadj3(T, a, b, c):
    """Count Ham paths containing consecutive triple (a, b, c)."""
    n = len(T)
    if not (T[a][b] and T[b][c]):
        return 0
    paths = ham_paths_list(T)
    count = 0
    for p in paths:
        for k in range(n-2):
            if p[k] == a and p[k+1] == b and p[k+2] == c:
                count += 1
                break
    return count


def main():
    print("=== Unmatched vs THM-013 Formula ===\n")

    # PART 1: Verify #U_T' - #U_T = delta_I at n=5 exhaustive
    print("--- Part 1: n=5 exhaustive ---")
    n = 5
    ok = 0
    total = 0
    for T in all_tournaments(n):
        for i in range(n):
            for j in range(i+1, n):
                if not T[i][j]:
                    continue
                Tp = flip_arc(T, i, j)
                _, _, uT, uTp = count_adj_and_unmatched(T, Tp, i, j)
                delta_unmatched = uTp - uT
                delta_I = compute_I(Tp) - compute_I(T)
                total += 1
                if delta_unmatched == delta_I:
                    ok += 1
    print(f"#U_T' - #U_T == delta_I: {ok}/{total}")

    # PART 2: Structural analysis of unmatched paths
    # At n=5: a path (..., x, i, j, ...) is blocked iff s_x=-1
    # nadj_T(x, i, j) = paths with consecutive x,i,j
    # nadj_T(i, j, x) = paths with consecutive i,j,x
    # #U_T = sum_{x: s=-1} [nadj(x,i,j) + nadj(i,j,x)] - sum_{x!=y, s=-1} nadj(x,i,j,y)
    print("\n--- Part 2: nadj structure at n=5 ---")
    n = 5
    # For each (T, i->j), decompose #U_T using nadj
    rng = random.Random(42)
    for trial in range(5):
        T = random_tournament(n, rng)
        i, j = rng.sample(range(n), 2)
        if not T[i][j]:
            i, j = j, i

        Tp = flip_arc(T, i, j)
        others = [x for x in range(n) if x != i and x != j]
        s = {x: 1 - T[x][i] - T[j][x] for x in others}

        print(f"\n  Trial {trial}: n={n}, flip {i}->{j}")
        print(f"  s values: {s}")

        _, _, uT, uTp = count_adj_and_unmatched(T, Tp, i, j)
        print(f"  #U_T={uT}, #U_T'={uTp}, delta={uTp-uT}")

        # nadj values
        for x in others:
            if s[x] == -1:
                n1 = nadj3(T, x, i, j)
                n2 = nadj3(T, i, j, x)
                print(f"  x={x} (s=-1): nadj(x,i,j)={n1}, nadj(i,j,x)={n2}, sum={n1+n2}")
            elif s[x] == 1:
                n1 = nadj3(Tp, x, j, i)
                n2 = nadj3(Tp, j, i, x)
                print(f"  x={x} (s=+1): nadj'(x,j,i)={n1}, nadj'(j,i,x)={n2}, sum={n1+n2}")

        # What is H(B_x)?
        for x in others:
            B_x = [v for v in range(n) if v != i and v != j and v != x]
            k = len(B_x)
            sub = [[0]*k for _ in range(k)]
            for a in range(k):
                for b in range(k):
                    sub[a][b] = T[B_x[a]][B_x[b]]
            h = hamiltonian_path_count(sub)
            print(f"  x={x}: H(B_x)={h}, B_x={B_x}")

    # PART 3: Key question - what is nadj(x,i,j) + nadj(i,j,x) as a function of B_x?
    # It should be related to H(B_x) but with endpoint-dependent corrections.
    print("\n\n--- Part 3: nadj(x,i,j) + nadj(i,j,x) vs H(B_x) at n=6 ---")
    n = 6
    rng = random.Random(123)
    data = []
    for _ in range(300):
        T = random_tournament(n, rng)
        i, j = rng.sample(range(n), 2)
        if not T[i][j]:
            i, j = j, i
        others = [x for x in range(n) if x != i and x != j]
        s = {x: 1 - T[x][i] - T[j][x] for x in others}

        for x in others:
            if s[x] != -1:
                continue
            n1 = nadj3(T, x, i, j)
            n2 = nadj3(T, i, j, x)

            B_x = [v for v in range(n) if v != i and v != j and v != x]
            k = len(B_x)
            sub = [[0]*k for _ in range(k)]
            for a in range(k):
                for b in range(k):
                    sub[a][b] = T[B_x[a]][B_x[b]]
            h = hamiltonian_path_count(sub)

            data.append({
                'nadj_sum': n1 + n2,
                'H_Bx': h,
                'nadj_xij': n1,
                'nadj_ijx': n2,
            })

    print(f"Collected {len(data)} (T, i->j, x with s=-1) triples")
    ratios = Counter()
    for d in data:
        ratios[(d['nadj_sum'], d['H_Bx'])] += 1
    print("\n(nadj_sum, H_Bx) -> count:")
    for key in sorted(ratios.keys()):
        print(f"  {key}: {ratios[key]}")

    # Check if nadj_sum / H_Bx is constant
    ratio_vals = Counter()
    for d in data:
        if d['H_Bx'] > 0:
            r = d['nadj_sum'] / d['H_Bx']
            ratio_vals[r] += 1
    print(f"\nnadj_sum / H_Bx ratios: {sorted(ratio_vals.items())}")

    # PART 4: What about endpoint-restricted H counts?
    # nadj(x,i,j) = #{paths of T[B_x] starting from q where T[j][q]=1,
    #                 split as prefix + suffix where prefix ends at p with T[p][x]=1}
    # Actually nadj(x,i,j) = #{full Ham paths of T with x,i,j consecutive}
    # = #{ways to split B_x into (prefix_set, suffix_set) with
    #     Ham path of prefix ending at p with T[p][x]=1
    #     Ham path of suffix starting at q with T[j][q]=1}

    # Simpler to think: nadj(x,i,j) = sum over q in B_x where T[j][q]=1:
    #   #{Ham paths of T starting x,i,j,q,...}
    # Actually no, that's only when j is at position 2 (i.e., x is at 0).

    # Let me try: define h_pairs(T[B_x]) = sum over (start s, end e) of
    #   #{paths from s to e} * {T[j][s]} * {T[e][x]}
    # This would be nadj(x,i,j) only if x is not at position 0 or n-1.
    # If x is at position 0: it's the start, so there's no predecessor, so
    # nadj(x,i,j) includes paths starting with x.

    # This is getting complicated. Let me just check the KEY AGGREGATE identity.
    print("\n\n--- Part 4: Aggregate identity at n=6 ---")
    # Test: #U_T' - #U_T = delta_I for n=6
    n = 6
    rng = random.Random(42)
    ok_agg = 0
    total_agg = 0
    for _ in range(200):
        T = random_tournament(n, rng)
        i, j = rng.sample(range(n), 2)
        if not T[i][j]:
            i, j = j, i
        Tp = flip_arc(T, i, j)
        _, _, uT, uTp = count_adj_and_unmatched(T, Tp, i, j)
        delta_unmatched = uTp - uT
        delta_I = compute_I(Tp) - compute_I(T)
        total_agg += 1
        if delta_unmatched == delta_I:
            ok_agg += 1
    print(f"#U_T' - #U_T == delta_I: {ok_agg}/{total_agg}")

    # PART 5: The involution idea.
    # Instead of decomposing by blocking vertex, can we find an involution
    # on U_T ∪ U_T' with exactly |delta_I| fixed points?
    # Or a weight-preserving map?
    print("\n\n--- Part 5: Involution structure ---")
    # For each unmatched T-path (...,x,i,j,...) where s_x=-1:
    # The 3-cycle {x,i,j} with x->i->j->x is a "certificate" for why the path is blocked.
    # The complement T[B_x] contains the rest of the path.
    # Can we pair unmatched T-paths with unmatched T'-paths via complement paths?

    # Key insight to test: for a vertex x with s_x = -1 in T,
    # in T' we have s_x = ?  Let's check.
    # s_x(T) = 1 - T[x][i] - T[j][x]
    # s_x(T') = 1 - T'[x][i] - T'[j][x] = 1 - T[x][i] - T[j][x] = s_x(T)
    # WAIT: T' only changes arcs between i and j. Arcs x->i and j->x are UNCHANGED.
    # So s_x is the SAME in T and T'!

    # This means: s_x = -1 blocks T-paths, and s_x = +1 blocks T'-paths,
    # and these are DIFFERENT vertices. The blocking structure is "orthogonal."

    # So #U_T involves vertices x with s=-1 (3-cycles using i->j)
    # and #U_T' involves vertices x with s=+1 (3-cycles using j->i)
    # and these are complementary sets!

    n = 5
    print(f"s_x distribution at n={n}:")
    counts = Counter()
    for T in all_tournaments(n):
        for i in range(n):
            for j in range(i+1, n):
                if not T[i][j]:
                    continue
                for x in range(n):
                    if x == i or x == j:
                        continue
                    sv = 1 - T[x][i] - T[j][x]
                    counts[sv] += 1
    for sv in sorted(counts):
        print(f"  s={sv:+d}: {counts[sv]}")

    # PART 6: Direct proof approach for n=6
    # At n=6: delta_I = -2*sum(s_x*H(B_x)) + 2*(D5-C5)
    # We need: #U_T' - #U_T = this formula.
    #
    # #U_T = sum_{x: s=-1} [nadj(x,i,j) + nadj(i,j,x)] - double_blocked
    # #U_T' = sum_{x: s=+1} [nadj'(x,j,i) + nadj'(j,i,x)] - double_blocked'
    #
    # Key idea: nadj(x,i,j) + nadj(i,j,x) counts paths through the 3-CYCLE {x,i,j}.
    # This is exactly 2 * #{Ham paths of T that "use" the 3-cycle as a consecutive triple}
    # DIVIDED by... no. nadj(x,i,j) has x before i, nadj(i,j,x) has x after j.
    #
    # Alternative: nadj(x,i,j) + nadj(i,j,x) = #{T-paths using arc i->j with x adjacent to {i,j}}
    #
    # At n=6, B_x has 3 vertices. Let's verify that:
    # nadj(x,i,j) + nadj(i,j,x) = h_startable(B_x, j) * ... no, this is path-dependent.
    #
    # Let me try a COMPLETELY different angle.
    # Define f(x) = nadj_T(x,i,j) + nadj_T(i,j,x) for s_x=-1
    # and g(x) = nadj_{T'}(x,j,i) + nadj_{T'}(j,i,x) for s_x=+1
    #
    # Then #U_T <= sum_{s=-1} f(x) (with double-count correction)
    # and #U_T' <= sum_{s=+1} g(x)
    #
    # What if we compute sum_{s=-1} f(x) - sum_{s=+1} g(x) directly?
    print("\n\n--- Part 6: Direct sum at n=6 ---")
    n = 6
    rng = random.Random(999)
    ok_direct = 0
    total_direct = 0

    for _ in range(200):
        T = random_tournament(n, rng)
        i, j = rng.sample(range(n), 2)
        if not T[i][j]:
            i, j = j, i
        Tp = flip_arc(T, i, j)
        others = [x for x in range(n) if x != i and x != j]
        s = {x: 1 - T[x][i] - T[j][x] for x in others}

        # Compute sum_{s=-1} f(x) and sum_{s=+1} g(x)
        sum_f = 0
        for x in others:
            if s[x] == -1:
                sum_f += nadj3(T, x, i, j) + nadj3(T, i, j, x)

        sum_g = 0
        for x in others:
            if s[x] == +1:
                sum_g += nadj3(Tp, x, j, i) + nadj3(Tp, j, i, x)

        # Compute double-blocking corrections
        _, _, uT, uTp = count_adj_and_unmatched(T, Tp, i, j)
        db_T = sum_f - uT  # double-blocked T-paths
        db_Tp = sum_g - uTp  # double-blocked T'-paths

        delta = uTp - uT
        delta_I = compute_I(Tp) - compute_I(T)

        total_direct += 1
        if delta == delta_I:
            ok_direct += 1

        if total_direct <= 5:
            print(f"  Trial {total_direct}: sum_f={sum_f}, sum_g={sum_g}, "
                  f"uT={uT}, uTp={uTp}, db_T={db_T}, db_Tp={db_Tp}, "
                  f"delta={delta}, delta_I={delta_I}")

    print(f"\n#U_T' - #U_T == delta_I: {ok_direct}/{total_direct}")

    # PART 7: Can we relate nadj sums to H(B_x) via endpoint counts?
    # For 3-vertex B_x: H(B_x) in {1, 3}. H=3 iff B_x is a 3-cycle.
    # nadj(x,i,j) = #{paths (...,a,x,i,j,b,...)} where a->x and j->b in T
    # The B_x portion: a path from some vertex in B_x to another, with constraints.
    print("\n\n--- Part 7: Endpoint analysis at n=6 ---")
    n = 6
    rng = random.Random(42)
    # For each x with s=-1, compute:
    # nadj(x,i,j) in terms of "paths of B_x starting at q where T[j][q]=1
    #              and ending at p where T[p][x]=1" — NO, it's more complex.
    # Let's just compute and compare.
    for trial in range(10):
        T = random_tournament(n, rng)
        i, j = rng.sample(range(n), 2)
        if not T[i][j]:
            i, j = j, i
        others = [x for x in range(n) if x != i and x != j]
        s = {x: 1 - T[x][i] - T[j][x] for x in others}

        neg_verts = [x for x in others if s[x] == -1]
        if not neg_verts:
            continue

        x = neg_verts[0]
        B_x = [v for v in range(n) if v != i and v != j and v != x]
        bk = len(B_x)
        sub = [[0]*bk for _ in range(bk)]
        for a in range(bk):
            for b in range(bk):
                sub[a][b] = T[B_x[a]][B_x[b]]
        h = hamiltonian_path_count(sub)

        # Endpoint-restricted counts
        # paths of B_x from start s to end e
        paths_Bx = ham_paths_list(sub)
        # Map local indices back to global
        h_se = defaultdict(int)  # (start_global, end_global) -> count
        for p in paths_Bx:
            s_g = B_x[p[0]]
            e_g = B_x[p[-1]]
            h_se[(s_g, e_g)] += 1

        # nadj(x,i,j): paths with ...a, x, i, j, b...
        # OR x at position 0: x, i, j, b... (no a)
        # OR j at last position: ...a, x, i, j (no b)
        n1 = nadj3(T, x, i, j)
        n2 = nadj3(T, i, j, x)

        # For nadj(x,i,j):
        # When x is not at start: need a in B_x with T[a][x]=1, and rest of B_x split
        # When x is at start: need b in B_x with T[j][b]=1, and path through rest

        # Actually, let's think of it as:
        # nadj(x,i,j) = sum over paths P of B_x:
        #   P can be placed as prefix or suffix or split
        #   such that x,i,j are consecutive in the full path
        # This equals: #{completions of (?,x,i,j,?) to full Ham path}

        # = h_end(B_x, a) * T[a][x] * [not at start] + ...
        # Too complex. Let me just verify numerically.

        # Compute "restricted H": h_restricted(B_x, j, x)
        #   = #{paths of B_x from q (where T[j][q]=1) to p (where T[p][x]=1)}
        #   + special terms for x at start/j at end
        h_jq_px = 0
        for p in paths_Bx:
            start_g = B_x[p[0]]
            end_g = B_x[p[-1]]
            if T[j][start_g] and T[end_g][x]:
                h_jq_px += 1
        # Also: when x is at position 0 (start of full path):
        # path = (x, i, j, ...) with ... a Ham path of B_x starting at q with T[j][q]=1
        h_j_start = sum(1 for p in paths_Bx if T[j][B_x[p[0]]])
        # When j is at last position:
        # path = (..., x, i, j) with ... a Ham path of B_x ending at a with T[a][x]=1
        h_x_end = sum(1 for p in paths_Bx if T[B_x[p[-1]]][x])

        if trial < 10:
            print(f"\n  Trial {trial}: x={x}, B_x={B_x}, H(B_x)={h}")
            print(f"    nadj(x,i,j)={n1}, nadj(i,j,x)={n2}")
            print(f"    h(j->start, end->x)={h_jq_px}")
            print(f"    h(j->start)={h_j_start}, h(end->x)={h_x_end}")
            # nadj(x,i,j) should equal h_jq_px (middle case) + contributions from endpoints
            # More precisely: for the full path (prefix, x, i, j, suffix):
            # prefix is a Ham path on some S subset B_x, ending at a with T[a][x]=1
            # suffix is a Ham path on B_x\S, starting at q with T[j][q]=1
            # So nadj(x,i,j) = sum_{S} h_end(T[S], a:T[a][x]=1) * h_start(T[B_x\S], q:T[j][q]=1)
            # When S = empty: prefix is empty (x starts path), suffix = full B_x starting at q
            # When S = B_x: prefix = full B_x ending at a, suffix = empty (j ends path)


if __name__ == "__main__":
    main()
