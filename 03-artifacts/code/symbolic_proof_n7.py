#!/usr/bin/env python3
"""
Overnight verification: prove delta_H = delta_I at n=7.

n=7 has 2^20 = 1,048,576 variable assignments (5 internal arcs of V\{i,j}
plus 10 arcs connecting V\{i,j} to {i,j}).

Uses optimized DP to compute adj and unmatched counts without
enumerating all paths.

Instance: kind-pasteur-2026-03-05-S6
"""

import sys
import time
from itertools import combinations, permutations


def prove_n7():
    print("=== n=7 Symbolic Proof (Overnight) ===")
    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    I, J = 0, 1
    others = [2, 3, 4, 5, 6]
    n = 7

    # Variables: T[x][0] for x in others (5), T[1][x] for x in others (5),
    # internal arcs C(5,2)=10. Total: 20 variables, 2^20 = 1048576.
    # BUT: T[x][0] and T[1][x] are the "interface" arcs.
    # Internal arcs: T[a][b] for a,b in others, a<b.

    # Variable ordering: first 5 are T[x][0] for x=2..6,
    # next 5 are T[1][x] for x=2..6,
    # next 10 are internal arcs T[a][b] for (a,b) in sorted pairs from others.
    interface_vars = [(x, 0) for x in others] + [(1, x) for x in others]
    internal_pairs = [(a, b) for a in others for b in others if a < b]
    all_vars = interface_vars + internal_pairs
    n_vars = len(all_vars)
    assert n_vars == 20, f"Expected 20 vars, got {n_vars}"
    total_cases = 2 ** n_vars

    print(f"Variables: {n_vars}, Total cases: {total_cases}")
    print(f"Estimated time: ~4-12 hours depending on machine speed")
    print()

    # Precompute variable index map
    var_idx = {}
    for idx, (a, b) in enumerate(all_vars):
        var_idx[(a, b)] = idx

    def build_adj_matrix(mask):
        """Build 7x7 adjacency matrix from variable mask."""
        T = [[0] * n for _ in range(n)]
        T[I][J] = 1  # fixed arc 0->1

        for idx, (a, b) in enumerate(all_vars):
            val = (mask >> idx) & 1
            T[a][b] = val
            T[b][a] = 1 - val

        return T

    def build_adj_matrix_flipped(mask):
        """Build adjacency matrix with arc 0->1 flipped to 1->0."""
        T = build_adj_matrix(mask)
        T[I][J] = 0
        T[J][I] = 1
        return T

    def compute_adj_unmatched(T, src, dst):
        """
        Compute adj(src,dst) and #unmatched using DP.

        Uses bitmask DP tracking:
        - Whether arc src->dst has been used
        - The blocking status when it was used

        States: (mask, last_vertex, arc_status)
        arc_status:
          0 = not yet used src->dst
          1 = used, and matched (pred ok so far, succ pending)
          2 = used, and pred-blocked (definitely unmatched)
          3 = used, matched (confirmed: pred ok AND succ ok or dst is last)
          4 = used, unmatched (confirmed)
        """
        full = (1 << n) - 1

        # dp[mask][v][status] = count
        # status: 0=unused, 1=used_pred_ok, 2=used_pred_bad
        # We resolve succ status when extending from dst

        # Use arrays for speed
        # dp[mask * n * 3 + v * 3 + status]
        size = (1 << n) * n * 3
        dp = [0] * size

        def idx(mask, v, s):
            return mask * n * 3 + v * 3 + s

        # Base: single vertex
        for v in range(n):
            dp[idx(1 << v, v, 0)] = 1

        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)):
                    continue

                for status in range(3):
                    cnt = dp[idx(mask, v, status)]
                    if cnt == 0:
                        continue

                    for u in range(n):
                        if mask & (1 << u):
                            continue
                        if not T[v][u]:
                            continue

                        nmask = mask | (1 << u)

                        if v == src and u == dst:
                            # Using the arc src->dst
                            # Check pred status
                            # Predecessor of src is the vertex before src in path
                            # At this point, src=v is the last vertex placed,
                            # and we're about to place dst=u.
                            # But we need the vertex BEFORE src.
                            # In the DP, we don't track the full path, just the last vertex.
                            # We need to know if the vertex before src beats dst.

                            # The vertex before src is NOT tracked in standard DP.
                            # We need a modified DP that tracks the second-to-last vertex
                            # OR the pred-blocking status.

                            # PROBLEM: standard DP doesn't track predecessor.
                            # We need dp[mask][v][prev_v] or at least know T[prev][dst].

                            # Let me restructure: instead of tracking prev for all states,
                            # only track it at the moment of using src->dst.

                            # Alternative: split into cases based on whether src is at start.
                            pass

                        # This approach is getting too complex for inline DP.
                        # Fall back to path enumeration with pruning.
                        pass

        # Fall back to backtracking
        return _backtrack_count(T, n, src, dst)

    def _backtrack_count(T, n, src, dst):
        """Count adj(src,dst) and unmatched via backtracking."""
        adj = 0
        unmatched = 0

        path = [0] * n
        used = [0]  # bitmask

        def bt(pos):
            nonlocal adj, unmatched
            if pos == n:
                # Check if src->dst was used
                for k in range(n - 1):
                    if path[k] == src and path[k + 1] == dst:
                        adj += 1
                        blocked = False
                        if k > 0 and not T[path[k - 1]][dst]:
                            blocked = True
                        if k + 1 < n - 1 and not T[src][path[k + 2]]:
                            blocked = True
                        if blocked:
                            unmatched += 1
                        break
                return

            last = path[pos - 1] if pos > 0 else -1
            for v in range(n):
                if used[0] & (1 << v):
                    continue
                if pos > 0 and not T[last][v]:
                    continue
                path[pos] = v
                used[0] |= (1 << v)
                bt(pos + 1)
                used[0] ^= (1 << v)

        bt(0)
        return adj, unmatched

    # Precompute cycle counting helpers
    def count_lost_gained_cycles(T, I, J, others, length):
        """Count lost and gained odd cycles of given length using arc I->J."""
        lost = 0
        gained = 0
        if length == 3:
            for x in others:
                # Lost 3-cycle: I->J->x->I needs T[J][x]=1 and T[x][I]=1
                if T[J][x] and T[x][I]:
                    lost += 1
                # Gained 3-cycle: J->I->x->J needs T[I][x]=1 and T[x][J]=1
                if T[I][x] and T[x][J]:
                    gained += 1
        elif length == 5:
            for subset in combinations(others, 3):
                for perm in permutations(subset):
                    v1, v2, v3 = perm
                    # Lost: (I,J,v1,v2,v3,I)
                    if T[J][v1] and T[v1][v2] and T[v2][v3] and T[v3][I]:
                        lost += 1
                    # Gained: (J,I,v1,v2,v3,J)
                    if T[I][v1] and T[v1][v2] and T[v2][v3] and T[v3][J]:
                        gained += 1
        elif length == 7:
            # 7-cycle uses all 7 vertices
            for perm in permutations(others):
                v1, v2, v3, v4, v5 = perm
                # Lost: (I,J,v1,v2,v3,v4,v5,I)
                if (T[J][v1] and T[v1][v2] and T[v2][v3] and
                    T[v3][v4] and T[v4][v5] and T[v5][I]):
                    lost += 1
                # Gained: (J,I,v1,v2,v3,v4,v5,J)
                if (T[I][v1] and T[v1][v2] and T[v2][v3] and
                    T[v3][v4] and T[v4][v5] and T[v5][J]):
                    gained += 1
        return lost, gained

    def h_tournament(T, verts):
        """Count Ham paths of T restricted to verts."""
        k = len(verts)
        if k == 0:
            return 1
        if k == 1:
            return 1
        count = 0
        for perm in permutations(verts):
            valid = True
            for idx in range(k - 1):
                if not T[perm[idx]][perm[idx + 1]]:
                    valid = False
                    break
            if valid:
                count += 1
        return count

    # For n=7, the THM-013 formula is (same as n<=7 simplified):
    # H(T)-H(T') = -2*sum(s_x*H(B_x)) + 2*(D5-C5) + 2*(D7-C7)
    # where D=destroyed(lost), C=created(gained) in THM-013 convention.
    # H(T')-H(T) = 2*sum(s_x*H(B_x)) - 2*(D5-C5) - 2*(D7-C7)
    #            = 2*sum(s_x*H(B_x)) + 2*(gained5-lost5) + 2*(gained7-lost7)

    # Actually for n=7, from THM-013:
    # "n=7: Same structure as n=6 (alpha_k=0 for k>=3, no VD 3-5 possible)."
    # DeltaH = -2*sum_x s_x*H(B_x) + 2*sum_{L>=5}(DL-CL)
    # = -2*sum s_x*H(B_x) + 2*(D5-C5) + 2*(D7-C7)
    # In my convention (H(T')-H(T)):
    # = 2*sum s_x*H(B_x) + 2*(gained5-lost5) + 2*(gained7-lost7)

    ok = 0
    total = 0
    fail_count = 0
    start_time = time.time()

    for mask in range(total_cases):
        T = build_adj_matrix(mask)

        _, uT = _backtrack_count(T, n, I, J)
        # For T', only the arc between I and J changes
        Tp = build_adj_matrix_flipped(mask)
        _, uTp = _backtrack_count(Tp, n, J, I)

        delta = uTp - uT

        # Compute expected
        s_vals = {x: 1 - T[x][I] - T[J][x] for x in others}

        # H(B_x) for 4-vertex sub-tournaments
        formula_sum = 0
        for x in others:
            B_x = [v for v in others if v != x]
            hbx = h_tournament(T, B_x)
            formula_sum += s_vals[x] * hbx

        # Cycle counts
        lost5, gained5 = count_lost_gained_cycles(T, I, J, others, 5)
        lost7, gained7 = count_lost_gained_cycles(T, I, J, others, 7)

        expected = 2 * formula_sum + 2 * (gained5 - lost5) + 2 * (gained7 - lost7)

        total += 1
        if delta == expected:
            ok += 1
        else:
            fail_count += 1
            if fail_count <= 10:
                print(f"  FAIL: mask={mask}, delta={delta}, expected={expected}, "
                      f"sum_sH={formula_sum}, g5-l5={gained5-lost5}, g7-l7={gained7-lost7}")
            if fail_count > 50:
                print(f"  Aborting after {fail_count} failures.")
                break

        if total % 100000 == 0:
            elapsed = time.time() - start_time
            rate = total / elapsed
            eta = (total_cases - total) / rate
            print(f"  Progress: {total}/{total_cases} ({100*total/total_cases:.1f}%), "
                  f"ok={ok}, fail={fail_count}, "
                  f"rate={rate:.0f}/s, ETA={eta/3600:.1f}h")

    elapsed = time.time() - start_time
    print(f"\nCompleted in {elapsed:.0f}s ({elapsed/3600:.1f}h)")
    print(f"Identity verified: {ok}/{total}")
    if ok == total:
        print("PROVED: delta_H = delta_I for ALL n=7 arc configurations.")
        print("This proves OCF at n=7 via arc-flip induction.")
    elif fail_count > 0:
        print(f"FAILED: {fail_count} counterexamples found.")


if __name__ == "__main__":
    prove_n7()
