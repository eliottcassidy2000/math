#!/usr/bin/env python3
"""
Fast n=7 verification of delta_H = delta_I (polynomial identity, THM-015).

Uses bitmask DP instead of backtracking. For each variable assignment:
1. Compute adj(0,1) and adj'(1,0) via DP
2. Compute U_T and U_{T'} (unmatched path counts)
3. Verify U_{T'} - U_T = delta_I from THM-013

DP approach: O(2^n * n^2) per tournament, vs O(n!) backtracking.

For unmatched counts, we need nadj(x,i,j) = #{paths with x,i,j consecutive}.
This requires a modified DP tracking the last TWO vertices.

Instance: kind-pasteur-2026-03-05-S7
"""

import time
from itertools import combinations, permutations

N = 7
FULL = (1 << N) - 1
I_V, J_V = 0, 1
OTHERS = [2, 3, 4, 5, 6]


def run_verification():
    print("=== n=7 Fast Symbolic Proof (DP-based) ===")
    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Variable ordering: T[x][0] for x in others (5 vars),
    # T[1][x] for x in others (5 vars),
    # T[a][b] for a<b in others (10 vars). Total: 20.
    interface_vars = [(x, 0) for x in OTHERS] + [(1, x) for x in OTHERS]
    internal_pairs = [(a, b) for a in OTHERS for b in OTHERS if a < b]
    all_vars = interface_vars + internal_pairs
    n_vars = len(all_vars)
    assert n_vars == 20
    total_cases = 1 << n_vars

    print(f"Variables: {n_vars}, Total cases: {total_cases:,}")

    # Precompute variable index
    var_idx = {}
    for idx, (a, b) in enumerate(all_vars):
        var_idx[(a, b)] = idx

    def build_T_array(mask):
        """Build adjacency as flat array T[a*N+b]."""
        T = [0] * (N * N)
        T[I_V * N + J_V] = 1  # 0->1
        # T[J_V * N + I_V] = 0  # already 0
        for idx, (a, b) in enumerate(all_vars):
            val = (mask >> idx) & 1
            T[a * N + b] = val
            T[b * N + a] = 1 - val
        return T

    def compute_unmatched_dp(T, src, dst):
        """
        Count unmatched paths using arc src->dst via DP.

        A path (..., x, src, dst, y, ...) is unmatched iff
        T[x*N+dst]=0 (pred-blocked) or T[src*N+y]=0 (succ-blocked).

        We use a DP tracking (mask, last, second_to_last) to know the
        context around the src->dst arc when it appears.

        But that's O(2^n * n^3) which is expensive. Instead, use the
        factored DP approach:

        adj(src,dst) = sum_S dp_end[S][src] * dp_start[complement][dst]
        where S contains src but not dst.

        For unmatched, we need the predecessor of src (last vertex in
        S-path before src) and successor of dst (first vertex in
        complement-path after dst).
        """
        if T[src * N + dst] == 0:
            return 0, 0

        # dp_end[mask][v] = # Ham paths on mask ending at v
        size = 1 << N
        dp_end = [[0] * N for _ in range(size)]
        for v in range(N):
            dp_end[1 << v][v] = 1

        for mask in range(1, size):
            for v in range(N):
                if not (mask & (1 << v)):
                    continue
                c = dp_end[mask][v]
                if c == 0:
                    continue
                for u in range(N):
                    if mask & (1 << u):
                        continue
                    if T[v * N + u]:
                        dp_end[mask | (1 << u)][u] += c

        # dp_end_prev[mask][v][prev] = # paths on mask ending at ...,prev,v
        # This is O(2^n * n^2) states which is 128*49 = 6272, very manageable
        dp_end_prev = [[[0] * N for _ in range(N)] for _ in range(size)]
        # Base: single vertex, no predecessor
        # We'll use prev=-1 for "no predecessor", but store separately
        # Actually, let's track it properly. For 2-vertex paths:
        for v in range(N):
            for u in range(N):
                if u == v:
                    continue
                if T[v * N + u]:
                    dp_end_prev[1 << v | 1 << u][u][v] = 1

        for mask in range(1, size):
            for v in range(N):
                if not (mask & (1 << v)):
                    continue
                for prev in range(N):
                    if prev == v or not (mask & (1 << prev)):
                        continue
                    c = dp_end_prev[mask][v][prev]
                    if c == 0:
                        continue
                    for u in range(N):
                        if mask & (1 << u):
                            continue
                        if T[v * N + u]:
                            dp_end_prev[mask | (1 << u)][u][v] += c

        # dp_start_next[mask][v][nxt] = # paths on mask starting at v
        # with second vertex = nxt
        # Build this by the REVERSE DP: dp on reversed tournament
        # Actually easier: dp_start[mask][v] = paths starting at v on mask
        # We need the second vertex. Let's build dp_start_next similarly.

        # dp_start[mask][v] = # Ham paths on mask starting at v
        # We build this by DP from the END:
        # dp_from_end[mask][v] = # paths on mask ending at v (same as dp_end)
        # Then dp_start[mask][v] = sum over all endpoints e of paths starting
        # at v and ending at e. But dp_end doesn't track start.

        # Alternative: dp_start_next[mask][v][nxt] where the path is v->nxt->...
        # Base: two vertices v->nxt
        dp_sn = [[[0] * N for _ in range(N)] for _ in range(size)]
        for v in range(N):
            for nxt in range(N):
                if nxt == v:
                    continue
                if T[v * N + nxt]:
                    dp_sn[(1 << v) | (1 << nxt)][v][nxt] = 1

        # Extend at the END of paths (adding vertices after the last one)
        # dp_sn_ext[mask][start][second][last] -- too many dimensions
        # Instead, let's combine: dp_sn[mask][start][second] = # paths on mask
        # starting with start->second->...

        # Build by extending dp_end with the constraint that first two vertices
        # are fixed. This is equivalent to:
        # dp_sn[mask][v][nxt] = dp_end[mask \ {v}][nxt -> ... -> last] starting from nxt
        # = sum over last endpoint: dp_end restricted paths starting at nxt on mask\{v}
        # Hmm, dp_end doesn't track start vertex.

        # Let me just use a different approach. For each pair (v, nxt) with T[v][nxt]=1:
        # dp_sn[(1<<v)|(1<<nxt)][v][nxt] = 1
        # Then extend by adding a vertex at the end:
        # Need dp with start vertex fixed. dp_fixed_start[mask][start][last]
        # = # paths on mask from start to last.
        # This is O(2^n * n^2) = 128*49 states. Build it.

        dp_fs = [[[0] * N for _ in range(N)] for _ in range(size)]
        for v in range(N):
            dp_fs[1 << v][v][v] = 1  # single vertex: start=v, last=v

        for mask in range(1, size):
            for s in range(N):
                if not (mask & (1 << s)):
                    continue
                for last in range(N):
                    if not (mask & (1 << last)):
                        continue
                    c = dp_fs[mask][s][last]
                    if c == 0:
                        continue
                    for u in range(N):
                        if mask & (1 << u):
                            continue
                        if T[last * N + u]:
                            dp_fs[mask | (1 << u)][s][u] += c

        # Now:
        # adj(src,dst) = sum over S containing src, not dst:
        #   (sum_prev dp_end_prev[S][src][prev]) * dp_fs[V\S][dst][some_end]
        #   where dp_fs is summed over all endpoints
        # = sum_S dp_end[S][src] * (sum_e dp_fs[V\S][dst][e])

        # For unmatched: path ...,prev,src,dst,nxt,...
        # blocked iff T[prev][dst]=0 OR T[src][nxt]=0
        # pred-blocked: T[prev*N+dst]=0
        # succ-blocked: T[src*N+nxt]=0

        # Case 1: src is first in path (no prev). Path = src,dst,...
        # No pred-block possible. Succ: T[src*N+nxt]=0 where nxt is second of complement path
        # For unmatched: need T[src*N+nxt]=0 for the vertex right after dst
        # Wait, if path is (src, dst, nxt, ...), the swap gives (dst, src, nxt, ...)
        # This needs T[dst][src] (which is 0 in T, since src->dst), so it's always blocked!
        # No wait... the swap exchanges POSITIONS of src and dst, not the arc.
        # The swap of (src, dst, v3, ..., vn) is (dst, src, v3, ..., vn).
        # For this to be valid in T': T'[dst][src]=1 (yes, T' has dst->src... wait no)
        # T has src->dst. T' flips to dst->src. So T'[dst][src]=1.
        # Also need T'[src][v3] = T[src][v3] (unchanged since v3 != dst).
        # So swap of (src,dst,v3,...,vn) = (dst,src,v3,...,vn) needs T[src][v3].
        # It's blocked iff T[src][v3]=0.
        # No predecessor for src (position 1), so no pred-block.
        # => unmatched iff n>2 and T[src][v3]=0

        # Case 2: dst is last in path (no succ). Path = ...,prev,src,dst.
        # Swap: ...,prev,dst,src. Need T[prev][dst]. Blocked iff T[prev][dst]=0.
        # No succ-block.

        # General case: path = ...,prev,src,dst,nxt,...
        # Blocked iff T[prev][dst]=0 or T[src][nxt]=0.

        src_bit = 1 << src
        dst_bit = 1 << dst

        adj_count = 0
        unmatched = 0

        # Enumerate all subsets S that contain src, not dst
        # S is the set of vertices BEFORE and including src in the path
        # The complement R = FULL ^ S contains dst and all vertices after src

        # For efficiency, iterate over subsets containing src but not dst
        # Use the "enumerate subsets of a mask" trick

        # All subsets of FULL that include src but not dst
        others_mask = FULL ^ src_bit ^ dst_bit  # vertices other than src, dst

        # S = src_bit | subset_of_others
        # R = dst_bit | (others_mask ^ subset_of_others)
        for sub_others in range(1 << (N - 2)):
            # Map sub_others bits to actual vertex positions
            # This is tricky. Let's enumerate differently.
            pass

        # Actually, let me just iterate over all masks containing src but not dst
        # For n=7, there are 2^5 = 32 possible S sizes, but we need to iterate
        # over all subsets. Total: 2^5 = 32 subsets of others.

        # Better approach: directly iterate
        others_list = [v for v in range(N) if v != src and v != dst]
        n_others = len(others_list)  # = 5

        for sub_mask in range(1 << n_others):
            # Build S = {src} union selected others
            S = src_bit
            for bit_idx in range(n_others):
                if sub_mask & (1 << bit_idx):
                    S |= (1 << others_list[bit_idx])

            R = (FULL ^ S) # contains dst and remaining others
            assert R & dst_bit

            # Paths ending at src on S
            # If S = {src} only (sub_mask=0): dp_end[S][src] = 1, no predecessor
            # Paths starting at dst on R

            cnt_end_at_src = dp_end[S][src]
            if cnt_end_at_src == 0:
                continue

            # Total paths on R starting at dst
            total_start_dst = sum(dp_fs[R][dst])
            if total_start_dst == 0:
                continue

            adj_count += cnt_end_at_src * total_start_dst

            # Now compute unmatched from these paths
            # Need to know prev (last before src in S-path) and nxt (first after dst in R-path)

            if S == src_bit:
                # No predecessor. Only succ-blocking possible.
                if R == dst_bit:
                    # Path is just (src, dst), n=2. No blocking.
                    pass
                else:
                    # nxt = second vertex in R-path starting at dst
                    for nxt in range(N):
                        if nxt == dst or not (R & (1 << nxt)):
                            continue
                        # Paths on R: dst -> nxt -> ... (dp_fs with start=dst, first step to nxt)
                        paths_dst_nxt = dp_fs[R][dst][nxt] if popcount(R) == 2 else 0
                        # Wait, dp_fs[R][dst][last] counts paths from dst to last.
                        # I need paths from dst where the SECOND vertex is nxt.
                        # That's dp_fs[R][dst][*] where first step is dst->nxt.
                        # = sum_e dp_fs[R\{dst}? no...
                        # Hmm, dp_fs doesn't track the second vertex.
                        pass
                    # This approach needs dp_fs to also track the second vertex,
                    # which is too expensive (O(2^n * n^3)).
                    # Let me use a different method.
                    pass

            # OK, the DP with fixed start doesn't give us the second vertex easily.
            # Let me switch to using dp_end_prev and a separate dp for the
            # complement that tracks the "next after start" vertex.

        # This is getting complicated. Let me use a cleaner decomposition.
        # Reset and use the direct nadj approach.

        adj_count = 0
        unmatched = 0

        # Iterate over all full-mask paths via dp_fs
        for s in range(N):
            for e in range(N):
                if dp_fs[FULL][s][e] == 0:
                    continue
                # dp_fs[FULL][s][e] counts paths from s to e.
                # But we need to find which of these use src->dst consecutively.
                # dp_fs doesn't track consecutive pairs.
                pass

        # I think the cleanest approach for n=7 is to enumerate paths directly
        # using the DP structure, but track the specific (prev, next) around
        # the src->dst arc. Let me build a specialized DP.

        # dp_adj[S][prev] = # Ham paths on S ending at ...,prev,src
        # Then for each such path, extend with dst and then paths on R starting at dst
        # dp_adj2[R][nxt] = # Ham paths on R starting at dst,nxt,...
        # Combined: adj with predecessor prev and successor nxt

        # dp_adj[S][prev]: paths on S ending at prev->src
        # S must contain both prev and src
        # dp_adj[S][prev] = dp_end_prev[S][src][prev]  (paths on S ending at src with prev before)
        # This is exactly dp_end_prev!

        # dp_adj2[R][nxt]: paths on R starting at dst->nxt->...
        # R must contain both dst and nxt
        # Build: for each R containing dst, dp_start_after[R][nxt] = # paths starting dst->nxt->...->end
        # dp_start_after[{dst,nxt}][nxt] = T[dst*N+nxt]  (1 if arc exists, 0 otherwise)
        # Extend by adding vertices at the end.

        # Actually this is just: dp_fs[R][dst][last] summed over last, but restricted to
        # paths where second vertex is nxt. Let me build dp_start_second:

        # dp_ss[R][nxt] = #{paths on R: dst->nxt->...}
        # = dp_fs[R\{dst}][nxt][*] where all paths on R\{dst} starting at nxt
        # ... only if T[dst][nxt]=1
        # = T[dst*N+nxt] * sum_e dp_fs[R ^ dst_bit][nxt][e]

        # That works! dp_ss[R][nxt] = T[dst*N+nxt] * sum_e dp_fs[R ^ dst_bit][nxt][e]

        for sub_mask in range(1 << n_others):
            S = src_bit
            for bit_idx in range(n_others):
                if sub_mask & (1 << bit_idx):
                    S |= (1 << others_list[bit_idx])
            R = FULL ^ S  # contains dst

            cnt_end_at_src = dp_end[S][src]
            if cnt_end_at_src == 0:
                continue

            if S == src_bit:
                # Path starts with src. No predecessor.
                # Succ-blocking only.
                if R == dst_bit:
                    # n=2 path (src, dst). No blocking.
                    total_r = 1
                    adj_count += total_r
                    # unmatched: no pred, no succ => not blocked
                    continue

                # Paths on R starting with dst
                total_r = sum(dp_fs[R][dst])
                if total_r == 0:
                    continue
                adj_count += total_r  # cnt_end_at_src = 1 for single vertex

                # Succ-blocked paths: those where T[src*N+nxt]=0
                for nxt in range(N):
                    if nxt == dst or not (R & (1 << nxt)):
                        continue
                    if not T[dst * N + nxt]:
                        continue
                    # Paths dst->nxt->... on R
                    paths_nxt = sum(dp_fs[R ^ dst_bit][nxt]) if T[dst * N + nxt] else 0
                    if T[src * N + nxt] == 0:
                        unmatched += paths_nxt  # succ-blocked

                continue

            # S has more than just src. We have a predecessor.
            R_no_dst = R ^ dst_bit  # R without dst

            for prev in range(N):
                if prev == src or not (S & (1 << prev)):
                    continue
                c_prev = dp_end_prev[S][src][prev]
                if c_prev == 0:
                    continue

                if R == dst_bit:
                    # dst is last. No successor.
                    # Pred-blocked: T[prev*N+dst]=0
                    adj_count += c_prev
                    if T[prev * N + dst] == 0:
                        unmatched += c_prev
                    continue

                # General case: both prev and nxt exist
                for nxt in range(N):
                    if not (R_no_dst & (1 << nxt)):
                        continue
                    if not T[dst * N + nxt]:
                        continue
                    # Paths on R starting dst->nxt->...
                    paths_nxt = sum(dp_fs[R_no_dst][nxt])
                    if paths_nxt == 0:
                        continue

                    total_here = c_prev * paths_nxt
                    adj_count += total_here

                    # Blocked iff T[prev][dst]=0 or T[src][nxt]=0
                    if T[prev * N + dst] == 0 or T[src * N + nxt] == 0:
                        unmatched += total_here

        return adj_count, unmatched

    # THM-013 formula components
    def compute_formula(T):
        """Compute delta_I from THM-013 for the flip of arc 0->1."""
        # s_x values
        s = {}
        for x in OTHERS:
            s[x] = 1 - T[x * N + I_V] - T[J_V * N + x]

        # H(B_x) for 4-vertex sub-tournaments B_x = OTHERS \ {x}
        def h_sub(verts):
            """Count Ham paths on sub-tournament."""
            k = len(verts)
            if k <= 1:
                return 1
            count = 0
            for p in permutations(verts):
                ok = True
                for idx in range(k - 1):
                    if not T[p[idx] * N + p[idx + 1]]:
                        ok = False
                        break
                if ok:
                    count += 1
            return count

        formula_sum = 0
        for x in OTHERS:
            Bx = [v for v in OTHERS if v != x]
            formula_sum += s[x] * h_sub(Bx)

        # 5-cycle counts
        lost5 = 0
        gained5 = 0
        for subset in combinations(OTHERS, 3):
            for perm in permutations(subset):
                v1, v2, v3 = perm
                if T[J_V * N + v1] and T[v1 * N + v2] and T[v2 * N + v3] and T[v3 * N + I_V]:
                    lost5 += 1
                if T[I_V * N + v1] and T[v1 * N + v2] and T[v2 * N + v3] and T[v3 * N + J_V]:
                    gained5 += 1

        # 7-cycle counts (all 7 vertices)
        lost7 = 0
        gained7 = 0
        for perm in permutations(OTHERS):
            v1, v2, v3, v4, v5 = perm
            if (T[J_V * N + v1] and T[v1 * N + v2] and T[v2 * N + v3] and
                    T[v3 * N + v4] and T[v4 * N + v5] and T[v5 * N + I_V]):
                lost7 += 1
            if (T[I_V * N + v1] and T[v1 * N + v2] and T[v2 * N + v3] and
                    T[v3 * N + v4] and T[v4 * N + v5] and T[v5 * N + J_V]):
                gained7 += 1

        expected = 2 * formula_sum + 2 * (gained5 - lost5) + 2 * (gained7 - lost7)
        return expected

    # Main verification loop
    ok = 0
    fail_count = 0
    start_time = time.time()

    for mask in range(total_cases):
        T = build_T_array(mask)

        _, uT = compute_unmatched_dp(T, I_V, J_V)

        # Build T' (flip 0->1 to 1->0)
        Tp = T[:]
        Tp[I_V * N + J_V] = 0
        Tp[J_V * N + I_V] = 1

        _, uTp = compute_unmatched_dp(Tp, J_V, I_V)

        delta = uTp - uT
        expected = compute_formula(T)

        if delta == expected:
            ok += 1
        else:
            fail_count += 1
            if fail_count <= 5:
                print(f"  FAIL mask={mask}: delta={delta}, expected={expected}")
            if fail_count > 20:
                print(f"  Aborting after {fail_count} failures.")
                break

        total = mask + 1
        if total % 50000 == 0:
            elapsed = time.time() - start_time
            rate = total / elapsed
            eta = (total_cases - total) / rate
            print(f"  Progress: {total:,}/{total_cases:,} ({100*total/total_cases:.1f}%), "
                  f"ok={ok}, fail={fail_count}, rate={rate:.0f}/s, ETA={eta:.0f}s")

    elapsed = time.time() - start_time
    print(f"\nCompleted in {elapsed:.1f}s")
    print(f"Identity verified: {ok}/{total_cases}")
    if ok == total_cases and fail_count == 0:
        print("PROVED: delta_H = delta_I for ALL n=7 arc configs (1,048,576 cases)")
    elif fail_count > 0:
        print(f"FAILED: {fail_count} counterexamples.")


def popcount(x):
    return bin(x).count('1')


if __name__ == "__main__":
    run_verification()
