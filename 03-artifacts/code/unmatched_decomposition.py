#!/usr/bin/env python3
"""
Decompose #U_T by blocking vertex and compare to H(B_x) terms.

THM-014 says: delta_H = #U_T - #U_T'
THM-013 says: delta_I = -2*sum_x s_x*H(B_x) + ... (for n<=7)

Question: does #U_T decompose per blocking vertex x as:
  #U_T = sum_{x: s_x=-1} f(x)
  #U_T' = sum_{x: s_x=+1} f'(x)
where f(x) relates to H(B_x)?

If so, the identity reduces to showing f(x) = 2*H(B_x) per vertex.

Instance: kind-pasteur-2026-03-05-S6
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
from collections import defaultdict
import random


def flip_arc(T, i, j):
    n = len(T)
    Tp = [row[:] for row in T]
    Tp[i][j] = 1 - T[i][j]
    Tp[j][i] = 1 - T[j][i]
    return Tp


def hamiltonian_paths(T):
    """Return all Hamiltonian paths as lists of vertices."""
    n = len(T)
    full = (1 << n) - 1
    # dp[mask][v] = list of paths ending at v using vertices in mask
    # Too expensive to store all paths for large n. Use count + specific queries.
    paths = []

    def backtrack(path, used):
        if len(path) == n:
            paths.append(list(path))
            return
        last = path[-1]
        for u in range(n):
            if not (used & (1 << u)) and T[last][u]:
                path.append(u)
                backtrack(path, used | (1 << u))
                path.pop()

    for v in range(n):
        backtrack([v], 1 << v)
    return paths


def analyze_flip(T, i, j):
    """Analyze the swap involution for flip i->j to j->i."""
    n = len(T)
    assert T[i][j] == 1
    Tp = flip_arc(T, i, j)

    # Compute s_x for all x != i,j
    others = [x for x in range(n) if x != i and x != j]
    s = {}
    for x in others:
        s[x] = 1 - T[x][i] - T[j][x]  # s_x in {-1, 0, +1}

    # Find all Ham paths of T and T'
    paths_T = hamiltonian_paths(T)
    paths_Tp = hamiltonian_paths(Tp)

    # adj(i,j) in T: paths using arc i->j consecutively
    adj_T = []
    for p in paths_T:
        for k in range(n-1):
            if p[k] == i and p[k+1] == j:
                adj_T.append((p, k))
                break

    # adj'(j,i) in T': paths using arc j->i consecutively
    adj_Tp = []
    for p in paths_Tp:
        for k in range(n-1):
            if p[k] == j and p[k+1] == i:
                adj_Tp.append((p, k))
                break

    # Try swap on each adj_T path
    matched_T = []
    unmatched_T = []
    for p, k in adj_T:
        # swap: (..., i, j, ...) -> (..., j, i, ...)
        can_swap = True
        blockers = []
        if k > 0 and not T[p[k-1]][j]:  # pred must beat j
            can_swap = False
            blockers.append(('pred', p[k-1]))
        if k+1 < n-1 and not T[i][p[k+2]]:  # i must beat succ
            can_swap = False
            blockers.append(('succ', p[k+2]))

        if can_swap:
            matched_T.append((p, k))
        else:
            unmatched_T.append((p, k, blockers))

    # Try swap on each adj_Tp path
    matched_Tp = []
    unmatched_Tp = []
    for p, k in adj_Tp:
        can_swap = True
        blockers = []
        if k > 0 and not Tp[p[k-1]][i]:  # pred must beat i (in T')
            can_swap = False
            blockers.append(('pred', p[k-1]))
        if k+1 < n-1 and not Tp[j][p[k+2]]:  # j must beat succ (in T')
            can_swap = False
            blockers.append(('succ', p[k+2]))

        if can_swap:
            matched_Tp.append((p, k))
        else:
            unmatched_Tp.append((p, k, blockers))

    return {
        's': s,
        'adj_T': len(adj_T),
        'adj_Tp': len(adj_Tp),
        'matched_T': len(matched_T),
        'matched_Tp': len(matched_Tp),
        'unmatched_T': unmatched_T,
        'unmatched_Tp': unmatched_Tp,
        'delta_H': hamiltonian_path_count(Tp) - hamiltonian_path_count(T),
    }


def main():
    print("=== Unmatched Path Decomposition by Blocking Vertex ===\n")

    # n=5 exhaustive
    print("--- n=5 exhaustive ---")
    n = 5
    per_vertex_data = []

    for T in all_tournaments(n):
        for i in range(n):
            for j in range(i+1, n):
                if not T[i][j]:
                    continue

                res = analyze_flip(T, i, j)
                s = res['s']
                uT = res['unmatched_T']
                uTp = res['unmatched_Tp']

                # Decompose unmatched by blocking vertex
                # A path can be blocked by pred, succ, or both
                # Count per blocking vertex x
                uT_by_x = defaultdict(int)
                uT_both = 0
                for p, k, blockers in uT:
                    if len(blockers) == 2:
                        uT_both += 1
                    for role, x in blockers:
                        uT_by_x[x] += 1

                uTp_by_x = defaultdict(int)
                uTp_both = 0
                for p, k, blockers in uTp:
                    if len(blockers) == 2:
                        uTp_both += 1
                    for role, x in blockers:
                        uTp_by_x[x] += 1

                others = [x for x in range(n) if x != i and x != j]

                for x in others:
                    B_x = [v for v in range(n) if v != i and v != j and v != x]
                    # Build sub-tournament on B_x
                    k = len(B_x)
                    sub = [[0]*k for _ in range(k)]
                    for a in range(k):
                        for b in range(k):
                            sub[a][b] = T[B_x[a]][B_x[b]]
                    h_bx = hamiltonian_path_count(sub)

                    per_vertex_data.append({
                        's_x': s[x],
                        'uT_x': uT_by_x.get(x, 0),
                        'uTp_x': uTp_by_x.get(x, 0),
                        'H_Bx': h_bx,
                    })

    # Analyze: for s_x = -1, is uT_x = 2*H(B_x)?
    # For s_x = +1, is uTp_x = 2*H(B_x)?
    print(f"Total (T, i->j, x) triples: {len(per_vertex_data)}")

    for sv in [-1, 0, 1]:
        subset = [d for d in per_vertex_data if d['s_x'] == sv]
        if not subset:
            continue
        print(f"\n  s_x = {sv:+d}: {len(subset)} cases")

        if sv == -1:
            matches = sum(1 for d in subset if d['uT_x'] == 2 * d['H_Bx'])
            print(f"    uT_x == 2*H(B_x): {matches}/{len(subset)}")
            ratios = set()
            for d in subset:
                if d['H_Bx'] > 0:
                    ratios.add((d['uT_x'], d['H_Bx']))
            print(f"    (uT_x, H_Bx) values: {sorted(ratios)}")

        if sv == 1:
            matches = sum(1 for d in subset if d['uTp_x'] == 2 * d['H_Bx'])
            print(f"    uTp_x == 2*H(B_x): {matches}/{len(subset)}")
            ratios = set()
            for d in subset:
                if d['H_Bx'] > 0:
                    ratios.add((d['uTp_x'], d['H_Bx']))
            print(f"    (uTp_x, H_Bx) values: {sorted(ratios)}")

        if sv == 0:
            # s_x=0: should have no contribution
            ut_nonzero = sum(1 for d in subset if d['uT_x'] > 0)
            utp_nonzero = sum(1 for d in subset if d['uTp_x'] > 0)
            print(f"    uT_x > 0: {ut_nonzero}/{len(subset)}")
            print(f"    uTp_x > 0: {utp_nonzero}/{len(subset)}")

    # n=6 random sample
    print("\n\n--- n=6 random sample (200 tournaments) ---")
    n = 6
    rng = random.Random(42)
    per_vertex_data_6 = []
    double_blocked_count = 0
    total_unmatched = 0

    for _ in range(200):
        T = random_tournament(n, rng)
        i, j = rng.sample(range(n), 2)
        if not T[i][j]:
            i, j = j, i

        res = analyze_flip(T, i, j)
        s = res['s']
        uT = res['unmatched_T']
        uTp = res['unmatched_Tp']

        # Count double-blocked
        for p, k, blockers in uT:
            total_unmatched += 1
            if len(blockers) == 2:
                double_blocked_count += 1
        for p, k, blockers in uTp:
            total_unmatched += 1
            if len(blockers) == 2:
                double_blocked_count += 1

        others = [x for x in range(n) if x != i and x != j]

        # Inclusion-exclusion: #U_T = sum_x uT_x(single) but need to handle double-blocking
        uT_by_x = defaultdict(int)
        for p, k, blockers in uT:
            for role, x in blockers:
                uT_by_x[x] += 1

        uTp_by_x = defaultdict(int)
        for p, k, blockers in uTp:
            for role, x in blockers:
                uTp_by_x[x] += 1

        for x in others:
            B_x = [v for v in range(n) if v != i and v != j and v != x]
            k = len(B_x)
            sub = [[0]*k for _ in range(k)]
            for a in range(k):
                for b in range(k):
                    sub[a][b] = T[B_x[a]][B_x[b]]
            h_bx = hamiltonian_path_count(sub)

            per_vertex_data_6.append({
                's_x': s[x],
                'uT_x': uT_by_x.get(x, 0),
                'uTp_x': uTp_by_x.get(x, 0),
                'H_Bx': h_bx,
            })

    print(f"Total unmatched paths: {total_unmatched}")
    print(f"Double-blocked: {double_blocked_count} ({100*double_blocked_count/max(1,total_unmatched):.1f}%)")

    for sv in [-1, 0, 1]:
        subset = [d for d in per_vertex_data_6 if d['s_x'] == sv]
        if not subset:
            continue
        print(f"\n  s_x = {sv:+d}: {len(subset)} cases")

        if sv == -1:
            matches = sum(1 for d in subset if d['uT_x'] == 2 * d['H_Bx'])
            print(f"    uT_x == 2*H(B_x): {matches}/{len(subset)}")
            # Show distribution of uT_x - 2*H_Bx
            diffs = [d['uT_x'] - 2*d['H_Bx'] for d in subset]
            print(f"    uT_x - 2*H(B_x) range: [{min(diffs)}, {max(diffs)}]")
            print(f"    mean: {sum(diffs)/len(diffs):.3f}")

        if sv == 1:
            matches = sum(1 for d in subset if d['uTp_x'] == 2 * d['H_Bx'])
            print(f"    uTp_x == 2*H(B_x): {matches}/{len(subset)}")
            diffs = [d['uTp_x'] - 2*d['H_Bx'] for d in subset]
            print(f"    uTp_x - 2*H(B_x) range: [{min(diffs)}, {max(diffs)}]")
            print(f"    mean: {sum(diffs)/len(diffs):.3f}")

        if sv == 0:
            ut_nonzero = sum(1 for d in subset if d['uT_x'] > 0)
            utp_nonzero = sum(1 for d in subset if d['uTp_x'] > 0)
            print(f"    uT_x > 0: {ut_nonzero}/{len(subset)}")
            print(f"    uTp_x > 0: {utp_nonzero}/{len(subset)}")

    # Key test: does #U_T - #U_T' = -2*sum(s_x*H(B_x))?
    # (ignoring 5-cycle terms for now)
    print("\n\n--- Aggregate identity test (n=6) ---")
    n = 6
    ok = 0
    total = 0
    for _ in range(200):
        T = random_tournament(n, rng)
        i, j = rng.sample(range(n), 2)
        if not T[i][j]:
            i, j = j, i

        res = analyze_flip(T, i, j)
        nU_T = len(res['unmatched_T'])
        nU_Tp = len(res['unmatched_Tp'])
        delta_actual = nU_T - nU_Tp  # should be -delta_H (since delta_H = adj_Tp - adj_T but adj = matched + unmatched...)

        # Wait, let's be careful:
        # adj(i,j) = matched_T + #U_T  (paths in T using i->j)
        # adj'(j,i) = matched_Tp + #U_Tp  (paths in T' using j->i)
        # matched_T = matched_Tp (by involution)
        # So delta_H = H(T') - H(T) = ... hmm
        # Actually: delta_H = adj'(j,i) - adj(i,j) + (common paths)
        # No: H(T) = adj(i,j) + non_adj(i,j) where non_adj = paths NOT using i->j
        # H(T') = adj'(j,i) + non_adj'(j,i) where non_adj' = paths of T' NOT using j->i
        # non_adj paths are the same in T and T' (they don't use the flipped arc at all!)
        # Wait, that's not right either. A T-path might use j->i (which doesn't exist in T),
        # but T-paths only use arcs in T. So non_adj_T = T-paths not using i->j.
        # Some of these may fail in T' (if they relied on i->j... no they don't use i->j).
        # But they might use j->... arcs differently. Hmm, T and T' differ only in arcs i->j vs j->i.
        # A T-path not using i->j: does not have i immediately before j.
        # Is it a T'-path? It must only use arcs of T'. The only arc that changed is i->j -> j->i.
        # So if a T-path uses i->j somewhere (i immediately before j), it's NOT in non_adj.
        # If it doesn't use i->j, then it uses the same arcs as T' (since all other arcs are same).
        # WAIT: a T-path might have j immediately before i. That arc exists in T (since T[j][i]=0 would mean no).
        # Actually T[j][i] = 0 in T (since T[i][j]=1). So no T-path has j immediately before i.
        # But a T-path that doesn't use i->j: it can still have arcs like i->k or k->j.
        # Those arcs are the same in T and T'. So yes, non_adj_T paths are all valid T'-paths!
        # Similarly non_adj_{T'} paths (not using j->i) are valid T-paths.

        # So: non_adj_T subset of T'-paths, non_adj_{T'} subset of T-paths.
        # Moreover non_adj_T = non_adj_{T'} (same set of paths).
        # Proof: P is non_adj_T iff P is a Ham path using only arcs from T and not i->j.
        #        Arcs of T minus {i->j} = arcs of T' minus {j->i}. So these are the same.
        # Therefore H(T) - H(T') = adj_T(i,j) - adj_{T'}(j,i)
        # and delta_H = H(T') - H(T) = adj_{T'}(j,i) - adj_T(i,j)
        #            = (matched_{T'} + #U_{T'}) - (matched_T + #U_T)
        #            = #U_{T'} - #U_T    (since matched are equal)

        delta_H = res['delta_H']
        computed = nU_Tp - nU_T  # should equal delta_H

        s = res['s']
        others = [x for x in range(n) if x != i and x != j]
        formula = 0
        for x in others:
            B_x = [v for v in range(n) if v != i and v != j and v != x]
            k = len(B_x)
            sub = [[0]*k for _ in range(k)]
            for a in range(k):
                for b in range(k):
                    sub[a][b] = T[B_x[a]][B_x[b]]
            h_bx = hamiltonian_path_count(sub)
            formula += -2 * s[x] * h_bx

        total += 1
        if computed == delta_H:
            ok += 1
        else:
            print(f"  BUG: computed={computed}, delta_H={delta_H}")

    print(f"#U_T' - #U_T == delta_H: {ok}/{total}")

    # Now test: #U_T' - #U_T == -2*sum(s_x*H(B_x)) + 5-cycle_correction
    # For n=5, no 5-cycles matter. Let's do n=5 exhaustive.
    print("\n\n--- n=5 exhaustive: #U_T' - #U_T vs -2*sum(s_x*H(B_x)) ---")
    n = 5
    ok_formula = 0
    total_5 = 0

    for T in all_tournaments(n):
        for i in range(n):
            for j in range(i+1, n):
                if not T[i][j]:
                    continue

                Tp = flip_arc(T, i, j)
                res = analyze_flip(T, i, j)
                nU_T = len(res['unmatched_T'])
                nU_Tp = len(res['unmatched_Tp'])
                delta = nU_Tp - nU_T

                s = res['s']
                others = [x for x in range(n) if x != i and x != j]
                formula = 0
                for x in others:
                    B_x = [v for v in range(n) if v != i and v != j and v != x]
                    k = len(B_x)
                    sub = [[0]*k for _ in range(k)]
                    for a in range(k):
                        for b in range(k):
                            sub[a][b] = T[B_x[a]][B_x[b]]
                    h_bx = hamiltonian_path_count(sub)
                    formula += -2 * s[x] * h_bx

                total_5 += 1
                if delta == formula:
                    ok_formula += 1

    print(f"Match: {ok_formula}/{total_5}")
    if ok_formula == total_5:
        print("PERFECT! At n=5: delta_H = -2*sum(s_x*H(B_x)) exactly.")
        print("This means the 5-cycle correction is zero (all 5-cycles cancel).")
        print("And the per-vertex decomposition #U_T' - #U_T = -2*sum s_x*H(B_x) holds.")

    # Per-vertex test at n=5: is each vertex's contribution exact?
    print("\n\n--- Per-vertex test at n=5 ---")
    print("For s_x=-1: is #(unmatched T-paths blocked by x) = 2*H(B_x)?")
    print("For s_x=+1: is #(unmatched T'-paths blocked by x) = 2*H(B_x)?")

    # Need to be more careful about double-blocking
    # Let's count: for each x with s_x=-1,
    # f_T(x) = #{T-paths (...,x,i,j,...)} + #{T-paths (...,i,j,x,...) where x blocks}
    # More precisely: paths where x is the predecessor of i (and x doesn't beat j)
    #                 or x is the successor of j (and i doesn't beat x)

    # Actually at n=5, B_x has 2 vertices. H(B_x) = 1 (transitive) or 1 (only 2 vertices).
    # Wait: 2-vertex tournament has exactly 1 Ham path. So H(B_x) = 1 always at n=5.
    # So the formula becomes: delta_H = -2*sum(s_x) = 2*(#gained - #lost) for 3-cycles
    # which is the simple formula (T028).

    print("At n=5, H(B_x)=1 always (2-vertex tournaments), so formula = -2*sum(s_x).")
    print("This is equivalent to the simple formula delta=2*(gained-lost).")


if __name__ == "__main__":
    main()
