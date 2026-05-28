#!/usr/bin/env python3
"""
king_tight_bound_analysis.py
opus-2026-05-27-S1

Deep analysis of the king H-increment bound tightness,
the "tight iff non-SC" conjecture at n>=5,
and related rigorous proofs.

Key finding from king_vertex_analysis.py:
  n=5: ALL tight cases (delta=2*|rivals|) are NON-strongly-connected.
  This suggests: T strongly connected, n>=5, Q=max-degree => H(T)-H(T-Q) > 2|N^-(Q)|.
"""

import itertools
from collections import defaultdict

# ---------- core utilities ----------

def all_tournaments(n):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(1 << len(arcs)):
        T = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arcs):
            if bits & (1 << k):
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

def outdegree(T, v):
    return sum(T[v])

def is_strongly_connected(T):
    n = len(T)
    def reach(src):
        vis = {src}
        q = [src]
        while q:
            u = q.pop()
            for w in range(n):
                if T[u][w] and w not in vis:
                    vis.add(w)
                    q.append(w)
        return vis
    return len(reach(0)) == n and all(0 in reach(v) for v in range(n))

def hamilton_paths(T):
    n = len(T)
    count = 0
    for perm in itertools.permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            count += 1
    return count

def tournament_minus_v(T, v):
    n = len(T)
    verts = [u for u in range(n) if u != v]
    m = len(verts)
    Tv = [[0]*m for _ in range(m)]
    for i, u in enumerate(verts):
        for j, w in enumerate(verts):
            Tv[i][j] = T[u][w]
    return Tv

def count_odd_cycles_through(T, v, length):
    """Count directed cycles of given length passing through vertex v."""
    n = len(T)
    count = 0
    others = [u for u in range(n) if u != v]
    for seq in itertools.permutations(others, length-1):
        # Cycle: v -> seq[0] -> seq[1] -> ... -> seq[-1] -> v
        path = [v] + list(seq) + [v]
        if all(T[path[i]][path[i+1]] for i in range(len(path)-1)):
            count += 1
    return count // length  # each cycle counted 'length' times

def all_odd_cycles_through(T, v):
    """Collect all odd cycles through v (by length)."""
    n = len(T)
    by_length = {}
    for L in range(3, n+1, 2):
        c = count_odd_cycles_through(T, v, L)
        if c > 0:
            by_length[L] = c
    return by_length

def three_cycles_through(T, v):
    """Count 3-cycles through v."""
    n = len(T)
    count = 0
    for a in range(n):
        if a == v or not T[v][a]:
            continue
        for b in range(n):
            if b == v or b == a or not T[b][v] or not T[a][b]:
                continue
            count += 1
    return count

def score_seq(T):
    return tuple(sorted(outdegree(T, v) for v in range(len(T))))

# ---------- analysis ----------

def analyze_tight_bound_deeply(n_max=7):
    """
    For each (n, T, Q):
      - compute H(T), H(T-Q), delta, bound = 2|rivals|
      - classify: tight vs not, SC vs not
      - when tight: investigate WHY (which mu(C) contributions)
    """
    print(f"\n{'='*70}")
    print("DEEP ANALYSIS: Tight Bound Characterization")
    print(f"{'='*70}")

    for n in range(3, n_max+1):
        tight_sc = []
        tight_nsc = []
        nontight_sc = []
        nontight_nsc = []

        for T in all_tournaments(n):
            degs = [outdegree(T, v) for v in range(n)]
            max_deg = max(degs)
            sc = is_strongly_connected(T)
            H_T = hamilton_paths(T)

            for Q in [v for v in range(n) if degs[v] == max_deg]:
                Tv = tournament_minus_v(T, Q)
                H_Tv = hamilton_paths(Tv)
                rivals = sum(1 for v in range(n) if v != Q and T[v][Q])
                delta = H_T - H_Tv
                bound = 2 * rivals
                excess = delta - bound

                entry = {
                    'sc': sc,
                    'score': score_seq(T),
                    'deg_Q': max_deg,
                    'n_rivals': rivals,
                    'H_T': H_T,
                    'delta': delta,
                    'bound': bound,
                    'excess': excess
                }

                if delta == bound:
                    if sc:
                        tight_sc.append(entry)
                    else:
                        tight_nsc.append(entry)
                else:
                    if sc:
                        nontight_sc.append(entry)
                    else:
                        nontight_nsc.append(entry)

        total = len(tight_sc) + len(tight_nsc) + len(nontight_sc) + len(nontight_nsc)
        print(f"\nn={n}: {total} (T,Q) pairs")
        print(f"  Tight+SC:     {len(tight_sc):5d}  |  Non-tight+SC:  {len(nontight_sc):6d}")
        print(f"  Tight+nonSC:  {len(tight_nsc):5d}  |  Non-tight+nSC: {len(nontight_nsc):6d}")
        print(f"  Conjecture 'tight => non-SC': {'HOLDS' if len(tight_sc)==0 else 'FAILS'}")

        # For tight+SC cases (if any), show what they look like
        if tight_sc:
            print(f"  TIGHT+SC cases:")
            for c in tight_sc[:5]:
                print(f"    score={c['score']}, deg_Q={c['deg_Q']}, n_rivals={c['n_rivals']}, H={c['H_T']}")

        # What scores achieve tight+SC?
        if tight_sc and n >= 4:
            sc_tight_scores = defaultdict(int)
            for c in tight_sc:
                sc_tight_scores[c['score']] += 1
            for ss, cnt in sorted(sc_tight_scores.items()):
                print(f"    tight+SC score {ss}: {cnt}")

        # For non-tight cases, show excess distribution
        all_excesses = [c['excess'] for c in nontight_sc + nontight_nsc]
        if all_excesses:
            print(f"  Non-tight excess: min={min(all_excesses)}, max={max(all_excesses)}, avg={sum(all_excesses)/len(all_excesses):.2f}")

        # What is the minimum excess for strongly connected tournaments?
        sc_excesses = [c['excess'] for c in nontight_sc]
        if sc_excesses:
            print(f"  SC excess: min={min(sc_excesses)}, avg={sum(sc_excesses)/len(sc_excesses):.2f}")

def analyze_tight_mechanism(n=5):
    """
    For tight cases at n=5: explain WHY delta = 2*|rivals| exactly.
    The delta = H(T) - H(T-Q) = 2 * sum_C mu(C) where C runs over odd cycles through Q.
    Tight means: sum_C mu(C) = |rivals| exactly.
    Since each rival must be in >= 1 triangle with Q, and each triangle gives mu >= 1:
    Tight iff: (1) only 3-cycles through Q (no longer cycles), and
               (2) each rival in EXACTLY 1 triangle with Q, and mu(C) = 1 for each.
    """
    print(f"\n{'='*70}")
    print(f"TIGHT MECHANISM ANALYSIS (n={n})")
    print(f"{'='*70}")

    # For n=5, we can compute mu(C) directly from the conflict graph
    # mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)
    # For a 3-cycle C=(v,a,b): C\{v}={a,b}; T-v restricted to cycles
    # avoiding {a,b}. At n=5, T-v has 4 vertices, so cycles avoiding {a,b}
    # live in the 2-vertex sub-tournament of the remaining 2 vertices.
    # With only 2 vertices, there are no odd cycles. So mu(C) = I({},2) = 1 for ALL 3-cycles!
    # This means: at n=5, all 3-cycles through Q have mu=1.
    # Also: 5-cycles through Q are possible. For 5-cycle C through v in n=5 tournament:
    # C uses all 5 vertices. So C\{v} = 4 vertices = all of T-v. Avoiding all 4 vertices
    # means no vertices left, so Omega(T-v)|_{avoid C\{v}} = empty graph. mu(5-cycle) = 1 too.
    # Wait: mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2). If avoiding C\{v} means restricting
    # to odd cycles in T-v that are vertex-disjoint from C\{v}, then for a 5-cycle at n=5:
    # C\{v} has 4 vertices = all of T-v. So no odd cycle in T-v can avoid all 4 vertices.
    # mu(5-cycle) = I(empty, 2) = 1.
    #
    # So at n=5: every odd cycle C through v has mu(C) = 1.
    # Hence delta = 2 * (# odd cycles through Q) = 2 * c_odd(Q).
    # Bound = 2 * |rivals|.
    # Tight iff c_odd(Q) = |rivals|.
    # Since c_odd(Q) >= |rivals| (each rival in >=1 triangle), tight iff
    # c_odd(Q) = |rivals|, i.e., each rival in EXACTLY 1 odd cycle with Q.

    print("At n=5: mu(C) = 1 for ALL odd cycles C through any vertex Q.")
    print("  Reason: T-v has 4 vertices. Any 3-cycle C's {C\\{v}} avoids 2 vertices,")
    print("  leaving 2 vertices in T-v -- too few for odd cycles. So mu(3-cycle) = 1.")
    print("  For 5-cycle: C\\{v} = 4 vertices = all of T-v. mu(5-cycle) = 1.")
    print("  For 7-cycle (not possible at n=5 -- max cycle length = n = 5).")
    print()
    print("Therefore at n=5: delta = 2 * c_odd(Q) where c_odd(Q) = #{odd cycles through Q}.")
    print("Tight (delta = 2|rivals|) iff c_odd(Q) = |rivals|.")
    print("Since c_odd(Q) >= |rivals| by King theorem, tight iff each rival in exactly 1 odd cycle with Q.")
    print()

    # Verify this at n=5
    for T in all_tournaments(n):
        degs = [outdegree(T, v) for v in range(n)]
        max_deg = max(degs)
        sc = is_strongly_connected(T)
        H_T = hamilton_paths(T)

        for Q in [v for v in range(n) if degs[v] == max_deg]:
            Tv = tournament_minus_v(T, Q)
            H_Tv = hamilton_paths(Tv)
            rivals = [v for v in range(n) if v != Q and T[v][Q]]
            delta = H_T - H_Tv
            bound = 2 * len(rivals)

            # Count all odd cycles through Q
            c_odd = 0
            for L in range(3, n+1, 2):
                c_odd += count_odd_cycles_through(T, Q, L)

            expected_delta = 2 * c_odd
            if expected_delta != delta:
                print(f"  DISCREPANCY: delta={delta}, 2*c_odd={expected_delta}, score={score_seq(T)}")

    print("Verification at n=5: delta = 2*c_odd(Q) for all (T,Q) pairs? Checking...")
    all_ok = True
    for T in all_tournaments(n):
        degs = [outdegree(T, v) for v in range(n)]
        max_deg = max(degs)
        for Q in [v for v in range(n) if degs[v] == max_deg]:
            Tv = tournament_minus_v(T, Q)
            H_Tv = hamilton_paths(Tv)
            delta = hamilton_paths(T) - H_Tv
            c_odd = sum(count_odd_cycles_through(T, Q, L) for L in range(3, n+1, 2))
            if delta != 2*c_odd:
                all_ok = False
                print(f"  FAIL: delta={delta}, 2c_odd={2*c_odd}")
    print(f"  Result: {'ALL OK - delta = 2*c_odd(Q) for all n=5 cases' if all_ok else 'FAILURES FOUND'}")
    print()

    # For larger n: when does mu(C) > 1?
    print("When does mu(C) > 1 first appear?")
    print("  mu(C) > 1 requires: in T-v, exists an odd cycle avoiding C\\{v}.")
    print("  For 3-cycle C=(v,a,b): need odd cycle in T-v avoiding {a,b}.")
    print("  T-v has n-1 vertices. Odd cycle avoiding {a,b} needs >= 3 remaining vertices.")
    print("  So n-1-2 = n-3 remaining vertices, need n-3 >= 3, i.e., n >= 6.")
    print("  At n=6: first possible mu(3-cycle) > 1.")
    print()

    # Verify for n=6: what is the range of mu(C) for 3-cycles?
    print("For n=6: checking range of mu(3-cycles)...")
    mu_vals_seen = set()
    for T in all_tournaments(6):
        v = 0
        for a in range(6):
            if a == v or not T[v][a]:
                continue
            for b in range(6):
                if b == v or b == a or not T[b][v] or not T[a][b]:
                    continue
                # 3-cycle (v, a, b): mu = I(Omega(T-v)|_{avoid {a,b}}, 2)
                # T-v has vertices {1,2,3,4,5} (excluding v=0)
                # Find odd cycles in T-v avoiding {a,b}
                Tv = tournament_minus_v(T, v)
                verts_tv = [u for u in range(6) if u != v]
                # Reindex: a and b in T map to positions in Tv
                a_idx = verts_tv.index(a)
                b_idx = verts_tv.index(b)
                avoid = {a_idx, b_idx}
                remaining = [i for i in range(5) if i not in avoid]
                # Find odd cycles in Tv restricted to 'remaining' vertices
                if len(remaining) < 3:
                    mu = 1
                else:
                    sub_n = len(remaining)
                    sub_T = [[0]*sub_n for _ in range(sub_n)]
                    for i, u in enumerate(remaining):
                        for j, w in enumerate(remaining):
                            sub_T[i][j] = Tv[u][w]
                    # Count independent sets weighted by 2^k in conflict graph
                    # of odd cycles in sub_T
                    sub_cycles = find_odd_cycles(sub_T)
                    mu = independence_poly_at_2_simple(sub_cycles)
                mu_vals_seen.add(mu)
        break  # just do first tournament to check
    print(f"  mu values for first tournament: {sorted(mu_vals_seen)}")

    # Now check all n=6 tournaments
    mu_vals_all = set()
    for T in all_tournaments(6):
        v = 0
        for a in range(6):
            if a == v or not T[v][a]:
                continue
            for b in range(6):
                if b == v or b == a or not T[b][v] or not T[a][b]:
                    continue
                Tv = tournament_minus_v(T, v)
                verts_tv = [u for u in range(6) if u != v]
                a_idx = verts_tv.index(a)
                b_idx = verts_tv.index(b)
                avoid = {a_idx, b_idx}
                remaining = [i for i in range(5) if i not in avoid]
                if len(remaining) < 3:
                    mu = 1
                else:
                    sub_n = len(remaining)
                    sub_T = [[0]*sub_n for _ in range(sub_n)]
                    for i, u in enumerate(remaining):
                        for j, w in enumerate(remaining):
                            sub_T[i][j] = Tv[u][w]
                    sub_cycles = find_odd_cycles(sub_T)
                    mu = independence_poly_at_2_simple(sub_cycles)
                mu_vals_all.add(mu)
    print(f"  mu values for all n=6 tournaments: {sorted(mu_vals_all)}")

def find_odd_cycles(T):
    """Find all directed odd cycles in T (as frozensets)."""
    n = len(T)
    cycles = set()
    def dfs(start, v, path, pathset):
        for w in range(n):
            if T[v][w]:
                if w == start and len(path) >= 3 and len(path) % 2 == 1:
                    cycles.add(frozenset(path))
                elif w not in pathset:
                    pathset.add(w)
                    path.append(w)
                    dfs(start, w, path, pathset)
                    path.pop()
                    pathset.remove(w)
    for s in range(n):
        dfs(s, s, [s], {s})
    return list(cycles)

def independence_poly_at_2_simple(cycles):
    """Evaluate I(conflict_graph(cycles), 2)."""
    if not cycles:
        return 1
    m = len(cycles)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True
    total = 0
    for mask in range(1 << m):
        subset = [i for i in range(m) if mask & (1 << i)]
        ok = True
        for a in range(len(subset)):
            for b in range(a+1, len(subset)):
                if adj[subset[a]][subset[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += 2**len(subset)
    return total

def analyze_sc_strong_excess(n_max=6):
    """
    For strongly connected T with n>=5:
    What is the minimum excess delta - 2*|rivals|?
    Is it always >= 2? >= 4? Characterize the minimum.
    """
    print(f"\n{'='*70}")
    print("SC EXCESS ANALYSIS: minimum of (H(T)-H(T-Q)) - 2|rivals| for SC tournaments")
    print(f"{'='*70}")

    for n in range(3, n_max+1):
        min_excess = None
        min_excess_cases = []

        for T in all_tournaments(n):
            if not is_strongly_connected(T):
                continue
            degs = [outdegree(T, v) for v in range(n)]
            max_deg = max(degs)
            H_T = hamilton_paths(T)

            for Q in [v for v in range(n) if degs[v] == max_deg]:
                Tv = tournament_minus_v(T, Q)
                H_Tv = hamilton_paths(Tv)
                rivals = sum(1 for v in range(n) if v != Q and T[v][Q])
                delta = H_T - H_Tv
                excess = delta - 2*rivals

                if min_excess is None or excess < min_excess:
                    min_excess = excess
                    min_excess_cases = [{
                        'score': score_seq(T),
                        'deg_Q': max_deg,
                        'n_rivals': rivals,
                        'H_T': H_T,
                        'delta': delta,
                        'excess': excess
                    }]
                elif excess == min_excess:
                    min_excess_cases.append({
                        'score': score_seq(T),
                        'deg_Q': max_deg,
                        'n_rivals': rivals,
                        'H_T': H_T,
                        'delta': delta,
                        'excess': excess
                    })

        print(f"\n  n={n}: min excess for SC = {min_excess}")
        for c in min_excess_cases[:5]:
            print(f"    score={c['score']}, deg_Q={c['deg_Q']}, n_rivals={c['n_rivals']}, "
                  f"H={c['H_T']}, delta={c['delta']}")
        if len(min_excess_cases) > 5:
            print(f"    ... ({len(min_excess_cases)} total cases)")

def analyze_qp_gap_exact(n_max=6):
    """
    Q-P gap and H: is H a function of Q-P gap?
    At n=3,4: perfect correlation (gap determines H).
    At n>=5: check if within each gap class, H is constant.
    """
    print(f"\n{'='*70}")
    print("Q-P GAP vs H: Does gap determine H?")
    print(f"{'='*70}")

    for n in range(3, n_max+1):
        by_gap = defaultdict(lambda: defaultdict(int))

        for T in all_tournaments(n):
            degs = [outdegree(T, v) for v in range(n)]
            gap = max(degs) - min(degs)
            H = hamilton_paths(T)
            by_gap[gap][H] += 1

        print(f"\nn={n}:")
        for gap in sorted(by_gap.keys()):
            h_dist = by_gap[gap]
            h_vals = sorted(h_dist.keys())
            if len(h_vals) == 1:
                print(f"  gap={gap}: H={h_vals[0]} (unique, {sum(h_dist.values())} tournaments)")
            else:
                print(f"  gap={gap}: H ∈ {h_vals}, counts={[h_dist[h] for h in h_vals]}")

def analyze_kings_and_h(n_max=7):
    """
    How many king vertices does T have, and how does it relate to H?
    Claim: #kings = n iff regular (odd n) or... special (even n)?
    Also: is there a monotone relationship between #kings and H?
    """
    print(f"\n{'='*70}")
    print("NUMBER OF KING VERTICES vs H")
    print(f"{'='*70}")

    for n in range(3, n_max+1):
        by_kings = defaultdict(list)

        for T in all_tournaments(n):
            num_kings = sum(1 for v in range(n) if is_king_v(T, v))
            H = hamilton_paths(T)
            sc = is_strongly_connected(T)
            by_kings[num_kings].append((H, sc))

        print(f"\nn={n}:")
        for k in sorted(by_kings.keys()):
            entries = by_kings[k]
            hs = [e[0] for e in entries]
            scs = sum(1 for e in entries if e[1])
            print(f"  #kings={k}: {len(entries)} tourn, H∈[{min(hs)},{max(hs)}], avg_H={sum(hs)/len(hs):.2f}, "
                  f"SC={scs}")

def is_king_v(T, v):
    n = len(T)
    reachable = {v}
    for w in range(n):
        if T[v][w]:
            reachable.add(w)
    for w in list(reachable):
        if w != v:
            for u in range(n):
                if T[w][u]:
                    reachable.add(u)
    return len(reachable) == n

def analyze_sc_tiling_count(n_max=8):
    """
    Count SC tilings and compute the non-SC fraction precisely.
    Check if there's a pattern / closed form.
    """
    print(f"\n{'='*70}")
    print("SC TILING COUNTS AND PATTERNS")
    print(f"{'='*70}")

    def tiles_for_cut_n(n, k):
        tiles = []
        idx = 0
        for y in range(n-2):
            for x in range(n-1, y+1, -1):
                if x >= k and y < k:
                    tiles.append(idx)
                idx += 1
        return tiles

    def is_sc_by_cuts(n, tile_bits):
        for k in range(1, n):
            cut_tiles = tiles_for_cut_n(n, k)
            if not any(tile_bits[i] == 1 for i in cut_tiles):
                return False
        return True

    print("n | m    | total | SC    | non-SC | frac_nonSC")
    print("-" * 55)
    for n in range(3, n_max+1):
        m = (n-1)*(n-2)//2
        total = 2**m
        sc_count = 0
        for bits in range(total):
            tile_bits = [(bits >> k) & 1 for k in range(m)]
            if is_sc_by_cuts(n, tile_bits):
                sc_count += 1
        nsc = total - sc_count
        print(f"{n} | {m:4d} | {total:5d} | {sc_count:5d} | {nsc:6d} | {nsc/total:.6f}")

    print()
    print("Non-SC counts: ", end="")
    non_sc_seq = []
    for n in range(3, n_max+1):
        m = (n-1)*(n-2)//2
        total = 2**m
        sc_count = sum(1 for bits in range(total)
                       if is_sc_by_cuts(n, [(bits>>k)&1 for k in range(m)]))
        non_sc_seq.append(total - sc_count)
    print(non_sc_seq)

    # Check the pattern in non-SC counts
    print("\nRatio of successive non-SC counts:")
    for i in range(1, len(non_sc_seq)):
        print(f"  n={i+3}/n={i+2}: {non_sc_seq[i]/non_sc_seq[i-1]:.4f}")

def analyze_apex_decrease_pattern(n_max=6):
    """
    When does flipping the apex tile DECREASE H?
    Expected: only when T is 'near-regular' (high H), so the apex
    tile's default orientation is 'correct' for high H.
    """
    print(f"\n{'='*70}")
    print("APEX TILE H-DECREASE PATTERN")
    print(f"{'='*70}")

    def tiling_to_T(n, tile_bits):
        tiles = []
        for y in range(n-2):
            for x in range(n-1, y+1, -1):
                tiles.append((x, y))
        T = [[0]*n for _ in range(n)]
        for k in range(1, n):
            T[k][k-1] = 1
        for idx, (x, y) in enumerate(tiles):
            if tile_bits[idx] == 0:
                T[x][y] = 1
            else:
                T[y][x] = 1
        return T

    for n in range(3, n_max+1):
        m = (n-1)*(n-2)//2
        decreases = []
        same = []
        increases = []

        for bits in range(1 << m):
            tile_bits = [(bits >> k) & 1 for k in range(m)]
            if tile_bits[0] != 0:  # apex already flipped
                continue

            T0 = tiling_to_T(n, tile_bits)
            flipped = [1 if i == 0 else tile_bits[i] for i in range(m)]
            T1 = tiling_to_T(n, flipped)

            H0 = hamilton_paths(T0)
            H1 = hamilton_paths(T1)
            delta = H1 - H0

            entry = {
                'H0': H0,
                'H1': H1,
                'delta': delta,
                'score0': score_seq(T0),
                'score1': score_seq(T1),
                'sc0': is_strongly_connected(T0)
            }

            if delta < 0:
                decreases.append(entry)
            elif delta == 0:
                same.append(entry)
            else:
                increases.append(entry)

        print(f"\nn={n}: apex flip from default (n-1->0) to upward (0->n-1)")
        print(f"  Increase: {len(increases)}, Same: {len(same)}, Decrease: {len(decreases)}")

        if decreases:
            print(f"  DECREASE cases:")
            for c in decreases[:4]:
                print(f"    H0={c['H0']}, H1={c['H1']}, delta={c['delta']}, "
                      f"score0={c['score0']}, SC0={c['sc0']}")
            # Pattern: when does flip decrease?
            h0_vals = [c['H0'] for c in decreases]
            print(f"  H0 values in decrease cases: min={min(h0_vals)}, max={max(h0_vals)}")
            # Are all decrease cases from SC tournaments?
            sc_decrease = sum(1 for c in decreases if c['sc0'])
            print(f"  SC among decrease cases: {sc_decrease}/{len(decreases)}")

        if same:
            print(f"  SAME cases: (H0=H1, no change)")
            for c in same[:3]:
                print(f"    H={c['H0']}, score={c['score0']}, SC={c['sc0']}")

if __name__ == '__main__':
    print("KING TIGHT BOUND — DEEP ANALYSIS")
    print("opus-2026-05-27-S1")
    print("=" * 70)

    analyze_tight_bound_deeply(n_max=6)
    analyze_tight_mechanism(n=5)
    analyze_sc_strong_excess(n_max=6)
    analyze_qp_gap_exact(n_max=6)
    analyze_kings_and_h(n_max=6)
    analyze_sc_tiling_count(n_max=8)
    analyze_apex_decrease_pattern(n_max=6)
