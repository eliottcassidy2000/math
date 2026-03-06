#!/usr/bin/env python3
"""
BlackSelf(8) v3: Survey all SC+Aut>1 tournaments comprehensively.

Check multiple interpretations of "H(T)/|Fix(beta)| is even":
  (A) Fix | H and H/Fix is even
  (B) (H - Fix)/2 is even (number of beta-orbit pairs is even)
  (C) (H + Fix)/2 is even (total orbit count is even)

Also: check if T has multiple conjugacy classes of anti-auts with different Fix.

Instance: opus-2026-03-05-S8
"""

from itertools import combinations, permutations
from collections import Counter


def count_ham_paths_dp(T, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1][v] for v in range(n))


def count_concordant_paths(T, n, alpha):
    half = n // 2
    orbits = []
    seen = set()
    for v in range(n):
        if v not in seen:
            orbits.append((v, alpha[v]))
            seen.add(v)
            seen.add(alpha[v])

    count = 0
    for perm in permutations(range(len(orbits))):
        for bits in range(1 << half):
            path = [0] * n
            for i in range(half):
                o = perm[i]
                path[i] = orbits[o][(bits >> i) & 1]
                path[n - 1 - i] = alpha[path[i]]
            valid = all(T[path[i]][path[i+1]] for i in range(n - 1))
            if valid:
                count += 1
    return count


def has_automorphism_beyond_identity(T, n):
    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    sg = {}
    for v, s in enumerate(scores):
        sg.setdefault(s, []).append(v)
    if all(len(g) == 1 for g in sg.values()):
        return False
    possible = [sg[scores[v]] for v in range(n)]

    def bt(v, perm, used):
        if v == n:
            return any(perm[i] != i for i in range(n))
        for t in possible[v]:
            if t in used:
                continue
            perm[v] = t
            ok = all(T[v][u] == T[t][perm[u]] and T[u][v] == T[perm[u]][t] for u in range(v))
            if ok:
                used.add(t)
                if bt(v + 1, perm, used):
                    return True
                used.remove(t)
            perm[v] = -1
        return False

    return bt(0, [-1]*n, set())


def count_automorphisms(T, n):
    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    sg = {}
    for v, s in enumerate(scores):
        sg.setdefault(s, []).append(v)
    possible = [sg[scores[v]] for v in range(n)]
    count = [0]

    def bt(v, perm, used):
        if v == n:
            count[0] += 1
            return
        for t in possible[v]:
            if t in used:
                continue
            perm[v] = t
            ok = all(T[v][u] == T[t][perm[u]] and T[u][v] == T[perm[u]][t] for u in range(v))
            if ok:
                used.add(t)
                bt(v + 1, perm, used)
                used.remove(t)
            perm[v] = -1

    bt(0, [-1]*n, set())
    return count[0]


def find_all_involutive_anti_auts(T, n):
    """Find all involutive anti-automorphisms of T."""
    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    # Anti-aut maps out-deg d to in-deg d = out-deg (n-1-d)
    score_map = {}
    for v, s in enumerate(scores):
        score_map.setdefault(n - 1 - s, []).append(v)

    results = []

    def bt(v, perm, used):
        if v == n:
            if all(perm[perm[i]] == i for i in range(n)):
                results.append(tuple(perm))
            return
        target_score = n - 1 - scores[v]
        candidates = score_map.get(scores[v], [])  # Wait, need to think about this

        # Anti-aut: T(alpha(u), alpha(v)) = T(v, u)
        # So alpha maps vertex with score s to vertex with score (n-1-s)
        # But wait: T^op(u,v) = T(v,u). In T^op, vertex u has out-deg = (n-1) - score_T(u).
        # So alpha must map vertex u (score s in T) to vertex alpha(u) that has score s in T^op,
        # i.e., score (n-1-s) in T. Wait no...
        # Actually T(alpha(u), alpha(v)) = T(v, u) means alpha is an isomorphism from T to T^op.
        # The out-degree of u in T is scores[u].
        # The out-degree of alpha(u) in T^op is scores[u] (since alpha preserves adjacency from T to T^op).
        # Out-degree of alpha(u) in T^op = in-degree of alpha(u) in T = (n-1) - scores_T[alpha(u)].
        # So (n-1) - scores_T[alpha(u)] = scores_T[u].
        # => scores_T[alpha(u)] = (n-1) - scores_T[u].

        target_score_val = (n - 1) - scores[v]
        for u in range(n):
            if u in used or scores[u] != target_score_val:
                continue
            perm[v] = u
            ok = True
            for w in range(v):
                if T[u][perm[w]] != T[w][v] or T[perm[w]][u] != T[v][w]:
                    ok = False
                    break
            if ok:
                # Check involution constraint: if perm[u] is set, must be v
                if u < v and perm[u] != -1 and perm[u] != v:
                    perm[v] = -1
                    continue
                # If u == v, it's a fixed point: need T(v,v)=T(v,v), trivially OK
                # but for tournament: if u has fixed point, T(v,w) = T(w,v) for
                # another fixed point w, impossible.
                used.add(u)
                bt(v + 1, perm, used)
                used.remove(u)
            perm[v] = -1

    bt(0, [-1]*n, set())
    return results


def build_self_converse_tournaments(n, alpha):
    free_arcs = []
    determined = set()
    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in determined:
                continue
            ai, aj = alpha[i], alpha[j]
            u, v = min(aj, ai), max(aj, ai)
            if (u, v) == (i, j):
                free_arcs.append(((i, j), None))
                determined.add((i, j))
            else:
                free_arcs.append(((i, j), (aj, ai)))
                determined.add((i, j))
                determined.add((u, v))

    num_free = len(free_arcs)
    tournaments = []
    for bits in range(1 << num_free):
        T = [[0]*n for _ in range(n)]
        for k, (primary, linked) in enumerate(free_arcs):
            val = (bits >> k) & 1
            i, j = primary
            T[i][j] = val
            T[j][i] = 1 - val
            if linked is not None:
                a, b = linked
                T[a][b] = val
                T[b][a] = 1 - val
        tournaments.append(T)
    return tournaments, num_free


def score_sequence(T, n):
    return tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)], reverse=True))


def find_all_3cycles(T, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if T[a][b] and T[b][c] and T[c][a]:
            cycles.append((a, b, c))
        elif T[a][c] and T[c][b] and T[b][a]:
            cycles.append((a, c, b))
    return cycles


def has_induced_C5(adj_sets, m):
    for v0 in range(m):
        for v1 in adj_sets[v0]:
            if v1 <= v0:
                continue
            for v2 in adj_sets[v1]:
                if v2 <= v0 or v2 in adj_sets[v0]:
                    continue
                for v3 in adj_sets[v2]:
                    if v3 <= v0 or v3 == v1 or v3 in adj_sets[v0] or v3 in adj_sets[v1]:
                        continue
                    for v4 in adj_sets[v3]:
                        if v4 <= v0 or v4 == v1 or v4 == v2:
                            continue
                        if v4 in adj_sets[v0] and v4 not in adj_sets[v1] and v4 not in adj_sets[v2]:
                            return True
    return False


def independence_poly(adj, m):
    if m > 22:
        return None
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        verts = [i for i in range(m) if mask & (1 << i)]
        independent = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            coeffs[len(verts)] += 1
    return coeffs


def mu_value(T, n, cycle_verts):
    remaining = [v for v in range(n) if v not in cycle_verts]
    if len(remaining) < 3:
        return 1
    sub_cycles = []
    for a, b, c in combinations(remaining, 3):
        if T[a][b] and T[b][c] and T[c][a]:
            sub_cycles.append((a, b, c))
        elif T[a][c] and T[c][b] and T[b][a]:
            sub_cycles.append((a, c, b))
    if not sub_cycles:
        return 1
    sm = len(sub_cycles)
    cs = [set(c) for c in sub_cycles]
    adj = [[0]*sm for _ in range(sm)]
    for i in range(sm):
        for j in range(i+1, sm):
            if cs[i] & cs[j]:
                adj[i][j] = adj[j][i] = 1
    coeffs = independence_poly(adj, sm)
    return sum(c * (2**k) for k, c in enumerate(coeffs))


def main():
    n = 8
    alpha = [7, 6, 5, 4, 3, 2, 1, 0]

    print("=" * 70)
    print("BLACKSELF(8) COMPREHENSIVE SURVEY")
    print("=" * 70)

    tournaments, num_free = build_self_converse_tournaments(n, alpha)
    print(f"Tournaments: {len(tournaments)}")

    # Collect ALL SC+Aut>1 data
    data = []
    for tidx, T in enumerate(tournaments):
        if tidx % 10000 == 0:
            print(f"  Phase 1: {tidx}/{len(tournaments)}")
        if not has_automorphism_beyond_identity(T, n):
            continue

        H = count_ham_paths_dp(T, n)
        fix = count_concordant_paths(T, n, alpha)
        aut = count_automorphisms(T, n)
        data.append((T, H, fix, aut, tidx))

    print(f"\nSC+Aut>1 tournaments: {len(data)}")
    print()

    # Survey all data points
    print("Complete survey (unique (H, Fix, Aut) combinations):")
    combos = Counter()
    for T, H, fix, aut, tidx in data:
        combos[(H, fix, aut)] += 1

    print(f"{'H':>6} {'Fix':>5} {'|Aut|':>5} {'(H-F)/2':>8} {'pairs%2':>7} {'F|H':>4} {'H/F':>5} {'Count':>5}")
    print("-" * 60)

    blackself_A = []  # Fix|H and H/Fix even
    blackself_B = []  # (H-Fix)/2 even
    blackself_C = []  # (H+Fix)/2 even

    for (H, fix, aut), cnt in sorted(combos.items()):
        pairs = (H - fix) // 2
        divides = "Y" if fix > 0 and H % fix == 0 else "N"
        ratio = H // fix if fix > 0 and H % fix == 0 else "-"
        print(f"{H:>6} {fix:>5} {aut:>5} {pairs:>8} {pairs%2:>7} {divides:>4} {str(ratio):>5} {cnt:>5}")

        if fix > 0 and H % fix == 0 and (H // fix) % 2 == 0:
            blackself_A.append((H, fix, aut, cnt))
        if pairs % 2 == 0:
            blackself_B.append((H, fix, aut, cnt))
        if (H + fix) % 4 == 0:
            blackself_C.append((H, fix, aut, cnt))

    print(f"\nInterpretation A (Fix|H, H/Fix even): {len(blackself_A)} types")
    for x in blackself_A:
        print(f"  H={x[0]}, Fix={x[1]}, |Aut|={x[2]}, count={x[3]}")

    print(f"\nInterpretation B ((H-Fix)/2 even): {len(blackself_B)} types")
    for x in blackself_B:
        print(f"  H={x[0]}, Fix={x[1]}, |Aut|={x[2]}, count={x[3]}, pairs={(x[0]-x[1])//2}")

    print(f"\nInterpretation C ((H+Fix)/2 even): {len(blackself_C)} types")
    for x in blackself_C:
        print(f"  H={x[0]}, Fix={x[1]}, |Aut|={x[2]}, count={x[3]}")

    # Now: for the most interesting candidates, check if they also have other anti-auts
    # giving different Fix values
    print(f"\n{'='*70}")
    print("Checking for alternative anti-auts on select candidates...")

    # Pick one tournament from each interesting class
    shown = set()
    for T, H, fix, aut, tidx in data:
        sig = (H, fix, aut)
        if sig in shown:
            continue
        shown.add(sig)

        if aut >= 4 or H > 100:  # Focus on interesting ones
            print(f"\n  H={H}, Fix(alpha0)={fix}, |Aut|={aut}")
            # Find all involutive anti-auts
            anti_auts = find_all_involutive_anti_auts(T, n)
            print(f"  Number of involutive anti-auts: {len(anti_auts)}")
            if len(anti_auts) > 1:
                fix_values = set()
                for aa in anti_auts[:20]:  # Cap at 20 to avoid slowness
                    f = count_concordant_paths(T, n, list(aa))
                    fix_values.add(f)
                print(f"  Distinct Fix values: {fix_values}")

            # Quick Omega check
            cycles_3 = find_all_3cycles(T, n)
            m = len(cycles_3)
            cs = [set(c) for c in cycles_3]
            adj_s = [set() for _ in range(m)]
            for i in range(m):
                for j in range(i+1, m):
                    if cs[i] & cs[j]:
                        adj_s[i].add(j)
                        adj_s[j].add(i)
            c5 = has_induced_C5(adj_s, m)
            print(f"  3-cycles: {m}, Omega has C5: {c5}")

            # mu distribution
            mu_dist = Counter()
            for c in cycles_3:
                mu_dist[mu_value(T, n, set(c))] += 1
            print(f"  mu dist: {dict(mu_dist)}")

            # Score sequence
            print(f"  Scores: {score_sequence(T, n)}")

            # Check if regular
            out_deg = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
            print(f"  Out-degrees: {out_deg}")


if __name__ == "__main__":
    main()
