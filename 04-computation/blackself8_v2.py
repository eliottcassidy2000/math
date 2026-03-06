#!/usr/bin/env python3
"""
BlackSelf(8) investigation v2.

Fix: compute Fix(beta) by directly counting concordant paths,
not via a naive quotient tournament.

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
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def count_concordant_paths(T, n, alpha):
    """Count Hamiltonian paths P where P[i] = alpha(P[n-1-i]) for all i.

    For n=8 even: first n/2 vertices determine the rest.
    Enumerate all orderings of n/2 orbit reps, check arc conditions.
    """
    half = n // 2
    # Orbits: each is {v, alpha(v)}
    orbits = []
    seen = set()
    for v in range(n):
        if v not in seen:
            orbits.append((v, alpha[v]))
            seen.add(v)
            seen.add(alpha[v])

    # A concordant path: choose one rep from each orbit for positions 0..half-1
    # Position i gets v_i, position n-1-i gets alpha(v_i)
    count = 0
    for perm in permutations(range(len(orbits))):
        # Try all 2^half choices of which representative goes in first half
        for bits in range(1 << half):
            # Build the full path
            path = [0] * n
            for i in range(half):
                orbit_idx = perm[i]
                if (bits >> i) & 1:
                    path[i] = orbits[orbit_idx][1]
                else:
                    path[i] = orbits[orbit_idx][0]
                path[n - 1 - i] = alpha[path[i]]

            # Check all consecutive arcs
            valid = True
            for i in range(n - 1):
                if not T[path[i]][path[i+1]]:
                    valid = False
                    break
            if valid:
                count += 1

    return count


def has_automorphism_beyond_identity(T, n):
    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    score_groups = {}
    for v, s in enumerate(scores):
        score_groups.setdefault(s, []).append(v)
    if all(len(g) == 1 for g in score_groups.values()):
        return False
    possible = [score_groups[scores[v]] for v in range(n)]

    def backtrack(v, perm, used):
        if v == n:
            return any(perm[i] != i for i in range(n))
        for target in possible[v]:
            if target in used:
                continue
            perm[v] = target
            ok = True
            for u in range(v):
                if T[v][u] != T[target][perm[u]] or T[u][v] != T[perm[u]][target]:
                    ok = False
                    break
            if ok:
                used.add(target)
                if backtrack(v + 1, perm, used):
                    return True
                used.remove(target)
            perm[v] = -1
        return False

    return backtrack(0, [-1]*n, set())


def count_automorphisms(T, n):
    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    score_groups = {}
    for v, s in enumerate(scores):
        score_groups.setdefault(s, []).append(v)
    possible = [score_groups[scores[v]] for v in range(n)]
    count = [0]

    def backtrack(v, perm, used):
        if v == n:
            count[0] += 1
            return
        for target in possible[v]:
            if target in used:
                continue
            perm[v] = target
            ok = True
            for u in range(v):
                if T[v][u] != T[target][perm[u]] or T[u][v] != T[perm[u]][target]:
                    ok = False
                    break
            if ok:
                used.add(target)
                backtrack(v + 1, perm, used)
                used.remove(target)
            perm[v] = -1

    backtrack(0, [-1]*n, set())
    return count[0]


def find_all_anti_automorphisms_involutive(T, n):
    """Find all involutive anti-automorphisms of T.
    alpha is involutive: alpha^2 = id.
    Anti-aut: T(alpha(u), alpha(v)) = T(v, u).
    """
    Top = [[T[j][i] for j in range(n)] for i in range(n)]
    scores_T = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    scores_Top = [sum(Top[i][j] for j in range(n) if j != i) for i in range(n)]

    # Anti-aut maps T to T^op, so alpha(v) must have same out-degree in T^op as v in T
    # out-degree in T^op = in-degree in T = n-1 - out-degree in T
    # So alpha maps vertex with out-deg d to vertex with out-deg n-1-d

    results = []

    def backtrack(v, perm, used):
        if v == n:
            # Must be involution
            if all(perm[perm[i]] == i for i in range(n)):
                results.append(tuple(perm))
            return
        target_score = n - 1 - scores_T[v]
        for u in range(n):
            if u in used:
                continue
            if scores_T[u] != target_score:
                continue
            perm[v] = u
            ok = True
            for w in range(v):
                # Check: T(perm[v], perm[w]) = T(w, v)
                if T[u][perm[w]] != T[w][v]:
                    ok = False
                    break
                if T[perm[w]][u] != T[v][w]:
                    ok = False
                    break
            if ok:
                used.add(u)
                backtrack(v + 1, perm, used)
                used.remove(u)
            perm[v] = -1

    backtrack(0, [-1]*n, set())
    return results


def build_self_converse_tournaments(n, alpha):
    """Build all tournaments with alpha as anti-automorphism."""
    # Find free arc orbits under the constraint T(alpha(u), alpha(v)) = T(v, u)
    # Equivalently T(i,j) = T(alpha(j), alpha(i))
    free_arcs = []
    determined = set()

    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in determined:
                continue
            # T(i,j) determines T(alpha(j), alpha(i))
            ai, aj = alpha[i], alpha[j]
            # T(aj, ai) = T(i, j) [since T(alpha(i), alpha(j)) = T(j, i) = 1-T(i,j),
            # so T(alpha(j), alpha(i)) = T(i, j)]
            u, v = min(aj, ai), max(aj, ai)
            if (u, v) == (i, j):
                # Self-paired: no extra constraint
                free_arcs.append(((i, j), None))
                determined.add((i, j))
            else:
                free_arcs.append(((i, j), (aj, ai)))  # Store actual direction
                determined.add((i, j))
                determined.add((u, v))

    num_free = len(free_arcs)
    tournaments = []

    for bits in range(1 << num_free):
        T = [[0]*n for _ in range(n)]
        for k, (primary, linked) in enumerate(free_arcs):
            val = (bits >> k) & 1
            i, j = primary
            if val:
                T[i][j] = 1
            else:
                T[j][i] = 1

            if linked is not None:
                # T(linked[0], linked[1]) = val (same as T(i,j))
                a, b = linked
                if val:
                    T[a][b] = 1
                else:
                    T[b][a] = 1

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


def independence_poly(adj, m):
    if m > 25:
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
                            return True, (v0, v1, v2, v3, v4)
    return False, None


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
    cycle_sets = [set(c) for c in sub_cycles]
    adj = [[0]*sm for _ in range(sm)]
    for i in range(sm):
        for j in range(i+1, sm):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = 1
    coeffs = independence_poly(adj, sm)
    return sum(c * (2**k) for k, c in enumerate(coeffs))


def main():
    n = 8
    print("=" * 70)
    print("BLACKSELF(8) INVESTIGATION v2")
    print("=" * 70)

    # Strategy: try multiple involutive anti-automorphisms
    alphas_to_try = [
        [7, 6, 5, 4, 3, 2, 1, 0],     # reversal
        [1, 0, 3, 2, 5, 4, 7, 6],     # (01)(23)(45)(67)
        [3, 2, 1, 0, 7, 6, 5, 4],     # (03)(12)(47)(56)
        [4, 5, 6, 7, 0, 1, 2, 3],     # (04)(15)(26)(37)
    ]

    all_blackself = []

    for alpha_idx, alpha in enumerate(alphas_to_try):
        # Verify it's an involution
        assert all(alpha[alpha[i]] == i for i in range(n))

        print(f"\n--- Alpha #{alpha_idx}: {alpha} ---")
        tournaments, num_free = build_self_converse_tournaments(n, alpha)
        print(f"Free bits: {num_free}, tournaments: {len(tournaments)}")

        # Verify anti-automorphism for first tournament
        T0 = tournaments[0]
        ok = all(T0[alpha[u]][alpha[v]] == T0[v][u]
                 for u in range(n) for v in range(n) if u != v)
        print(f"Anti-aut check on T0: {ok}")
        if not ok:
            print("  SKIPPING (anti-aut construction is buggy for this alpha)")
            continue

        found = 0
        for tidx, T in enumerate(tournaments):
            if tidx % 5000 == 0 and tidx > 0:
                print(f"  Progress: {tidx}/{len(tournaments)}, found {found}")

            # Quick filter: need |Aut| > 1
            if not has_automorphism_beyond_identity(T, n):
                continue

            H = count_ham_paths_dp(T, n)

            # Count concordant paths (Fix(beta))
            fix_beta = count_concordant_paths(T, n, alpha)

            # Must have Fix(beta) odd, and H/Fix(beta) even
            if fix_beta == 0 or fix_beta % 2 == 0:
                continue
            if H % fix_beta != 0:
                # This shouldn't happen if alpha is truly an anti-aut
                # But beta-orbits should partition non-fixed paths into pairs
                # so (H - fix_beta) should be even
                if (H - fix_beta) % 2 != 0:
                    print(f"  BUG: H={H}, Fix(beta)={fix_beta}, H-Fix not even!")
                continue

            ratio = H // fix_beta
            if ratio % 2 == 0:
                aut_size = count_automorphisms(T, n)
                scores = score_sequence(T, n)
                found += 1
                print(f"  *** BLACKSELF: H={H}, |Aut|={aut_size}, Fix(beta)={fix_beta}, "
                      f"ratio={ratio}, scores={scores}")
                all_blackself.append((T, H, aut_size, fix_beta, ratio, scores, alpha))

        print(f"  Found {found} with this alpha")

    print(f"\n{'='*70}")
    print(f"TOTAL BLACKSELF(8) CANDIDATES: {len(all_blackself)}")

    if not all_blackself:
        print("\nNo candidates found. Let's also scan by finding anti-auts of random SC tournaments...")
        # Alternative: generate random tournaments and check all conditions
        import random
        random.seed(42)
        print("\nBroadening search: random tournaments with |Aut| > 1...")

        found_random = 0
        for trial in range(5000):
            # Random tournament
            T = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        T[i][j] = 1
                    else:
                        T[j][i] = 1

            if not has_automorphism_beyond_identity(T, n):
                continue

            # Check if self-converse: find any anti-automorphism
            anti_auts = find_all_anti_automorphisms_involutive(T, n)
            if not anti_auts:
                continue

            H = count_ham_paths_dp(T, n)

            for alpha_tuple in anti_auts:
                alpha = list(alpha_tuple)
                fix_beta = count_concordant_paths(T, n, alpha)
                if fix_beta == 0 or fix_beta % 2 == 0:
                    continue
                if H % fix_beta != 0:
                    continue
                ratio = H // fix_beta
                if ratio % 2 == 0:
                    aut_size = count_automorphisms(T, n)
                    scores = score_sequence(T, n)
                    found_random += 1
                    print(f"  *** RANDOM BLACKSELF: H={H}, |Aut|={aut_size}, "
                          f"Fix(beta)={fix_beta}, ratio={ratio}, scores={scores}")
                    all_blackself.append((T, H, aut_size, fix_beta, ratio, scores, alpha))
                    break  # One alpha is enough

            if found_random > 0 and found_random % 5 == 0:
                print(f"  ({trial} trials, {found_random} found)")

    # Deep analysis
    if all_blackself:
        # Take first candidate for deep dive
        T, H, aut_size, fix_beta, ratio, scores, alpha = all_blackself[0]
        print(f"\n{'='*70}")
        print(f"DEEP ANALYSIS of first BlackSelf(8)")
        print(f"  H(T) = {H}")
        print(f"  |Aut(T)| = {aut_size}")
        print(f"  Fix(beta) = {fix_beta}")
        print(f"  H/Fix(beta) = {ratio}")
        print(f"  Score: {scores}")

        print(f"\n  Adjacency matrix:")
        for row in T:
            print(f"    {''.join(str(x) for x in row)}")

        out_deg = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
        print(f"  Out-degrees: {out_deg}")
        print(f"  Regular: {len(set(out_deg)) == 1}")

        cycles_3 = find_all_3cycles(T, n)
        print(f"  3-cycles: {len(cycles_3)}")

        m = len(cycles_3)
        cycle_sets = [set(c) for c in cycles_3]
        adj = [[0]*m for _ in range(m)]
        for i in range(m):
            for j in range(i+1, m):
                if cycle_sets[i] & cycle_sets[j]:
                    adj[i][j] = adj[j][i] = 1
        adj_sets = [set(j for j in range(m) if adj[i][j]) for i in range(m)]

        c5_result = has_induced_C5(adj_sets, m)
        print(f"  Omega has C5: {c5_result[0]}")

        if m <= 20:
            coeffs = independence_poly(adj, m)
            I2 = sum(c * (2**k) for k, c in enumerate(coeffs))
            print(f"  I(Omega, 2) = {I2}, H = {H}, OCF: {I2 == H}")

        mu_dist = Counter()
        for c in cycles_3:
            mu_dist[mu_value(T, n, set(c))] += 1
        print(f"  mu distribution (3-cycles): {dict(mu_dist)}")
    else:
        # Even broader: just look for SC tournaments with |Aut|>1 and analyze their beta-orbits
        print("\nLet's just survey all SC+Aut>1 tournaments for this alpha and show H, Fix(beta)...")
        alpha = [7, 6, 5, 4, 3, 2, 1, 0]
        tournaments, _ = build_self_converse_tournaments(n, alpha)

        survey = []
        for tidx, T in enumerate(tournaments):
            if not has_automorphism_beyond_identity(T, n):
                continue
            H = count_ham_paths_dp(T, n)
            fb = count_concordant_paths(T, n, alpha)
            survey.append((H, fb, H - fb, (H - fb) % 2, fb % 2))
            if len(survey) <= 20:
                print(f"  H={H}, Fix(beta)={fb}, H-Fix={H-fb}, "
                      f"(H-Fix)%2={survey[-1][3]}, Fix%2={survey[-1][4]}")

        print(f"\nSurvey of {len(survey)} SC+Aut>1 tournaments:")
        h_fix_odd = sum(1 for s in survey if s[4] == 1)
        h_fix_even = sum(1 for s in survey if s[4] == 0)
        print(f"  Fix(beta) odd: {h_fix_odd}, even: {h_fix_even}")

        # Among Fix(beta) odd: check H/Fix(beta) parity
        for s in survey:
            H, fb = s[0], s[1]
            if fb % 2 == 1 and fb > 0 and H % fb == 0 and (H // fb) % 2 == 0:
                print(f"  *** BLACKSELF: H={H}, Fix={fb}, ratio={H//fb}")


if __name__ == "__main__":
    main()
