#!/usr/bin/env python3
"""
Fast investigation of the n=8 anomaly: BlackSelf(8).

Uses efficient automorphism detection and avoids full permutation enumeration.
Instance: opus-2026-03-05-S8
"""

import sys
from itertools import permutations, combinations
from collections import Counter


def count_ham_paths(T, n):
    """Count Hamiltonian paths using DP over bitmasks."""
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def score_sequence(T, n):
    return tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)], reverse=True))


def canonical_form(T, n):
    """Simple canonical form: flatten upper triangle as tuple."""
    bits = []
    for i in range(n):
        for j in range(i+1, n):
            bits.append(T[i][j])
    return tuple(bits)


def apply_perm_and_check(T, perm, n, target):
    """Apply perm to T and check if it equals target. Returns early on mismatch."""
    for i in range(n):
        for j in range(i+1, n):
            if T[perm[i]][perm[j]] != target[i][j]:
                return False
    return True


def has_automorphism_beyond_identity(T, n):
    """Check if T has any non-trivial automorphism. Uses score-based pruning."""
    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
    score_groups = {}
    for v, s in enumerate(scores):
        score_groups.setdefault(s, []).append(v)

    # An automorphism must map vertices to same-score vertices
    # Quick check: if all scores distinct, only identity is possible
    if all(len(g) == 1 for g in score_groups.values()):
        return False

    # Build possible images for each vertex
    possible = [score_groups[scores[v]] for v in range(n)]

    # Try to find a non-trivial automorphism via backtracking
    def backtrack(v, perm, used):
        if v == n:
            # Check: is this not the identity?
            if any(perm[i] != i for i in range(n)):
                return True
            return False
        for target in possible[v]:
            if target in used:
                continue
            # Check consistency with arcs from v to all assigned vertices
            ok = True
            for u in range(v):
                if T[v][u] != T[perm[v] if v < n else target][perm[u]]:
                    # Wait, perm[v] not set yet. Let me fix.
                    pass
            # Actually, let me just check directly
            perm[v] = target
            ok = True
            for u in range(v):
                if T[v][u] != T[target][perm[u]]:
                    ok = False
                    break
                if T[u][v] != T[perm[u]][target]:
                    ok = False
                    break
            if ok:
                used.add(target)
                if backtrack(v + 1, perm, used):
                    return True
                used.remove(target)
            perm[v] = -1
        return False

    perm = [-1] * n
    return backtrack(0, perm, set())


def count_automorphisms(T, n):
    """Count all automorphisms. Score-pruned backtracking."""
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
                if T[v][u] != T[target][perm[u]]:
                    ok = False
                    break
                if T[u][v] != T[perm[u]][target]:
                    ok = False
                    break
            if ok:
                used.add(target)
                backtrack(v + 1, perm, used)
                used.remove(target)
            perm[v] = -1

    perm = [-1] * n
    backtrack(0, perm, set())
    return count[0]


def enumerate_ham_paths(T, n):
    """Enumerate all Hamiltonian paths."""
    paths = []
    def backtrack(path, visited):
        if len(path) == n:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v not in visited and T[last][v]:
                visited.add(v)
                path.append(v)
                backtrack(path, visited)
                path.pop()
                visited.remove(v)
    for start in range(n):
        backtrack([start], {start})
    return paths


def build_self_converse_tournaments(n, alpha):
    """Build all tournaments with alpha as anti-automorphism.

    Constraint: T(alpha(u), alpha(v)) = T(v, u) for all u != v.
    """
    # Find free arcs: pairs (i,j) with i<j, grouped by orbits under the constraint
    pairs_list = []
    seen = set()
    self_paired_arcs = []

    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in seen:
                continue
            # The linked arc: T(alpha(j), alpha(i)) must equal T(i, j)
            # since T(alpha(i), alpha(j)) = T(j, i) = 1 - T(i, j)
            # and T(alpha(j), alpha(i)) = T(i, j)
            ai, aj = alpha[i], alpha[j]
            # Linked pair in canonical order
            u, v = (min(aj, ai), max(aj, ai))
            linked = (u, v)

            if (i, j) == linked:
                self_paired_arcs.append((i, j))
                seen.add((i, j))
            else:
                pairs_list.append(((i, j), linked))
                seen.add((i, j))
                seen.add(linked)

    free_arcs = self_paired_arcs + [p[0] for p in pairs_list]
    num_free = len(free_arcs)

    tournaments = []
    for bits in range(2**num_free):
        T = [[0]*n for _ in range(n)]
        for k in range(num_free):
            i, j = free_arcs[k]
            val = (bits >> k) & 1
            if val:
                T[i][j] = 1
            else:
                T[j][i] = 1

            if k >= len(self_paired_arcs):
                li, lj = pairs_list[k - len(self_paired_arcs)][1]
                # T(alpha(j), alpha(i)) = T(i, j) = val
                # So T(lj_small, li_big) or T(li, lj) depending on ordering
                # Actually, linked = (min(aj,ai), max(aj,ai))
                # We need T(aj, ai) = T(i, j) = val
                ai_orig, aj_orig = alpha[free_arcs[k][0]], alpha[free_arcs[k][1]]
                if val:
                    T[aj_orig][ai_orig] = 1
                else:
                    T[ai_orig][aj_orig] = 1

        tournaments.append(T)

    return tournaments, num_free


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
    """Fast induced C5 check."""
    for v0 in range(m):
        for v1 in adj_sets[v0]:
            if v1 <= v0:
                continue
            for v2 in adj_sets[v1]:
                if v2 <= v0 or v2 == v0 or v2 in adj_sets[v0]:
                    continue
                for v3 in adj_sets[v2]:
                    if v3 <= v0 or v3 == v0 or v3 == v1 or v3 in adj_sets[v0] or v3 in adj_sets[v1]:
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


def find_5cycles_through(T, n, v):
    """Find all directed 5-cycles through vertex v."""
    cycles = []
    others = [u for u in range(n) if u != v]
    for quad in combinations(others, 4):
        verts = (v,) + quad
        for perm in permutations(verts):
            if perm[0] != v:
                continue
            if all(T[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                cycles.append(perm)
                break  # one orientation per vertex set
    return cycles


def main():
    n = 8
    print("=" * 70)
    print("FAST INVESTIGATION: BlackSelf(8)")
    print("=" * 70)

    # Try multiple involutive anti-automorphisms
    # Standard reversal: alpha = (0 7)(1 6)(2 5)(3 4)
    alphas = [
        [7, 6, 5, 4, 3, 2, 1, 0],  # reversal
        [1, 0, 3, 2, 5, 4, 7, 6],  # (01)(23)(45)(67)
        [7, 6, 5, 4, 3, 2, 1, 0],  # same as first, but we'll also try others
    ]
    # Actually just use the reversal since it's symmetric enough
    alpha = [7, 6, 5, 4, 3, 2, 1, 0]

    print(f"Anti-automorphism alpha: {alpha}")
    print(f"Building all self-converse tournaments with this alpha...")

    tournaments, num_free = build_self_converse_tournaments(n, alpha)
    print(f"Free parameters: {num_free}, total tournaments: {len(tournaments)}")

    # Phase 1: Filter to those with |Aut| > 1
    print(f"\nPhase 1: Filter by |Aut| > 1...")
    candidates = []
    for idx, T in enumerate(tournaments):
        if idx % 1000 == 0 and idx > 0:
            print(f"  Checked {idx}/{len(tournaments)}, found {len(candidates)} with |Aut|>1")
        if has_automorphism_beyond_identity(T, n):
            candidates.append((idx, T))

    print(f"Found {len(candidates)} tournaments with |Aut| > 1")

    # Phase 2: Compute H(T) and Fix(beta) for candidates
    print(f"\nPhase 2: Check H(T) and Fix(beta)...")
    blackself = []
    for idx, T in candidates:
        H = count_ham_paths(T, n)

        # Compute H(Q) where Q is the decisive quotient
        # alpha pairs: (0,7), (1,6), (2,5), (3,4)
        # Representatives: 0, 1, 2, 3
        reps = [0, 1, 2, 3]
        Q = [[0]*4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                if i != j:
                    Q[i][j] = T[reps[i]][reps[j]]
        HQ = count_ham_paths(Q, 4)

        if HQ % 2 == 0:
            continue  # Need Fix(beta) odd
        if HQ == 0:
            continue
        if H % HQ != 0:
            print(f"  WARNING: H={H} not divisible by HQ={HQ}")
            continue

        ratio = H // HQ
        if ratio % 2 == 0:
            scores = score_sequence(T, n)
            aut_size = count_automorphisms(T, n)
            blackself.append((T, H, aut_size, HQ, ratio, scores))
            print(f"  *** BLACKSELF CANDIDATE: H={H}, |Aut|={aut_size}, "
                  f"Fix(beta)={HQ}, ratio={ratio}, scores={scores}")

    print(f"\nTotal BlackSelf(8) candidates: {len(blackself)}")

    if not blackself:
        print("No candidates found with reversal alpha. Try (01)(23)(45)(67)...")
        # Try second alpha
        alpha2 = [1, 0, 3, 2, 5, 4, 7, 6]
        print(f"\nAnti-automorphism alpha2: {alpha2}")
        tournaments2, nf2 = build_self_converse_tournaments(n, alpha2)
        print(f"Free parameters: {nf2}, tournaments: {len(tournaments2)}")

        candidates2 = []
        for idx, T in enumerate(tournaments2):
            if idx % 1000 == 0 and idx > 0:
                print(f"  Checked {idx}/{len(tournaments2)}")
            if has_automorphism_beyond_identity(T, n):
                candidates2.append((idx, T))

        print(f"Found {len(candidates2)} with |Aut| > 1")

        for idx, T in candidates2:
            H = count_ham_paths(T, n)
            reps = [0, 2, 4, 6]
            Q = [[0]*4 for _ in range(4)]
            for i in range(4):
                for j in range(4):
                    if i != j:
                        Q[i][j] = T[reps[i]][reps[j]]
            HQ = count_ham_paths(Q, 4)

            if HQ % 2 == 0 or HQ == 0:
                continue
            if H % HQ != 0:
                continue
            ratio = H // HQ
            if ratio % 2 == 0:
                scores = score_sequence(T, n)
                aut_size = count_automorphisms(T, n)
                blackself.append((T, H, aut_size, HQ, ratio, scores))
                print(f"  *** BLACKSELF CANDIDATE: H={H}, |Aut|={aut_size}, "
                      f"Fix(beta)={HQ}, ratio={ratio}, scores={scores}")

        if not blackself:
            print("Still no candidates. The BlackSelf(8) class may use a different alpha.")
            return

    # Phase 3: Deep analysis
    # Deduplicate by isomorphism (use canonical form heuristic: score + H + |Aut|)
    seen_sigs = set()
    unique = []
    for entry in blackself:
        T, H, aut_size, HQ, ratio, scores = entry
        sig = (H, aut_size, HQ, scores)
        if sig not in seen_sigs:
            seen_sigs.add(sig)
            unique.append(entry)

    print(f"\nDistinct candidates (by signature): {len(unique)}")

    for ci, (T, H, aut_size, HQ, ratio, scores) in enumerate(unique):
        print(f"\n{'='*60}")
        print(f"CANDIDATE {ci}")
        print(f"  H(T) = {H}")
        print(f"  |Aut(T)| = {aut_size}")
        print(f"  |Fix(beta)| = H(Q) = {HQ}")
        print(f"  H(T)/Fix(beta) = {ratio} (EVEN)")
        print(f"  H(T)/|Aut| = {H//aut_size}")
        print(f"  Score sequence: {scores}")

        # Adjacency matrix
        print(f"\n  Adjacency matrix:")
        for row in T:
            print(f"    {''.join(str(x) for x in row)}")

        # Out-degrees
        out_deg = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
        print(f"  Out-degrees: {out_deg}")

        # 3-cycles
        cycles_3 = find_all_3cycles(T, n)
        print(f"  Number of 3-cycles: {len(cycles_3)}")

        # Build Omega and check for C5
        m = len(cycles_3)
        cycle_sets = [set(c) for c in cycles_3]
        adj = [[0]*m for _ in range(m)]
        for i in range(m):
            for j in range(i+1, m):
                if cycle_sets[i] & cycle_sets[j]:
                    adj[i][j] = adj[j][i] = 1
        adj_sets = [set(j for j in range(m) if adj[i][j]) for i in range(m)]

        has_c5, c5_wit = has_induced_C5(adj_sets, m)
        print(f"  Omega(T) has induced C5: {has_c5}")
        if has_c5:
            print(f"    C5 cycles: {[cycles_3[v] for v in c5_wit]}")

        # Independence polynomial
        if m <= 20:
            coeffs = independence_poly(adj, m)
            I2 = sum(c * (2**k) for k, c in enumerate(coeffs))
            print(f"  I(Omega, 2) = {I2}")
            print(f"  OCF holds: {I2 == H}")
            print(f"  Coeffs: {coeffs}")

        # mu values for 3-cycles
        print(f"\n  mu values for 3-cycles:")
        mu_dist = Counter()
        mu_examples = {}
        for c in cycles_3:
            mu = mu_value(T, n, set(c))
            mu_dist[mu] += 1
            if mu > 1 and mu not in mu_examples:
                mu_examples[mu] = c
        print(f"    Distribution: {dict(mu_dist)}")
        for mu_val, ex in mu_examples.items():
            print(f"    Example with mu={mu_val}: {ex}")

        # 5-cycles through vertex 0
        print(f"\n  5-cycles through vertex 0:")
        c5_v0 = find_5cycles_through(T, n, 0)
        print(f"    Count: {len(c5_v0)}")
        if c5_v0:
            mu5_dist = Counter()
            for c in c5_v0:
                mu = mu_value(T, n, set(c))
                mu5_dist[mu] += 1
            print(f"    mu distribution: {dict(mu5_dist)}")

        # Check Paley connection
        print(f"\n  Regularity check:")
        is_regular = len(set(out_deg)) == 1
        print(f"    Regular: {is_regular}")
        if is_regular:
            # Doubly regular check
            common_out = {}
            for i in range(n):
                for j in range(i+1, n):
                    ct = sum(1 for k in range(n) if k != i and k != j and T[i][k] and T[j][k])
                    key = "adj" if T[i][j] else "non-adj"
                    common_out.setdefault(key, set()).add(ct)
            print(f"    Common out-neighbors: {common_out}")
            doubly_reg = all(len(v) == 1 for v in common_out.values())
            print(f"    Doubly regular: {doubly_reg}")

        # Signed position identity check (sample a few pairs)
        print(f"\n  Signed position identity (sampling 3 pairs):")
        paths = enumerate_ham_paths(T, n)
        print(f"    Total Ham paths: {len(paths)}")
        checked = 0
        for i in range(n):
            for j in range(i+1, n):
                if T[i][j] and checked < 3:
                    lhs = sum((-1)**P.index(i) for P in paths if P.index(i) < P.index(j))
                    rhs = sum((-1)**P.index(j) for P in paths if P.index(j) < P.index(i))
                    print(f"    ({i},{j}): LHS={lhs}, RHS={rhs}, match={lhs==rhs}")
                    checked += 1

        # OCF per-vertex analysis
        print(f"\n  OCF per-vertex (D_v = sum_C mu(C) for 3-cycles through v):")
        for v in range(min(4, n)):
            v_cycles = [c for c in cycles_3 if v in c]
            D_v = sum(mu_value(T, n, set(c)) for c in v_cycles)
            # S_v = Type I count, R_v = Type II count
            # For now just report D_v
            print(f"    v={v}: {len(v_cycles)} 3-cycles, D_v = {D_v}")


if __name__ == "__main__":
    main()
