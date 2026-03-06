#!/usr/bin/env python3
"""
Deep investigation of the n=8 anomaly: BlackSelf(8).

Find the unique tournament T* at n=8 that is:
  - Self-converse (T ~ T^op)
  - |Aut(T)| > 1
  - |Fix(beta)| odd (= H(Q) where Q is decisive quotient)
  - H(T)/|Fix(beta)| is even

Then examine:
  1. Is T* the Paley tournament on GF(7) + one vertex? Or Hadamard-related?
  2. What does Omega(T*) look like? Does it have a C5?
  3. What are the mu values for 3-cycles and 5-cycles of T*?
  4. How does the pos/signed-position identity behave on T*?

Instance: opus-2026-03-05-S8
"""

import sys
import random
from itertools import permutations, combinations
from collections import Counter

sys.path.insert(0, '03-artifacts/code')


def make_adj(T, n):
    """Convert tournament matrix to adjacency dict."""
    return {(i, j): T[i][j] for i in range(n) for j in range(n) if i != j}


def tournament_from_bits(bits, n):
    """Create tournament from upper-triangle bits."""
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]:
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T


def transpose_tournament(T, n):
    """T^op: reverse all arcs."""
    Top = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            Top[i][j] = T[j][i]
    return Top


def apply_perm(T, perm, n):
    """Apply vertex permutation to tournament."""
    T2 = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            T2[perm[i]][perm[j]] = T[i][j]
    return T2


def tournaments_equal(T1, T2, n):
    for i in range(n):
        for j in range(n):
            if T1[i][j] != T2[i][j]:
                return False
    return True


def is_isomorphic(T1, T2, n):
    """Check if T1 ~ T2. Returns witnessing permutation or None."""
    # For small n, brute force
    for perm in permutations(range(n)):
        if tournaments_equal(apply_perm(T1, perm, n), T2, n):
            return perm
    return None


def find_anti_automorphisms(T, n):
    """Find all anti-automorphisms: perm alpha s.t. T(alpha(u),alpha(v)) = T(v,u)."""
    Top = transpose_tournament(T, n)
    auts = []
    for perm in permutations(range(n)):
        if tournaments_equal(apply_perm(T, perm, n), Top, n):
            auts.append(perm)
    return auts


def find_automorphisms(T, n):
    """Find all automorphisms."""
    auts = []
    for perm in permutations(range(n)):
        if tournaments_equal(apply_perm(T, perm, n), T, n):
            auts.append(perm)
    return auts


def is_involution(perm):
    """Check if permutation is an involution (order divides 2)."""
    return all(perm[perm[i]] == i for i in range(len(perm)))


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


def decisive_quotient_H(T, n, alpha):
    """Compute H(Q) where Q is the decisive quotient tournament.

    alpha is an involutive anti-automorphism.
    Pair up vertices: {v, alpha(v)}.
    For n even: all pairs are 2-element.
    Q has n/2 vertices; Q(pair_i, pair_j) = T(v_i, v_j) where v_i, v_j are reps.
    """
    # Find orbits of alpha
    orbits = []
    seen = set()
    fixed = []
    for v in range(n):
        if v in seen:
            continue
        if alpha[v] == v:
            fixed.append(v)
            seen.add(v)
        else:
            orbits.append((v, alpha[v]))
            seen.add(v)
            seen.add(alpha[v])

    # For even n, should have n/2 orbits of size 2 (no fixed points)
    # For odd n, one fixed point + (n-1)/2 orbits

    m = len(orbits) + len(fixed)
    # Build quotient tournament Q on orbit representatives
    reps = [orb[0] for orb in orbits] + fixed
    Q = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j:
                Q[i][j] = T[reps[i]][reps[j]]

    return count_ham_paths(Q, m), orbits, fixed


def count_concordant_paths(T, n, alpha):
    """Count paths P such that beta(P) = P.

    beta acts on paths: if P = (v1,...,vn), then beta(P) = (alpha(vn),...,alpha(v1)).
    P is concordant (fixed by beta) iff v_i = alpha(v_{n+1-i}) for all i.
    """
    paths = enumerate_ham_paths(T, n)
    concordant = 0
    for P in paths:
        is_conc = True
        for i in range(n):
            if P[i] != alpha[P[n-1-i]]:
                is_conc = False
                break
        if is_conc:
            concordant += 1
    return concordant


def find_all_3cycles(T, n):
    """Find all directed 3-cycles."""
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if T[a][b] and T[b][c] and T[c][a]:
            cycles.append((a, b, c))
        elif T[a][c] and T[c][b] and T[b][a]:
            cycles.append((a, c, b))
    return cycles


def find_all_5cycles(T, n):
    """Find all directed 5-cycles."""
    cycles = []
    for verts in combinations(range(n), 5):
        for perm in permutations(verts):
            if all(T[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                # Canonical form: start with smallest vertex
                min_idx = perm.index(min(perm))
                canon = tuple(perm[min_idx:] + perm[:min_idx])
                if canon not in cycles:
                    cycles.append(canon)
                break  # Only need one orientation per vertex set if it exists
    # Deduplicate
    return list(set(cycles))


def build_omega(T, n):
    """Build the odd-cycle conflict graph Omega(T).

    Vertices of Omega = directed 3-cycles of T.
    Two 3-cycles are adjacent in Omega iff they share a vertex.
    """
    cycles_3 = find_all_3cycles(T, n)
    m = len(cycles_3)
    cycle_sets = [set(c) for c in cycles_3]

    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = 1

    return cycles_3, adj


def has_induced_C5(adj, m):
    """Check if graph has an induced 5-cycle."""
    adj_sets = [set(j for j in range(m) if adj[i][j]) for i in range(m)]
    for verts in combinations(range(m), 5):
        # Check if these 5 vertices form an induced C5
        v = list(verts)
        for perm in permutations(v):
            edges_ok = True
            for i in range(5):
                # Adjacent in cycle
                if perm[(i+1)%5] not in adj_sets[perm[i]]:
                    edges_ok = False
                    break
                # Non-adjacent (skip)
                if perm[(i+2)%5] in adj_sets[perm[i]]:
                    edges_ok = False
                    break
            if edges_ok:
                return True, perm
    return False, None


def independence_poly(adj, m):
    """Compute independence polynomial of a graph.

    I(G, x) = sum over independent sets S of x^|S|.
    Returns list of coefficients [alpha_0, alpha_1, ...].
    """
    # Use inclusion via bitmask for small m
    if m > 20:
        return None  # Too large

    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        # Check if mask is an independent set
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
    """Compute mu(C) = I(Omega(T)|_{avoid cycle_verts}, 2).

    The "avoid" subgraph: odd cycles on V \ cycle_verts.
    """
    remaining = [v for v in range(n) if v not in cycle_verts]
    if len(remaining) < 3:
        return 1  # No odd cycles possible

    # Find 3-cycles among remaining vertices
    sub_cycles = []
    for a, b, c in combinations(remaining, 3):
        if T[a][b] and T[b][c] and T[c][a]:
            sub_cycles.append((a, b, c))
        elif T[a][c] and T[c][b] and T[b][a]:
            sub_cycles.append((a, c, b))

    if not sub_cycles:
        return 1

    # Build conflict graph among these sub-cycles
    m = len(sub_cycles)
    cycle_sets = [set(c) for c in sub_cycles]
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = 1

    coeffs = independence_poly(adj, m)
    return sum(c * (2**k) for k, c in enumerate(coeffs))


def score_sequence(T, n):
    return tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)], reverse=True))


def signed_position_identity_check(T, n, i, j):
    """Check sum_{P: i->j} (-1)^{pos(i)} = sum_{P': j->i} (-1)^{pos(j)}.

    pos(v) = position of v in the path (0-indexed).
    """
    paths = enumerate_ham_paths(T, n)

    lhs = 0  # paths where i appears before j
    rhs = 0  # paths where j appears before i

    for P in paths:
        pi = P.index(i)
        pj = P.index(j)
        if pi < pj:
            lhs += (-1)**pi
        else:
            rhs += (-1)**pj

    return lhs, rhs, lhs == rhs


# =============================================================================
# MAIN INVESTIGATION
# =============================================================================

def investigate_blackself8():
    n = 8
    num_arcs = n * (n - 1) // 2  # 28

    print("=" * 70)
    print("INVESTIGATION: BlackSelf(8) — the n=8 exceptional class")
    print("=" * 70)
    print(f"n={n}, {num_arcs} arcs, 2^{num_arcs} = {2**num_arcs} tournaments")
    print()

    # Strategy: sample random self-converse tournaments at n=8
    # A self-converse tournament has T ~ T^op.
    # For n=8 (even), an involutive anti-automorphism alpha is a
    # fixed-point-free involution (4 two-cycles) such that
    # T(alpha(u), alpha(v)) = T(v, u).

    # More efficient: enumerate by fixing an involution alpha and
    # constructing tournaments compatible with it.

    # Fix alpha = (0 7)(1 6)(2 5)(3 4) — the "reversal" involution
    alpha = [7, 6, 5, 4, 3, 2, 1, 0]

    print(f"Fixed anti-automorphism alpha = {alpha}")
    print(f"Orbits: (0,7) (1,6) (2,5) (3,4)")
    print()

    # For T to have alpha as anti-automorphism:
    # T(alpha(u), alpha(v)) = T(v, u)
    # i.e., T(7-u, 7-v) = T(v, u)
    #
    # This constrains the 28 arc variables.
    # Free arcs: those (i,j) with i<j that are NOT determined by another arc.
    # T(i,j) determines T(7-j, 7-i).
    # So arcs (i,j) and (7-j, 7-i) are linked.
    # If (i,j) = (7-j, 7-i), then i = 7-j, j = 7-i => i+j = 7.
    # These "diagonal" arcs satisfy T(i, 7-i) = T(7-i, i) = 1 - T(i, 7-i),
    # which is impossible! So no arc is self-paired.
    # Wait: T(alpha(i), alpha(j)) = T(j, i) = 1 - T(i, j) (for i != j).
    # T(7-i, 7-j) = 1 - T(i, j).
    # For i < j, consider pair (i, j). The linked pair is (7-j, 7-i).
    # If 7-j < 7-i (i.e., i < j), the linked arc is (7-j, 7-i) with 7-j < 7-i.
    # Are (i,j) and (7-j, 7-i) the same? Only if i = 7-j and j = 7-i, i.e., i+j=7.
    # For such pairs: T(i, 7-i) and T(7-(7-i), 7-i) = T(i, 7-i).
    # Constraint: T(7-i, 7-(7-i)) = T(7-i, i) = 1 - T(i, 7-i).
    # But from anti-aut: T(alpha(i), alpha(7-i)) = T(7-i, i).
    # alpha(i) = 7-i, alpha(7-i) = i. So T(7-i, i) = T(7-i, i). Tautology!
    # Hmm, let me redo this.

    # Anti-automorphism condition: T(alpha(u), alpha(v)) = T(v, u) for all u != v.
    # With alpha(x) = 7-x:
    # T(7-u, 7-v) = T(v, u) for all u != v.
    #
    # For arc (i, j) with i < j:
    #   T(i, j) is free, and it determines T(7-j, 7-i) = T(j, i) = 1 - T(i, j).
    #   But (7-j, 7-i): since i < j, 7-j < 7-i, so this is really T(7-j, 7-i).
    #
    #   Also: T(7-i, 7-j) = T(j, i) = 1 - T(i, j).
    #   So T(7-j, 7-i) = T(i, j) and T(7-i, 7-j) = 1 - T(i, j).
    #   Since 7-j < 7-i, the arc (7-j, 7-i) has T(7-j, 7-i) = T(i, j).
    #
    # So arcs (i,j) and (7-j, 7-i) [with 7-j < 7-i since i < j] have the SAME value.
    # They're in the same orbit if (i,j) != (7-j, 7-i).
    # (i,j) = (7-j, 7-i) iff i = 7-j, j = 7-i, i.e., i+j = 7.

    # Count free arcs:
    pairs = []
    seen = set()
    self_paired = []
    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in seen:
                continue
            linked = (7-j, 7-i) if 7-j < 7-i else (7-i, 7-j)
            # Actually let me be more careful
            u, v = 7-j, 7-i
            if u > v:
                u, v = v, u
            linked = (u, v)
            if (i, j) == linked:
                self_paired.append((i, j))
                seen.add((i, j))
            else:
                pairs.append(((i, j), linked))
                seen.add((i, j))
                seen.add(linked)

    print(f"Self-paired arcs (i+j=7): {self_paired}")
    print(f"Paired arcs: {len(pairs)} pairs")
    print(f"Total free parameters: {len(self_paired) + len(pairs)}")
    print()

    # For self-paired arcs (i, 7-i): T(i, 7-i) = T(i, 7-i). No constraint from alpha.
    # For paired arcs: one bit determines both.
    # Total: len(self_paired) + len(pairs) free bits.

    num_free = len(self_paired) + len(pairs)
    print(f"Enumerating all 2^{num_free} = {2**num_free} self-converse tournaments with alpha=(0 7)(1 6)(2 5)(3 4)")
    print()

    # Enumerate
    free_arcs = self_paired + [p[0] for p in pairs]

    results = []

    for bits in range(2**num_free):
        T = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(free_arcs):
            val = (bits >> k) & 1
            if val:
                T[i][j] = 1
            else:
                T[j][i] = 1
            # Set linked arc
            if k < len(self_paired):
                pass  # Self-paired, no linked arc
            else:
                li, lj = pairs[k - len(self_paired)][1]
                if val:
                    T[li][lj] = 1
                else:
                    T[lj][li] = 1

        # Verify anti-automorphism
        ok = True
        for u in range(n):
            for v in range(n):
                if u != v:
                    if T[7-u][7-v] != T[v][u]:
                        ok = False
                        break
            if not ok:
                break
        if not ok:
            continue

        # Check |Aut(T)| > 1
        auts = find_automorphisms(T, n)
        if len(auts) <= 1:
            continue

        # Compute H(T)
        H = count_ham_paths(T, n)

        # Compute |Fix(beta)| = H(Q)
        HQ, orbits, fixed = decisive_quotient_H(T, n, alpha)

        # Check: |Fix(beta)| odd
        if HQ % 2 == 0:
            continue

        # Check: H(T)/|Fix(beta)| is even
        if HQ == 0:
            continue
        ratio = H // HQ
        if H % HQ != 0:
            # H should be divisible by |Fix(beta)| (by orbit counting)
            print(f"  WARNING: H={H} not divisible by Fix(beta)={HQ}")
            continue

        if ratio % 2 == 0:
            scores = score_sequence(T, n)
            print(f"*** BLACKSELF(8) CANDIDATE ***")
            print(f"  H(T) = {H}")
            print(f"  |Aut(T)| = {len(auts)}")
            print(f"  |Fix(beta)| = H(Q) = {HQ}")
            print(f"  H(T)/Fix(beta) = {ratio} (EVEN)")
            print(f"  Score sequence: {scores}")
            results.append((T, H, len(auts), HQ, ratio, scores, bits))

    print(f"\nTotal BlackSelf(8) candidates found: {len(results)}")

    if not results:
        print("No candidates found with this alpha. Try other involutions.")
        return

    # Group by isomorphism class
    classes = []
    for entry in results:
        T = entry[0]
        found = False
        for cls in classes:
            if is_isomorphic(T, cls[0][0], n) is not None:
                cls.append(entry)
                found = True
                break
        if not found:
            classes.append([entry])

    print(f"Distinct isomorphism classes: {len(classes)}")

    for ci, cls in enumerate(classes):
        T, H, aut_size, HQ, ratio, scores, bits = cls[0]
        print(f"\n{'='*60}")
        print(f"CLASS {ci}: {len(cls)} representatives")
        print(f"  H(T) = {H}, |Aut| = {aut_size}, Fix(beta) = {HQ}")
        print(f"  H/Fix(beta) = {ratio}, Score = {scores}")

        # Deep analysis of this class
        print(f"\n  --- Omega(T) analysis ---")
        cycles_3, omega_adj = build_omega(T, n)
        print(f"  Number of 3-cycles: {len(cycles_3)}")

        if len(cycles_3) <= 20:
            # Check for C5 in Omega
            has_c5, c5_witness = has_induced_C5(omega_adj, len(cycles_3))
            print(f"  Induced C5 in Omega(T): {has_c5}")
            if has_c5:
                print(f"    C5 witness: {[cycles_3[v] for v in c5_witness[:5]]}")

            # Independence polynomial
            coeffs = independence_poly(omega_adj, len(cycles_3))
            I_at_2 = sum(c * (2**k) for k, c in enumerate(coeffs))
            print(f"  I(Omega, 2) = {I_at_2}")
            print(f"  H(T) = {H}")
            print(f"  OCF holds: {I_at_2 == H}")
            print(f"  Independence poly coeffs: {coeffs}")

        # mu values for 3-cycles
        print(f"\n  --- mu values for 3-cycles ---")
        mu_vals = Counter()
        for c in cycles_3:
            mu = mu_value(T, n, set(c))
            mu_vals[mu] += 1
        print(f"  mu value distribution: {dict(mu_vals)}")

        # 5-cycles
        print(f"\n  --- 5-cycles ---")
        cycles_5 = find_all_5cycles(T, n)
        print(f"  Number of 5-cycles: {len(cycles_5)}")
        if cycles_5:
            mu5_vals = Counter()
            for c in cycles_5:
                mu = mu_value(T, n, set(c))
                mu5_vals[mu] += 1
            print(f"  mu values for 5-cycles: {dict(mu5_vals)}")

        # Signed position identity
        print(f"\n  --- Signed position identity ---")
        fail_count = 0
        for i in range(min(n, 4)):
            for j in range(i+1, min(n, 5)):
                if T[i][j]:
                    lhs, rhs, ok = signed_position_identity_check(T, n, i, j)
                    if not ok:
                        fail_count += 1
                        print(f"    FAIL: i={i}, j={j}: LHS={lhs}, RHS={rhs}")
        if fail_count == 0:
            print(f"    All tested pairs PASS")

        # Print adjacency matrix
        print(f"\n  --- Adjacency matrix ---")
        for row in T:
            print(f"    {''.join(str(x) for x in row)}")

        # Check if doubly regular
        print(f"\n  --- Regularity check ---")
        out_degrees = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
        print(f"  Out-degrees: {out_degrees}")
        is_regular = len(set(out_degrees)) == 1
        print(f"  Regular: {is_regular}")

        if is_regular:
            # Check doubly regular: for each pair (i,j), count common dominatees
            common_counts = set()
            for i in range(n):
                for j in range(i+1, n):
                    common = sum(1 for k in range(n) if k != i and k != j and T[i][k] and T[j][k])
                    common_counts.add(common)
            print(f"  Common out-neighbor counts: {common_counts}")
            print(f"  Doubly regular: {len(common_counts) <= 2}")  # At most 2 values (adj vs non-adj)

        # Check: is this the Paley tournament P(7)?
        # Paley tournament on GF(7): i -> j iff j - i is a QR mod 7
        # QRs mod 7: {1, 2, 4}
        print(f"\n  --- Paley tournament P(7) comparison ---")
        P7 = [[0]*7 for _ in range(7)]
        qr = {1, 2, 4}
        for i in range(7):
            for j in range(7):
                if i != j and (j - i) % 7 in qr:
                    P7[i][j] = 1
        H_P7 = count_ham_paths(P7, 7)
        print(f"  H(P7) = {H_P7}")
        print(f"  P7 score sequence: {score_sequence(P7, 7)}")

        # Check if T minus any vertex is isomorphic to P7
        for v in range(n):
            subT = [[T[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
            H_sub = count_ham_paths(subT, 7)
            sub_score = score_sequence(subT, 7)
            # Quick check: P7 is regular with degree 3
            if sub_score == score_sequence(P7, 7):
                print(f"  T - vertex {v}: score matches P7! H={H_sub}")


if __name__ == "__main__":
    investigate_blackself8()
