#!/usr/bin/env python3
"""
Check self-converse status of ALL 22 non-circulant vertex-transitive tournaments
at n=21, using McKay's digraph6 data from users.cecs.anu.edu.au/~bdm/data/.

Also checks the 88 circulant ones for comparison.

kind-pasteur-2026-03-06-S25e
"""

import sys
from itertools import permutations

def decode_digraph6(s):
    """Decode a digraph6 string to adjacency matrix."""
    s = s.strip()
    if s.startswith('&'):
        s = s[1:]

    data = [ord(c) - 63 for c in s]
    if data[0] <= 62:
        n = data[0]
        idx = 1
    else:
        n = (data[1] << 12) | (data[2] << 6) | data[3]
        idx = 4

    bits = []
    for d in data[idx:]:
        for k in range(5, -1, -1):
            bits.append((d >> k) & 1)

    adj = [[0]*n for _ in range(n)]
    bit_idx = 0
    for i in range(n):
        for j in range(n):
            if i == j:
                bit_idx += 1  # digraph6 includes diagonal entries
                continue
            adj[i][j] = bits[bit_idx]
            bit_idx += 1
    return n, adj

def canonical_invariant(adj, n):
    """Compute a detailed invariant for isomorphism testing."""
    # Out-degrees
    out_deg = tuple(sorted(sum(adj[i]) for i in range(n)))

    # 3-cycle count
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                s = adj[i][j] + adj[j][i] + adj[i][k] + adj[k][i] + adj[j][k] + adj[k][j]
                # 3-cycle if all edges go same direction around triangle
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1

    return (out_deg, c3)

def is_anti_automorphism(adj, n, sigma):
    """Check if sigma is an anti-automorphism: adj[sigma[i]][sigma[j]] = adj[j][i] for all i,j."""
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if adj[sigma[i]][sigma[j]] != adj[j][i]:
                return False
    return True

def check_self_converse_backtrack(adj, n, max_attempts=10000000):
    """
    Check if tournament is self-converse by searching for an anti-automorphism.
    Uses backtracking with pruning.
    """
    # Build reverse adjacency
    rev = [[adj[j][i] for j in range(n)] for i in range(n)]

    # For each vertex v, compute its "out-neighborhood profile"
    out_sets = [frozenset(j for j in range(n) if adj[i][j]) for i in range(n)]
    in_sets = [frozenset(j for j in range(n) if adj[j][i]) for i in range(n)]

    # Rev profiles
    rev_out_sets = [frozenset(j for j in range(n) if rev[i][j]) for i in range(n)]

    # Group vertices by (out_degree, in_degree, 3-cycle count)
    def vertex_profile(a, mat):
        od = sum(mat[a])
        # Count 3-cycles through a
        c3 = 0
        for j in range(n):
            if j == a:
                continue
            for k in range(j+1, n):
                if k == a:
                    continue
                if mat[a][j] and mat[j][k] and mat[k][a]:
                    c3 += 1
                if mat[a][k] and mat[k][j] and mat[j][a]:
                    c3 += 1
        return (od, c3)

    profiles_adj = [vertex_profile(v, adj) for v in range(n)]
    profiles_rev = [vertex_profile(v, rev) for v in range(n)]

    # sigma maps adj -> rev, so profile of sigma(v) in rev must match profile of v in adj
    # But rev profile = (in_degree in adj, same 3-cycles)
    # Actually for anti-aut: adj[sigma(i)][sigma(j)] = rev[i][j] = adj[j][i]
    # So sigma is an isomorphism from adj to rev
    # Profile of v in adj must match profile of sigma(v) in rev

    candidates = {}
    for v in range(n):
        p = profiles_adj[v]
        candidates[v] = [w for w in range(n) if profiles_rev[w] == p]

    # Check if any vertex has no candidates
    for v in range(n):
        if not candidates[v]:
            return False, None

    # Backtracking search
    sigma = [None] * n
    used = [False] * n
    attempts = [0]

    def backtrack(pos):
        if attempts[0] >= max_attempts:
            return None  # inconclusive
        if pos == n:
            return True

        for w in candidates[pos]:
            if used[w]:
                continue

            # Check consistency: for all already-assigned i < pos,
            # adj[sigma[i]][w] must equal adj[pos][i] (= rev[i][pos] = adj[pos][i])
            ok = True
            for i in range(pos):
                if sigma[i] is None:
                    continue
                # Need adj[sigma[i]][w] == adj[pos][i]
                # and adj[w][sigma[i]] == adj[i][pos]
                if adj[sigma[i]][w] != adj[pos][i]:
                    ok = False
                    break
                if adj[w][sigma[i]] != adj[i][pos]:
                    ok = False
                    break

            if ok:
                sigma[pos] = w
                used[w] = True
                attempts[0] += 1
                result = backtrack(pos + 1)
                if result is True:
                    return True
                if result is None:
                    return None
                sigma[pos] = None
                used[w] = False

        return False

    result = backtrack(0)
    if result is True:
        return True, list(sigma)
    elif result is None:
        return None, None  # inconclusive
    else:
        return False, None

def main():
    # Read non-circulant VT tournaments
    with open('othervttourn21.d6', 'r') as f:
        nc_lines = [l.strip() for l in f if l.strip()]

    print(f"=== Non-circulant VT tournaments at n=21: {len(nc_lines)} ===")

    nc_sc_count = 0
    nc_nsc_count = 0
    nc_unknown = 0

    for idx, line in enumerate(nc_lines):
        n, adj = decode_digraph6(line)
        out_degs = sorted(sum(adj[i]) for i in range(n))
        regular = all(d == 10 for d in out_degs)

        is_sc, sigma = check_self_converse_backtrack(adj, n)

        if is_sc is True:
            nc_sc_count += 1
            status = "SC"
        elif is_sc is False:
            nc_nsc_count += 1
            status = "NOT SC"
        else:
            nc_unknown += 1
            status = "UNKNOWN (timeout)"

        reg_str = "regular" if regular else f"degs={out_degs[0]}-{out_degs[-1]}"
        print(f"  T{idx+1:2d}: {status:12s} | {reg_str}")
        sys.stdout.flush()

    print(f"\nNon-circulant summary: {nc_sc_count} SC, {nc_nsc_count} NOT SC, {nc_unknown} unknown")

    # Read circulant tournaments
    with open('circtourn21.d6', 'r') as f:
        circ_lines = [l.strip() for l in f if l.strip()]

    print(f"\n=== Circulant tournaments at n=21: {len(circ_lines)} ===")

    circ_sc_count = 0
    circ_nsc_count = 0
    circ_unknown = 0

    for idx, line in enumerate(circ_lines):
        n, adj = decode_digraph6(line)

        # For circulants, self-converse via sigma: i -> -i mod n
        sigma_neg = [(n - i) % n for i in range(n)]
        if is_anti_automorphism(adj, n, sigma_neg):
            circ_sc_count += 1
        else:
            # Try backtracking
            is_sc, sigma = check_self_converse_backtrack(adj, n, max_attempts=1000000)
            if is_sc is True:
                circ_sc_count += 1
            elif is_sc is False:
                circ_nsc_count += 1
                print(f"  Circ T{idx+1}: NOT SC!")
            else:
                circ_unknown += 1

    print(f"Circulant summary: {circ_sc_count} SC, {circ_nsc_count} NOT SC, {circ_unknown} unknown")

    print(f"\n=== OVERALL at n=21 ===")
    print(f"Total VT: {len(nc_lines) + len(circ_lines)}")
    print(f"Total SC: {nc_sc_count + circ_sc_count}")
    print(f"Total NOT SC: {nc_nsc_count + circ_nsc_count}")
    print(f"Total unknown: {nc_unknown + circ_unknown}")

    if nc_nsc_count > 0:
        print(f"\nCONFIRMED: {nc_nsc_count} non-circulant VT tournaments at n=21 are NOT self-converse!")
        print("This validates our F_21 counterexample and shows the phenomenon is not isolated.")

    print("\nDONE")

if __name__ == '__main__':
    main()
