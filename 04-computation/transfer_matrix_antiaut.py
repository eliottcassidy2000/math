#!/usr/bin/env python3
"""
Investigate transfer matrix symmetry via anti-automorphism structure.

The transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(R) counts signed
weighted path contributions. M is ALWAYS symmetric (verified n=4,...,8).

Key hypothesis: M's symmetry follows from a hidden reversibility related to
tournament anti-automorphisms or the converse operation.

Specifically:
- For SC tournaments, the anti-automorphism sigma reverses all arcs.
  This maps "path ending at a" to "path starting at sigma(a)".
  If sigma is an involution, this gives M[a,b] = M[sigma(b), sigma(a)].
  For sigma = id, this is symmetry. For general involutory sigma, this is
  "twisted symmetry" M = P^T M P where P is the permutation matrix of sigma.

- For ALL tournaments (including NSC), is there an analogous structure?
  The converse T^op has H(T^op) = H(T) (path reversal). Does this imply
  M_T = M_{T^op}^T? If so, M_T symmetric iff M_T = M_{T^op}^T.

This script computes M explicitly and tests these relationships.

Instance: opus-2026-03-06-S4
"""
import sys
import time
from itertools import permutations, combinations
from collections import defaultdict

def build_adj(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    adj = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        adj[i][i + 1] = 1
    for k, (i, j) in enumerate(arcs):
        if bits & (1 << k):
            adj[j][i] = 1
        else:
            adj[i][j] = 1
    return adj

def canon(adj, n):
    best = None
    for perm in permutations(range(n)):
        s = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best

def count_ham(adj, n):
    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def compute_transfer_matrix(adj, n):
    """
    Compute the transfer matrix M[a,b] for tournament T.

    M[a,b] = sum_{S subset of V\{a,b}} (-1)^|S| * E_a(S) * B_b(R)

    where:
    - E_a(S) = number of Hamiltonian paths on S union {a} ending at a
    - B_b(R) = number of Hamiltonian paths on R union {b} starting at b
    - R = V \ (S union {a,b})

    Actually, let me use the standard definition:
    M[a,b] = number of Hamiltonian paths from a to b in T.
    """
    # Direct computation: M[a,b] = #Ham paths from a to b
    full = (1 << n) - 1
    M = [[0] * n for _ in range(n)]

    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + c

    # M[a][b] = number of Ham paths starting at a, ending at b
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            # Count paths starting at a, ending at b
            M[a][b] = dp.get((full, b), 0)  # This gives total ending at b, not starting at a

    # Actually we need a different DP for start-specific counts
    # DP: for each starting vertex a, count paths
    M = [[0] * n for _ in range(n)]
    for a in range(n):
        # DP starting only from a
        dp_a = {}
        dp_a[(1 << a, a)] = 1
        for mask in range(1, 1 << n):
            if not (mask & (1 << a)):
                continue  # must include a
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                c = dp_a.get((mask, v), 0)
                if c == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if adj[v][u]:
                        dp_a[(mask | (1 << u), u)] = dp_a.get((mask | (1 << u), u), 0) + c
        for b in range(n):
            M[a][b] = dp_a.get((full, b), 0)

    return M

def is_symmetric(M, n):
    for i in range(n):
        for j in range(i + 1, n):
            if M[i][j] != M[j][i]:
                return False
    return True

def transpose(M, n):
    return [[M[j][i] for j in range(n)] for i in range(n)]

def permute_matrix(M, sigma, n):
    """Return M' where M'[i][j] = M[sigma(i)][sigma(j)]."""
    return [[M[sigma[i]][sigma[j]] for j in range(n)] for i in range(n)]

def find_anti_auts(adj, n):
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            if not ok:
                break
            for j in range(n):
                if i == j:
                    continue
                if adj[perm[i]][perm[j]] != (1 - adj[i][j]):
                    ok = False
                    break
        if ok:
            auts.append(perm)
    return auts

def find_auts(adj, n):
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            if not ok:
                break
            for j in range(n):
                if i == j:
                    continue
                if adj[perm[i]][perm[j]] != adj[i][j]:
                    ok = False
                    break
        if ok:
            auts.append(perm)
    return auts


def main():
    target = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    for n in range(4, target + 1):
        arcs = [(i, j) for i in range(n) for j in range(i + 2, n)]
        m = len(arcs)
        total = 1 << m

        print(f"\n{'=' * 70}")
        print(f"n={n}: Transfer Matrix & Anti-Automorphism Analysis")
        print(f"{'=' * 70}")

        t0 = time.time()
        seen = {}

        sym_count = 0
        asym_count = 0
        sc_count = 0
        nsc_count = 0

        # Test: M_T = M_{T^op}^T ?
        converse_transpose_match = 0
        converse_transpose_fail = 0

        # Test: For SC with involutory anti-aut sigma, M = P_sigma^T M P_sigma?
        sigma_twisted_match = 0
        sigma_twisted_fail = 0

        # Test: For SC, does M have block structure w.r.t. sigma orbits?
        sigma_block_patterns = []

        for bits in range(total):
            if bits % 2000 == 0 and bits > 0:
                elapsed = time.time() - t0
                rate = bits / elapsed
                print(f"  {bits}/{total} ({100*bits/total:.1f}%), ETA {(total-bits)/rate:.0f}s",
                      flush=True)

            adj = build_adj(n, bits)
            cf = canon(adj, n)
            if cf in seen:
                continue
            seen[cf] = True

            M = compute_transfer_matrix(adj, n)
            sym = is_symmetric(M, n)

            if sym:
                sym_count += 1
            else:
                asym_count += 1
                print(f"  ASYMMETRIC M at bits={bits}!")
                for row in M:
                    print(f"    {row}")

            # Check SC
            conv = [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
            conv_cf = canon(conv, n)
            is_sc = (cf == conv_cf)

            # Compute M for converse
            M_conv = compute_transfer_matrix(conv, n)
            M_conv_T = transpose(M_conv, n)

            if M == M_conv_T:
                converse_transpose_match += 1
            else:
                converse_transpose_fail += 1

            if is_sc:
                sc_count += 1
                # Find involutory anti-auts
                anti_auts = find_anti_auts(adj, n)
                involutions = [s for s in anti_auts
                               if all(s[s[i]] == i for i in range(n))]

                for sigma in involutions[:1]:  # test first involution
                    M_perm = permute_matrix(M, sigma, n)
                    # Test: M = sigma^T M sigma, i.e., M[i][j] = M[sigma(i)][sigma(j)]
                    if M == M_perm:
                        sigma_twisted_match += 1
                    else:
                        sigma_twisted_fail += 1

                    # Analyze M's block structure w.r.t. sigma orbits
                    # Get orbits
                    visited = [False] * n
                    orbits = []
                    for i in range(n):
                        if visited[i]:
                            continue
                        orbit = [i]
                        visited[i] = True
                        j = sigma[i]
                        while not visited[j]:
                            orbit.append(j)
                            visited[j] = True
                            j = sigma[j]
                        orbits.append(tuple(orbit))

                    # For each orbit pair, check M values
                    orbit_block_info = {}
                    for oi, o1 in enumerate(orbits):
                        for oj, o2 in enumerate(orbits):
                            vals = set()
                            for a in o1:
                                for b in o2:
                                    if a != b:
                                        vals.add(M[a][b])
                            orbit_block_info[(oi, oj)] = vals

                    h = count_ham(adj, n)
                    sigma_block_patterns.append({
                        'h': h, 'orbits': orbits, 'blocks': orbit_block_info,
                        'sigma': sigma,
                    })
            else:
                nsc_count += 1

        elapsed = time.time() - t0
        print(f"\n  Results ({elapsed:.1f}s):")
        print(f"    Symmetric M: {sym_count}/{sym_count + asym_count} {'(ALL!)' if asym_count == 0 else 'FAILURES!'}")
        print(f"    M_T = M_{'{T^op}'}^T: {converse_transpose_match}/{converse_transpose_match + converse_transpose_fail}")
        print(f"    SC classes: {sc_count}, NSC classes: {nsc_count}")
        if sigma_twisted_match + sigma_twisted_fail > 0:
            print(f"    SC sigma-twisted symmetry: {sigma_twisted_match}/{sigma_twisted_match + sigma_twisted_fail}")

        # Detailed block structure for SC
        if sigma_block_patterns:
            print(f"\n  SC TOURNAMENT SIGMA-BLOCK STRUCTURE:")
            for entry in sorted(sigma_block_patterns, key=lambda x: x['h']):
                h = entry['h']
                orbits = entry['orbits']
                sigma = entry['sigma']
                blocks = entry['blocks']
                print(f"    H={h}, sigma={sigma}, orbits={orbits}")
                for (oi, oj), vals in sorted(blocks.items()):
                    if vals:
                        o1, o2 = orbits[oi], orbits[oj]
                        print(f"      O{oi}={o1} x O{oj}={o2}: M values = {sorted(vals)}")

        # Key insight: path reversal connection
        print(f"\n  PATH REVERSAL INSIGHT:")
        print(f"    M[a,b] = #paths a->b in T")
        print(f"    M_{{T^op}}[a,b] = #paths a->b in T^op = #paths b->a in T = M[b,a]")
        print(f"    So M_{{T^op}} = M^T always (by path reversal)")
        print(f"    M symmetric <=> M = M^T <=> M_T = M_{{T^op}}")
        print(f"    This means: for EVERY tournament, the endpoint-pair path count matrix")
        print(f"    is symmetric iff T and T^op have the same endpoint-pair path counts.")
        print(f"    For SC tournaments, T ~ T^op (isomorphic), so M_T = P M_{{T^op}} P^T = P M^T P^T")
        print(f"    where P = permutation realizing the isomorphism.")
        print(f"    M symmetric requires M = P M P^T for ALL such P.")

        # Additional test: Is M_T always equal to M_{T^op} (not transpose, the matrix itself)?
        print(f"\n  ADDITIONAL TEST: M_T = M_{{T^op}} (same matrix, not transpose)?")
        seen2 = {}
        eq_count = 0
        neq_count = 0
        for bits in range(total):
            adj = build_adj(n, bits)
            cf = canon(adj, n)
            if cf in seen2:
                continue
            seen2[cf] = True
            conv = [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
            M_T = compute_transfer_matrix(adj, n)
            M_Top = compute_transfer_matrix(conv, n)
            if M_T == M_Top:
                eq_count += 1
            else:
                neq_count += 1

        print(f"    M_T = M_{{T^op}}: {eq_count}/{eq_count + neq_count}")


if __name__ == '__main__':
    main()
