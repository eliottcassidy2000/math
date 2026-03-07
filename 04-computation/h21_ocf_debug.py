#!/usr/bin/env python3
"""
Debug: why is H != I(Omega, 2) for some n=8 tournaments?
Must be undercounting cycles.

At n=8 with 10 three-cycles, no 5-cycles detected, H=81 but I=21.
H=81 via OCF means I(Omega_TRUE, 2) = 81.
So there are MANY more cycles than 10.

Check: does my 5-cycle finder work correctly at n=8?

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations, permutations
from collections import Counter
import random

def find_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

def find_5cycles_bruteforce(adj, n):
    """Find all VERTEX SETS that support at least one directed 5-cycle."""
    cycle_sets = set()
    for verts in combinations(range(n), 5):
        v = list(verts)
        # Check all 5!/5 = 24 cyclic permutations
        found = False
        for perm in permutations(range(5)):
            ok = True
            for idx in range(5):
                if not adj[v[perm[idx]]][v[perm[(idx+1)%5]]]:
                    ok = False
                    break
            if ok:
                found = True
                break
        if found:
            cycle_sets.add(frozenset(verts))
    return list(cycle_sets)

def find_5cycles_count(adj, n):
    """Count DIRECTED 5-cycles per vertex set."""
    result = []
    for verts in combinations(range(n), 5):
        v = list(verts)
        count = 0
        for perm in permutations(range(5)):
            ok = True
            for idx in range(5):
                if not adj[v[perm[idx]]][v[perm[(idx+1)%5]]]:
                    ok = False
                    break
            if ok:
                count += 1
        num_cycles = count // 5  # divide by cyclic rotations
        if num_cycles > 0:
            result.append((frozenset(verts), num_cycles))
    return result

def find_7cycles_count(adj, n):
    """Count DIRECTED 7-cycles per vertex set."""
    result = []
    for verts in combinations(range(n), 7):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        for S in range(1, 128):
            for i in range(7):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(7):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, 7):
            if (127, j) in dp and adj[v[j]][v[0]]:
                count += dp[(127, j)]
        num = count // 7
        if num > 0:
            result.append((frozenset(verts), num))
    return result

def independence_polynomial_from_sets(cycle_vertex_sets, x):
    """Compute I(Omega, x) where Omega vertices are VERTEX SETS of cycles.
    Two cycle vertex sets are adjacent iff they share a vertex.
    Multiple directed cycles on same vertex set are one vertex in Omega.

    WAIT: Actually in Omega(T), each DIRECTED cycle is a separate vertex.
    But directed cycles on the same vertex set are always adjacent (share all vertices).
    So grouping by vertex set doesn't affect alpha_k.
    Each vertex set with d directed cycles contributes d vertices to Omega,
    all forming a clique. An independent set picks at most 1 from each clique.
    """
    n = len(cycle_vertex_sets)
    if n == 0:
        return 1

    # Build conflict graph on vertex sets
    adj_omega = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if cycle_vertex_sets[i] & cycle_vertex_sets[j]:
                adj_omega[i][j] = True
                adj_omega[j][i] = True

    # Enumerate independent sets of vertex SETS
    total = 0
    for mask in range(1 << n):
        independent = True
        bits = []
        for i in range(n):
            if mask & (1 << i):
                bits.append(i)
        for a in range(len(bits)):
            for b in range(a+1, len(bits)):
                if adj_omega[bits[a]][bits[b]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            k = len(bits)
            total += x**k
    return total


def held_karp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full])


def main():
    n = 8
    random.seed(42)

    # Find a tournament with t3=10 (3-cycles only from previous run)
    for trial in range(100000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        c3 = find_3cycles(adj, n)
        if len(c3) != 10:
            continue

        # CAREFULLY count ALL 5-cycles
        c5_sets = find_5cycles_bruteforce(adj, n)
        c5_count = find_5cycles_count(adj, n)

        # Count 7-cycles
        c7_count = find_7cycles_count(adj, n)

        H = held_karp(adj, n)

        print(f"Tournament with t3=10 found at trial {trial}")
        print(f"  H = {H}")
        print(f"  3-cycles: {len(c3)} vertex sets")
        print(f"  5-cycles: {len(c5_sets)} vertex sets, details: {c5_count}")
        print(f"  7-cycles: {len(c7_count)} vertex sets, details: {c7_count}")

        # Total cycles as vertex sets
        all_sets = list(set(c for c in c3))
        for vs, _ in c5_count:
            all_sets.append(vs)
        for vs, _ in c7_count:
            all_sets.append(vs)

        # Deduplicate vertex sets
        unique_sets = list(set(all_sets))
        print(f"  Total unique cycle vertex sets: {len(unique_sets)}")

        # But Omega has DIRECTED cycles as vertices, not vertex sets.
        # Multiple directed cycles on same vertex set = multiple vertices in Omega,
        # but they form a clique. For I(Omega,2), we need to account for this.

        # Actually: I(Omega,2) with directed cycles:
        # An independent set picks directed cycles, at most one from each vertex set,
        # and no two sharing a vertex.
        # Within a vertex set having d directed cycles: contributes d choices (not 1).
        # So the independence polynomial is:
        # sum over independent sets S of vertex sets: product_{vs in S} (d_vs * x)
        # Wait, no. Let me think again.

        # Omega(T) has V = {all directed odd cycles}, E = {sharing vertex}.
        # Multiple cycles on same vertex set form a CLIQUE (all share all vertices).
        # An independent set picks at most 1 cycle from each clique.
        # So alpha_k counts independent sets of size k.
        # alpha_1 = total number of directed cycles (counting multiplicity per vertex set).
        # alpha_2 = sum over disjoint vertex set pairs (i,j) of d_i * d_j.

        # So I(Omega, 2) = sum_{independent sets S of vertex-sets} product_{vs in S} (d_vs * 2)
        # where d_vs = number of directed cycles on vertex set vs.
        # Wait, no. Let's be precise.

        # If vertex set vs has d directed cycles, they form a K_d clique in Omega.
        # An independent set picks at most 1 from this clique.
        # The contribution of this clique to I(G,x) where G = whole Omega:
        # For clique K_d: I(K_d, x) = 1 + d*x.

        # For Omega = clique1 * clique2 * ... (where * means join on shared vertices):
        # It's NOT a simple product because different vertex sets might share vertices.

        # The correct approach: list ALL directed cycles, build full Omega, compute I.

        # Let me build the FULL Omega:
        all_directed_cycles = []
        for vs in c3:
            all_directed_cycles.append(vs)  # Each 3-set has exactly 1 directed 3-cycle
        for vs, d in c5_count:
            for _ in range(d):
                all_directed_cycles.append(vs)
        for vs, d in c7_count:
            for _ in range(d):
                all_directed_cycles.append(vs)

        print(f"  Total directed cycles: {len(all_directed_cycles)}")
        print(f"    3-cycles: {len(c3)}")
        print(f"    5-cycles: {sum(d for _, d in c5_count)}")
        print(f"    7-cycles: {sum(d for _, d in c7_count)}")

        # Build Omega and compute I(Omega, 2)
        nc = len(all_directed_cycles)
        if nc <= 20:  # can enumerate
            I_val = 0
            for mask in range(1 << nc):
                bits = [i for i in range(nc) if mask & (1 << i)]
                independent = True
                for a in range(len(bits)):
                    for b in range(a+1, len(bits)):
                        if all_directed_cycles[bits[a]] & all_directed_cycles[bits[b]]:
                            independent = False
                            break
                    if not independent:
                        break
                if independent:
                    I_val += 2**len(bits)
            print(f"  I(Omega, 2) = {I_val}")
            print(f"  H = I(Omega, 2)? {H == I_val}")
        else:
            print(f"  Too many cycles ({nc}) to enumerate I(Omega, 2)")

        # Score sequence
        scores = sorted([sum(adj[i]) for i in range(n)])
        print(f"  Score: {scores}")
        print()

        break  # Just check first example


if __name__ == "__main__":
    main()
