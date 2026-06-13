#!/usr/bin/env python3
"""
Verify critical claims from the n=8 deep dive.
1. Is T_657 minus vertex 3 actually isomorphic to P(7)?
2. Is the signed position identity really failing, or is my code wrong?
3. Build the FULL Omega (including 5-cycles and 7-cycles) and check OCF.

Instance: opus-2026-03-05-S8
"""

from itertools import combinations, permutations
from collections import Counter


def count_ham_dp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
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


def build_sc(n, alpha):
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
            T[i][j] = val; T[j][i] = 1 - val
            if linked:
                a, b = linked
                T[a][b] = val; T[b][a] = 1 - val
        tournaments.append(T)
    return tournaments


def is_isomorphic(T1, T2, n):
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if T1[i][j] != T2[perm[i]][perm[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return perm
    return None


def find_odd_cycles(T, n, max_len=None):
    """Find all directed odd cycles of length 3, 5, 7, ..."""
    if max_len is None:
        max_len = n
    cycles = []
    for length in range(3, max_len + 1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(T[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    # Canonicalize: start with min vertex
                    min_idx = perm.index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    cycles.append(canon)
                    break
    # Deduplicate
    return list(set(cycles))


def independence_poly(adj, m):
    if m > 25:
        return None
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        verts = [i for i in range(m) if mask & (1 << i)]
        ok = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            coeffs[len(verts)] += 1
    return coeffs


def main():
    n = 8
    alpha = [7, 6, 5, 4, 3, 2, 1, 0]

    # Find T_657
    tournaments = build_sc(n, alpha)
    T_657 = None
    T_A = None

    for T in tournaments:
        out_deg = sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)])
        if out_deg != [3, 3, 3, 3, 4, 4, 4, 4]:
            continue
        H = count_ham_dp(T, n)
        if H == 657:
            # Check |Aut|>1 quickly
            scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
            if len(set(scores)) < n:
                T_657 = T
                break

    if T_657 is None:
        # Just find ANY T with H=657
        for T in tournaments:
            H = count_ham_dp(T, n)
            if H == 657:
                T_657 = T
                break

    print("=" * 60)
    print("VERIFICATION 1: T_657 minus vertex 3 vs P(7)")
    print("=" * 60)

    # Build P(7)
    P7 = [[0]*7 for _ in range(7)]
    qr = {1, 2, 4}  # Quadratic residues mod 7
    for i in range(7):
        for j in range(7):
            if i != j and (j - i) % 7 in qr:
                P7[i][j] = 1

    print("P(7) adjacency matrix:")
    for row in P7:
        print(f"  {''.join(str(x) for x in row)}")
    print(f"P(7) H = {count_ham_dp(P7, 7)}")
    print(f"P(7) out-degrees: {[sum(P7[i][j] for j in range(7) if j!=i) for i in range(7)]}")

    if T_657:
        print(f"\nT_657 adjacency matrix:")
        for row in T_657:
            print(f"  {''.join(str(x) for x in row)}")

        for v in range(n):
            verts = [u for u in range(n) if u != v]
            subT = [[T_657[verts[i]][verts[j]] for j in range(7)] for i in range(7)]
            sub_out = [sum(subT[i][j] for j in range(7) if j!=i) for i in range(7)]
            if all(d == 3 for d in sub_out):
                sub_H = count_ham_dp(subT, 7)
                iso = is_isomorphic(subT, P7, 7)
                print(f"\n  T_657 - vertex {v}:")
                print(f"    H = {sub_H}")
                print(f"    Isomorphic to P(7): {iso is not None}")
                if iso:
                    print(f"    Witnessing permutation: {iso}")

    print(f"\n{'='*60}")
    print("VERIFICATION 2: Signed position identity")
    print("=" * 60)

    # Test on SMALL tournament where we KNOW it should hold
    # Use n=5 transitive + one flip
    n5 = 5
    T5 = [[0]*5 for _ in range(5)]
    for i in range(5):
        for j in range(i+1, 5):
            T5[i][j] = 1
    # Flip 0->1 to 1->0
    T5[0][1] = 0; T5[1][0] = 1

    paths5 = []
    def enum(path, visited, TT, nn):
        if len(path) == nn:
            paths5.append(tuple(path))
            return
        last = path[-1]
        for v in range(nn):
            if v not in visited and TT[last][v]:
                visited.add(v)
                path.append(v)
                enum(path, visited, TT, nn)
                path.pop()
                visited.remove(v)
    for s in range(n5):
        enum([s], {s}, T5, n5)

    print(f"n=5 test tournament, H={len(paths5)}")
    for i in range(n5):
        for j in range(i+1, n5):
            if T5[i][j]:
                # "i before j" means i appears at a smaller index than j in the path
                lhs = sum((-1)**P.index(i) for P in paths5 if P.index(i) < P.index(j))
                rhs = sum((-1)**P.index(j) for P in paths5 if P.index(j) < P.index(i))
                status = "OK" if lhs == rhs else "FAIL"
                print(f"  arc {i}->{j}: LHS={lhs}, RHS={rhs} [{status}]")

    # Now test on P(7) itself
    print(f"\nP(7) test:")
    paths7 = []
    for s in range(7):
        enum([s], {s}, P7, 7)
    print(f"  H(P7) = {len(paths7)}")
    fails7 = 0
    for i in range(7):
        for j in range(i+1, 7):
            if P7[i][j]:
                lhs = sum((-1)**P.index(i) for P in paths7 if P.index(i) < P.index(j))
                rhs = sum((-1)**P.index(j) for P in paths7 if P.index(j) < P.index(i))
                if lhs != rhs:
                    fails7 += 1
                    if fails7 <= 3:
                        print(f"  FAIL: arc {i}->{j}: LHS={lhs}, RHS={rhs}")
    print(f"  Failures: {fails7} / {sum(sum(row) for row in P7)//2}")

    print(f"\n{'='*60}")
    print("VERIFICATION 3: Full Omega OCF check on T_657")
    print("=" * 60)

    if T_657:
        # Find ALL odd cycles (3, 5, 7)
        print("Finding all odd cycles of T_657...")
        all_cycles = find_odd_cycles(T_657, n)
        len_dist = Counter(len(c) for c in all_cycles)
        print(f"  Cycle length distribution: {dict(len_dist)}")

        m = len(all_cycles)
        print(f"  Total odd cycles: {m}")

        # Build full Omega
        cycle_sets = [set(c) for c in all_cycles]
        adj = [[0]*m for _ in range(m)]
        for i in range(m):
            for j in range(i+1, m):
                if cycle_sets[i] & cycle_sets[j]:
                    adj[i][j] = adj[j][i] = 1

        # Compute I(Omega, 2)
        if m <= 22:
            coeffs = independence_poly(adj, m)
            I2 = sum(c * (2**k) for k, c in enumerate(coeffs))
            print(f"  I(full Omega, 2) = {I2}")
            print(f"  H(T_657) = 657")
            print(f"  OCF: {I2 == 657}")
            print(f"  Coefficients: {[c for c in coeffs if c > 0]}")
        else:
            print(f"  Too many cycles ({m}) for exact independence polynomial")
            print(f"  (Would need specialized algorithm)")


if __name__ == "__main__":
    main()
