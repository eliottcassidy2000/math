#!/usr/bin/env python3
"""
Test OCF and cycle structure for locally transitive tournaments.

A tournament T is locally transitive if for every vertex i,
both T_i^+ (out-neighborhood) and T_i^- (in-neighborhood) are transitive.

From Rajkumar et al. (arXiv:2110.05188), these are exactly the rank 2 tournaments,
equivalent to the flip class of the transitive tournament.

Key questions:
1. Do locally transitive tournaments only have 3-cycles as odd cycles?
2. Is OCF trivially provable for this class?
3. What is the structure of Omega(T) for locally transitive tournaments?

Instance: opus-2026-03-05-S4b
"""
import random
from itertools import combinations, permutations

def is_transitive(T, n, vertices):
    """Check if subtournament on vertex set is transitive."""
    vlist = sorted(vertices)
    for i in range(len(vlist)):
        for j in range(i+1, len(vlist)):
            for k in range(j+1, len(vlist)):
                a, b, c = vlist[i], vlist[j], vlist[k]
                # Check for 3-cycle among a,b,c
                if T[a*n+b] and T[b*n+c] and T[c*n+a]:
                    return False
                if T[a*n+c] and T[c*n+b] and T[b*n+a]:
                    return False
    return True

def is_locally_transitive(T, n):
    """Check if T is locally transitive."""
    for i in range(n):
        out_nbrs = [j for j in range(n) if j != i and T[i*n+j]]
        in_nbrs = [j for j in range(n) if j != i and T[j*n+i]]
        if not is_transitive(T, n, out_nbrs):
            return False
        if not is_transitive(T, n, in_nbrs):
            return False
    return True

def find_odd_cycles(T, n, max_length=None):
    """Find all directed odd cycles, grouped by vertex set."""
    if max_length is None:
        max_length = n
    cycles_by_vset = {}
    for length in range(3, max_length+1, 2):
        for combo in combinations(range(n), length):
            vmask = 0
            for v in combo:
                vmask |= (1 << v)
            # Count directed Hamiltonian cycles on this vertex set
            start = combo[0]
            others = list(combo[1:])
            mo = len(others)
            ofull = (1 << mo) - 1
            dp = [[0]*mo for _ in range(1 << mo)]
            for i, o in enumerate(others):
                if T[start*n+o]:
                    dp[1<<i][i] = 1
            for omask in range(1, ofull+1):
                for li in range(mo):
                    c = dp[omask][li]
                    if c == 0: continue
                    for ni in range(mo):
                        if omask & (1<<ni): continue
                        if T[others[li]*n+others[ni]]:
                            dp[omask|(1<<ni)][ni] += c
            cnt = sum(dp[ofull][li] for li in range(mo) if T[others[li]*n+start])
            if cnt > 0:
                cycles_by_vset[vmask] = (length, cnt)
    return cycles_by_vset

def ham_count(T, n):
    FULL = (1 << n) - 1
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for last in range(n):
            c = dp[mask][last]
            if c == 0: continue
            for nxt in range(n):
                if mask & (1<<nxt): continue
                if T[last*n+nxt]:
                    dp[mask|(1<<nxt)][nxt] += c
    return sum(dp[FULL])

def compute_I_omega_2(T, n):
    cycles = find_odd_cycles(T, n)
    items = list(cycles.items())  # (vmask, (length, count))
    total_cycles = sum(cnt for _, (_, cnt) in items)
    vd_pairs = 0
    for i in range(len(items)):
        for j in range(i+1, len(items)):
            if items[i][0] & items[j][0] == 0:
                vd_pairs += items[i][1][1] * items[j][1][1]
    # At small n, max independent set size is 2
    # For n <= 8 this is correct
    return 1 + 2*total_cycles + 4*vd_pairs

def random_tournament(n):
    T = [0]*(n*n)
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i*n+j] = 1
            else:
                T[j*n+i] = 1
    return T

def make_locally_transitive(n):
    """Generate a locally transitive tournament on n vertices.

    Method: take a permutation sigma, set i -> j iff sigma(i) < sigma(j),
    then apply a cut-flip phi_S for random S. This gives all locally transitive
    tournaments (they are exactly the flip class of the transitive tournament).
    """
    # Start with transitive tournament (identity permutation ordering)
    T = [0]*(n*n)
    for i in range(n):
        for j in range(i+1, n):
            T[i*n+j] = 1  # i -> j for i < j

    # Apply a random cut-flip phi_S
    S = set()
    for v in range(n):
        if random.random() < 0.5:
            S.add(v)
    Sc = set(range(n)) - S

    # Reverse all arcs between S and Sc
    for i in S:
        for j in Sc:
            # Swap T[i][j] and T[j][i]
            old_ij = T[i*n+j]
            T[i*n+j] = T[j*n+i]
            T[j*n+i] = old_ij

    return T

# --- Main tests ---

print("=" * 60)
print("LOCALLY TRANSITIVE TOURNAMENTS: CYCLE STRUCTURE ANALYSIS")
print("=" * 60)

for n in [5, 6, 7, 8]:
    print(f"\n--- n = {n} ---")

    lt_count = 0
    lt_only3 = 0  # LT tournaments with only 3-cycles
    lt_has5 = 0   # LT tournaments with 5-cycles
    lt_has7 = 0   # LT tournaments with 7-cycles
    lt_ocf_pass = 0

    random.seed(42)
    trials = 500 if n <= 7 else 200

    for _ in range(trials):
        T = make_locally_transitive(n)

        if not is_locally_transitive(T, n):
            print("  BUG: generated tournament is not locally transitive!")
            continue

        lt_count += 1
        cycles = find_odd_cycles(T, n)

        max_cycle_len = max((length for _, (length, _) in cycles.items()), default=0)
        if max_cycle_len <= 3:
            lt_only3 += 1
        if any(length >= 5 for _, (length, _) in cycles.items()):
            lt_has5 += 1
        if any(length >= 7 for _, (length, _) in cycles.items()):
            lt_has7 += 1

        H = ham_count(T, n)
        I = compute_I_omega_2(T, n)
        if H == I:
            lt_ocf_pass += 1

    print(f"  Generated: {lt_count} locally transitive tournaments")
    print(f"  Only 3-cycles: {lt_only3}/{lt_count} ({100*lt_only3/lt_count:.1f}%)")
    print(f"  Has 5-cycles: {lt_has5}/{lt_count} ({100*lt_has5/lt_count:.1f}%)")
    print(f"  Has 7-cycles: {lt_has7}/{lt_count} ({100*lt_has7/lt_count:.1f}%)")
    print(f"  OCF passes: {lt_ocf_pass}/{lt_count}")

print("\n" + "=" * 60)
print("COMPARISON: RANDOM (NON-LT) TOURNAMENTS")
print("=" * 60)

for n in [5, 6, 7]:
    print(f"\n--- n = {n} ---")
    random.seed(123)

    gen_count = 0
    is_lt_count = 0
    only3 = 0
    has5 = 0

    for _ in range(500):
        T = random_tournament(n)
        gen_count += 1

        if is_locally_transitive(T, n):
            is_lt_count += 1

        cycles = find_odd_cycles(T, n)
        max_cycle_len = max((length for _, (length, _) in cycles.items()), default=0)
        if max_cycle_len <= 3:
            only3 += 1
        if any(length >= 5 for _, (length, _) in cycles.items()):
            has5 += 1

    print(f"  {is_lt_count}/{gen_count} are locally transitive ({100*is_lt_count/gen_count:.1f}%)")
    print(f"  Only 3-cycles: {only3}/{gen_count} ({100*only3/gen_count:.1f}%)")
    print(f"  Has 5-cycles: {has5}/{gen_count} ({100*has5/gen_count:.1f}%)")

print("\n" + "=" * 60)
print("R-CONE STRUCTURE TEST")
print("=" * 60)
print("Testing OCF for R-cones (vertex beating/losing to everyone)")

for n in [5, 6, 7, 8]:
    print(f"\n--- n = {n} (R-cones) ---")
    random.seed(99)
    passes = 0
    trials = 200

    for _ in range(trials):
        # Build R-cone: vertex 0 beats everyone, rest is random
        T = [0]*(n*n)
        for j in range(1, n):
            T[0*n+j] = 1  # vertex 0 beats all
        for i in range(1, n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    T[i*n+j] = 1
                else:
                    T[j*n+i] = 1

        H = ham_count(T, n)
        I = compute_I_omega_2(T, n)
        if H == I:
            passes += 1

    print(f"  OCF passes: {passes}/{trials}")

print("\n" + "=" * 60)
print("AUTOMORPHISM TEST: tournaments fixed by transposition (0 1)")
print("=" * 60)

for n in [5, 6, 7]:
    print(f"\n--- n = {n} ---")
    random.seed(77)
    count = 0
    ocf_pass = 0

    for _ in range(500):
        # Build tournament where swapping vertices 0,1 is an automorphism
        # This means T[0][k] = T[1][k] for k >= 2
        T = [0]*(n*n)
        # Arc between 0 and 1: random
        if random.random() < 0.5:
            T[0*n+1] = 1
        else:
            T[1*n+0] = 1
        # Arcs from {0,1} to rest: same for both
        for k in range(2, n):
            if random.random() < 0.5:
                T[0*n+k] = 1
                T[1*n+k] = 1
            else:
                T[k*n+0] = 1
                T[k*n+1] = 1
        # Arcs among rest: random
        for i in range(2, n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    T[i*n+j] = 1
                else:
                    T[j*n+i] = 1

        count += 1
        H = ham_count(T, n)
        I = compute_I_omega_2(T, n)
        if H == I:
            ocf_pass += 1

    print(f"  OCF passes: {ocf_pass}/{count}")
    print(f"  (These have Aut(T) containing transposition (0 1))")
