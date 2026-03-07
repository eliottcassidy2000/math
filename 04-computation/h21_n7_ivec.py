#!/usr/bin/env python3
"""
At n=7, sample tournaments with i_1 values that could give H=21
and check what i_2 they produce.
"""
import random, time, sys
from itertools import combinations, permutations
from collections import Counter, defaultdict

def H_matrix(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not (mask & (1 << v)) or c == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += c
    return sum(dp[full])

def find_directed_odd_cycles(A, n):
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            first = verts[0]
            for perm in permutations(verts[1:]):
                path = (first,) + perm
                valid = True
                for i in range(length):
                    if not A[path[i]][path[(i+1) % length]]:
                        valid = False
                        break
                if valid:
                    cycles.append(path)
    return cycles

def ocf_i1_i2(A, n):
    """Compute just i_1 and i_2 (faster than full enumeration)."""
    cycles = find_directed_odd_cycles(A, n)
    m = len(cycles)
    i1 = m

    if m <= 1:
        return i1, 0

    vsets = [frozenset(c) for c in cycles]
    i2 = 0
    for a in range(m):
        for b in range(a+1, m):
            if not (vsets[a] & vsets[b]):  # disjoint = independent pair
                i2 += 1
    return i1, i2


# Sample n=7 tournaments and record (i_1, i_2, H)
n = 7
print(f"Sampling n={n} tournaments for (i_1, i_2, H) triples...")
print("Looking for combinations that could give H=21")
print()

# Target: find tournaments with H near 21 or i_1 values in [4,6,8,10]
random.seed(42)
i1_i2_h = defaultdict(list)  # (i1, i2) -> list of H values
h_near21 = []

t0 = time.time()
for trial in range(5000):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    h = H_matrix(A, n)

    # Only compute OCF for interesting cases
    if h <= 50 or (15 <= h <= 27):
        i1, i2 = ocf_i1_i2(A, n)
        i1_i2_h[(i1, i2)].append(h)

        # Verify OCF
        expected_h = 1 + 2*i1 + 4*i2  # only if i_3=0
        if h != expected_h and h <= 50:
            # Has higher-order terms
            pass

        if abs(h - 21) <= 4:
            h_near21.append((h, i1, i2))

    if (trial+1) % 1000 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/5000 ({elapsed:.0f}s)", flush=True)

elapsed = time.time() - t0
print(f"Total: {elapsed:.0f}s\n")

# Show all (i_1, i_2) -> H mappings
print("(i_1, i_2) -> H values (distinct):")
for key in sorted(i1_i2_h.keys()):
    h_vals = sorted(set(i1_i2_h[key]))
    i1, i2 = key
    s = 2*i1 + 4*i2
    print(f"  (i_1={i1:2d}, i_2={i2:2d}): S=2*i_1+4*i_2={s:3d}, H values={h_vals}")

# Check if any gives H=21
print(f"\nTournaments near H=21:")
for h, i1, i2 in sorted(h_near21):
    s = 2*i1 + 4*i2
    remainder = h - 1 - s
    print(f"  H={h}: i_1={i1}, i_2={i2}, S={s}, higher_order_contrib={remainder}")

# Check: for each i_1 that could target H=21, what i_2 appears?
print(f"\nFor H=21, need 2*i_1 + 4*i_2 + higher = 20:")
target_i1 = [4, 6, 8, 10]
for ti1 in target_i1:
    matching = [(i2, h_vals) for (i1, i2), h_vals in i1_i2_h.items() if i1 == ti1]
    if matching:
        for i2, h_vals in sorted(matching):
            needed_remainder = 20 - 2*ti1 - 4*i2
            print(f"  i_1={ti1}, i_2={i2}: H values={sorted(set(h_vals))}, "
                  f"need remainder={needed_remainder} from i_3+")
    else:
        print(f"  i_1={ti1}: NOT FOUND in sample")

sys.stdout.flush()
