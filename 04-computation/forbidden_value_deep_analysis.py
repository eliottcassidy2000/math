"""
forbidden_value_deep_analysis.py — opus-S74

Deep investigation of which H values are truly forbidden, and whether
the k-nacci trace / Mersenne number connection is real or numerological.

Questions:
1. Is 31 = 2^5 - 1 achievable? (Mersenne, but claimed in writeups to be "forbidden")
2. Is 63 = 2^6 - 1 achievable?
3. What are ALL permanent gaps up to some bound?
4. Does the k-nacci trace sequence [1,3,7,11,21,39,71,...] predict forbidden values?
5. What is the actual mechanism that makes EXACTLY 7 and 21 forbidden?
"""
from itertools import permutations
import numpy as np
from collections import defaultdict

def tournament_from_bits(bits, n):
    """Create adjacency matrix from bit encoding."""
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hamiltonian_paths(A, n):
    """Count Hamiltonian paths using Held-Karp DP."""
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    count += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    count += 1
    return count

# Part 1: Exhaustive H-spectrum at n=5,6,7
print("=" * 60)
print("PART 1: Exhaustive H-spectrum analysis")
print("=" * 60)

for n in [5, 6, 7]:
    num_edges = n * (n - 1) // 2
    total = 1 << num_edges
    h_values = set()
    h_examples = {}  # Store one example for interesting H values

    for bits in range(total):
        A = tournament_from_bits(bits, n)
        h = count_hamiltonian_paths(A, n)
        h_values.add(h)
        if h in [7, 21, 31, 63, 127]:
            if h not in h_examples:
                h_examples[h] = bits

    odd_in_range = set(range(1, max(h_values)+1, 2))
    gaps = sorted(odd_in_range - h_values)

    print(f"\nn={n}: {total} tournaments, {len(h_values)} distinct H values")
    print(f"  Range: [{min(h_values)}, {max(h_values)}]")
    print(f"  Gaps (odd values NOT achieved): {gaps[:30]}")

    # Check specific values
    for target in [7, 21, 31, 63, 127]:
        if target <= max(h_values):
            status = "ACHIEVED" if target in h_values else "GAP"
            extra = ""
            if target in h_examples:
                extra = f" (example: bits={h_examples[target]})"
            print(f"  H={target} ({target}=2^{int(np.log2(target+1))}-1 Mersenne): {status}{extra}")

print("\n" + "=" * 60)
print("PART 2: k-nacci trace sequence vs forbidden values")
print("=" * 60)

def knacci_traces(k, nmax):
    """Compute Tr(M_k^n) for n=1..nmax using Newton's identities."""
    # Characteristic polynomial of k-nacci: x^k - x^{k-1} - ... - x - 1
    # Elementary symmetric polynomials: e_1=1, e_j = -1 for j=2..k
    e = [0] * (k + 1)
    e[1] = 1
    for j in range(2, k + 1):
        e[j] = -1

    p = [0] * (nmax + 1)
    for n in range(1, nmax + 1):
        if n <= k:
            # Newton's identity: p_n = sum_{j=1}^{n-1} (-1)^{j-1} e_j p_{n-j} + (-1)^{n-1} n e_n
            s = 0
            for j in range(1, n):
                s += ((-1) ** (j - 1)) * e[j] * p[n - j]
            s += ((-1) ** (n - 1)) * n * e[n]
            p[n] = s
        else:
            # Recurrence: p_n = sum_{j=1}^{k} (-1)^{j-1} e_j p_{n-j}
            s = 0
            for j in range(1, k + 1):
                s += ((-1) ** (j - 1)) * e[j] * p[n - j]
            p[n] = s
    return p[1:]

# All achievable H values at n=7
n = 7
num_edges = n * (n - 1) // 2
total = 1 << num_edges
all_h_n7 = set()
for bits in range(total):
    A = tournament_from_bits(bits, n)
    h = count_hamiltonian_paths(A, n)
    all_h_n7.add(h)

print("\nk-nacci traces and their achievability as H values (checked at n≤7):")
for k in [2, 3, 4, 5]:
    traces = knacci_traces(k, 12)
    print(f"\n  k={k}-nacci: {traces}")
    for i, t in enumerate(traces):
        if t <= max(all_h_n7) and t > 0:
            status = "ACHIEVED" if t in all_h_n7 else "FORBIDDEN"
            marker = " *** FORBIDDEN ***" if status == "FORBIDDEN" else ""
            print(f"    p_{i+1} = {t}: {status}{marker}")

print("\n" + "=" * 60)
print("PART 3: Deeper OCF analysis — why 7 and 21 specifically?")
print("=" * 60)

# For each impossible H, show what (alpha_1, alpha_2, ...) would be needed
# H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
# So (H-1)/2 = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...

def ocf_decompositions(h, max_terms=6):
    """Find all decompositions of (h-1)/2 into alpha_1 + 2*alpha_2 + 4*alpha_3 + ..."""
    target = (h - 1) // 2
    results = []

    def backtrack(remaining, depth, current):
        if remaining == 0:
            results.append(tuple(current))
            return
        if depth > max_terms:
            return
        coeff = 2 ** (depth - 1) if depth >= 2 else 1
        for val in range(remaining // coeff + 1):
            current.append(val)
            backtrack(remaining - val * coeff, depth + 1, current)
            current.pop()

    backtrack(target, 1, [])
    return results

for h in [3, 5, 7, 9, 11, 13, 15, 21, 31]:
    decomps = ocf_decompositions(h, max_terms=4)
    print(f"\nH={h}: (H-1)/2={(h-1)//2}")
    print(f"  Need alpha_1 + 2*i_2 + 4*i_3 + ... = {(h-1)//2}")
    print(f"  Number of decompositions: {len(decomps)}")
    if len(decomps) <= 20:
        for d in decomps[:10]:
            alpha1 = d[0] if len(d) > 0 else 0
            i2 = d[1] if len(d) > 1 else 0
            i3 = d[2] if len(d) > 2 else 0
            print(f"    (alpha_1={alpha1}, i_2={i2}, i_3={i3})")

print("\n" + "=" * 60)
print("PART 4: Achievable (alpha_1, i_2) pairs at n=5,6,7")
print("=" * 60)

for n in [5, 6, 7]:
    num_edges = n * (n - 1) // 2
    total = 1 << num_edges

    alpha_i2_pairs = set()
    alpha_h_map = defaultdict(set)

    for bits in range(total):
        A = tournament_from_bits(bits, n)

        # Find all 3-cycles
        cycles_3 = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if A[i][j] and A[j][k] and A[k][i]:
                        cycles_3.append((i,j,k))
                    elif A[i][k] and A[k][j] and A[j][i]:
                        cycles_3.append((i,k,j))

        # Find 5-cycles (directed)
        cycles_5 = []
        if n >= 5:
            for combo in __import__('itertools').combinations(range(n), 5):
                verts = list(combo)
                for perm in permutations(verts):
                    if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                        canon = min(perm[i:] + perm[:i] for i in range(5))
                        if tuple(canon) not in [tuple(c) for c in cycles_5]:
                            cycles_5.append(list(canon))

        all_cycles = [(set(c), '3') for c in cycles_3] + [(set(c), '5') for c in cycles_5]
        alpha_1 = len(all_cycles)

        # Count i_2 (independent pairs)
        i_2 = 0
        for a in range(len(all_cycles)):
            for b in range(a+1, len(all_cycles)):
                if all_cycles[a][0].isdisjoint(all_cycles[b][0]):
                    i_2 += 1

        h = count_hamiltonian_paths(A, n)
        alpha_i2_pairs.add((alpha_1, i_2))
        alpha_h_map[(alpha_1, i_2)].add(h)

    print(f"\nn={n}: Achievable (alpha_1, i_2) pairs:")
    for pair in sorted(alpha_i2_pairs):
        h_vals = sorted(alpha_h_map[pair])
        print(f"  ({pair[0]:2d}, {pair[1]:2d}) -> H in {h_vals}")

print("\n" + "=" * 60)
print("PART 5: Is 31 achievable? Finding tournament with H=31")
print("=" * 60)

# We already know from n=7 exhaustive. Let's find the actual tournament.
for n in [6, 7]:
    num_edges = n * (n - 1) // 2
    total = 1 << num_edges
    found = False
    for bits in range(total):
        A = tournament_from_bits(bits, n)
        h = count_hamiltonian_paths(A, n)
        if h == 31:
            print(f"\nH=31 ACHIEVED at n={n}!")
            print(f"  bits = {bits}")
            print(f"  Adjacency matrix:")
            print(A)
            t3 = count_3cycles(A, n)
            print(f"  t3 = {t3}")
            found = True
            break
    if found:
        break
    else:
        if n == 7:
            print(f"H=31 not found at n={n}")

print("\n" + "=" * 60)
print("PART 6: Summary — which Mersenne/tribonacci values are forbidden?")
print("=" * 60)

tribonacci = knacci_traces(3, 10)
mersenne = [2**k - 1 for k in range(2, 11)]

print(f"\nTribonacci traces: {tribonacci}")
print(f"Mersenne numbers: {mersenne}")

all_achieved = set()
for n in [5, 6, 7]:
    num_edges = n * (n - 1) // 2
    total = 1 << num_edges
    for bits in range(total):
        A = tournament_from_bits(bits, n)
        h = count_hamiltonian_paths(A, n)
        all_achieved.add(h)

max_h = max(all_achieved)
print(f"\nAll odd values achieved at n≤7: range [1, {max_h}], {len(all_achieved)} values")

permanent_gaps = []
for h in range(1, max_h + 1, 2):
    if h not in all_achieved:
        permanent_gaps.append(h)

print(f"Permanent gaps (odd values NOT achieved at ANY n≤7): {permanent_gaps[:30]}")

print(f"\nCheck tribonacci traces:")
for i, t in enumerate(tribonacci):
    if 0 < t <= max_h:
        status = "FORBIDDEN" if t not in all_achieved else "achievable"
        print(f"  p_{i+1} = {t}: {status}")

print(f"\nCheck Mersenne numbers:")
for k in range(2, 11):
    m = 2**k - 1
    if m <= max_h:
        status = "FORBIDDEN" if m not in all_achieved else "achievable"
        print(f"  2^{k}-1 = {m}: {status}")
