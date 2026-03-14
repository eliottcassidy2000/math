#!/usr/bin/env python3
"""
WILD IDEAS — FINAL BATCH
opus-2026-03-14-S71m

More wild explorations:
1. The "H-code" — can tournament be reconstructed from H and local info?
2. H and Ramsey numbers
3. Tournament algebra — direct sum and tensor product
4. H under random walk on S_n (relabeling vertices)
5. The H-gap structure: which odd numbers are NOT achievable?
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math

def adj_matrix(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hp(n, A):
    dp = [[0]*n for _ in range(1 << n)]
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
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

print("=" * 70)
print("WILD IDEAS — FINAL BATCH")
print("opus-2026-03-14-S71m")
print("=" * 70)

# ======================================================================
# IDEA 1: THE H-GAP STRUCTURE
# ======================================================================
print("\n" + "=" * 70)
print("IDEA 1: THE H-GAP STRUCTURE")
print("=" * 70)

print("""
  Which odd numbers are NOT achievable as H values?
  Let gaps(n) = {odd numbers in [1, max_H] that are not in the H-spectrum}.
""")

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    total = 1 << m

    H_set = set()
    for bits in range(total):
        A = adj_matrix(n, bits)
        H_set.add(count_hp(n, A))

    max_H = max(H_set)
    all_odd = set(range(1, max_H + 1, 2))
    gaps = sorted(all_odd - H_set)

    print(f"\n  n={n}: H-spectrum = {sorted(H_set)}")
    print(f"    max H = {max_H}")
    print(f"    Gaps (missing odd values): {gaps}")
    print(f"    Number of gaps: {len(gaps)}")
    print(f"    Gap fraction: {len(gaps)}/{len(all_odd)} = {len(gaps)/len(all_odd):.4f}")

    # Factor the gaps
    if gaps:
        print(f"    Gap factorizations:")
        for g in gaps:
            factors = []
            temp = g
            for p in [3, 5, 7, 11, 13, 17, 19, 23]:
                while temp % p == 0:
                    factors.append(p)
                    temp //= p
            if temp > 1:
                factors.append(temp)
            print(f"      {g} = {' × '.join(map(str, factors))}")

# ======================================================================
# IDEA 2: THE H SPECTRUM AS AN ADDITIVE STRUCTURE
# ======================================================================
print("\n" + "=" * 70)
print("IDEA 2: ADDITIVE STRUCTURE OF H-SPECTRUM")
print("=" * 70)

print("""
  Is the H-spectrum closed under any algebraic operation?
  For example, if h1 and h2 are achievable, is h1 + h2 - 1 achievable?
  (This would make the set {(H-1)/2} closed under addition.)
""")

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    total = 1 << m

    H_set = set()
    for bits in range(total):
        A = adj_matrix(n, bits)
        H_set.add(count_hp(n, A))

    max_H = max(H_set)
    H_sorted = sorted(H_set)

    # Check: is {(H-1)/2 : H in spectrum} a numerical semigroup?
    k_set = set((h-1)//2 for h in H_set)
    max_k = max(k_set)

    # Check closure: a+b in k_set for all a,b in k_set with a+b <= max_k?
    closed = True
    violations = []
    for a in k_set:
        for b in k_set:
            if a + b <= max_k and a + b not in k_set:
                closed = False
                violations.append((a, b, a+b))

    print(f"\n  n={n}: k = (H-1)/2 spectrum = {sorted(k_set)}")
    if closed:
        print(f"    k-set is CLOSED under addition up to {max_k}")
    else:
        print(f"    k-set NOT closed: {len(violations)} violations")
        for a, b, s in violations[:5]:
            print(f"      k={a} + k={b} = {s} NOT in spectrum (H={2*s+1} missing)")

    # Check: is {H mod 2^k} eventually periodic?
    for k in [2, 3, 4]:
        residues = sorted(set(h % (2**k) for h in H_set))
        print(f"    H mod {2**k}: {residues}")

# ======================================================================
# IDEA 3: H AND THE PERMUTOHEDRON
# ======================================================================
print("\n" + "=" * 70)
print("IDEA 3: H AND THE PERMUTOHEDRON")
print("=" * 70)

print("""
  The permutohedron Pi_n is the convex hull of all permutations of (1,2,...,n).
  Its vertices are the n! permutations.
  Its edges connect permutations that differ by an adjacent transposition.

  For a tournament T, the set of Hamiltonian paths forms a subset of
  the permutohedron vertices (those permutations that are valid HPs).

  H(T) = |{vertices of Pi_n that are valid HPs of T}|.

  QUESTION: Does the CONVEX HULL of the valid HP permutations have
  interesting properties? Is its volume related to H?
""")

# For n=4, compute the "HP polytope" for each tournament class
n = 4
m = n * (n-1) // 2

H_to_perms = defaultdict(set)
for bits in range(1 << m):
    A = adj_matrix(n, bits)
    H = count_hp(n, A)
    for perm in permutations(range(n)):
        valid = True
        for k in range(n-1):
            if not A[perm[k]][perm[k+1]]:
                valid = False
                break
        if valid:
            H_to_perms[(bits, H)].add(perm)

# Aggregate by H value
HP_by_H = defaultdict(list)
for (bits, H), perms in H_to_perms.items():
    HP_by_H[H].append(perms)

for H_val in sorted(set(h for _, h in H_to_perms.keys())):
    sample_perms = HP_by_H[H_val][0]
    print(f"\n  H={H_val}: {len(sample_perms)} HPs")
    # Show the permutations
    for p in sorted(sample_perms)[:5]:
        print(f"    {p}")
    if len(sample_perms) > 5:
        print(f"    ... ({len(sample_perms)} total)")

    # Inversion count distribution
    inv_counts = []
    for p in sample_perms:
        inv = sum(1 for i in range(n) for j in range(i+1, n) if p[i] > p[j])
        inv_counts.append(inv)
    print(f"    Inversion counts: {sorted(Counter(inv_counts).items())}")

# ======================================================================
# IDEA 4: H AND GRAPH EIGENVALUES
# ======================================================================
print("\n" + "=" * 70)
print("IDEA 4: SKEW-ADJACENCY EIGENVALUES AND H")
print("=" * 70)

print("""
  For tournament T, the SKEW-ADJACENCY matrix S has S[i][j] = 1 if i→j, -1 if j→i, 0 on diagonal.
  S is skew-symmetric: S^T = -S.
  Eigenvalues of S are purely imaginary: ±iλ_1, ±iλ_2, ..., ±iλ_{n/2} (plus 0 for odd n).

  The characteristic polynomial det(xI - S) encodes these eigenvalues.
  For n=3: the only tournament up to iso gives det(xI - S).

  QUESTION: Does H relate to the eigenvalues of S?
  Specifically, det(S) (for even n) = Pfaffian^2 = product of λ_k.
""")

# Compute the characteristic polynomial of S for small tournaments
for n in [3, 4]:
    m = n * (n-1) // 2
    H_seen = {}

    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        ss = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

        if (H, ss) in H_seen:
            continue
        H_seen[(H, ss)] = bits

        # Build skew-adjacency matrix
        S = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    S[i][j] = 1 if A[i][j] else -1

        # Compute det(S) for small n
        if n == 3:
            det_S = (S[0][0]*(S[1][1]*S[2][2]-S[1][2]*S[2][1])
                    -S[0][1]*(S[1][0]*S[2][2]-S[1][2]*S[2][0])
                    +S[0][2]*(S[1][0]*S[2][1]-S[1][1]*S[2][0]))
            print(f"\n  n=3, H={H}, ss={ss}: det(S) = {det_S}")
        elif n == 4:
            # 4x4 determinant
            def det4(M):
                result = 0
                for j in range(4):
                    minor = [[M[i][k] for k in range(4) if k != j] for i in range(1, 4)]
                    det3 = (minor[0][0]*(minor[1][1]*minor[2][2]-minor[1][2]*minor[2][1])
                           -minor[0][1]*(minor[1][0]*minor[2][2]-minor[1][2]*minor[2][0])
                           +minor[0][2]*(minor[1][0]*minor[2][1]-minor[1][1]*minor[2][0]))
                    result += ((-1)**j) * M[0][j] * det3
                return result

            det_S = det4(S)
            # Also S^2
            S2 = [[sum(S[i][k]*S[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
            tr_S2 = sum(S2[i][i] for i in range(n))
            print(f"\n  n=4, H={H}, ss={ss}: det(S) = {det_S}, tr(S^2) = {tr_S2}")

            # S^2 = -S*S^T for skew-symmetric S
            # tr(S^2) = -tr(SS^T) = -sum |S[i][j]|^2 = -(n^2 - n)
            # since |S[i][j]| = 1 for i≠j and S[i][i] = 0
            print(f"    Expected tr(S^2) = -{n*(n-1)} = {-n*(n-1)}")

# ======================================================================
# IDEA 5: THE H-SPECTRUM AND FIBONACCI NUMBERS
# ======================================================================
print("\n" + "=" * 70)
print("IDEA 5: H-SPECTRUM AND FIBONACCI")
print("=" * 70)

print("""
  Fibonacci sequence: 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, ...
  Tribonacci: 1, 1, 2, 4, 7, 13, 24, 44, 81, ...

  At n=5: H-spectrum = {1, 3, 5, 9, 11, 13, 15}
  Gap at H=7: the 7th value is missing!
  7 = Φ_3(2) = tribonacci(5) from our tribonacci sequence!

  At n=6: gaps = {7, 21, 35, 39}
  7 = Φ_3(2), 21 = 3*7 = Φ_3(2)*3, 35 = 5*7 = Φ_3(2)*5
  39 = 3*13

  Are all gaps multiples of 7 for n≤6?
""")

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    total = 1 << m
    H_set = set()
    for bits in range(total):
        A = adj_matrix(n, bits)
        H_set.add(count_hp(n, A))

    max_H = max(H_set)
    gaps = sorted(set(range(1, max_H + 1, 2)) - H_set)

    div7 = [g for g in gaps if g % 7 == 0]
    not_div7 = [g for g in gaps if g % 7 != 0]

    print(f"\n  n={n}: gaps = {gaps}")
    print(f"    Divisible by 7: {div7}")
    print(f"    NOT divisible by 7: {not_div7}")

# ======================================================================
# IDEA 6: HOW MANY TOURNAMENTS DETERMINE H?
# ======================================================================
print("\n" + "=" * 70)
print("IDEA 6: H-DETERMINING CAPACITY")
print("=" * 70)

print("""
  H partitions the 2^m tournaments into |H-spectrum| classes.
  How does the number of classes grow with n?

  |H-spectrum|: the number of distinct H values.
  This is related to the MAXIMUM H minus minimum H divided by 2, minus gaps.
""")

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    total = 1 << m
    H_set = set()
    H_counter = Counter()
    for bits in range(total):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        H_set.add(H)
        H_counter[H] += 1

    k = len(H_set)
    max_H = max(H_set)

    print(f"\n  n={n}: m={m}")
    print(f"    2^m = {total}")
    print(f"    |H-spectrum| = {k}")
    print(f"    max H = {max_H}")
    print(f"    Average class size = {total/k:.1f}")
    print(f"    log2(|spectrum|) = {math.log2(k):.4f}")
    print(f"    log2(|spectrum|) / m = {math.log2(k)/m:.4f}")

    # Entropy of H
    entropy = sum(-p/total * math.log2(p/total) for p in H_counter.values() if p > 0)
    print(f"    H(H) = {entropy:.4f} bits")
    print(f"    H(H)/m = {entropy/m:.4f}")

# ======================================================================
# IDEA 7: ITERATED ARC FLIPS AND MIXING TIME
# ======================================================================
print("\n" + "=" * 70)
print("IDEA 7: RANDOM WALK MIXING")
print("=" * 70)

print("""
  Random walk on Q_m (flip a random arc each step).
  Starting from transitive tournament (H=1),
  how quickly does H approach its stationary mean?

  Since the walk is on Q_m (symmetric), the stationary distribution is uniform.
  E[H] at stationarity = mean H = n!/2^{n-1}.

  But H might mix FASTER or SLOWER than the overall position.
  The mixing time for H (the FUNCTION) is determined by the
  eigenvalues of Q_m that appear in H's Walsh decomposition.

  For H with purely degree-2 Walsh components (n=3,4):
  The relevant eigenvalue is m - 2*2 = m - 4.
  Mixing time ~ 1/(1 - (m-4)/m) = m/4.

  For n=3: mixing time ~ 3/4. Very fast!
  For n=4: mixing time ~ 6/4 = 1.5. Still fast.
  For n=5: m=10, dominant eigenvalue = 10-4=6, so ρ = 6/10 = 0.6.
  Mixing time ~ 1/(-ln(0.6)) ≈ 1.96.
""")

for n in [3, 4, 5]:
    m = n * (n-1) // 2

    # Simulate random walk
    import random
    random.seed(42)

    num_walks = 1000
    steps = 20
    H_by_step = [[] for _ in range(steps + 1)]

    for walk in range(num_walks):
        # Start from transitive tournament (bits = 0)
        bits = 0
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        H_by_step[0].append(H)

        for step in range(1, steps + 1):
            arc = random.randint(0, m - 1)
            bits ^= (1 << arc)
            A = adj_matrix(n, bits)
            H = count_hp(n, A)
            H_by_step[step].append(H)

    # E[H] at each step
    print(f"\n  n={n}: Random walk from transitive tournament (H=1)")
    mean_H_stationary = math.factorial(n) / 2**(n-1)
    print(f"    Stationary E[H] = {mean_H_stationary:.4f}")
    for step in range(min(steps + 1, 10)):
        mean_H = sum(H_by_step[step]) / num_walks
        print(f"    Step {step:2d}: E[H] = {mean_H:.4f} "
              f"(gap = {abs(mean_H - mean_H_stationary):.4f})")

# ======================================================================
# IDEA 8: ISOMORPHISM CLASS STRUCTURE
# ======================================================================
print("\n" + "=" * 70)
print("IDEA 8: ISOMORPHISM CLASSES AND H MULTIPLICITY")
print("=" * 70)

for n in [3, 4, 5, 6]:
    m = n * (n-1) // 2
    total = 1 << m

    # Group by (H, score_seq)
    classes = defaultdict(int)
    for bits in range(total):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        ss = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        classes[(H, ss)] += 1

    print(f"\n  n={n}: {len(classes)} (H, score-seq) classes")

    # Orbit sizes
    orbit_sizes = Counter(classes.values())
    print(f"    Orbit size distribution:")
    for size, count in sorted(orbit_sizes.items()):
        # Orbit size should divide n!
        aut_size = math.factorial(n) // size
        print(f"      Size {size:5d} (|Aut|={aut_size:3d}): {count} classes")

    if n == 6:
        break

print("\n" + "=" * 70)
print("DONE — WILD IDEAS FINAL BATCH")
print("=" * 70)
