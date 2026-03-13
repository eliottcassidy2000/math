"""
two_three_signs.py -- kind-pasteur-2026-03-13-S62

INVESTIGATION: b0 <= 0 and b1 > 0 — are these universal?

b0 = I(Omega, -1) = 1 - a1 + a2
b1 = a1 - 2*a2

b0 <= 0  <=>  a1 >= a2 + 1  <=>  a1 - a2 >= 1
b1 > 0   <=>  a1 > 2*a2

Combined: a1 > 2*a2 >= 2*0 = 0, so a1 >= 1 (at least one directed odd cycle).
Also: a1 > 2*a2 and a1 >= a2 + 1.

If a2 = 0: b0 = 1 - a1 <= 0 iff a1 >= 1. b1 = a1 > 0 iff a1 >= 1.
So both conditions reduce to: every tournament has at least 1 directed 3-cycle
(for n >= 3). This is FALSE: transitive tournaments have c3 = 0.

Wait: at n=3, the transitive tournament has a1=0, H=1, b0=1, b1=0.
So b0 > 0 and b1 = 0 at n=3 for transitive!

And at n=5, bits=0 is transitive: a1=0, H=1, b0=1, b1=0.
So b0 = 1 > 0 for transitive tournaments.

But at n=7, all 400 samples had b0 <= 0. Let me check: can n=7
have a1 = 0? A transitive tournament on 7 vertices has 0 directed
3-cycles. So a1 = 0, a2 = 0, H = 1, b0 = 1 > 0!

I must have gotten unlucky: random tournaments on 7 vertices
almost never have a1 = 0. Let me check with specific examples.
"""

import numpy as np
from itertools import combinations
from collections import Counter

def bits_to_adj(bits, n):
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

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_ham_cycles_exact(A_sub, k):
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask_size in range(2, k+1):
        for mask in range(1 << k):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):
                continue
            for v in range(k):
                if not (mask & (1 << v)):
                    continue
                if v == 0:
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(k):
                    if (prev_mask & (1 << u)) and A_sub[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << k) - 1
    total_cycles = 0
    for v in range(1, k):
        if A_sub[v][0] and dp.get((full, v), 0):
            total_cycles += dp[(full, v)]
    return total_cycles

def get_alpha_1_2(A, n):
    cycles = []
    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            verts = list(combo)
            sub = np.zeros((size, size), dtype=int)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            if size <= 5:
                c = int(np.trace(np.linalg.matrix_power(sub, size))) // size
            else:
                c = count_ham_cycles_exact(sub, size)
            for _ in range(c):
                cycles.append(frozenset(combo))
    alpha_1 = len(cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1
    return alpha_1, alpha_2

# ============================================================
# PART 1: Check transitive tournament
# ============================================================
print("=" * 70)
print("PART 1: Transitive tournament (a1=0)")
print("=" * 70)

for n in [3, 5, 7]:
    # Transitive: i beats j iff i > j. bits = 0 gives j beats i for i<j.
    # Actually bits=0 means for each pair (i,j) with i<j, bit=0 means A[j][i]=1.
    # So j beats i for all j > i. Score seq = (n-1, n-2, ..., 1, 0). Transitive.
    A = bits_to_adj(0, n)
    H = count_ham_paths(A, n)
    c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3
    a1, a2 = get_alpha_1_2(A, n)
    b0 = 1 - a1 + a2
    b1 = a1 - 2*a2

    print(f"\n  n={n} transitive: H={H}, c3={c3}, a1={a1}, a2={a2}")
    print(f"    b0={b0}, b1={b1}")
    print(f"    b0 > 0: {b0 > 0}")

# ============================================================
# PART 2: Exhaustive at n=5 — check b0 sign
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Exhaustive b0 sign at n=5")
print("=" * 70)

n = 5
total_bits = n*(n-1)//2

b0_neg = 0
b0_zero = 0
b0_pos = 0
b0_examples = {'neg': [], 'zero': [], 'pos': []}

for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    a1 = c3 + c5
    a2 = 0  # always 0 at n=5

    b0 = 1 - a1 + a2

    if b0 < 0:
        b0_neg += 1
        if len(b0_examples['neg']) < 3:
            b0_examples['neg'].append((bits, a1, b0))
    elif b0 == 0:
        b0_zero += 1
        if len(b0_examples['zero']) < 3:
            b0_examples['zero'].append((bits, a1, b0))
    else:
        b0_pos += 1
        if len(b0_examples['pos']) < 3:
            b0_examples['pos'].append((bits, a1, b0))

print(f"  n={n}: b0 < 0: {b0_neg}, b0 = 0: {b0_zero}, b0 > 0: {b0_pos}")
print(f"    total = {b0_neg + b0_zero + b0_pos}")

for sign_name in ['neg', 'zero', 'pos']:
    for bits, a1, b0 in b0_examples[sign_name]:
        print(f"    b0={b0}: bits={bits}, a1={a1}")

# ============================================================
# PART 3: n=7 — find b0 > 0 cases
# ============================================================
print("\n" + "=" * 70)
print("PART 3: n=7 — searching for b0 > 0")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2

# Try transitive first
A = bits_to_adj(0, n)
a1, a2 = get_alpha_1_2(A, n)
b0 = 1 - a1 + a2
print(f"  Transitive (bits=0): a1={a1}, a2={a2}, b0={b0}")

# Also try the all-1 bits (reverse transitive)
A = bits_to_adj((1 << total_bits) - 1, n)
a1, a2 = get_alpha_1_2(A, n)
b0 = 1 - a1 + a2
print(f"  Reverse transitive (bits=max): a1={a1}, a2={a2}, b0={b0}")

# Try near-transitive
for shift in [1, 2, 4, 8, 16]:
    A = bits_to_adj(shift, n)
    a1, a2 = get_alpha_1_2(A, n)
    b0 = 1 - a1 + a2
    H = count_ham_paths(A, n)
    print(f"  bits={shift}: a1={a1}, a2={a2}, b0={b0}, H={H}")

# Random search for b0 > 0
np.random.seed(42)
b0_pos_found = 0
for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)
    b0 = 1 - a1 + a2
    if b0 > 0:
        b0_pos_found += 1
        H = count_ham_paths(A, n)
        print(f"  FOUND b0 > 0: bits={bits}, a1={a1}, a2={a2}, b0={b0}, H={H}")

if b0_pos_found == 0:
    print(f"  No b0 > 0 found in 2000 random samples")

# ============================================================
# PART 4: Is b0 <= 0 equivalent to a1 >= 1?
# ============================================================
print("\n" + "=" * 70)
print("PART 4: b0 sign at n=5 vs a1")
print("=" * 70)

n = 5
total_bits = n*(n-1)//2

print(f"  b0 = 1 - a1 (at n=5, since a2=0)")
print(f"  b0 <= 0 iff a1 >= 1 iff tournament has at least one directed 3-cycle")
print(f"  b0 > 0 iff a1 = 0 iff tournament is transitive")

# Count transitive tournaments at n=5
trans = sum(1 for bits in range(2**total_bits)
            if int(np.trace(np.linalg.matrix_power(bits_to_adj(bits, n), 3))) == 0)
print(f"  # transitive n=5: {trans}/{2**total_bits}")

# At n=7: transitive has a1=0. Are there non-transitive with a1=0?
# No: a1=0 means NO directed odd cycles at all. For a tournament this means
# no 3-cycles (acyclic). Acyclic tournament = transitive (unique ordering).
# At n=7, there's exactly 1 transitive tournament per vertex ordering, so
# 7!/(number of automorphisms) = ... well, there are n!/n = (n-1)! = 720
# distinct transitive tournaments up to labeling... no.
# Actually for labeled tournaments: bits=0 is one, and any permutation
# gives a different bits value. But the representation depends on the
# encoding in bits_to_adj.

print(f"\n  Transitive iff a1=0 iff acyclic (for tournaments)")
print(f"  At n=7: transitive tournament exists, has b0=1 > 0")
print(f"  Random sampling missed it because P(transitive) = 7!/2^21 = 5040/2097152 ~ 0.0024")

# Verify: how many labeled transitive tournaments at n=7?
# A transitive tournament = a total order = a permutation of vertices
# There are n! = 5040 total orders, each giving a unique tournament
# But wait: different total orders can give the same tournament only if
# the tournament has non-trivial automorphism group.
# For transitive tournament: Aut = trivial (identity only, since the
# score sequence is (0,1,2,...,n-1), all different).
# So # labeled transitive tournaments = n! = 5040... no that's wrong.
# Each total order sigma gives a tournament where sigma(i) beats sigma(j)
# iff i > j. Different orders give different tournaments? Not quite:
# The tournament is DETERMINED by the total order (vertex ordering).
# Two orderings give the same tournament iff they are the same ordering.
# So there are exactly n! labeled transitive tournaments.
# Wait, that's not right either. A total order on n vertices gives a
# unique tournament. And conversely, an acyclic tournament gives a unique
# total order (topological sort). So there's a bijection between total
# orders and transitive tournaments.
# # total orders on n labeled vertices = n! (permutations).
# So # labeled transitive tournaments on n vertices = n!.

# But in our encoding, bits represents one specific labeling.
# Hmm, actually each total order gives a unique SET of arc directions.
# A total order v_{sigma(1)} < v_{sigma(2)} < ... < v_{sigma(n)} gives
# the tournament where v_{sigma(i)} beats v_{sigma(j)} iff i > j.
# Different permutations sigma give different tournaments.
# So there are n! transitive tournaments on labeled vertices.

# But in our bits encoding, each tournament has a unique bits value.
# How many bits values correspond to transitive tournaments?
n = 7
total_bits = n*(n-1)//2

# To find: # transitive tournaments = # bits values with c3 = 0
# This is expensive for n=7 (2^21 = 2M), let me estimate
from math import factorial
print(f"\n  Expected # transitive tournaments at n={n}:")
print(f"    n! = {factorial(n)}")
print(f"    2^C(n,2) = {2**total_bits}")
print(f"    P(transitive) = {factorial(n)}/{2**total_bits} = {factorial(n)/2**total_bits:.6f}")

# ============================================================
# PART 5: b1 > 0 — when can it fail?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Can b1 <= 0?")
print("=" * 70)

print(f"  b1 = a1 - 2*a2")
print(f"  b1 <= 0 iff a1 <= 2*a2")
print(f"  This requires: # disjoint pairs >= (# total cycles)/2")
print(f"  At n=5: a2=0, so b1 = a1 >= 0. b1=0 only for transitive.")
print(f"  At n=7: need a1 <= 2*a2. Since a2 <= C(alpha_1 choose 2)...")

# At n=7, max a2 is achieved when cycles are spread out.
# Two disjoint 3-cycles need 6 vertices, leaving 1.
# With 7 vertices, max disjoint 3-cycle pairs from C(7,3)=35 vertex sets:
# Each disjoint pair uses 6 vertices, and for each pair of disjoint 3-vertex
# sets, we need both to be cyclic.

# When does a1 <= 2*a2?
# If all cycles are 3-cycles (c5=c7=0), and there are many disjoint pairs:
# a1 = c3 (total directed 3-cycles)
# a2 = #{disjoint pairs of directed 3-cycles}

# For n=7 regular tournament (Paley):
# c3 = 14 (7 * 2 / 1 = 14? no)
# Actually c3 = C(7,3) * (# 3-cycles per triple) / ... Let's compute
# For Paley T_7: S = {1,2,4}, each triple has probability determined by QR structure

# Let me just check n=9 where alpha_3 > 0 is possible
# (three disjoint 3-cycles need 9 vertices)

print(f"\n  Checking b1 sign at n=9 (small sample):")
n = 9
total_bits = n*(n-1)//2
np.random.seed(42)

b1_neg_count = 0
for trial in range(50):
    bits = int(np.random.randint(0, 2**31)) | (int(np.random.randint(0, 1 << (total_bits - 31))) << 31) if total_bits > 31 else np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    # Only count 3-cycles and 5-cycles (7 and 9 too expensive for full computation)
    c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3

    # For a rough estimate: a1 ~ c3 + c5, a2 ~ disjoint 3-cycle pairs
    c3_cycles = []
    for combo in combinations(range(n), 3):
        sub = np.zeros((3, 3), dtype=int)
        verts = list(combo)
        for a in range(3):
            for b in range(3):
                sub[a][b] = A[verts[a]][verts[b]]
        c = int(np.trace(sub @ sub @ sub)) // 3
        if c > 0:
            c3_cycles.append(frozenset(combo))

    a1_c3 = len(c3_cycles)
    a2_c3 = 0
    for i in range(len(c3_cycles)):
        for j in range(i+1, len(c3_cycles)):
            if not (c3_cycles[i] & c3_cycles[j]):
                a2_c3 += 1

    b1_est = a1_c3 - 2*a2_c3
    if b1_est <= 0:
        b1_neg_count += 1
        print(f"    trial {trial}: a1_c3={a1_c3}, a2_c3={a2_c3}, b1_est={b1_est}")

if b1_neg_count == 0:
    print(f"    b1 (3-cycle only estimate) always > 0 in 50 samples")

# ============================================================
# PART 6: The inequality a1 > 2*a2 — a graph theory statement
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Why a1 > 2*a2? A graph theory perspective")
print("=" * 70)

print("""
Claim: For the conflict graph Omega(T) of any non-transitive tournament T,
  alpha_1 > 2*alpha_2

(alpha_1 = # vertices of Omega, alpha_2 = # independent edges of Omega-bar)

Actually, alpha_2 = # pairs of NON-adjacent vertices in Omega
(= # edges of Omega-bar-complement... no).

Wait: alpha_2 = # independent sets of size 2 in Omega
     = # pairs of NON-adjacent vertices in Omega
     = # edges of the COMPLEMENT of Omega
     = C(alpha_1, 2) - (# edges of Omega)

So alpha_2 = C(a1, 2) - e(Omega)
where e(Omega) = # edges of Omega (conflict edges).

Then: b1 = a1 - 2*alpha_2 = a1 - 2*(C(a1,2) - e(Omega))
     = a1 - a1*(a1-1) + 2*e(Omega)
     = a1 - a1^2 + a1 + 2*e(Omega)
     = 2*a1 - a1^2 + 2*e(Omega)
     = -a1*(a1-2) + 2*e(Omega)

b1 > 0 iff 2*e(Omega) > a1*(a1-2)
       iff e(Omega) > a1*(a1-2)/2

The max possible edges is C(a1,2) = a1*(a1-1)/2.
So the condition is: e(Omega) > a1*(a1-2)/2 = C(a1,2) - a1/2 + 1/2

This means: the complement has FEWER than a1/2 edges.
In other words: Omega is "almost complete" — its complement is sparse.

For tournament conflict graphs at n=7:
  alpha_1 ~ 30, e(Omega) ~ C(30,2) - alpha_2 ~ 435 - 5 = 430
  That's e/C(a1,2) ~ 430/435 ~ 0.99 — extremely dense!

The conflict graph Omega is ALMOST complete because at n=7,
most pairs of directed odd cycles share at least one vertex
(since all cycles live in a 7-vertex universe).
""")

# Verify density of Omega at n=7
n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

print(f"  Omega density at n=7:")
for trial in range(10):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)

    max_edges = a1*(a1-1)//2 if a1 >= 2 else 0
    e_omega = max_edges - a2
    density = e_omega / max_edges if max_edges > 0 else 0

    print(f"    a1={a1:>3}, a2={a2:>3}, e(Omega)={e_omega:>5}, "
          f"C(a1,2)={max_edges:>5}, density={density:.4f}, "
          f"b1={a1-2*a2:>3}")

# ============================================================
# PART 7: Threshold for b0 <= 0
# ============================================================
print("\n" + "=" * 70)
print("PART 7: At what n does b0 <= 0 become 'almost universal'?")
print("=" * 70)

for n in [3, 4, 5, 6, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    n_samples = min(2**total_bits, 2000)

    b0_neg = 0
    b0_zero = 0
    b0_pos = 0

    for trial in range(n_samples):
        if trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)

        # alpha_1 only (fast)
        a1 = 0
        for size in range(3, n+1, 2):
            for combo in combinations(range(n), size):
                sub = np.zeros((size, size), dtype=int)
                verts = list(combo)
                for a_idx in range(size):
                    for b_idx in range(size):
                        sub[a_idx][b_idx] = A[verts[a_idx]][verts[b_idx]]
                if size <= 5:
                    c = int(np.trace(np.linalg.matrix_power(sub, size))) // size
                else:
                    c = count_ham_cycles_exact(sub, size)
                a1 += c

        # At these n values, a2 is small relative to a1
        # b0 = 1 - a1 + a2 ~ 1 - a1 (lower bound)
        # For b0 > 0 we need a1 < 1 + a2
        # Without computing a2, b0 > 0 requires a1 = 0 (transitive)

        # But actually for n <= 5, a2 = 0 (no disjoint pairs fit in n=5)
        # For n=6: can have disjoint 3-cycle pair (uses 6 vertices)
        # For n=7: same

        # For exact check at small n, let's just use a1 vs 0
        if n <= 5:
            b0 = 1 - a1
        else:
            # Need a2 too — skip for speed, just estimate
            b0 = 1 - a1  # lower bound (since a2 >= 0)

        if b0 < 0:
            b0_neg += 1
        elif b0 == 0:
            b0_zero += 1
        else:
            b0_pos += 1

    total = b0_neg + b0_zero + b0_pos
    print(f"  n={n}: b0<0: {b0_neg} ({100*b0_neg/total:.1f}%), "
          f"b0=0: {b0_zero} ({100*b0_zero/total:.1f}%), "
          f"b0>0: {b0_pos} ({100*b0_pos/total:.1f}%)")

    if n <= 5:
        # b0 > 0 iff a1 = 0 iff transitive
        print(f"    (b0 > 0 iff transitive: {b0_pos} tournaments)")
    else:
        print(f"    (b0 estimate uses lower bound 1-a1, actual b0 >= this)")

print("\n\nDone.")
