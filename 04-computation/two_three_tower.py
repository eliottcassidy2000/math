"""
two_three_tower.py -- kind-pasteur-2026-03-13-S62

The 2-adic tower and the 2-3 duality.

Now that we understand:
  H = I(Omega, 2) where Omega uses directed cycles as vertices
  2 = number of orientations per odd cycle with same parity

Explore:
1. The 2-adic tower: I(Omega, 2^k) for k = 1, 2, 3
2. Why the factor 2 in Claim A: H(T) - H(T\v) = 2 * sum mu(C)
3. The Worpitzky decomposition and x=2 vs x=3
4. H mod 2*3 = H mod 6 and the CRT structure
5. The 48 = 2 * 24 in tr(Sigma^3) gap
6. The "2 orientations" interpretation at higher n
7. What I(Omega, 3) counts (3 "phantom orientations"?)
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

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

# ============================================================
# PART 1: The derivative of I(Omega, x) at x=2
# ============================================================
print("=" * 70)
print("PART 1: I(Omega, x) derivative structure")
print("=" * 70)

print("""
I(Omega, x) = sum alpha_k * x^k

At x=2: H = I(Omega, 2) = sum alpha_k * 2^k
At x=3: I(Omega, 3) = sum alpha_k * 3^k

Difference: I(3) - I(2) = sum alpha_k * (3^k - 2^k)
  = alpha_1 * 1 + alpha_2 * 5 + alpha_3 * 19 + ...

Derivative: I'(Omega, 2) = sum k * alpha_k * 2^{k-1}
  = alpha_1 + 2*alpha_2*2 + 3*alpha_3*4 + ...
  = alpha_1 + 4*alpha_2 + 12*alpha_3

At x=2: dH/dx = I'(Omega, 2)
This measures how sensitive H is to the evaluation point.

Logarithmic derivative: I'(Omega,2)/I(Omega,2) = d/dx[log I(Omega,x)] at x=2
This is the "average number of cycles" in a random independent set
weighted by x^|S| / I(Omega, x).
""")

# ============================================================
# PART 2: The Worpitzky transform: H = sum W_k * C(x-1, k) at x=1
# ============================================================
print("=" * 70)
print("PART 2: Worpitzky decomposition F(T,x) and the role of 2,3")
print("=" * 70)

for n in [4, 5, 6]:
    total_bits = n*(n-1)//2

    print(f"\n--- n={n} ---")
    print(f"F(T,x) = sum F_k x^k where F_k = #Ham paths with k ascents")
    print(f"Worpitzky: F(T,x) = sum W_k * C(x-1,k)")
    print(f"  W_0 = F(T,1) = H(T)")
    print(f"  W_k = sum (-1)^{'{j}'} C(k,j) F(T, k-j+1) = Euler transform")

    n_samples = min(2**total_bits, 200)

    w_stats = defaultdict(list)

    for trial in range(n_samples):
        if trial < 2**total_bits:
            bits = trial
        else:
            np.random.seed(42 + trial)
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)

        # Compute F_k
        F_coeffs = [0] * n
        for perm in permutations(range(n)):
            valid = True
            asc = 0
            for i in range(n-1):
                if not A[perm[i]][perm[i+1]]:
                    valid = False
                    break
                if perm[i] < perm[i+1]:
                    asc += 1
            if valid:
                F_coeffs[asc] += 1

        H = sum(F_coeffs)

        # Compute F(T, x) at x = 1, 2, 3
        FT = lambda x: sum(c * x**k for k, c in enumerate(F_coeffs))

        # Worpitzky coefficients via difference operator
        # W_0 = F(T, 1) = H
        # W_1 = F(T, 2) - F(T, 1)
        # W_2 = F(T, 3) - 2*F(T, 2) + F(T, 1)
        W = [0] * n
        from math import comb as mcomb
        for k in range(n):
            W[k] = sum((-1)**j * mcomb(k, j) * FT(k - j + 1)
                       for j in range(k+1))

        for k in range(n):
            w_stats[k].append(W[k])

    print(f"\n  Worpitzky coefficient W_k statistics:")
    for k in range(n):
        vals = w_stats[k]
        print(f"    W_{k}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}, "
              f"always_even={all(v%2==0 for v in vals)}, "
              f"always_div3={all(v%3==0 for v in vals)}")

# ============================================================
# PART 3: The 48 = 2 * 24 = 2 * 4! in tr(Sigma^3) gap
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Why 48 = 2 * 4! in the tr(Sigma^3) gap?")
print("=" * 70)

print("""
At n=7, ambiguous simplex profiles (same (sigma,lambda,delta) spectrum
but different c7) always differ by |delta(c7)| <= 3.

The tr(Sigma^3) values differ by EXACTLY 48 = 2 * 24 = 2 * 4!
between tournaments in the same ambiguous group.

Why 48?
  n=7, n-2=5
  4! = 24 = (n-3)!
  48 = 2 * (n-3)!

The factor 2 is the "OCF weight" (2 orientations per cycle).
The factor (n-3)! = 4! = 24 counts... what?

At n=7, a Vitali atom reversal changes 4 arcs.
The 4 arcs involve 4 vertices. 4! is the number of permutations
of these 4 vertices. But the actual number of arc flips is 4.

Alternative: 48 = 2 * C(5,2) * C(3,1) * ... hmm.
Let me compute this more carefully.
""")

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

# Find ambiguous simplex profiles and measure tr(Sigma^3) gaps
from collections import defaultdict

simplex_profiles = defaultdict(list)

for trial in range(3000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    # Build sigma matrix
    sigma = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            s = 0
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][w] and A[v][w]) or (A[w][u] and A[w][v]):
                    s += 1
            sigma[u][v] = sigma[v][u] = s

    # Simplex profile = sorted upper triangle of sigma
    profile = tuple(sorted(sigma[u][v] for u in range(n) for v in range(u+1, n)))

    tr2 = int(np.trace(sigma @ sigma))
    tr3 = int(np.trace(sigma @ sigma @ sigma))
    tr4 = int(np.trace(np.linalg.matrix_power(sigma, 4)))

    c7 = int(np.trace(np.linalg.matrix_power(A, 7))) // 7

    simplex_profiles[profile].append({
        'c7': c7, 'tr2': tr2, 'tr3': tr3, 'tr4': tr4, 'bits': bits
    })

# Find ambiguous profiles
ambig = [(p, entries) for p, entries in simplex_profiles.items()
         if len(set(e['c7'] for e in entries)) > 1]

print(f"\n  Ambiguous simplex profiles: {len(ambig)}")

for profile, entries in ambig[:10]:
    c7_vals = sorted(set(e['c7'] for e in entries))
    tr3_vals = sorted(set(e['tr3'] for e in entries))
    tr2_val = entries[0]['tr2']
    tr4_val = entries[0]['tr4']

    print(f"\n    Profile (first 5): {profile[:5]}...")
    print(f"      c7 values: {c7_vals}, gap = {max(c7_vals)-min(c7_vals)}")
    print(f"      tr(S^2) = {tr2_val} (constant)")
    print(f"      tr(S^3) values: {tr3_vals}, gap = {max(tr3_vals)-min(tr3_vals)}")
    print(f"      tr(S^4) = {tr4_val}")

    # Check if tr3 gap is exactly 48
    if len(tr3_vals) == 2:
        gap = tr3_vals[1] - tr3_vals[0]
        print(f"      tr3 gap = {gap} = 2 * {gap//2}")
        if gap == 48:
            print(f"      = 2 * 24 = 2 * 4! = 2 * (n-3)!  YES!")

# ============================================================
# PART 4: Check tr3 gap at other n values
# ============================================================
print("\n" + "=" * 70)
print("PART 4: tr(Sigma^3) gap at different n")
print("=" * 70)

for n in [5, 6, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    simplex_profiles = defaultdict(list)

    n_samples = min(2**total_bits, 5000) if n <= 6 else 3000

    for trial in range(n_samples):
        if n <= 5 and trial < 2**total_bits:
            bits = trial
        elif n == 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)

        sigma = np.zeros((n, n), dtype=int)
        for u in range(n):
            for v in range(u+1, n):
                s = 0
                for w in range(n):
                    if w == u or w == v: continue
                    if (A[u][w] and A[v][w]) or (A[w][u] and A[w][v]):
                        s += 1
                sigma[u][v] = sigma[v][u] = s

        profile = tuple(sorted(sigma[u][v] for u in range(n) for v in range(u+1, n)))

        tr3 = int(np.trace(sigma @ sigma @ sigma))

        ck = int(np.trace(np.linalg.matrix_power(A, n))) // n if n % 2 == 1 else 0

        simplex_profiles[profile].append({
            'ck': ck, 'tr3': tr3
        })

    ambig = [(p, entries) for p, entries in simplex_profiles.items()
             if len(set(e['ck'] for e in entries)) > 1]

    print(f"\n  n={n}: {len(ambig)} ambiguous profiles")

    tr3_gaps = []
    for profile, entries in ambig:
        ck_vals = sorted(set(e['ck'] for e in entries))
        tr3_vals = sorted(set(e['tr3'] for e in entries))

        if len(tr3_vals) == 2:
            gap = tr3_vals[1] - tr3_vals[0]
            tr3_gaps.append(gap)

    if tr3_gaps:
        print(f"    tr3 gap values: {Counter(tr3_gaps)}")
        for gap in sorted(set(tr3_gaps)):
            import math
            print(f"      gap={gap}: {gap}=2*{gap//2}, "
                  f"(n-3)!={(math.factorial(n-3))}, "
                  f"2*(n-3)!={2*math.factorial(n-3)}, "
                  f"match={gap==2*math.factorial(n-3)}")

# ============================================================
# PART 5: The 2-adic tower I(Omega, 2^k)
# ============================================================
print("\n" + "=" * 70)
print("PART 5: 2-adic tower I(Omega, 2^k)")
print("=" * 70)

print("""
I(Omega, 2) = H (the OCF — Hamiltonian path count)
I(Omega, 4) = sum alpha_k * 4^k = sum alpha_k * (2^2)^k
I(Omega, 8) = sum alpha_k * 8^k

If x = 2^m:
  I(Omega, 2^m) = sum alpha_k * 2^{mk}

This is a "2-adic tower" where each level refines the previous.

I(Omega, 2^m) mod 2^m gives information about alpha_0 = 1.
I(Omega, 2^m) mod 2^{2m} gives information about alpha_0 + alpha_1*2^m.
etc.

The point x=2 is the "first level" of this tower.
""")

n = 5
total_bits = n*(n-1)//2

print(f"\n  n={n}:")
print(f"  {'bits':>5} {'H':>5} {'I(O,4)':>8} {'I(O,8)':>8} {'I(O,2)/I(O,4)':>14}")

for bits in [0, 10, 40, 100, 200, 500, 1023]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Count odd cycles with multiplicity
    c3_list = []
    c5 = 0

    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                sub = np.zeros((3,3), dtype=int)
                verts = [i,j,k]
                for a in range(3):
                    for b in range(3):
                        sub[a][b] = A[verts[a]][verts[b]]
                c = int(np.trace(np.linalg.matrix_power(sub, 3))) // 3
                if c > 0:
                    c3_list.append(frozenset(verts))

    # 5-cycle count (directed, with multiplicity)
    sub = A.copy()
    c5 = int(np.trace(np.linalg.matrix_power(sub, 5))) // 5

    # At n=5, all cycles share vertices (5 vertices only)
    # So Omega has alpha_1 = len(c3_list) + c5, alpha_k = 0 for k >= 2
    alpha_1 = len(c3_list) + c5

    io2 = 1 + 2 * alpha_1
    io4 = 1 + 4 * alpha_1
    io8 = 1 + 8 * alpha_1

    print(f"  {bits:>5} {H:>5} {io4:>8} {io8:>8} {'N/A' if io4 == 0 else f'{io2/io4:.4f}':>14}")

# ============================================================
# PART 6: H mod 3 and c3 mod 3 correlation
# ============================================================
print("\n" + "=" * 70)
print("PART 6: H mod 3 vs c3 mod 3 — is there a pattern?")
print("=" * 70)

for n in [5, 6, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    cross_tab = Counter()

    n_samples = min(2**total_bits, 5000)

    for trial in range(n_samples):
        if n <= 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3

        cross_tab[(H % 3, c3 % 3)] += 1

    print(f"\n  n={n}: H mod 3 vs c3 mod 3 crosstab:")
    print(f"  {'':>12} c3=0(3)  c3=1(3)  c3=2(3)")
    for h_r in [0, 1, 2]:
        row = [cross_tab.get((h_r, c_r), 0) for c_r in [0, 1, 2]]
        print(f"    H={h_r}(3):  {row[0]:>7}  {row[1]:>7}  {row[2]:>7}")

    # Check if H = 1 + 2*c3 mod 3 (i.e. H mod 3 = (1 + 2*c3) mod 3)
    # Since H = 1 + 2*alpha_1 + 4*alpha_2 + ...
    # mod 3: H = 1 + 2*alpha_1 + alpha_2 + 2*alpha_3 + alpha_4 + ...
    # The alpha_k involve c5, c7 too, not just c3

# ============================================================
# PART 7: The "3-cycle number" and "2-orientation" coupling
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Coupling between 2 and 3 in H(T)")
print("=" * 70)

print("""
H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

KEY IDENTITY: H = 1 (mod 2) always (Redei)
This follows from I(G, 2) = I(G, 0) = 1 (mod 2) for ANY graph G.

DEEPER: H - 1 = 2*alpha_1 + 4*alpha_2 + ...
So (H-1)/2 = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...
           = I(Omega, 2) evaluated at the "shifted" polynomial!

Wait: (H-1)/2 = (I(Omega,2) - 1)/2 = sum_{k>=1} alpha_k * 2^{k-1}

This is the number of "odd cycle choices with halved weight."

And (H-1)/2 mod 2:
  = alpha_1 mod 2  (since higher terms are divisible by 2)
  = (# directed odd cycles in T) mod 2

So the second bit of H is determined by the PARITY of the number
of directed odd cycles!

At n=5: odd cycles = 3-cycles + 5-cycles
  alpha_1 = c3 + c5
  (H-1)/2 = alpha_1 (since alpha_2 = 0 at n=5)
""")

n = 5
total_bits = n*(n-1)//2

# Verify (H-1)/2 = alpha_1 at n=5
print(f"\n  Verifying (H-1)/2 mod 2 = alpha_1 mod 2 at n={n}:")
matches = 0
for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    alpha_1 = c3 + c5

    if ((H-1)//2) % 2 == alpha_1 % 2:
        matches += 1

print(f"    Matches: {matches}/{2**total_bits}")

# At n=7, alpha_1 = c3 + c5_sets + c7_directed_cycles
# More complex but same principle
n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

print(f"\n  Verifying (H-1)/2 mod 2 = alpha_1 mod 2 at n={n}:")
matches = 0
n_samples = 1000
for trial in range(n_samples):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # alpha_1 = total number of directed odd cycles (with multiplicity)
    alpha_1 = 0
    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            sub = np.zeros((size, size), dtype=int)
            verts = list(combo)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            c = int(np.trace(np.linalg.matrix_power(sub, size))) // size
            alpha_1 += c

    if ((H-1)//2) % 2 == alpha_1 % 2:
        matches += 1

print(f"    Matches: {matches}/{n_samples}")

# ============================================================
# PART 8: What does H mod 3 correspond to?
# ============================================================
print("\n" + "=" * 70)
print("PART 8: H mod 3 = I(Omega, 2) mod 3")
print("=" * 70)

print("""
H mod 3 = I(Omega, 2) mod 3
  = (1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...) mod 3
  = (1 + 2*alpha_1 + alpha_2 + 2*alpha_3 + alpha_4 + ...) mod 3
  = (1 + 2*(alpha_1 + alpha_3 + ...) + (alpha_2 + alpha_4 + ...)) mod 3

Since 2 = -1 mod 3 and 4 = 1 mod 3:
  H mod 3 = 1 - alpha_1 + alpha_2 - alpha_3 + alpha_4 - ... mod 3
          = I(Omega, -1) mod 3 !!

REMARKABLE: H mod 3 = I(Omega, -1) mod 3.

I(Omega, -1) is the evaluation at x=-1, which counts independent sets
with alternating signs: sum (-1)^k * alpha_k.

For bipartite graphs, I(G, -1) = 0 always. For triangle-free graphs,
I(G, -1) = 1 - alpha_1 + alpha_2.
""")

# Verify: H mod 3 = I(Omega, -1) mod 3 at n=5
n = 5
total_bits = n*(n-1)//2

print(f"\n  Verifying H mod 3 = I(Omega, -1) mod 3 at n={n}:")
matches = 0
for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # At n=5: alpha_1 = c3 + c5, alpha_2 = 0, alpha_k = 0 for k >= 2
    c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    alpha_1 = c3 + c5

    I_neg1 = 1 - alpha_1  # since alpha_k = 0 for k >= 2

    if H % 3 == I_neg1 % 3:
        matches += 1

print(f"    Matches: {matches}/{2**total_bits}")

# Verify at n=6 (need full alpha computation)
n = 6
total_bits = n*(n-1)//2
np.random.seed(42)

print(f"\n  Verifying H mod 3 = I(Omega, -1) mod 3 at n={n} (sampled):")
matches = 0
n_samples = 200
for trial in range(n_samples):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Compute alpha_1 and alpha_2 (sufficient for I(Omega,-1) mod 3 at small n)
    # alpha_1 = total directed odd cycles with multiplicity
    alpha_1 = 0
    cycle_sets = []
    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            sub = np.zeros((size, size), dtype=int)
            verts = list(combo)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            c = int(np.trace(np.linalg.matrix_power(sub, size))) // size
            if c > 0:
                alpha_1 += c
                for _ in range(c):
                    cycle_sets.append(frozenset(combo))

    # alpha_2 = vertex-disjoint pairs
    alpha_2 = 0
    for a in range(len(cycle_sets)):
        for b in range(a+1, len(cycle_sets)):
            if not (cycle_sets[a] & cycle_sets[b]):
                alpha_2 += 1

    I_neg1 = 1 - alpha_1 + alpha_2
    if H % 3 == I_neg1 % 3:
        matches += 1

print(f"    Matches: {matches}/{n_samples}")

print("\n\nDone.")
