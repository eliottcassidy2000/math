"""
hmod3_identity.py -- kind-pasteur-2026-03-13-S62

KEY IDENTITY: H(T) mod 3 = I(Omega(T), -1) mod 3

Since 2 = -1 (mod 3), evaluating I(Omega, x) at x=2 gives H,
and reducing mod 3 gives I(Omega, -1) mod 3.

This means: H mod 3 = sum_{k=0}^{alpha(Omega)} (-1)^k * alpha_k (mod 3)

I(G, -1) has a name: it's the ALTERNATING INDEPENDENT SET COUNT.
For a graph G: I(G, -1) = sum (-1)^k * alpha_k

Connections:
  - For bipartite G: I(G, -1) = 0 (matched by even/odd partition)
  - For clique K_n: I(K_n, -1) = 1 - n (since alpha_k = 0 for k >= 2)
  - I(G, -1) = 0 iff chromatic number chi(G) > 1 AND G is bipartite? No...
  - Actually I(G, -1) relates to the EULER CHARACTERISTIC of Omega!

KEY: I(G, -1) = chi(Ind(G)) where Ind(G) is the INDEPENDENCE COMPLEX
     (simplicial complex of independent sets of G).
     This is the reduced Euler characteristic (plus 1).

So: H mod 3 = chi(Ind(Omega(T))) mod 3 !!

The Euler characteristic of the independence complex of the
conflict graph controls H modulo 3.
"""

import numpy as np
from itertools import combinations
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
# PART 1: Verify the identity H mod 3 = I(Omega, -1) mod 3
# ============================================================
print("=" * 70)
print("PART 1: H mod 3 = I(Omega, -1) mod 3 = chi(Ind(Omega)) mod 3")
print("=" * 70)

print("""
PROOF:
  H = I(Omega, 2) = sum alpha_k * 2^k
  2 = -1 (mod 3)
  Therefore: H mod 3 = sum alpha_k * (-1)^k (mod 3) = I(Omega, -1) mod 3
  QED.

This is trivially true for ANY graph G and ANY x with x = -1 mod 3.
But the specific interpretation is profound:

I(G, -1) = 1 - alpha_1 + alpha_2 - alpha_3 + ...
         = sum_{S: indep set} (-1)^|S|
         = chi(Ind(G)) + 1

where chi is the reduced Euler characteristic and Ind(G) is the
independence complex (simplicial complex of independent sets of G).

So: H mod 3 = (1 + chi_tilde(Ind(Omega))) mod 3

where chi_tilde = reduced Euler characteristic.

CONSEQUENCE: H = 0 mod 3 iff chi_tilde(Ind(Omega)) = -1 mod 3
             H = 1 mod 3 iff chi_tilde(Ind(Omega)) = 0 mod 3
             H = 2 mod 3 iff chi_tilde(Ind(Omega)) = 1 mod 3
""")

# ============================================================
# PART 2: Worpitzky W_k divisibility by 2 and 3
# ============================================================
print("=" * 70)
print("PART 2: Worpitzky W_k divisibility structure")
print("=" * 70)

from math import comb as mcomb, factorial

print("""
F(T, x) = sum_{k=0}^{n-1} W_k * C(x-1, k) (Worpitzky/Eulerian decomposition)

where W_k = sum_{j=0}^k (-1)^j C(k,j) F(T, k-j+1)

Key Eulerian identity: sum_T W_k(T) / 2^{C(n,2)} = A(n,k)
where A(n,k) = Eulerian numbers.

OBSERVATION: W_2 is ALWAYS EVEN, W_3 is ALWAYS DIV BY 6.

Why?

W_2 = F(T,3) - 2*F(T,2) + F(T,1)
    = F(T,3) - 2*F(T,2) + H

F(T,x) = sum F_k x^k where F_k = # paths with k ascents.

W_2 = sum_k F_k (3^k - 2*2^k + 1)
For k=0: 1 - 2 + 1 = 0
For k=1: 3 - 4 + 1 = 0
For k=2: 9 - 8 + 1 = 2
For k=3: 27 - 16 + 1 = 12
For k=4: 81 - 32 + 1 = 50
For k=5: 243 - 64 + 1 = 180

So W_2 = 2*F_2 + 12*F_3 + 50*F_4 + 180*F_5 + ...
All coefficients are even! Hence W_2 is always even.

Actually, the coefficient of F_k in W_2 is 3^k - 2*2^k + 1.
This equals (3-1)^k + 0 = 2^k for k >= 0... no.
(3^k - 2*2^k + 1) for k=0: 0, k=1: 0, k=2: 2, k=3: 12

CLAIM: 3^k - 2*2^k + 1 is always even for k >= 0.
Proof: 3^k is odd. 2*2^k is even. 1 is odd. Odd - even + odd = even. QED.

So W_2 is a sum of F_k times even coefficients, hence ALWAYS EVEN.
More: for k >= 2, 3^k - 2*2^k + 1 = 2 mod 4 when k=2,
                                     = 0 mod 4 when k >= 3.
So W_2 = 2*F_2 (mod 4).
""")

# Verify
for n in [4, 5, 6, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    w2_div = Counter()
    w3_div = Counter()

    n_samples = min(2**total_bits, 1000) if n <= 6 else 500

    for trial in range(n_samples):
        if n <= 5 and trial < 2**total_bits:
            bits = trial
        elif n == 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)

        # Compute F_k
        from itertools import permutations
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

        FT = lambda x: sum(c * x**k for k, c in enumerate(F_coeffs))

        W = [0] * n
        for k in range(n):
            W[k] = sum((-1)**j * mcomb(k, j) * FT(k - j + 1) for j in range(k+1))

        # Divisibility
        if n >= 3:
            v2_w2 = 0
            val = W[2]
            if val != 0:
                while val % 2 == 0:
                    v2_w2 += 1
                    val //= 2
            else:
                v2_w2 = 99
            w2_div[v2_w2] += 1

        if n >= 4:
            v_w3 = 0
            val = W[3]
            if val != 0:
                while val % 6 == 0:
                    v_w3 += 1
                    val //= 6
            else:
                v_w3 = 99
            w3_div[v_w3] += 1

    print(f"\n  n={n} ({n_samples} samples):")
    if n >= 3:
        print(f"    W_2 min v_2: {min(w2_div.keys())}")
        print(f"    W_2 v_2 distribution: {dict(sorted(w2_div.items()))}")
    if n >= 4:
        print(f"    W_3 v_6 distribution: {dict(sorted(w3_div.items()))}")

# ============================================================
# PART 3: W_k modular structure
# ============================================================
print("\n" + "=" * 70)
print("PART 3: W_k mod 2 and mod 3 structure")
print("=" * 70)

print("""
For the Worpitzky coefficient W_k:
  coefficient of F_j in W_k = sum_{i=0}^k (-1)^i C(k,i) (k-i+1)^j
  = S(j, k) * k!  (Stirling number times k!)

Actually: W_k = k! * S(?, k) related to Stirling numbers.

The key question: when is W_k divisible by 2? By 3? By 6?

Let me just compute the modular structure for all k.
""")

for n in [5, 6, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    mod_structure = defaultdict(lambda: Counter())

    n_samples = min(2**total_bits, 500) if n <= 6 else 200

    for trial in range(n_samples):
        if n <= 5 and trial < 2**total_bits:
            bits = trial
        elif n == 6:
            np.random.seed(42 + trial)
            bits = np.random.randint(0, 1 << total_bits)
        else:
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)

        from itertools import permutations
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

        FT = lambda x: sum(c * x**k for k, c in enumerate(F_coeffs))

        W = [0] * n
        for k in range(n):
            W[k] = sum((-1)**j * mcomb(k, j) * FT(k - j + 1) for j in range(k+1))

        for k in range(n):
            mod_structure[(n, k, 'mod2')][W[k] % 2] += 1
            mod_structure[(n, k, 'mod3')][W[k] % 3] += 1
            mod_structure[(n, k, 'mod6')][W[k] % 6] += 1

    print(f"\n  n={n}:")
    for k in range(n):
        m2 = dict(sorted(mod_structure[(n, k, 'mod2')].items()))
        m3 = dict(sorted(mod_structure[(n, k, 'mod3')].items()))
        m6 = dict(sorted(mod_structure[(n, k, 'mod6')].items()))
        print(f"    W_{k}: mod2={m2}, mod3={m3}, mod6={m6}")

# ============================================================
# PART 4: The 2-3 coupling in the OCF
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Synthesis - the fundamental 2-3 coupling")
print("=" * 70)

print("""
SYNTHESIS OF THE 2-3 DUALITY:

1. H = I(Omega, 2) — the OCF at x=2
2. H mod 2 = 1 always (Redei) — comes from I(G, 0) = 1 for any G
3. H mod 3 = I(Omega, -1) mod 3 — comes from 2 = -1 (mod 3)

So x=2 simultaneously:
  - Gives H (the full count) via direct evaluation
  - Gives H mod 2 = 1 via I(G, 0) = 1
  - Gives H mod 3 via I(Omega, -1) = Euler char of independence complex

The NUMBER 2 is both:
  (a) The fugacity of the hard-core lattice gas (weight per cycle)
  (b) The number of orientations of an odd permutation cycle
  (c) The generator of Z/3Z* (multiplicative group mod 3)

The NUMBER 3 is both:
  (a) The smallest odd cycle length
  (b) The modulus where 2 = -1 (linking OCF to topology)
  (c) The number of vertices per atom in the simplex (sigma, lambda, delta)

The PAIR (2, 3):
  - 2 and 3 are the only consecutive primes
  - 2 = -1 (mod 3): links even/odd to ternary structure
  - 3 = -1 (mod 4): links ternary to quaternary
  - 2*3 = 6: H mod 6 is determined by Redei (mod 2) + topology (mod 3)

THE DEEPEST CONNECTION:
  H(T) = I(Omega(T), 2)
  x=2 is the UNIQUE positive integer that:
    - Makes I(Omega, x) always odd (x = 0 mod 2)
    - Makes I(Omega, x) encode topology (x = -1 mod 3)
    - Is the SMALLEST such integer (next would be x=8: 0 mod 2, -1 mod 3)
    - Actually x=2 is the unique x with x mod 2 = 0 AND x mod 3 = 2
      By CRT: x = 2 mod 6. Solutions: 2, 8, 14, 20, ...
      x=2 is the SMALLEST.

  So x=2 is the unique smallest positive integer x such that:
    I(G, x) is odd for all G (Redei's theorem)
    AND I(G, x) mod 3 = chi(Ind(G)) + 1 (mod 3) (topological content)
""")

# ============================================================
# PART 5: What does I(Omega, 8) mean?
# ============================================================
print("=" * 70)
print("PART 5: I(Omega, 8) — the next level of the tower")
print("=" * 70)

print("""
x=8 = 2^3 is the NEXT element of the tower:
  8 mod 2 = 0 (still odd output)
  8 mod 3 = 2 = -1 (same topology as x=2)
  8 mod 4 = 0 (new: divisible by 4)
  8 mod 5 = 3 (new information)
  8 mod 7 = 1 (new information)

I(Omega, 8) = sum alpha_k * 8^k = sum alpha_k * 2^{3k}

This gives a "3-fold magnification" of the independence set structure.
The ratio I(Omega, 8) / I(Omega, 2) depends on the alpha_k distribution.
""")

n = 5
total_bits = n*(n-1)//2

# Compare I(Omega, 2), I(Omega, 8), and I(Omega, -1) for all n=5
print(f"\n  n={n}:")
print(f"  {'bits':>5} {'H=I(2)':>8} {'I(8)':>8} {'I(-1)':>6} "
      f"{'H%3':>4} {'I(-1)%3':>8}")

for bits in range(0, 2**total_bits, 100):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Count alpha_1 (all directed odd cycles)
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

    # At n=5, alpha_2 = 0 (all cycles share vertices)
    I_2 = 1 + 2 * alpha_1
    I_8 = 1 + 8 * alpha_1
    I_neg1 = 1 - alpha_1

    print(f"  {bits:>5} {I_2:>8} {I_8:>8} {I_neg1:>6} {H%3:>4} {I_neg1%3:>8}"
          f"  {'MATCH' if H%3 == I_neg1%3 else 'FAIL'}")

# ============================================================
# PART 6: The tr(Sigma^3) gap divisibility by 12
# ============================================================
print("\n" + "=" * 70)
print("PART 6: tr(Sigma^3) gaps always divisible by 12 = 4 * 3")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

simplex_profiles = defaultdict(list)

for trial in range(5000):
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
    c7 = int(np.trace(np.linalg.matrix_power(A, 7))) // 7

    simplex_profiles[profile].append({'c7': c7, 'tr3': tr3})

ambig = [(p, entries) for p, entries in simplex_profiles.items()
         if len(set(e['c7'] for e in entries)) > 1]

all_gaps = []
for profile, entries in ambig:
    tr3_vals = sorted(set(e['tr3'] for e in entries))
    if len(tr3_vals) >= 2:
        for i in range(len(tr3_vals)):
            for j in range(i+1, len(tr3_vals)):
                gap = tr3_vals[j] - tr3_vals[i]
                all_gaps.append(gap)

print(f"  All tr3 gaps: {sorted(set(all_gaps))}")
print(f"  GCD of all gaps: {np.gcd.reduce(all_gaps) if all_gaps else 'N/A'}")

if all_gaps:
    gcd_val = int(np.gcd.reduce(all_gaps))
    print(f"  GCD = {gcd_val}")
    print(f"  {gcd_val} = {gcd_val//2} * 2 = {gcd_val//3} * 3 = {gcd_val//4} * 4")
    print(f"  {gcd_val} = 2^{bin(gcd_val).count('0') if gcd_val > 0 else 0}...")

    # Factorize the GCD
    g = gcd_val
    factors = {}
    for p in [2, 3, 5, 7, 11]:
        while g % p == 0:
            factors[p] = factors.get(p, 0) + 1
            g //= p
    print(f"  GCD factorization: {factors}")

    # Check all gaps divisible by 12
    all_div_12 = all(g % 12 == 0 for g in all_gaps)
    print(f"  All gaps divisible by 12: {all_div_12}")

    # Distribution of gap/12
    print(f"  gap/12 distribution: {Counter(g//12 for g in all_gaps)}")

print("\n\nDone.")
