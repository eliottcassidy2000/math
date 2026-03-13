"""
two_three_deep.py -- kind-pasteur-2026-03-13-S62

Deep exploration of the numbers 2 and 3 in tournament parity theory.

WHY is x=2 the evaluation point in H(T) = I(Omega(T), 2)?
WHY do 3-cycles control everything?
What is the 2-adic and 3-adic structure of H(T)?

Key threads:
1. I(Omega(T), x) at x=1,2,3,4,...  — what's special about x=2?
2. F(T,x) mod 2 = (1+x)^{n-1} — the 2 is deeply embedded
3. F(T,omega) divisible by 9 for n>=6 — the 3 is deeply embedded
4. The factor of 2 in Claim A: H(T)-H(T-v) = 2*sum mu(C)
5. total_lambda = 3*c3 — why factor 3?
6. n-2 = sigma + lambda + delta — the triple decomposition
7. 2-adic valuation v_2(H(T)) — how does it vary?
8. 3-adic valuation v_3(H(T)) — how does it vary?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
from math import gcd, factorial

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
    """Count Hamiltonian paths using inclusion-exclusion / DP."""
    # Held-Karp DP
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

def build_omega(A, n):
    """Build conflict graph Omega(T)."""
    cycles_3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if {i,j,k} forms a directed 3-cycle
                verts = [i, j, k]
                for perm in permutations(verts):
                    if A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and A[perm[2]][perm[0]]:
                        cycles_3.append(frozenset(verts))
                        break
    cycles_3 = list(set(cycles_3))

    # Build adjacency: two cycles conflict if they share a vertex
    nc = len(cycles_3)
    omega_adj = np.zeros((nc, nc), dtype=int)
    for a in range(nc):
        for b in range(a+1, nc):
            if cycles_3[a] & cycles_3[b]:
                omega_adj[a][b] = omega_adj[b][a] = 1
    return cycles_3, omega_adj

def independence_poly(adj, nc):
    """Compute I(G, x) as polynomial coefficients [c0, c1, c2, ...]."""
    # Enumerate all independent sets
    coeffs = [0] * (nc + 1)
    for mask in range(1 << nc):
        # Check independence
        verts = [i for i in range(nc) if mask & (1 << i)]
        independent = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            coeffs[len(verts)] += 1
    return coeffs

def eval_poly(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))

def v2(n):
    """2-adic valuation"""
    if n == 0: return float('inf')
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v

def v3(n):
    """3-adic valuation"""
    if n == 0: return float('inf')
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v

def count_ck(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

# ============================================================
# PART 1: Why x=2? I(Omega, x) at various evaluation points
# ============================================================
print("=" * 70)
print("PART 1: WHY x=2? Independence polynomial at various points")
print("=" * 70)

for n in [5, 6, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    print(f"\n--- n={n} ---")

    # Collect I(Omega, x) for x = -1, 0, 1, 2, 3, 4
    results = defaultdict(list)

    n_samples = min(200, 2**total_bits) if n <= 5 else 500

    for trial in range(n_samples):
        if n <= 5 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        c3 = count_ck(A, n, 3)

        cycles, omega_adj = build_omega(A, n)
        nc = len(cycles)

        if nc <= 20:  # Can compute full IP
            ip_coeffs = independence_poly(omega_adj, nc)

            for x in [-1, 0, 1, 2, 3, 4]:
                val = eval_poly(ip_coeffs, x)
                results[(n, x)].append((val, H, c3))

    for x in [-1, 0, 1, 2, 3, 4]:
        key = (n, x)
        if key in results:
            vals = [r[0] for r in results[key]]
            Hs = [r[1] for r in results[key]]

            # Check if I(Omega, x) = H
            matches_H = sum(1 for v, h in zip(vals, Hs) if v == h)

            print(f"  I(Omega, {x}): range [{min(vals)}, {max(vals)}], "
                  f"mean={np.mean(vals):.1f}, "
                  f"matches H: {matches_H}/{len(vals)}, "
                  f"I mod 2: {Counter(v % 2 for v in vals)}")

# ============================================================
# PART 2: The 2-adic structure of H(T)
# ============================================================
print("\n" + "=" * 70)
print("PART 2: 2-ADIC VALUATION v_2(H(T))")
print("=" * 70)

for n in [3, 4, 5, 6, 7]:
    total_bits = n*(n-1)//2

    v2_counts = Counter()
    H_mod4 = Counter()
    H_mod8 = Counter()

    n_samples = min(2**total_bits, 5000)

    for trial in range(n_samples):
        if n <= 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)
            np.random.seed(42 + trial)

        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)

        v2_counts[v2(H)] += 1
        H_mod4[H % 4] += 1
        H_mod8[H % 8] += 1

    print(f"\n--- n={n} (from {n_samples} tournaments) ---")
    print(f"  v_2(H) distribution: {dict(sorted(v2_counts.items()))}")
    print(f"  H mod 4: {dict(sorted(H_mod4.items()))}")
    print(f"  H mod 8: {dict(sorted(H_mod8.items()))}")

    # v_2(n!) for reference
    nfact = factorial(n)
    print(f"  v_2({n}!) = {v2(nfact)}")

# ============================================================
# PART 3: The 3-adic structure of H(T)
# ============================================================
print("\n" + "=" * 70)
print("PART 3: 3-ADIC VALUATION v_3(H(T))")
print("=" * 70)

for n in [3, 4, 5, 6, 7]:
    total_bits = n*(n-1)//2

    v3_counts = Counter()
    H_mod3 = Counter()
    H_mod9 = Counter()
    H_mod27 = Counter()

    n_samples = min(2**total_bits, 5000)

    for trial in range(n_samples):
        if n <= 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)
            np.random.seed(42 + trial)

        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)

        v3_counts[v3(H)] += 1
        H_mod3[H % 3] += 1
        H_mod9[H % 9] += 1
        H_mod27[H % 27] += 1

    print(f"\n--- n={n} (from {n_samples} tournaments) ---")
    print(f"  v_3(H) distribution: {dict(sorted(v3_counts.items()))}")
    print(f"  H mod 3: {dict(sorted(H_mod3.items()))}")
    print(f"  H mod 9: {dict(sorted(H_mod9.items()))}")
    print(f"  H mod 27: {dict(sorted(H_mod27.items()))}")
    print(f"  v_3({n}!) = {v3(factorial(n))}")

# ============================================================
# PART 4: The factor of 2 in Claim A
# ============================================================
print("\n" + "=" * 70)
print("PART 4: H(T) - H(T-v) = 2 * sum mu(C)")
print("=" * 70)
print("The deletion ratio R(T) = sum_v [H(T) - H(T\\v)] / (n * H(T))")
print("H(T) - H(T\\v) is always EVEN (from Claim A, the factor of 2)")

for n in [5, 6, 7]:
    total_bits = n*(n-1)//2

    delta_parities = Counter()
    delta_mod4 = Counter()

    n_samples = min(2**total_bits, 1000) if n <= 6 else 300

    for trial in range(n_samples):
        if n <= 5:
            bits = trial
        elif n == 6 and trial < 2**total_bits:
            bits = trial
        else:
            np.random.seed(42 + trial)
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)

        for v in range(n):
            # Delete vertex v
            remaining = [u for u in range(n) if u != v]
            A_del = np.zeros((n-1, n-1), dtype=int)
            for i, u in enumerate(remaining):
                for j, w in enumerate(remaining):
                    A_del[i][j] = A[u][w]
            H_del = count_ham_paths(A_del, n-1)
            delta = H - H_del
            delta_parities[delta % 2] += 1
            delta_mod4[delta % 4] += 1

    print(f"\n--- n={n} ---")
    print(f"  H(T)-H(T\\v) mod 2: {dict(sorted(delta_parities.items()))}")
    print(f"  H(T)-H(T\\v) mod 4: {dict(sorted(delta_mod4.items()))}")
    even_count = delta_parities.get(0, 0)
    total = sum(delta_parities.values())
    print(f"  Always even: {even_count}/{total} ({100*even_count/total:.1f}%)")

# ============================================================
# PART 5: Why 3-cycles? The role of 3 in cycle structure
# ============================================================
print("\n" + "=" * 70)
print("PART 5: WHY 3-CYCLES? The counting structure")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

# Count total_lambda = 3*c3 verification
# AND: what other multiples appear?
# total_lambda / c3 = 3 always
# total cycles of length k per vertex set: always the same?

c_to_h = defaultdict(list)
for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    c3 = count_ck(A, n, 3)
    c5 = count_ck(A, n, 5)
    c7 = count_ck(A, n, 7)

    c_to_h[(c3, c5, c7)].append(H)

# H as function of cycle counts
print(f"  Distinct (c3,c5,c7) triples: {len(c_to_h)}")
h_determined = sum(1 for hs in c_to_h.values() if len(set(hs)) == 1)
print(f"  c7 determines H uniquely: checking via OCF")
print(f"  H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...")
print(f"  where alpha_k = #{'{'}independent sets of size k in Omega{'}'}")

# The binary structure: H = sum alpha_k * 2^k
print(f"\n  H in base 2 (binary expansion via independence sets):")
np.random.seed(42)
for trial in range(5):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    cycles, omega_adj = build_omega(A, n)
    nc = len(cycles)

    if nc <= 20:
        ip = independence_poly(omega_adj, nc)
        H_check = eval_poly(ip, 2)
        binary = bin(H)
        print(f"    H={H}, I(Omega,2)={H_check}, alpha={ip[:6]}, "
              f"binary={binary}")

# ============================================================
# PART 6: The 2-3 interplay: H mod 6
# ============================================================
print("\n" + "=" * 70)
print("PART 6: H mod 6 (Chinese Remainder: mod 2 and mod 3 combined)")
print("=" * 70)

for n in [3, 4, 5, 6, 7]:
    total_bits = n*(n-1)//2

    H_mod6 = Counter()
    n_samples = min(2**total_bits, 5000)

    for trial in range(n_samples):
        if n <= 6 and trial < 2**total_bits:
            bits = trial
        else:
            np.random.seed(42 + trial)
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        H_mod6[H % 6] += 1

    print(f"\n  n={n}: H mod 6 = {dict(sorted(H_mod6.items()))}")
    print(f"    Redei: H is odd, so H mod 6 in {{1,3,5}}")
    if all(r in [1,3,5] for r in H_mod6.keys()):
        print(f"    CONFIRMED: only odd residues")
    # For n>=6, H mod 3 might be constrained
    H_mod3_from6 = Counter()
    for r, c in H_mod6.items():
        H_mod3_from6[r % 3] += c
    print(f"    H mod 3: {dict(sorted(H_mod3_from6.items()))}")

# ============================================================
# PART 7: The 2-expansion of I(Omega, x)
# ============================================================
print("\n" + "=" * 70)
print("PART 7: I(Omega, x) = sum alpha_k * x^k — coefficient structure")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

alpha_distributions = defaultdict(list)

for trial in range(500):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    cycles, omega_adj = build_omega(A, n)
    nc = len(cycles)

    if nc <= 20:
        ip = independence_poly(omega_adj, nc)
        for k, a in enumerate(ip):
            if a > 0:
                alpha_distributions[k].append(a)

print(f"  alpha_k statistics at n={n}:")
for k in sorted(alpha_distributions.keys()):
    vals = alpha_distributions[k]
    print(f"    alpha_{k}: min={min(vals)}, max={max(vals)}, "
          f"mean={np.mean(vals):.1f}, "
          f"always_odd={all(v%2==1 for v in vals)}, "
          f"always_div3={all(v%3==0 for v in vals)}")

# ============================================================
# PART 8: The n-2 = sigma + lambda + delta decomposition
#          and the role of 2 and 3
# ============================================================
print("\n" + "=" * 70)
print("PART 8: n-2 = sigma + lambda + delta — the triple budget")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

# For each pair, n-2 = 5 is split into (sigma, lambda, delta)
# The possible splits of 5 into 3 non-negative integers
print(f"  At n={n}, n-2={n-2}")
print(f"  Possible (sigma, lambda, delta) triples summing to {n-2}:")

triple_counts = Counter()
for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    for u in range(n):
        for v in range(u+1, n):
            # lambda = # 3-cycles containing u,v
            lam = 0
            sig = 0
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    lam += 1
                # sigma = common successors + common predecessors
                if (A[u][w] and A[v][w]) or (A[w][u] and A[w][v]):
                    sig += 1
            delta = n - 2 - sig - lam
            triple_counts[(sig, lam, delta)] += 1

print(f"  Distribution of (sigma, lambda, delta) per pair:")
for triple in sorted(triple_counts.keys()):
    print(f"    {triple}: {triple_counts[triple]}")

# The mean coordinates
total = sum(triple_counts.values())
mean_sig = sum(t[0]*c for t, c in triple_counts.items()) / total
mean_lam = sum(t[1]*c for t, c in triple_counts.items()) / total
mean_del = sum(t[2]*c for t, c in triple_counts.items()) / total
print(f"\n  Mean: sigma={mean_sig:.3f}, lambda={mean_lam:.3f}, delta={mean_del:.3f}")
print(f"  Ratio sigma:lambda:delta = {mean_sig/mean_lam:.3f}:{1:.3f}:{mean_del/mean_lam:.3f}")
print(f"  sigma/lambda ~ 2:1 ? {mean_sig/mean_lam:.4f}")

# ============================================================
# PART 9: F(T, x) at x=2 vs x=3
# ============================================================
print("\n" + "=" * 70)
print("PART 9: F(T, x) at special points")
print("=" * 70)
print("F(T, x) = sum over Ham paths P of x^{asc(P)}")
print("F(T, 1) = n! / n = (n-1)! ... wait, F(T,1) = H(T)")
print("Actually F(T, x) = sum_k F_k * x^k where F_k = # Ham paths with k ascents")

for n in [4, 5, 6]:
    total_bits = n*(n-1)//2

    print(f"\n  --- n={n} ---")

    n_samples = min(2**total_bits, 500)

    fx_at = defaultdict(list)

    for trial in range(n_samples):
        if trial < 2**total_bits:
            bits = trial
        else:
            np.random.seed(42 + trial)
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)

        # Enumerate all Hamiltonian paths and count ascents
        F_coeffs = [0] * n

        for perm in permutations(range(n)):
            # Check if perm is a Ham path
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
        F_at_2 = sum(c * 2**k for k, c in enumerate(F_coeffs))
        F_at_3 = sum(c * 3**k for k, c in enumerate(F_coeffs))
        F_at_neg1 = sum(c * (-1)**k for k, c in enumerate(F_coeffs))

        fx_at['H'].append(H)
        fx_at['F(2)'].append(F_at_2)
        fx_at['F(3)'].append(F_at_3)
        fx_at['F(-1)'].append(F_at_neg1)

    for key in ['H', 'F(2)', 'F(3)', 'F(-1)']:
        vals = fx_at[key]
        print(f"    {key}: range [{min(vals)}, {max(vals)}], "
              f"mod2={Counter(v%2 for v in vals)}, "
              f"mod3={Counter(v%3 for v in vals)}")

    # F(T,1) should be H(T)
    # F(T,-1) should be +-1 (Redei-like?)
    print(f"    F(T,1)=H(T): {all(a==b for a,b in zip(fx_at['H'], fx_at['H']))}")
    print(f"    F(T,-1) values: {sorted(set(fx_at['F(-1)']))[:20]}")

    # What about F(T,2)/H(T) and F(T,3)/H(T)?
    ratios_2 = [f/h for f, h in zip(fx_at['F(2)'], fx_at['H']) if h > 0]
    ratios_3 = [f/h for f, h in zip(fx_at['F(3)'], fx_at['H']) if h > 0]

    print(f"    F(2)/H: range [{min(ratios_2):.3f}, {max(ratios_2):.3f}], mean={np.mean(ratios_2):.3f}")
    print(f"    F(3)/H: range [{min(ratios_3):.3f}, {max(ratios_3):.3f}], mean={np.mean(ratios_3):.3f}")

print("\n\nDone with Part 1-9.")
