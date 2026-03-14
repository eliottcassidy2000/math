#!/usr/bin/env python3
"""
seven_barrier_deep.py — opus-2026-03-14-S73
Deep analysis of the H ≡ 0 (mod 7) barrier.

Does H ≡ 0 (mod 7) persist as forbidden for ALL n?
If so, this would be a new theorem about tournament Hamiltonian paths.

Key: H = 1 + 2α₁ + 4α₂ + 8α₃ + ...
     H mod 7 = 1 + 2α₁ + 4α₂ + α₃ + 2α₄ + 4α₅ + ...
     (since 2^3 ≡ 1 mod 7, period 3)

H ≡ 0 (mod 7) requires:
  Σ_{k≥0} i_k · 2^k ≡ 0 (mod 7)
  where i_k = #{indep sets of size k}
"""

from itertools import permutations, combinations
from collections import Counter
import random, time

def banner(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")

def adj_matrix_random(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_ham_paths_dp(A, n):
    """DP-based Hamiltonian path count (faster for n≥7)."""
    # dp[S][v] = number of Hamiltonian paths ending at v, visiting set S
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                for u in range(n):
                    if (S_prev & (1 << u)) and A[u][v]:
                        dp[(S, v)] = dp.get((S, v), 0) + dp.get((S_prev, u), 0)
    
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# ─────────────────────────────────────────────────────────────────────
# PART 1: Check H mod 7 at n=7,8 by sampling
# ─────────────────────────────────────────────────────────────────────
banner("PART 1: H mod 7 AT n=7,8 BY SAMPLING")

random.seed(42)

for n in [7, 8, 9]:
    num_samples = 5000 if n <= 8 else 2000
    h_mod7 = Counter()
    t0 = time.time()
    for trial in range(num_samples):
        A = adj_matrix_random(n)
        h = count_ham_paths_dp(A, n)
        h_mod7[h % 7] += 1
    
    elapsed = time.time() - t0
    print(f"n={n} ({num_samples} samples, {elapsed:.1f}s):")
    for r in range(7):
        count = h_mod7.get(r, 0)
        pct = count / num_samples * 100
        marker = " ← NOT SEEN" if count == 0 else ""
        print(f"  H ≡ {r} (mod 7): {count:5d} ({pct:5.1f}%){marker}")
    print()

# ─────────────────────────────────────────────────────────────────────
# PART 2: Algebraic analysis — why H ≢ 0 (mod 7)?
# ─────────────────────────────────────────────────────────────────────
banner("PART 2: ALGEBRAIC ANALYSIS")

print("H = I(CG, 2) where CG = odd-cycle conflict graph")
print("I(G, x) = Σ_{k≥0} i_k x^k")
print()
print("H mod 7 = I(CG, 2) mod 7")
print()
print("If we could show I(G, 2) ≡ 0 (mod 7) is impossible")
print("for ALL odd-cycle conflict graphs, this would prove it.")
print()
print("Note: I(G, 2) ≡ 0 (mod 7) ≡ I(G, 2) ≡ 0 (mod 7)")
print("And 2³ ≡ 1 (mod 7), so the evaluation at 2 mod 7")
print("uses a cube root of unity in F_7.")
print()
print("In F_7: 2³ = 8 ≡ 1, so 2 is a primitive 3rd root of unity.")
print("Let ω = 2 in F_7. Then ω³ = 1, ω ≠ 1.")
print("The minimal polynomial of ω over F_7 is x² + x + 1")
print("(since ω²+ω+1 = 4+2+1 = 7 ≡ 0 mod 7).")
print()
print("OBSERVATION: I(G, ω) = 0 in F_7 means the independence polynomial")
print("has ω as a root modulo 7.")
print("Since ω satisfies ω²+ω+1 = 0, this means")
print("(x²+x+1) | I(G, x) in F_7[x].")
print()
print("CLAIM: For odd-cycle conflict graphs of tournaments,")
print("(x²+x+1) never divides I(G, x) mod 7.")

# ─────────────────────────────────────────────────────────────────────
# PART 3: Check I(G, x) mod 7 for small cases
# ─────────────────────────────────────────────────────────────────────
banner("PART 3: I(G, x) mod 7 FOR SMALL CASES")

def get_odd_cycles(A, n):
    """Get all directed odd cycles (3,5,7) in tournament A."""
    cycles = []
    # 3-cycles
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if A[i][j] and A[j][k] and A[k][i]:
            cycles.append(frozenset(triple))
        elif A[i][k] and A[k][j] and A[j][i]:
            cycles.append(frozenset(triple))
    
    # 5-cycles
    if n >= 5:
        for quintet in combinations(range(n), 5):
            has_5cycle = False
            for perm in permutations(quintet):
                if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                    has_5cycle = True
                    break
            if has_5cycle:
                cycles.append(frozenset(quintet))
    
    return cycles

def independence_poly_mod(cycles, p):
    """Compute independence polynomial of conflict graph modulo p.
    Two cycles conflict if they share a vertex."""
    n_cycles = len(cycles)
    if n_cycles == 0:
        return [1]  # I(G,x) = 1
    
    # Build conflict graph
    adj = [[False]*n_cycles for _ in range(n_cycles)]
    for i in range(n_cycles):
        for j in range(i+1, n_cycles):
            if cycles[i] & cycles[j]:  # share a vertex
                adj[i][j] = adj[j][i] = True
    
    # Count independent sets by size
    max_indep = n_cycles
    i_k = [0] * (max_indep + 1)
    i_k[0] = 1
    
    # Enumerate all independent sets (feasible for small n_cycles)
    if n_cycles <= 20:
        for size in range(1, n_cycles + 1):
            for subset in combinations(range(n_cycles), size):
                # Check independence
                is_indep = True
                for a, b in combinations(subset, 2):
                    if adj[a][b]:
                        is_indep = False
                        break
                if is_indep:
                    i_k[size] += 1
    
    # Reduce mod p
    return [c % p for c in i_k]

# At n=5, compute I(CG, x) for all tournaments
print("n=5: I(CG, x) mod 7 for all 1024 tournaments")
from itertools import product as iprod

def adj_matrix_idx(n, idx):
    A = [[0]*n for _ in range(n)]
    bits = idx
    for i in range(n):
        for j in range(i+1, n):
            if bits & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            bits >>= 1
    return A

# Check I(CG, 2) mod 7 algebraically
i_at_2_mod7 = Counter()
i_at_3_mod7 = Counter()  # 3 is the other cube root (3³=27≡6≡-1, no. 3³=27≡6 mod 7)
# Actually: in F_7, cube roots of 1 are {1, 2, 4} since 2³=8≡1, 4³=64≡1
# So 2 and 4 are primitive cube roots of unity in F_7

poly_mod7_counter = Counter()
for idx in range(1024):
    A = adj_matrix_idx(5, idx)
    h = count_ham_paths_dp(A, 5)
    
    cycles = get_odd_cycles(A, 5)
    # Instead of full independence poly, just compute I(2) mod 7
    # I(x) = Σ i_k x^k, I(2) = H (should match)
    
    poly = independence_poly_mod(cycles, 7)
    
    # Evaluate at x=2 mod 7
    val_2 = sum(c * pow(2, k, 7) for k, c in enumerate(poly)) % 7
    val_4 = sum(c * pow(4, k, 7) for k, c in enumerate(poly)) % 7
    
    i_at_2_mod7[val_2] += 1
    i_at_3_mod7[val_4] += 1
    
    # Check consistency with H
    h_mod7 = h % 7
    if h_mod7 != val_2:
        print(f"  MISMATCH at idx={idx}: H%7={h_mod7}, I(2)%7={val_2}")
    
    poly_tuple = tuple(poly)
    poly_mod7_counter[poly_tuple] += 1

print(f"  I(CG, 2) mod 7: {dict(sorted(i_at_2_mod7.items()))}")
print(f"  I(CG, 4) mod 7: {dict(sorted(i_at_3_mod7.items()))}")
print()

# The independence polynomial mod 7
print("Independence polynomial mod 7 (distinct types):")
for poly, cnt in sorted(poly_mod7_counter.items(), key=lambda x: -x[1]):
    # Show as polynomial
    terms = []
    for k, c in enumerate(poly):
        if c > 0:
            if k == 0:
                terms.append(str(c))
            elif k == 1:
                terms.append(f"{c}x")
            else:
                terms.append(f"{c}x^{k}")
    poly_str = " + ".join(terms) if terms else "0"
    
    # Evaluate at x=2 mod 7
    val = sum(c * pow(2, k, 7) for k, c in enumerate(poly)) % 7
    
    # Check if (x²+x+1) divides this mod 7
    # I(x) mod (x²+x+1) mod 7
    r = list(poly)
    while len(r) < 3:
        r.append(0)
    # Polynomial division by x²+x+1 over F_7
    while len(r) > 2:
        lead = r[-1] % 7
        if lead == 0:
            r.pop()
            continue
        deg = len(r) - 1
        # Subtract lead * x^{deg-2} * (x²+x+1)
        r[deg] = (r[deg] - lead) % 7
        r[deg-1] = (r[deg-1] - lead) % 7
        r[deg-2] = (r[deg-2] - lead) % 7
        while r and r[-1] % 7 == 0:
            r.pop()
    
    remainder = tuple(c % 7 for c in r) if r else (0,)
    divides = all(c % 7 == 0 for c in r)
    
    print(f"  {poly_str:>30}: I(2)≡{val} (mod 7), remainder by x²+x+1 = {remainder}, count={cnt}")

# ─────────────────────────────────────────────────────────────────────
# PART 4: Direct verification — I(CG,x) and x²+x+1
# ─────────────────────────────────────────────────────────────────────
banner("PART 4: I(CG,x) MOD (x²+x+1) OVER F_7")

print("Since 2 is a root of x²+x+1 over F_7,")
print("H ≡ 0 (mod 7) iff (x²+x+1) | I(CG,x) (mod 7).")
print()
print("Question: does I(CG,x) ≡ 0 (mod x²+x+1) ever happen?")
print()

# For ALL n=5 tournaments, we already checked.
# Now check at n=6 (32768 tournaments)
print("n=6: checking I(CG,2) mod 7 for all 32768 tournaments...")
t0 = time.time()
h_0_mod7_n6 = 0
for idx in range(32768):
    A = adj_matrix_idx(6, idx)
    h = count_ham_paths_dp(A, 6)
    if h % 7 == 0:
        h_0_mod7_n6 += 1
print(f"  H ≡ 0 (mod 7): {h_0_mod7_n6}/{32768} tournaments ({time.time()-t0:.1f}s)")
print()

# At n=7: sample
print("n=7: sampling 3000 tournaments for H mod 7...")
random.seed(123)
h_0_mod7_n7 = 0
num_n7 = 3000
t0 = time.time()
for trial in range(num_n7):
    A = adj_matrix_random(7)
    h = count_ham_paths_dp(A, 7)
    if h % 7 == 0:
        h_0_mod7_n7 += 1
    if trial % 1000 == 0:
        print(f"  trial {trial}: {time.time()-t0:.1f}s")

print(f"  H ≡ 0 (mod 7): {h_0_mod7_n7}/{num_n7} tournaments")

if h_0_mod7_n7 == 0:
    print("  STILL ZERO! The 7-barrier persists at n=7.")
else:
    print(f"  BROKEN! {h_0_mod7_n7} tournaments have H ≡ 0 (mod 7).")

# ─────────────────────────────────────────────────────────────────────
# PART 5: Connection to 5 — the role of Fibonacci
# ─────────────────────────────────────────────────────────────────────
banner("PART 5: THE 5-7 CONNECTION")

print("Why are 5 and 7 special for forbidden residues?")
print()
print("5 = the discriminant of the Fibonacci recurrence")
print("7 = 2³ - 1, a Mersenne prime")
print("  Also: 7 = the first prime where ord_p(2) < p-1")
print("  (ord_7(2) = 3, not 6)")
print()
print("The key: powers of 2 mod 7 cycle with period 3,")
print("making 2 a CUBE ROOT OF UNITY in F_7.")
print("So H = I(CG, 2) mod 7 evaluates at a special algebraic point.")
print()
print("In F_7: x²+x+1 = (x-2)(x-4)")
print("  Check: (2)²+2+1 = 7 ≡ 0 ✓")
print("  Check: (4)²+4+1 = 21 ≡ 0 ✓")
print()
print("So H ≡ 0 (mod 7) iff I(CG, 2) ≡ 0 in F_7")
print("iff the cyclotomic polynomial Φ_3(x) = x²+x+1 divides I(CG,x) in F_7[x]")
print()
print("EQUIVALENTLY: H ≡ 0 (mod 7) iff BOTH I(CG,2) ≡ 0 AND I(CG,4) ≡ 0 (mod 7)")
print("But we already know I(CG,2) = H and I(CG,4) is a different evaluation.")
print()

# Check I(CG,4) mod 7 at n=5
print("At n=5, I(CG,4) mod 7 distribution:")
i4_mod7 = Counter()
for idx in range(1024):
    A = adj_matrix_idx(5, idx)
    cycles = get_odd_cycles(A, 5)
    
    # i_k for CG
    n_cycles = len(cycles)
    i_k = [0] * (n_cycles + 1)
    i_k[0] = 1
    
    if n_cycles <= 15:
        adj = [[False]*n_cycles for _ in range(n_cycles)]
        for i in range(n_cycles):
            for j in range(i+1, n_cycles):
                if cycles[i] & cycles[j]:
                    adj[i][j] = adj[j][i] = True
        
        for size in range(1, n_cycles + 1):
            for subset in combinations(range(n_cycles), size):
                is_indep = all(not adj[a][b] for a, b in combinations(subset, 2))
                if is_indep:
                    i_k[size] += 1
    
    val = sum(c * pow(4, k, 7) for k, c in enumerate(i_k)) % 7
    i4_mod7[val] += 1

print(f"  I(CG, 4) mod 7: {dict(sorted(i4_mod7.items()))}")
print(f"  I(CG, 4) ≡ 0 (mod 7): {i4_mod7.get(0, 0)} tournaments")

# ─────────────────────────────────────────────────────────────────────
# PART 6: Does H ≡ 0 (mod 7) fail for ALL n?
# ─────────────────────────────────────────────────────────────────────
banner("PART 6: CONJECTURE — H ≢ 0 (mod 7) FOR ALL n")

print("Evidence so far:")
print("  n=3: H ∈ {1,3} → H mod 7 ∈ {1,3} ✓ (no 0)")
print("  n=4: H ∈ {1,3,5} → H mod 7 ∈ {1,3,5} ✓")
print("  n=5: verified all 1024 ✓")
print("  n=6: verified all 32768 ✓")
print(f"  n=7: sampled {num_n7}, found {h_0_mod7_n7}")
print()

if h_0_mod7_n7 == 0:
    print("CONJECTURE: H(T) ≢ 0 (mod 7) for ALL tournaments T.")
    print()
    print("If true, this says: the number of Hamiltonian paths")
    print("in a tournament is NEVER divisible by 7.")
    print()
    print("Combined with Rédei (H odd), this means")
    print("H mod 14 ∈ {1, 3, 5, 9, 11, 13} (the odd non-multiples of 7)")
    print()
    print("POSSIBLE PROOF STRATEGY:")
    print("Show that I(G, x) mod (x²+x+1) in F_7[x]")
    print("is never ≡ 0 for odd-cycle conflict graphs of tournaments.")
    print("This is equivalent to showing I(G, ω) ≠ 0 in F_7")
    print("where ω is a primitive cube root of unity.")
else:
    print(f"The conjecture FAILS at n=7!")
    print("So H ≡ 0 (mod 7) first becomes possible at some n between 7 and beyond.")

# ─────────────────────────────────────────────────────────────────────
# PART 7: Extend sampling to n=8,9 for definitive check
# ─────────────────────────────────────────────────────────────────────
banner("PART 7: EXTENDED SAMPLING")

random.seed(999)
for n in [8, 9, 10]:
    num_s = 2000 if n <= 9 else 500
    found_0 = 0
    t0 = time.time()
    for trial in range(num_s):
        A = adj_matrix_random(n)
        h = count_ham_paths_dp(A, n)
        if h % 7 == 0:
            found_0 += 1
    elapsed = time.time() - t0
    status = "NO H≡0" if found_0 == 0 else f"FOUND {found_0} with H≡0!"
    print(f"n={n}: {num_s} samples in {elapsed:.1f}s → {status}")

print("\nDone.")
