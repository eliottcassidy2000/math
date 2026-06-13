"""
two_three_deeper.py -- kind-pasteur-2026-03-13-S62

Deeper exploration of the 2:1:1 ratio and the algebraic structure.

From Part 1 findings:
- sigma:lambda:delta = 2:1:1 EXACTLY (proved below)
- I(Omega, x) odd iff x even (because I(G,0)=1 and 2|x => I(G,x) = I(G,0) = 1 mod 2)
- This means x=2 is the SMALLEST positive integer making I(Omega,x) always odd

Key questions:
A. Can we PROVE sigma/lambda = 2 exactly?
B. What is the connection between the 2:1:1 ratio and x=2?
C. Is there a "3-adic" analog?
D. The Worpitzky decomposition F(T,x) = sum F_k * C(x-1,k)
E. Why does 3 appear in total_lambda = 3*c3?
F. The "2*3 = 6" structure: H mod 6 is constrained
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
from math import factorial, comb

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

def count_ck(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

# ============================================================
# PART A: PROVING sigma/lambda = 2 exactly
# ============================================================
print("=" * 70)
print("PART A: PROOF that E[sigma]/E[lambda] = 2 exactly")
print("=" * 70)

print("""
For a RANDOM tournament (each arc independently 50/50):

For a pair (u,v), consider a witness vertex w != u,v.
The 3 arcs involving w are: u-w, v-w, u-v (but u-v is fixed for this pair).

Given arc direction u->v (WLOG by symmetry):
  lambda witness: w forms a 3-cycle with u,v.
    Either u->v->w->u or u<-v<-w<-u.
    Since u->v: need v->w->u (probability 1/4)
    Or v->u: need u->w->v, but u->v so this means v->u (contradiction).
    So lambda = Pr(v->w AND w->u) = 1/4.

  sigma witness: w has same relation to both u and v.
    common successor: u->w AND v->w, prob = 1/4
    common predecessor: w->u AND w->v, prob = 1/4
    So sigma = 1/4 + 1/4 = 1/2.

  delta witness: transitive triple through w
    u->v->w AND u->w: prob = 1/2 * 1/2 = 1/4 (given u->v)
    w->u->v AND w->v: prob = 1/2 * 1/2 = 1/4 (given u->v)
    Hmm but need to be more careful.

  Actually, for w, the 2 arcs u-w and v-w are independent.
  There are 4 cases (each prob 1/4):
    Case 1: u->w, v->w (common successor) => sigma witness
    Case 2: u->w, w->v (given u->v: u->v, u->w, w->v => transitive) => delta
    Case 3: w->u, v->w (given u->v: u->v, w->u, v->w => 3-cycle v->w->u->v) => lambda
    Case 4: w->u, w->v (common predecessor) => sigma witness

  So: sigma = 2/4 = 1/2, lambda = 1/4, delta = 1/4
  Ratio: sigma:lambda:delta = 2:1:1

  THIS IS EXACT, not just asymptotic!
""")

# Verify computationally
for n in [5, 6, 7, 8]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    total_sig = 0
    total_lam = 0
    total_del = 0
    n_pairs = 0

    n_samples = min(2**total_bits, 3000)

    for trial in range(n_samples):
        if n <= 5 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)

        for u in range(n):
            for v in range(u+1, n):
                lam = 0
                sig = 0
                for w in range(n):
                    if w == u or w == v:
                        continue
                    # 3-cycle?
                    if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                        lam += 1
                    # Common successor or predecessor?
                    if (A[u][w] and A[v][w]) or (A[w][u] and A[w][v]):
                        sig += 1
                delta = n - 2 - sig - lam
                total_sig += sig
                total_lam += lam
                total_del += delta
                n_pairs += 1

    avg_sig = total_sig / n_pairs
    avg_lam = total_lam / n_pairs
    avg_del = total_del / n_pairs

    print(f"  n={n}: avg sigma={avg_sig:.4f}, lambda={avg_lam:.4f}, delta={avg_del:.4f}")
    print(f"         sigma/lambda = {avg_sig/avg_lam:.6f}")
    print(f"         Expected: sigma = (n-2)/2 = {(n-2)/2}, lambda = (n-2)/4 = {(n-2)/4}")

# ============================================================
# PART B: The 2:1 ratio and x=2
# ============================================================
print("\n" + "=" * 70)
print("PART B: Connection between sigma:lambda = 2:1 and x=2")
print("=" * 70)

print("""
The OCF formula: H(T) = I(Omega(T), 2)

The independence polynomial at x=2 counts independent sets weighted by 2^|S|.
Each independent set of k 3-cycles contributes 2^k.

WHY x=2? Because:
1. H(T) = #{Hamiltonian paths} = sum over odd-cycle coverings of (-1)^...
2. The Rédei function at x=1 counts Ham paths
3. The OCF expands: each 3-cycle contributes a BINARY choice (include or not)
4. The weight 2 comes from: for each independent 3-cycle,
   there are 2 orientations of the cycle that don't conflict

Wait — let's think about this differently.

H(T) = I(Omega(T), 2).
But also H(T) = F(T, 1) where F is the Rédei-Berge function.

The Worpitzky decomposition: F(T, x) = sum_{k=0}^{n-1} F_k(T) * C(x,k)
where F_k counts Ham paths with exactly k ascents.

At x=1: F(T,1) = sum F_k C(1,k) = F_0 + F_1 = H(T)
At x=2: F(T,2) = sum F_k C(2,k) = F_0 + 2*F_1 + F_2
At x=-1: F(T,-1) = sum F_k C(-1,k) = F_0 - F_1 + F_2 - F_3 + ...

Actually with C(x-1, k) binomial: F(T,x) = sum W_k * C(x-1, k)
At x=1: F(T,1) = W_0 = H(T)

So H = W_0. And I(Omega, 2) = H. The 2 in I(Omega, 2) is NOT directly
the 2 from F(T,2). They're different evaluations of different polynomials.

The real connection:
  I(Omega, x) = sum_{k=0}^{alpha(Omega)} alpha_k * x^k

  H = I(Omega, 2) = alpha_0 + 2*alpha_1 + 4*alpha_2 + ...

  Each independent k-set contributes 2^k to H.
  The binary weights 1, 2, 4, 8, ... reflect that each independent
  3-cycle DOUBLES the contribution.

  WHY doubling? Because each independent 3-cycle has 2 directed orientations,
  and choosing orientations independently gives 2^k total orientations
  for k independent cycles.
""")

# Verify: alpha_1 = c3 (number of 3-cycles)
n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

print(f"\n  Checking alpha_1 = c3 at n={n}:")
for trial in range(10):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    c3 = count_ck(A, n, 3)

    # Count 3-cycles (vertex sets)
    cycle_sets = set()
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                verts = [i, j, k]
                for perm in permutations(verts):
                    if A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and A[perm[2]][perm[0]]:
                        cycle_sets.add(frozenset(verts))
                        break

    alpha_1 = len(cycle_sets)
    print(f"    c3={c3}, alpha_1={alpha_1}, match: {c3 == alpha_1}")

# ============================================================
# PART C: The 3-adic structure and total_lambda = 3*c3
# ============================================================
print("\n" + "=" * 70)
print("PART C: Why total_lambda = 3 * c3")
print("=" * 70)

print("""
total_lambda = sum_{u<v} lambda(u,v) = sum_{u<v} #{3-cycles containing u,v}

Each 3-cycle {a,b,c} contributes 1 to lambda(a,b), lambda(a,c), lambda(b,c).
That's 3 pairs. So total_lambda = 3 * c3.

The 3 comes from C(3,2) = 3: each 3-vertex cycle has 3 edges (pairs).

Similarly:
  total_sigma = sum_{u<v} sigma(u,v)
  total_delta = sum_{u<v} delta(u,v)

  total = C(n,2) * (n-2) (since each pair has n-2 witnesses)

  total_lambda = 3 * c3
  total_sigma = sum(s_i^2) - n(n-1)/2 = sum(s_i^2) - C(n,2) (THM-179)
  total_delta = C(n,2)(n-2) - total_sigma - total_lambda
""")

# The factor 3 is C(3,2). What about 5-cycles?
# Each 5-cycle has C(5,2)=10 pairs. But not all pairs are "lambda_5" witnesses.
print(f"\n  For k-cycles:")
for k in [3, 5, 7]:
    pairs = comb(k, 2)
    print(f"    {k}-cycle has {pairs} pairs = C({k},2)")

print("""
  For 3-cycles: each contributes to 3 pairs (factor C(3,2)=3)
  For 5-cycles: each spans 10 pairs but contributes DIFFERENTLY
  For 7-cycles: each spans 21 pairs

  The factor 3 is special because:
  - It's the SMALLEST odd prime
  - 3-cycles are COMPLETE (every pair is an edge of the cycle)
  - For k-cycles with k>3, not every pair of vertices is adjacent in the cycle
  - In a 3-cycle, ALL C(3,2)=3 pairs are edges; in a 5-cycle, only 5 of C(5,2)=10 are
""")

# ============================================================
# PART D: H = 1 + 2*c3 + ... (the leading terms)
# ============================================================
print("=" * 70)
print("PART D: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3")
print("        = 1 + 2*c3 + 4*i2 + 8*i3")
print("=" * 70)

# i2 = # independent pairs of 3-cycles
# i3 = # independent triples of 3-cycles

print(f"\n  At n=5:")
n = 5
total_bits = n*(n-1)//2

H_formula_check = 0
H_mismatch = 0

for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Find 3-cycles
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                for perm in permutations([i,j,k]):
                    if A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and A[perm[2]][perm[0]]:
                        cycles.append(frozenset([i,j,k]))
                        break

    alpha_1 = len(cycles)  # c3

    # Independent pairs
    alpha_2 = 0
    for a in range(len(cycles)):
        for b in range(a+1, len(cycles)):
            if not (cycles[a] & cycles[b]):
                alpha_2 += 1

    H_formula = 1 + 2*alpha_1 + 4*alpha_2

    if H != H_formula:
        H_mismatch += 1
    H_formula_check += 1

print(f"    H = 1 + 2*c3 + 4*i2 matches: {H_formula_check - H_mismatch}/{H_formula_check}")
print(f"    (At n=5, max independence number of Omega = 2, so alpha_3=0)")

print(f"\n  At n=6:")
n = 6
total_bits = n*(n-1)//2

H_formula_check = 0
H_mismatch = 0
alpha3_count = 0

for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                for perm in permutations([i,j,k]):
                    if A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and A[perm[2]][perm[0]]:
                        cycles.append(frozenset([i,j,k]))
                        break

    alpha_1 = len(cycles)

    # Build adjacency for Omega
    nc = len(cycles)

    # Independent sets of size 2
    alpha_2 = 0
    indep_pairs = []
    for a in range(nc):
        for b in range(a+1, nc):
            if not (cycles[a] & cycles[b]):
                alpha_2 += 1
                indep_pairs.append((a,b))

    # Independent sets of size 3
    alpha_3 = 0
    for a in range(nc):
        for b in range(a+1, nc):
            if cycles[a] & cycles[b]:
                continue
            for c_idx in range(b+1, nc):
                if not (cycles[a] & cycles[c_idx]) and not (cycles[b] & cycles[c_idx]):
                    alpha_3 += 1

    if alpha_3 > 0:
        alpha3_count += 1

    H_formula = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3

    if H != H_formula:
        H_mismatch += 1
    H_formula_check += 1

print(f"    H = 1 + 2*c3 + 4*i2 + 8*i3 matches: {H_formula_check - H_mismatch}/{H_formula_check}")
print(f"    Tournaments with alpha_3 > 0: {alpha3_count}/{H_formula_check}")

# ============================================================
# PART E: The "triangle" 2-3 connection
# ============================================================
print("\n" + "=" * 70)
print("PART E: The triangle connection — 2 and 3 in one formula")
print("=" * 70)

print("""
H(T) = I(Omega(T), 2) = 1 + 2*c3 + 4*i2 + 8*i3 + ...

The COMPLETE formula for H in terms of 3-cycles:
  H = sum_{k=0}^{alpha(Omega)} 2^k * alpha_k

where alpha_k = # independent sets of k 3-cycles.

The numbers 2 and 3 appear simultaneously:
  - 2 is the evaluation point (each cycle doubles contribution)
  - 3 is the cycle length (the atoms are 3-cycles)

The pair (2,3) is the unique pair (p, p+1) where:
  - p is the smallest prime
  - p+1 is the smallest odd prime
  - Actually (2,3) are the ONLY consecutive primes!

In the H formula:
  - The 2^k comes from 2 orientations per cycle (binary choice)
  - The 3-cycles come from the minimal odd cycle structure
  - Claim A has factor 2: H(T)-H(T\v) = 2 * sum mu(C)

  The deletion formula:
    H(T) - H(T\v) = 2 * sum_{C: 3-cycle through v} mu(C)

  where mu(C) = I(Omega(T\v) restricted to avoid C's vertices, 2)
""")

# ============================================================
# PART F: What is I(Omega, 3)?
# ============================================================
print("=" * 70)
print("PART F: I(Omega, 3) — the '3-adic analog'")
print("=" * 70)

# At x=3, each independent set of k cycles contributes 3^k
# This is like "3 orientations" — but cycles only have 2!
# So I(Omega, 3) overcounts by a factor related to 3/2

n = 5
total_bits = n*(n-1)//2

print(f"\n  At n={n}:")
io2_vals = []
io3_vals = []
h_vals = []

for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                for perm in permutations([i,j,k]):
                    if A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and A[perm[2]][perm[0]]:
                        cycles.append(frozenset([i,j,k]))
                        break

    c3 = len(cycles)

    # Independent pairs
    i2 = 0
    for a in range(len(cycles)):
        for b in range(a+1, len(cycles)):
            if not (cycles[a] & cycles[b]):
                i2 += 1

    io2 = 1 + 2*c3 + 4*i2
    io3 = 1 + 3*c3 + 9*i2

    io2_vals.append(io2)
    io3_vals.append(io3)
    h_vals.append(H)

# Check I(Omega, 2) = H
match2 = sum(1 for a, b in zip(io2_vals, h_vals) if a == b)
print(f"    I(Omega, 2) = H: {match2}/{len(h_vals)}")

# I(Omega, 3) / H ratio
ratios = [io3/h for io3, h in zip(io3_vals, h_vals)]
print(f"    I(Omega, 3)/H: min={min(ratios):.3f}, max={max(ratios):.3f}, mean={np.mean(ratios):.3f}")

# I(Omega, 3) - I(Omega, 2) = c3 + 5*i2
diffs = [io3 - io2 for io3, io2 in zip(io3_vals, io2_vals)]
print(f"    I(Omega,3) - I(Omega,2) = c3 + 5*i2:")
for bits in range(min(10, 2**total_bits)):
    A = bits_to_adj(bits, n)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                for perm in permutations([i,j,k]):
                    if A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and A[perm[2]][perm[0]]:
                        cycles.append(frozenset([i,j,k]))
                        break
    c3 = len(cycles)
    i2 = 0
    for a in range(len(cycles)):
        for b in range(a+1, len(cycles)):
            if not (cycles[a] & cycles[b]):
                i2 += 1
    diff = (1+3*c3+9*i2) - (1+2*c3+4*i2)
    print(f"      c3={c3}, i2={i2}, diff={diff}, c3+5*i2={c3+5*i2}, match={diff==c3+5*i2}")

# ============================================================
# PART G: The derivative I'(Omega, x) at x=2
# ============================================================
print("\n" + "=" * 70)
print("PART G: I'(Omega, 2) = c3 + 4*i2 + 12*i3")
print("=" * 70)

print("""
I(Omega, x) = 1 + c3*x + i2*x^2 + i3*x^3 + ...
I'(Omega, x) = c3 + 2*i2*x + 3*i3*x^2 + ...
I'(Omega, 2) = c3 + 4*i2 + 12*i3 + ...

This is the "sensitivity" of H to the evaluation point.
Note: I'(Omega, 2) = d/dx[H] — how fast H grows if we change x from 2.

The ratio I'(Omega,2) / I(Omega,2) = (c3 + 4*i2 + 12*i3)/(1 + 2*c3 + 4*i2 + 8*i3)
""")

n = 6
total_bits = n*(n-1)//2

deriv_ratios = []
for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                for perm in permutations([i,j,k]):
                    if A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and A[perm[2]][perm[0]]:
                        cycles.append(frozenset([i,j,k]))
                        break

    nc = len(cycles)
    c3 = nc

    # Independent pairs and triples
    i2 = 0
    for a in range(nc):
        for b in range(a+1, nc):
            if not (cycles[a] & cycles[b]):
                i2 += 1

    i3 = 0
    for a in range(nc):
        for b in range(a+1, nc):
            if cycles[a] & cycles[b]: continue
            for c_idx in range(b+1, nc):
                if not (cycles[a] & cycles[c_idx]) and not (cycles[b] & cycles[c_idx]):
                    i3 += 1

    Iprime = c3 + 4*i2 + 12*i3
    if H > 0:
        deriv_ratios.append(Iprime / H)

print(f"  At n={n}:")
print(f"    I'(Omega,2)/H: min={min(deriv_ratios):.4f}, max={max(deriv_ratios):.4f}, mean={np.mean(deriv_ratios):.4f}")

# ============================================================
# PART H: The "2-adic logarithm" of H
# ============================================================
print("\n" + "=" * 70)
print("PART H: H = 2^0 + 2^1*c3 + 2^2*i2 + 2^3*i3 — binary decomposition")
print("=" * 70)

print("""
H in binary: 1 + 2*c3 + 4*i2 + 8*i3 + ...

The binary representation of H is:
  bit 0 (units): always 1 (from the empty set — no cycles)
  bit 1 (twos):  determined by c3 mod 2 (parity of 3-cycle count)
  bit 2 (fours): determined by i2 mod 2 (parity of independent pairs)
  ...

Actually this isn't quite right because carries can propagate.
Let's check the actual binary structure.
""")

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

bit_positions = defaultdict(Counter)

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    for bit_pos in range(10):
        bit_val = (H >> bit_pos) & 1
        bit_positions[bit_pos][bit_val] += 1

print(f"  At n={n}, binary structure of H:")
for pos in range(10):
    c = bit_positions[pos]
    p0 = c[0] / (c[0] + c[1])
    p1 = c[1] / (c[0] + c[1])
    print(f"    bit {pos} (2^{pos}={2**pos}): P(0)={p0:.3f}, P(1)={p1:.3f}")

# ============================================================
# PART I: H mod 2^k for small k
# ============================================================
print("\n" + "=" * 70)
print("PART I: H mod 2^k — the 2-adic expansion")
print("=" * 70)

for n in [5, 6, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    print(f"\n  n={n}:")

    for k in [1, 2, 3, 4, 5, 6]:
        mod = 2**k
        H_mod = Counter()

        n_samples = min(2**total_bits, 3000)
        for trial in range(n_samples):
            if n <= 5 and trial < 2**total_bits:
                bits = trial
            else:
                bits = np.random.randint(0, 1 << total_bits)

            A = bits_to_adj(bits, n)
            H = count_ham_paths(A, n)
            H_mod[H % mod] += 1

        n_residues = len(H_mod)
        odd_only = all(r % 2 == 1 for r in H_mod.keys())
        print(f"    H mod {mod}: {n_residues} residues, odd-only: {odd_only}, "
              f"distribution: {dict(sorted(H_mod.items()))}")

# ============================================================
# PART J: The magical relationship: 2*3 = 6 and Rédei
# ============================================================
print("\n" + "=" * 70)
print("PART J: Rédei structure mod 6")
print("=" * 70)

print("""
Rédei: H(T) is always ODD. So H mod 6 in {1, 3, 5}.
This gives H mod 6 = H mod 3 (since H is odd, H mod 6 determines H mod 3 and vice versa).

H mod 2 = 1 (Rédei — the number 2)
H mod 3 varies (the number 3)

For n >= 6, THM-085 says 9 | F(T, omega) where omega = e^{2pi*i/3}.
But this constrains F(T, omega), not H(T) = F(T, 1).

However, F(T, x) mod 3 has special structure (THM-086):
  F(T, x) mod 3 = alpha * (x-1)^{n-1} mod 3

At x=1: F(T,1) = H = alpha * 0 = 0 mod 3... NO that gives H=0 mod 3.
But H is not always 0 mod 3! Let me check.

Actually, (x-1)^{n-1} at x=1 is 0 for n>=2. So we need the Taylor expansion.
F(T,x) = sum c_j (x-1)^j. THM-086 says c_j = 0 mod 3 for j < val(n).
The leading surviving coefficient is c_0 = F(T,1) = H(T).
Wait, c_0 = F(T,1) = H(T). THM-086 says c_j = 0 for j < val(n).
val(n) = 2*floor((n-1)/2). For n=7: val(7) = 6.

So c_0 = H(T) and THM-086 says c_j = 0 for j < 6.
But c_0 = H(T) and 0 < 6, so THM-086 says 3 | H(T) for n >= 7?!

Let me re-read THM-086...
""")

# Quick check: is 3 | H(T) for all n=7 tournaments?
n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

h_mod3 = Counter()
for trial in range(3000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    h_mod3[H % 3] += 1

print(f"\n  H mod 3 at n=7: {dict(sorted(h_mod3.items()))}")
print(f"  So 3 does NOT always divide H at n=7!")
print(f"  This means c_0 = H(T) is NOT always 0 mod 3")
print(f"  THM-086 says c_j = 0 mod 3 for j < val(n), but j starts at 1, not 0!")
print(f"  c_0 = H(T) is the SPECIAL coefficient not covered by the universal zeros")

# What IS H mod 3?
print(f"\n  H mod 3 distribution at various n:")
for n in [3,4,5,6,7,8]:
    total_bits = n*(n-1)//2
    np.random.seed(42)
    h_mod3 = Counter()
    n_samples = min(2**total_bits, 3000)
    for trial in range(n_samples):
        if n <= 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)
        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        h_mod3[H % 3] += 1
    print(f"    n={n}: {dict(sorted(h_mod3.items()))}")

print("\n\nDone.")
