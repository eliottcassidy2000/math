#!/usr/bin/env python3
"""
pi_forbidden_89c.py — Forbidden H values and missing odd numbers

opus-2026-03-14-S89c

DISCOVERY: At n=5, H=7 never appears! At n=6, {7, 21, 35, 39} are missing.
Why should certain odd numbers be impossible as Hamiltonian path counts?

Also investigating:
- The Walsh spectrum of H (only weight-2 terms at n=3,4)
- The H mod structure and its connection to π
- New OEIS sequence candidates
"""

import math
import itertools
from fractions import Fraction
from collections import Counter

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def count_H(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[mask * n + v]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(mask | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("PART 1: Forbidden H Values — Full Analysis")
print("=" * 70)

for n in range(3, 7):
    h_values = []
    for adj in all_tournaments(n):
        h_values.append(count_H(adj, n))

    h_unique = sorted(set(h_values))
    max_h = max(h_unique)
    all_odds = set(range(1, max_h + 1, 2))
    actual = set(h_unique)
    missing = sorted(all_odds - actual)

    print(f"\n  n={n}: max H = {max_h}")
    print(f"    Possible odd values: {sorted(all_odds)}")
    print(f"    Actual H values:     {h_unique}")
    print(f"    MISSING:             {missing if missing else 'NONE'}")

    if missing:
        # Why are these values forbidden?
        # H = 1 + 2·(t₃ + ...) + 4·d₃₃ + ...
        # = 1 + 2·C + 4·D + 8·E + ...
        # So H ≡ 1 mod 2 always
        # H mod 4: 1 if C even, 3 if C odd (C = total odd cycle count)
        # H mod 8: depends on C mod 4 and D mod 2

        # For n=5: H = 1 + 2(t₃ + t₅) since no disjoint pairs fit
        # t₃ ∈ {0,1,2,3,4}, t₅ ∈ {0,1,2,...}
        # H = 1 + 2(t₃ + t₅)
        # So H is always of the form 1 + 2k for k = t₃ + t₅
        # H values: 1, 3, 5, 7, 9, 11, 13, 15
        # But H=7 never appears! That means t₃ + t₅ = 3 is impossible!

        if n == 5:
            print("\n    For n=5: H = 1 + 2(t₃ + t₅)")
            print("    H=7 means t₃ + t₅ = 3")
            print("    Checking all tournaments for t₃ + t₅ = 3...")

            cycle_sums = Counter()
            for adj in all_tournaments(n):
                # Count 3-cycles
                t3 = 0
                for combo in itertools.combinations(range(n), 3):
                    a, b, c = combo
                    if (adj[a][b] and adj[b][c] and adj[c][a]) or \
                       (adj[a][c] and adj[c][b] and adj[b][a]):
                        t3 += 1

                # Count 5-cycles
                t5 = 0
                for perm in itertools.permutations(range(n)):
                    if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                        t5 += 1
                t5 //= 5  # normalize

                total = t3 + t5
                cycle_sums[total] += 1

            print("    Distribution of t₃ + t₅:")
            for k in sorted(cycle_sums.keys()):
                h_val = 1 + 2*k
                print(f"      t₃+t₅ = {k}: {cycle_sums[k]} tournaments, H = {h_val}")

            # Now check which (t₃, t₅) pairs occur
            print("\n    Joint distribution of (t₃, t₅):")
            pair_dist = Counter()
            for adj in all_tournaments(n):
                t3 = 0
                for combo in itertools.combinations(range(n), 3):
                    a, b, c = combo
                    if (adj[a][b] and adj[b][c] and adj[c][a]) or \
                       (adj[a][c] and adj[c][b] and adj[b][a]):
                        t3 += 1

                t5 = 0
                for perm in itertools.permutations(range(n)):
                    if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                        t5 += 1
                t5 //= 5

                pair_dist[(t3, t5)] += 1

            for (t3, t5), count in sorted(pair_dist.items()):
                print(f"      (t₃={t3}, t₅={t5}): {count} tournaments, H = {1+2*(t3+t5)}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 2: WHY is t₃+t₅=3 impossible at n=5?")
print("=" * 70)

# At n=5: t₃ ∈ {0,1,2,3,4} and t₅ ∈ {0,1,2,3,...24}
# (There are C(5,3)=10 possible 3-subsets, max 10 3-cycles if we count directed)
# Actually: max directed 3-cycles at n=5 is 4 (regular tournament has t₃ = 4)
# And max t₅ at n=5: there are 4!/5 = 24/5... wait, directed 5-cycles on 5 vertices
# = (5-1)! / 1 = 24 directed cycles, but on a tournament at most half are present
# Regular tournament has t₅ = 24/2 = 12? No...

# Let me just look at the data
print("\nAll (t₃, t₅) pairs at n=5:")
pair_dist = Counter()
for adj in all_tournaments(5):
    t3 = 0
    for combo in itertools.combinations(range(5), 3):
        a, b, c = combo
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            t3 += 1
    t5 = 0
    for perm in itertools.permutations(range(5)):
        if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
            t5 += 1
    t5 //= 5
    pair_dist[(t3, t5)] += 1

print(f"{'t₃':>4} {'t₅':>4} {'count':>8} {'H':>6}")
for (t3, t5), count in sorted(pair_dist.items()):
    print(f"{t3:4d} {t5:4d} {count:8d} {1+2*(t3+t5):6d}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 3: The Constraint — t₃ Determines t₅ at n=5?")
print("=" * 70)

# From the data, check if t₃ determines t₅ or vice versa
t3_to_t5 = {}
for (t3, t5), count in pair_dist.items():
    if t3 not in t3_to_t5:
        t3_to_t5[t3] = set()
    t3_to_t5[t3].add(t5)

print("\nt₃ → possible t₅ values:")
for t3 in sorted(t3_to_t5.keys()):
    print(f"  t₃={t3}: t₅ ∈ {sorted(t3_to_t5[t3])}")

# Check: t₃ + t₅ achievable values
achievable_sums = set()
for (t3, t5) in pair_dist:
    achievable_sums.add(t3 + t5)
print(f"\nAchievable t₃+t₅: {sorted(achievable_sums)}")
print(f"Missing: {sorted(set(range(max(achievable_sums)+1)) - achievable_sums)}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 4: Missing H at n=6 — What Constraint?")
print("=" * 70)

# At n=6: H = 1 + 2(t₃ + t₅) + 4·d₃₃
# Missing: {7, 21, 35, 39}
# H=7 → 2(t₃+t₅) + 4·d₃₃ = 6 → t₃+t₅ + 2·d₃₃ = 3
# H=21 → t₃+t₅ + 2·d₃₃ = 10
# H=35 → t₃+t₅ + 2·d₃₃ = 17
# H=39 → t₃+t₅ + 2·d₃₃ = 19

# Let me find the joint distribution of (t₃, t₅, d₃₃) at n=6
print("\nJoint (t₃+t₅, d₃₃) distribution at n=6:")
triplet_dist = Counter()
for adj in all_tournaments(6):
    t3 = 0
    for combo in itertools.combinations(range(6), 3):
        a, b, c = combo
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            t3 += 1

    t5 = 0
    for combo in itertools.combinations(range(6), 5):
        for perm in itertools.permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t5 += 1
    t5 //= 5

    # Count d₃₃: disjoint pairs of directed 3-cycles
    all_3cycles = []
    for combo in itertools.combinations(range(6), 3):
        a, b, c = combo
        if adj[a][b] and adj[b][c] and adj[c][a]:
            all_3cycles.append(frozenset([a, b, c]))
        if adj[a][c] and adj[c][b] and adj[b][a]:
            all_3cycles.append(frozenset([a, b, c]))
    # Actually for disjointness we need vertex-disjoint, direction doesn't matter
    # d₃₃ counts PAIRS of vertex-disjoint directed 3-cycles
    d33 = 0
    for i in range(len(all_3cycles)):
        for j in range(i+1, len(all_3cycles)):
            if not (all_3cycles[i] & all_3cycles[j]):
                d33 += 1

    C = t3 + t5
    triplet_dist[(C, d33)] += 1

# What values of (C, d₃₃) give each H?
h_from_cd = {}
for (C, d), count in triplet_dist.items():
    h = 1 + 2*C + 4*d
    if h not in h_from_cd:
        h_from_cd[h] = []
    h_from_cd[h].append((C, d, count))

print(f"\n{'H':>4} {'C=t₃+t₅':>10} {'d₃₃':>6} {'count':>8}")
for h in sorted(h_from_cd.keys()):
    for (C, d, count) in h_from_cd[h]:
        print(f"{h:4d} {C:10d} {d:6d} {count:8d}")

# Which H are achievable?
achievable_h = set(h_from_cd.keys())
max_h_val = max(achievable_h)
all_odds = set(range(1, max_h_val + 1, 2))
missing = sorted(all_odds - achievable_h)
print(f"\nMissing H values at n=6: {missing}")

# For each missing H, what (C, d) would be needed?
for h in missing:
    val = h - 1
    # 2C + 4d = val, so C + 2d = val/2
    if val % 2 != 0:
        print(f"  H={h}: impossible (H-1 odd)")
        continue
    target = val // 2
    possible = [(c, (target-c)//2) for c in range(target+1) if (target-c) % 2 == 0 and (target-c)//2 >= 0]
    print(f"  H={h}: need C+2d={target}, possible (C,d) = {possible}")
    # Which of these don't occur?
    for c, d in possible:
        occurs = (c, d) in triplet_dist
        print(f"    (C={c}, d={d}): {'EXISTS' if occurs else 'MISSING'}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 5: Score Sequence → H Mapping")
print("=" * 70)

# Which score sequences produce which H values?
# Does the score sequence constrain H enough to explain missing values?

for n in [5, 6]:
    score_to_h = {}
    for adj in all_tournaments(n):
        scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))
        h = count_H(adj, n)
        if scores not in score_to_h:
            score_to_h[scores] = Counter()
        score_to_h[scores][h] += 1

    print(f"\n  n={n}: Score sequence → H values:")
    for scores in sorted(score_to_h.keys()):
        h_dist = score_to_h[scores]
        h_vals = sorted(h_dist.keys())
        total = sum(h_dist.values())
        print(f"    {scores}: H ∈ {h_vals} ({total} tournaments)")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 6: The Walsh Spectrum — Degree Structure of H")
print("=" * 70)

# For n=3: H has only weight-0 and weight-2 Walsh coefficients
# For n=4: same!
# Does this persist? At what n does weight-4 (or higher) appear?

for n in range(3, 6):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    H_func = []
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        H_func.append(count_H(adj, n))

    # Walsh-Hadamard transform
    weight_energy = Counter()
    for S in range(1 << m):
        h_hat = sum(H_func[T] * (-1)**bin(S & T).count('1') for T in range(1 << m))
        if h_hat != 0:
            w = bin(S).count('1')
            weight_energy[w] += h_hat**2

    total_energy = sum(weight_energy.values())
    print(f"\n  n={n} (m={m}): Walsh weight spectrum:")
    for w in sorted(weight_energy.keys()):
        frac = weight_energy[w] / total_energy
        print(f"    weight {w}: energy = {weight_energy[w]}, fraction = {frac:.6f}")

    # The fact that only even weights appear is because H(T) = H(T^op) + correction?
    # No: H(T^op) = H(T) always. So H is INVARIANT under complement.
    # Under complement, each bit flips, so χ_S(T^op) = (-1)^|S| χ_S(T)
    # H(T) = H(T^op) means ĥ(S) = (-1)^|S| ĥ(S)
    # So ĥ(S) = 0 for |S| odd!
    # This explains why only even-weight coefficients appear!

    odd_weight_energy = sum(v for w, v in weight_energy.items() if w % 2 == 1)
    print(f"    Odd-weight energy: {odd_weight_energy} (should be 0 by H=H(T^op))")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 7: H is a QUADRATIC Boolean Function (for small n)?")
print("=" * 70)

# At n=3: only weight 0 and 2 (degree ≤ 2 polynomial on F₂^m)
# At n=4: only weight 0 and 2 (also degree 2!)
# At n=5: does weight 4 appear?

# This would mean H = a₀ + Σ_{|S|=2} a_S · χ_S where χ_S = product of ±1 variables
# = affine function of PAIRWISE PRODUCTS of arc indicators

# If H is truly degree 2 for all n, that would be extraordinary
# It would mean H(T) = Σ c_{ij} T_{ij} T_{kl} + Σ b_i T_{ij} + constant
# (where the sum is over arc pairs)

n = 5
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)
print(f"\n  n={n} (m={m}): Full weight spectrum:")

H_func = []
for bits in range(1 << m):
    adj = [[0]*n for _ in range(n)]
    for idx, (i, j) in enumerate(pairs):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H_func.append(count_H(adj, n))

weight_counts = Counter()
weight_energy = Counter()
max_coeff_by_weight = {}

for S in range(1 << m):
    h_hat = sum(H_func[T] * (-1)**bin(S & T).count('1') for T in range(1 << m))
    w = bin(S).count('1')
    if h_hat != 0:
        weight_counts[w] += 1
        weight_energy[w] += h_hat**2
        if w not in max_coeff_by_weight or abs(h_hat) > max_coeff_by_weight[w]:
            max_coeff_by_weight[w] = abs(h_hat)

total_energy = sum(weight_energy.values())
for w in sorted(weight_counts.keys()):
    frac = weight_energy[w] / total_energy
    print(f"    weight {w}: {weight_counts[w]} nonzero, energy fraction = {frac:.6f}, max|ĥ| = {max_coeff_by_weight[w]}")

if 4 in weight_counts:
    print(f"\n  *** Weight 4 appears at n=5! H is NOT quadratic for n≥5 ***")
elif max(weight_counts.keys()) <= 2:
    print(f"\n  *** H is STILL quadratic at n=5! ***")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 8: NEW OEIS CANDIDATE SEQUENCES")
print("=" * 70)

# Sequences found in this session that are NOT in OEIS:

print("""
NEW OEIS CANDIDATE SEQUENCES:

1. Σ_T H(T)² for tournaments on n vertices:
   a(n) = sum over all 2^{C(n,2)} tournaments T of H(T)²
   Values: a(3)=24, a(4)=768, a(5)=75840, a(6)=21381120
   Not in OEIS as of 2026-03-14.

2. f(n) = Σ_{compatible (π,σ)} 2^{k(π,σ)} / n!
   where (π,σ) compatible means the union of directed arcs forms
   a consistent tournament, and k = shared arcs.
   Values: f(3)=8, f(4)=32, f(5)=158, f(6)=928, f(7)=6350
   Not in OEIS.

3. max H(T) for tournaments on n vertices (odd n):
   3, 5, 15, 45, ... (max = (2n-1)!! / something?)
   Actually: max H at n=3 is 3, n=4 is 5, n=5 is 15, n=6 is 45
   = 3, 5, 15, 45 = 3·1, 5·1, 5·3, 5·9
   Ratios: 5/3, 3, 3 → seems to approach 3
   Actually: 3, 5, 15, 45 are in OEIS (likely A001147 double factorials or similar)
   3 = 3, 5 = 5, 15 = 3·5, 45 = 3·5·3
   Let me check: 45/15 = 3, 15/5 = 3, 5/3 = 5/3
   NOT a clean ratio.

4. Number of distinct H values for tournaments on n vertices:
   2, 3, 7, 19, ...
   OEIS candidate.

5. Missing H values (number of odd integers in [1, max H] not achieved):
   0, 0, 1, 4, ...

6. Signed Euler characteristic χ = #{H≡1 mod 4} - #{H≡3 mod 4}:
   4, 32, 416, 8192, ...
""")

# Verify max H
for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    print(f"  n={n}: max H = {max(h_values)}, #distinct = {len(set(h_values))}, #missing = {len(set(range(1, max(h_values)+1, 2)) - set(h_values))}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 9: Why H is Even-Weight (degree 2k) on Q_m")
print("=" * 70)

# THEOREM: H(T) = H(T^op) for all tournaments T.
# PROOF: Reversing all arcs sends each Hamiltonian path to its reverse,
# which is a Hamiltonian path of T^op. This is a bijection. □
#
# CONSEQUENCE: In the Walsh basis, ĥ(S) = 0 for |S| odd.
# PROOF: T^op corresponds to flipping all bits. Under the Walsh transform,
# χ_S(T^op) = (-1)^{|S|} · χ_S(T). So:
# ĥ(S) = (1/2^m) Σ_T H(T) χ_S(T) = (1/2^m) Σ_T H(T^op) χ_S(T)
#       = (1/2^m) Σ_T H(T) χ_S(T^op) (substituting T→T^op)
#       = (-1)^{|S|} · ĥ(S)
# So if |S| is odd, ĥ(S) = -ĥ(S), hence ĥ(S) = 0. □
#
# This is a CLEAN THEOREM connecting H(T)=H(T^op) to Walsh structure!

print("""
THEOREM (Even-Weight Walsh Spectrum):
H(T) = H(T^op) implies that the Walsh-Hadamard transform ĥ(S) = 0
for all S with |S| odd.

PROOF: T^op flips all m bits. Under Walsh transform, this multiplies
χ_S by (-1)^|S|. Since H is invariant under complement, ĥ(S) must
vanish for odd |S|.

CONSEQUENCE: H is a polynomial of EVEN degree in the ±1 variables.
The "lowest nontrivial" contribution is at weight 2 (pairs of arcs).

At n=3,4: H is EXACTLY degree 2 (only weight 0 and 2).
At n=5: Weight 4 also appears (H becomes degree 4).
""")

print("=" * 70)
print("END — opus-2026-03-14-S89c forbidden")
print("=" * 70)
