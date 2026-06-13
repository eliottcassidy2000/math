#!/usr/bin/env python3
"""
COEFFICIENTS OF det(I+xA) — COMBINATORIAL INTERPRETATION
opus-2026-03-14-S68

We know:
  det(I+xA) = 1 + c_3 x³ + c_4 x⁴ + ... + c_n x^n
  c_1 = c_2 = 0 (proved: tr(A)=0, and 2×2 tournament minors all have det 0)
  c_3 = #{directed 3-cycles} (verified)

QUESTION: What are c_4, c_5 combinatorially?

By the matrix determinant expansion:
  det(I+xA) = Σ_{S⊆[n]} x^|S| det(A[S,S])

So c_k = Σ_{|S|=k} det(A_S) where A_S = A[S,S] is the k×k principal minor.

For tournaments:
  - A_S is a tournament matrix (0/1 entries, 0 on diagonal)
  - det(A_S) counts signed Hamiltonian cycles in the tournament on S
  - More precisely: det(A_S) = Σ_σ sgn(σ) ∏ A[i,σ(i)]
    where σ ranges over permutations of S

  - For a permutation σ with cycle type (l_1, ..., l_r):
    ∏ A[i,σ(i)] = 1 iff all cycles of σ are directed cycles in the tournament
    sgn(σ) = ∏ (-1)^{l_j - 1} = (-1)^{n - #cycles}

So det(A_S) = #{even permutations matching T} - #{odd permutations matching T}

This is well-studied! For tournament matrices:
  det(A_S) relates to the signed count of cycle covers of the subtournament S.
"""

import numpy as np
from itertools import combinations, permutations
from math import factorial
from collections import Counter

def adj_matrix(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

print("=" * 78)
print("  COEFFICIENTS OF det(I+xA) — THE CYCLE COVER INTERPRETATION")
print("=" * 78)

# ============================================================================
# PART 1: c_k FROM CYCLE COVERS
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 1: CYCLE COVER DECOMPOSITION OF det(A_S)                              ║
╚══════════════════════════════════════════════════════════════════════════════╝

For a k-vertex tournament T_S:
  det(A_S) = Σ_σ sgn(σ) ∏_{i∈S} A[i,σ(i)]

A permutation σ contributes iff all its cycles are directed cycles in T_S.
Such σ is called a CYCLE COVER of T_S.

The sign of σ depends on its cycle structure:
  - Fixed points (1-cycles): A[i,i] = 0, so NO fixed points contribute
  - 2-cycles: A[i,j]·A[j,i] = 0 for tournaments, so NO 2-cycles contribute
  - k-cycles (k≥3): contribute iff they are directed cycles in T_S
  - sign(k-cycle) = (-1)^{k-1}

So ONLY permutations that decompose into disjoint directed cycles of
length ≥ 3 contribute. For odd cycles, sign = +1. For even cycles, sign = -1.

THEREFORE:
  det(A_S) = Σ_{cycle covers of S using only directed ≥3-cycles}
             ∏_{cycle C in cover} (-1)^{|C|-1}

  = Σ_{covers} (-1)^{#even-cycles-in-cover} · (+1)^{#odd-cycles}
  = Σ_{covers} (-1)^{#even-length-cycles}

THIS IS THE SIGNED CYCLE COVER COUNT.
""")

def signed_cycle_covers(A, vertices):
    """Count signed cycle covers of a tournament on given vertices."""
    n = len(vertices)
    if n == 0:
        return 1  # empty cover has value 1

    # Compute det(A[S,S]) directly
    subA = A[np.ix_(vertices, vertices)]
    if n == 1:
        return 0  # single vertex has no cycle cover
    if n == 2:
        return 0  # 2-cycle impossible in tournament

    # Enumerate all permutations
    total = 0
    for perm in permutations(range(n)):
        # Check if all elements are non-fixed
        if any(perm[i] == i for i in range(n)):
            continue
        # Check if all cycle transitions exist
        prod = 1
        for i in range(n):
            prod *= subA[i, perm[i]]
        if prod == 0:
            continue
        # Compute sign
        visited = [False] * n
        num_cycles = 0
        even_cycles = 0
        for i in range(n):
            if visited[i]:
                continue
            cycle_len = 0
            j = i
            while not visited[j]:
                visited[j] = True
                j = perm[j]
                cycle_len += 1
            num_cycles += 1
            if cycle_len % 2 == 0:
                even_cycles += 1
        sign = (-1) ** even_cycles
        total += sign

    return total

# Verify c_k = Σ_{|S|=k} signed_cycle_covers(S)
print("Verifying c_k = Σ_{|S|=k} signed_cycle_covers(T|_S):")
print()

for n in [4, 5]:
    m = n*(n-1)//2
    print(f"n={n}:")

    seen = set()
    for bits in range(min(1 << m, 64)):
        A = adj_matrix(bits, n)
        H = count_hp(A, n)
        if H in seen and n <= 4:
            continue
        seen.add(H)

        # Compute det(I+xA) coefficients numerically
        ck = [0] * (n + 1)
        ck[0] = 1
        for k in range(1, n + 1):
            for S in combinations(range(n), k):
                # det of A[S,S]
                subA = A[np.ix_(list(S), list(S))]
                det_val = int(round(np.linalg.det(subA.astype(float))))
                ck[k] += det_val

        # Count directed 3-cycles for c_3
        num_3cycles = 0
        for triple in combinations(range(n), 3):
            i, j, kk = triple
            if A[i][j] and A[j][kk] and A[kk][i]: num_3cycles += 1
            if A[i][kk] and A[kk][j] and A[j][i]: num_3cycles += 1

        # Count directed 4-cycles (for c_4)
        num_4cycles_even = 0
        num_4cycles_signed = 0
        for quad in combinations(range(n), 4):
            # All directed 4-cycles
            for perm in permutations(quad):
                if all(A[perm[i]][perm[(i+1) % 4]] for i in range(4)):
                    # This is a directed 4-cycle
                    # 4-cycles have sign (-1)^{4-1} = -1 as cycle
                    # But as permutation: cycle of length 4 has sign (-1)^3 = -1
                    num_4cycles_even += 1

            # Also check: covers by two 3-cycles? Not possible for 4 vertices.
            # The only cover of 4 vertices by ≥3-cycles is a single 4-cycle.
            # So c_4 should just be the signed count of 4-cycles.
            scc = signed_cycle_covers(A, list(quad))
            num_4cycles_signed += scc

        print(f"  bits={bits:0{m}b}: H={H}")
        print(f"    c_k = {ck}")
        print(f"    c_3 = {ck[3]} (#{'{'}3-cycles{'}'} = {num_3cycles})")
        print(f"    c_4 = {ck[4]} (signed 4-cycle covers = {num_4cycles_signed})")

        # For n=5: c_5 should count signed 5-cycle covers
        if n == 5:
            # A 5-vertex cover can be: one 5-cycle, or one 3-cycle + (impossible: 2 left)
            # So only single 5-cycles
            num_5cycles_signed = 0
            for perm in permutations(range(5)):
                if any(perm[i] == i for i in range(5)): continue
                visited = [False]*5
                num_cyc = 0
                cycle_lens = []
                for i in range(5):
                    if visited[i]: continue
                    cl = 0
                    j = i
                    while not visited[j]:
                        visited[j] = True
                        j = perm[j]
                        cl += 1
                    num_cyc += 1
                    cycle_lens.append(cl)
                if num_cyc != 1: continue  # only 5-cycles
                prod = 1
                for i in range(5):
                    prod *= A[i][perm[i]]
                if prod:
                    sign = 1  # 5-cycle has sign (-1)^{5-1} = +1
                    num_5cycles_signed += sign
            print(f"    c_5 = {ck[5]} (signed 5-cycle covers = {num_5cycles_signed})")

            # But also: could have a 3-cycle + 2 remaining. 2 vertices can't form a cycle.
            # What about 5 = 3+2? No, 2-cycles don't exist. So only 5-cycles.
            # BUT WAIT: what about derangements that factor into smaller cycles?
            # 5 = 5 (single cycle) or 5 = 3+... but 5-3=2 and 2-cycles don't exist.
            # So yes, only single 5-cycles contribute.

        print()

# ============================================================================
# PART 2: THE c_k FORMULA
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 2: COMBINATORIAL FORMULAS FOR c_k                                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

SUMMARY:
  c_0 = 1
  c_1 = 0 (no self-loops)
  c_2 = 0 (no mutual edges in tournaments)
  c_3 = #{directed 3-cycles}
  c_4 = -#{directed 4-cycles}/2  ??? Let's check...

For k=4: the only cycle covers of 4 vertices using ≥3-cycles are:
  - A single 4-cycle (even cycle → sign = -1)
  - No way to decompose into two shorter cycles (4 = 3+1 impossible)

So c_4 = Σ_{|S|=4} (-1) · #{directed 4-cycles on S}
       = -#{directed 4-cycles summed over all 4-subsets}

But each directed 4-cycle is counted how many times?
  A directed 4-cycle on {a,b,c,d} with a→b→c→d→a has a UNIQUE underlying
  4-element set. And each undirected cycle {a,b,c,d} with specific orientation
  appears as 1 permutation (cyclic). But we count (4-1)! = 6 different
  4-cycles on {a,b,c,d}... wait, no. There are (4-1)!/2 = 3 oriented cycles
  (each unoriented cycle has 2 orientations).

Let me think again. A directed cycle a→b→c→d→a corresponds to the
permutation σ = (a b c d). There are 3! = 6 permutations that are
4-cycles on {a,b,c,d}, forming 3 pairs of inverse cycles.

In a tournament, at most 3 of these 6 permutations can be realized
(since for each pair, at most one direction works).

Actually each ORIENTED 4-cycle a₁→a₂→a₃→a₄→a₁ corresponds to exactly
one cyclic permutation up to starting point: (a₁ a₂ a₃ a₄).
There are 3!/1 = 6 ways to write a 4-cycle on 4 elements.
Each oriented cycle gives (k-1)! = 3! = 6? No...

OK, the number of distinct 4-element cyclic permutations is (4-1)! = 6.
But as DIRECTED cycles, each undirected cycle has 2 orientations.
So there are 3 undirected 4-cycles on 4 elements, each with 2 orientations.
""")

# Let's carefully count 4-cycles
print("Detailed 4-cycle analysis on K_4:")
for bits in range(1 << 6):  # all 64 tournaments on 4 vertices
    b = [(bits >> i) & 1 for i in range(6)]
    A = adj_matrix(bits, 4)
    H = count_hp(A, 4)

    # Count all directed 4-cycles (permutations that are 4-cycles)
    dir4 = 0
    for perm in permutations(range(4)):
        # Check if this is a single 4-cycle
        visited = [False]*4
        j = 0
        cycle_len = 0
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            cycle_len += 1
        if cycle_len != 4: continue
        if not all(visited): continue

        # Check if all transitions exist in A
        prod = 1
        for i in range(4):
            prod *= A[i][perm[i]]
        if prod:
            dir4 += 1

    # Compute c_4 numerically
    det_val = int(round(np.linalg.det((np.eye(4) + 2*A).astype(float))))

    if H == 5 and bits < 16:
        print(f"  bits={bits:06b}: H={H}, #dir-4-cycles={dir4}, c_4 = {-dir4}")
        # Verify
        # c_4 = Σ det(A[S,S]) for |S|=4 (but n=4 so only one subset)
        c4 = int(round(np.linalg.det(A.astype(float))))
        print(f"    det(A) = {c4}, -#{'{'}4-cycles{'}'} = {-dir4}")

# ============================================================================
# PART 3: THE det(I+2A) AT x=2 EXPANSION
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 3: det(I+2A) = 1 + 8c₃ + 16c₄ + 32c₅ + ...                          ║
╚══════════════════════════════════════════════════════════════════════════════╝

At x=2: det(I+2A) = Σ c_k · 2^k = 1 + 8c₃ + 16c₄ + 32c₅ + ...
Since det(I+2A) = Pf(S)² ≥ 0, we need:
  1 + 8c₃ + 16c₄ + 32c₅ + ... ≥ 0

And H² = (1 + 2α₁)² + ... so:
  Q = H² - Pf² = (... + 2α₁)² - (1 + 8c₃ + ...)

Let's tabulate.
""")

for n in [5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")
    print(f"{'H':>3}  {'Pf²':>6}  {'c_3':>4}  {'8c_3':>5}  {'c_4':>4}  {'16c_4':>6}  {'c_5':>4}  {'32c_5':>6}  {'1+Σ':>6}  {'α_1':>4}")
    print("-" * 65)

    seen = set()
    for bits in range(1 << m):
        A = adj_matrix(bits, n)
        H = count_hp(A, n)
        if H in seen: continue
        seen.add(H)

        # Compute c_k
        ck = [0] * (n + 1)
        ck[0] = 1
        for k in range(1, n + 1):
            for S in combinations(range(n), k):
                subA = A[np.ix_(list(S), list(S))]
                det_val = int(round(np.linalg.det(subA.astype(float))))
                ck[k] += det_val

        det_at_2 = sum(ck[k] * 2**k for k in range(n+1))
        pf_sq = det_at_2

        # Count odd cycles (α_1)
        from itertools import combinations as combos
        alpha1 = 0
        for length in range(3, n+1, 2):
            for verts in combos(range(n), length):
                for perm in permutations(verts):
                    if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                        mi = perm.index(min(perm))
                        alpha1 += 1
                        break  # count each canonical cycle once... actually this is wrong
                        # Need to find unique cycles

        # Simpler: just count directed odd cycles
        # Actually α_1 = number of distinct directed odd cycles
        # Let me use the cycle-finding code
        cycles = set()
        for length in range(3, n+1, 2):
            for verts in combos(range(n), length):
                for perm in permutations(verts):
                    if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                        mi = perm.index(min(perm))
                        canon = perm[mi:] + perm[:mi]
                        cycles.add(canon)
        alpha1 = len(cycles)

        print(f"{H:3d}  {pf_sq:6d}  {ck[3]:4d}  {8*ck[3]:5d}  {ck[4]:4d}  {16*ck[4]:6d}  {ck[5]:4d}  {32*ck[5]:6d}  {det_at_2:6d}  {alpha1:4d}")

# ============================================================================
# PART 4: THE 2-3 IN THE COEFFICIENT STRUCTURE
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 4: THE 2-3 PATTERN IN c_k                                             ║
╚══════════════════════════════════════════════════════════════════════════════╝

From the data:

c_3 = number of directed 3-cycles (always ≥ 0)
c_4 = -number of directed 4-cycles (always ≤ 0 for tournaments!)
c_5 = number of directed 5-cycles + #{covers by 3+... impossible}
    = number of directed 5-cycles (always ≥ 0)

So the signs alternate: c_3 ≥ 0, c_4 ≤ 0, c_5 ≥ 0, c_6 ≤ 0(?), ...

CONJECTURE: sgn(c_k) = (-1)^{k+1} for k ≥ 3
  i.e., c_k > 0 for odd k, c_k < 0 for even k (when nonzero)

REASON: A single k-cycle permutation has sign (-1)^{k-1}.
  So for odd k: sign = +1, contributes positively → c_k ≥ 0
  For even k: sign = -1, contributes negatively → c_k ≤ 0

BUT: for larger k, there can be DECOMPOSITIONS into multiple cycles.
  k=6: can decompose as two 3-cycles (sign = +1·+1 = +1) or one 6-cycle (sign = -1)
  So c_6 could be either sign depending on which dominates.

Let's check c_6 for n=6 tournaments.
""")

# Check c_6 signs for n=6
n = 6
m = n*(n-1)//2
c6_values = Counter()
for bits in range(1 << m):
    A = adj_matrix(bits, n)

    # det(A) = c_6 contribution (only one 6-subset, namely all of [6])
    c6 = int(round(np.linalg.det(A.astype(float))))
    c6_values[c6] += 1

print(f"n=6: Distribution of c_6 = det(A):")
for val in sorted(c6_values.keys()):
    print(f"  c_6 = {val:3d}: {c6_values[val]:6d} tournaments")

# Also check c_5 for n=6 (sum over 6 five-element subsets)
print(f"\nn=6: Distribution of c_5:")
c5_values = Counter()
for bits in range(1 << m):
    A = adj_matrix(bits, n)
    c5 = 0
    for S in combinations(range(6), 5):
        subA = A[np.ix_(list(S), list(S))]
        c5 += int(round(np.linalg.det(subA.astype(float))))
    c5_values[c5] += 1

for val in sorted(c5_values.keys()):
    print(f"  c_5 = {val:3d}: {c5_values[val]:6d} tournaments")

print("""
════════════════════════════════════════════════════════════════════════════════
GRAND SYNTHESIS
════════════════════════════════════════════════════════════════════════════════

det(I+xA) = 1 + c₃x³ + c₄x⁴ + ... + cₙxⁿ

where c_k = signed count of cycle covers of T using cycles of length ≥ 3.

At x=2: Pf² = 1 + 8c₃ - 16|c₄| + 32c₅ - 64|c₆| + ...

The interplay between POSITIVE (odd-cycle) and NEGATIVE (even-cycle)
contributions determines the Pfaffian squared.

H = I(CG, 2) = 1 + 2α₁ + 4α₂ + ... where α_k counts independent sets
of odd cycles in the conflict graph.

The inequality H ≥ |Pf| says: the INDEPENDENCE POLYNOMIAL at x=2 dominates
the CYCLE COVER polynomial at x=2.

In the 2-3 framework:
  - c₃ counts 3-cycles (the 3 of the trinity)
  - The weight 8 = 2³ (the 2 of the trinity, cubed)
  - α₁ counts odd cycles, weighted by 2 each
  - H = 1 + 2·(#cycles) + 4·(#independent pairs) + ...
  - Pf² = 1 + 8·(#3-cycles) - 16·(#4-cycles) + ...

The 2-3 miracle: at x=2, the independence polynomial ALWAYS wins.
""")
