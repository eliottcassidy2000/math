#!/usr/bin/env python3
"""
PROJECTIVE & ALGEBRAIC GEOMETRY — PART 4: THE BINARY COEFFICIENT THEOREM
opus-2026-03-14-S71n

EXPLOSIVE FINDINGS FROM PART 3:
1. ALL multilinear coefficients of H are +-2^k
2. Score sequence determines H for n<=4 (constant on strata!)
3. H(T) = H(T^op) for ALL T when n is odd

This script investigates these three findings deeply.

QUESTIONS:
A. WHY are all coefficients powers of 2? Is this related to the
   Walsh basis being over F_2?
B. What exactly causes the score-H correspondence to break at n=5?
   Only stratum (1,2,2,2,3) has variance — WHAT INVARIANT splits it?
C. Is there a CLOSED FORM for the multilinear coefficients?
D. Does the +-2^k property extend to M[a,b] and other functions?
E. What is the 2-adic structure of H?
"""

from itertools import combinations, permutations
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

def score_seq(n, A):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

def get_all_tournaments(n):
    m = n*(n-1)//2
    results = {}
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        results[bits] = H
    return results

def count_3cycles(n, A):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] + A[j][k] + A[k][i] == 3:
                    count += 1
                if A[i][k] + A[k][j] + A[j][i] == 3:
                    count += 1
    return count

print("=" * 70)
print("THE BINARY COEFFICIENT THEOREM AND SCORE-H STRATIFICATION")
print("opus-2026-03-14-S71n")
print("=" * 70)

# ======================================================================
# PART A: MULTILINEAR COEFFICIENTS — STRUCTURE AND PATTERN
# ======================================================================
print(f"\n{'='*70}")
print("PART A: MULTILINEAR COEFFICIENTS — ALL +-2^k")
print(f"{'='*70}")

for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    # Multilinear coefficients via Mobius inversion
    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    # Analyze the coefficients
    print(f"\n  n={n} (m={m}): Multilinear coefficient analysis")

    # Group by degree
    by_degree = defaultdict(list)
    for S, v in ml_coeffs.items():
        deg = bin(S).count('1')
        by_degree[deg].append((S, v))

    for deg in sorted(by_degree.keys()):
        entries = by_degree[deg]
        vals = [v for _, v in entries]
        abs_vals = sorted(set(abs(v) for v in vals))
        # How many positive, negative?
        pos = sum(1 for v in vals if v > 0)
        neg = sum(1 for v in vals if v < 0)
        print(f"    Degree {deg}: {len(entries)} nonzero, values: {abs_vals}, +:{pos} -:{neg}")

        # Are the signs determined by something about the arc set S?
        if deg == 2 and n <= 5:
            # For degree-2 monomials (pair of arcs): what determines the sign?
            print(f"      Degree-2 sign analysis:")
            arc_list = []
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    arc_list.append((i, j))
                    idx += 1

            for S, v in entries[:15]:  # Show first 15
                bits = [i for i in range(m) if S & (1 << i)]
                arc1 = arc_list[bits[0]]
                arc2 = arc_list[bits[1]]
                # Do the arcs share a vertex?
                shared = set(arc1) & set(arc2)
                print(f"        arcs {arc1},{arc2}: coeff={v:+d}, shared_vertices={len(shared)}")

    # 2-adic valuation distribution
    print(f"    2-adic valuation of |coefficients|:")
    v2_dist = Counter()
    for S, v in ml_coeffs.items():
        abs_v = abs(v)
        v2 = 0
        while abs_v % 2 == 0 and abs_v > 0:
            v2 += 1
            abs_v //= 2
        v2_dist[v2] += 1
    for val, count in sorted(v2_dist.items()):
        print(f"      v_2 = {val}: {count} coefficients")

# ======================================================================
# PART B: SIGN PATTERN — WHAT DETERMINES THE SIGN?
# ======================================================================
print(f"\n{'='*70}")
print("PART B: SIGN PATTERN OF DEGREE-2 COEFFICIENTS")
print(f"{'='*70}")

# At n=3: coefficients of x_i*x_j are ±2
# Arc pairs: (01,02), (01,12), (02,12) with coefficients -2, +2, -2
# Pattern: the arcs form a triangle, and the sign alternates?

# Let's think about this more carefully.
# The arc encoding is:
# x_0 = arc (0,1), x_1 = arc (0,2), x_2 = arc (1,2)

# The degree-2 coefficients are:
# x_0*x_1: coeff = H(111) - H(110) - H(011) + H(010)
# where H(abc) means H with x_0=a, x_1=b, x_2=c

# Actually let's check: do adjacent arcs (sharing a vertex) always
# have the SAME sign pattern?

print("\n  Degree-2 coefficient sign vs arc adjacency:")

for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    arc_list = []
    for i in range(n):
        for j in range(i+1, n):
            arc_list.append((i, j))

    # For each pair of arcs, classify by adjacency type
    adj_signs = defaultdict(list)  # (adj_type) -> [signs]

    for S, v in ml_coeffs.items():
        bits = [i for i in range(m) if S & (1 << i)]
        if len(bits) != 2:
            continue
        arc1 = arc_list[bits[0]]
        arc2 = arc_list[bits[1]]
        shared = set(arc1) & set(arc2)
        n_shared = len(shared)

        # More specific: what is the shared vertex position?
        if n_shared == 1:
            sv = shared.pop()
            # Is the shared vertex an endpoint of both arcs?
            # Arc (a,b) with a<b: shared vertex could be a or b
            pos1 = 0 if arc1[0] == sv else 1
            pos2 = 0 if arc2[0] == sv else 1
            adj_type = f"adj({pos1},{pos2})"
        else:
            adj_type = "disjoint"

        adj_signs[adj_type].append(v)

    print(f"\n  n={n}:")
    for atype, vals in sorted(adj_signs.items()):
        pos = sum(1 for v in vals if v > 0)
        neg = sum(1 for v in vals if v < 0)
        print(f"    {atype}: {len(vals)} pairs, +:{pos} -:{neg}, values: {sorted(set(vals))}")

# ======================================================================
# PART C: SCORE-H STRATIFICATION AT n=5 — WHAT SPLITS (1,2,2,2,3)?
# ======================================================================
print(f"\n{'='*70}")
print("PART C: THE (1,2,2,2,3) STRATUM AT n=5")
print(f"{'='*70}")

n = 5
m = n*(n-1)//2
all_T5 = get_all_tournaments(5)

# The only stratum with variance is (1,2,2,2,3)
target_score = (1, 2, 2, 2, 3)
stratum = []
for bits, H in all_T5.items():
    A = adj_matrix(n, bits)
    ss = score_seq(n, A)
    if ss == target_score:
        t3 = count_3cycles(n, A)
        stratum.append((bits, H, t3))

H_counter = Counter(H for _, H, _ in stratum)
t3_counter = Counter((H, t3) for _, H, t3 in stratum)

print(f"\n  Stratum (1,2,2,2,3): {len(stratum)} tournaments")
print(f"  H distribution: {dict(sorted(H_counter.items()))}")
print(f"  (H, t3) distribution:")
for (H, t3), count in sorted(t3_counter.items()):
    print(f"    H={H:3d}, t3={t3}: {count} tournaments")

# So t3 distinguishes within the (1,2,2,2,3) stratum?
# Check: does score + t3 determine H at n=5?
print(f"\n  Score + t3 -> H determination at n=5:")
all_strata = defaultdict(list)
for bits, H in all_T5.items():
    A = adj_matrix(n, bits)
    ss = score_seq(n, A)
    t3 = count_3cycles(n, A)
    all_strata[(ss, t3)].append(H)

all_determined = True
for key, H_list in sorted(all_strata.items()):
    if len(set(H_list)) > 1:
        all_determined = False
        print(f"    {key}: H = {sorted(set(H_list))} (NOT DETERMINED)")

if all_determined:
    print(f"    YES! Score + t3 uniquely determines H at n=5!")

# ======================================================================
# PART D: SCORE-H AT n=6 — HOW BADLY DOES IT BREAK?
# ======================================================================
print(f"\n{'='*70}")
print("PART D: SCORE-H STRATIFICATION AT n=6")
print(f"{'='*70}")

n = 6
m = n*(n-1)//2
print(f"\n  Computing all 2^{m} = {1<<m} tournaments at n={n}...")
all_T6 = get_all_tournaments(6)
print(f"  Done.")

# Group by score sequence
score_strata = defaultdict(list)
for bits, H in all_T6.items():
    A = adj_matrix(n, bits)
    ss = score_seq(n, A)
    score_strata[ss].append(H)

print(f"\n  Score strata at n=6:")
non_constant = 0
for ss in sorted(score_strata.keys()):
    H_list = score_strata[ss]
    H_set = sorted(set(H_list))
    if len(H_set) > 1:
        non_constant += 1
        print(f"    {ss}: |stratum|={len(H_list)}, H_range=[{min(H_list)},{max(H_list)}], "
              f"#distinct={len(H_set)}")
    else:
        print(f"    {ss}: |stratum|={len(H_list)}, H={H_set[0]} (constant)")

print(f"\n  Non-constant strata: {non_constant}/{len(score_strata)}")

# Does score + t3 determine H at n=6?
print(f"\n  Score + t3 -> H determination at n=6:")
score_t3_strata = defaultdict(list)
for bits, H in all_T6.items():
    A = adj_matrix(n, bits)
    ss = score_seq(n, A)
    t3 = count_3cycles(n, A)
    score_t3_strata[(ss, t3)].append(H)

not_determined_count = 0
for key, H_list in sorted(score_t3_strata.items()):
    if len(set(H_list)) > 1:
        not_determined_count += 1
        if not_determined_count <= 5:
            print(f"    {key}: H = {sorted(set(H_list))}")

print(f"  Strata where score+t3 doesn't determine H: {not_determined_count}/{len(score_t3_strata)}")

# ======================================================================
# PART E: THE 2-ADIC STRUCTURE OF H
# ======================================================================
print(f"\n{'='*70}")
print("PART E: 2-ADIC STRUCTURE OF H")
print(f"{'='*70}")

# H is always odd (PROVED: H ≡ 1 mod 2)
# What about mod 4, mod 8?
for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n) if n <= 6 else None

    if all_T:
        print(f"\n  n={n}: H modular structure")
        for p in [2, 4, 8, 16]:
            mod_counts = Counter(H % p for H in all_T.values())
            print(f"    mod {p:2d}: {dict(sorted(mod_counts.items()))}")

        # 2-adic valuation of H-1
        v2_dist = Counter()
        for H in all_T.values():
            val = H - 1
            if val == 0:
                v2_dist['inf'] += 1
                continue
            v2 = 0
            while val % 2 == 0:
                v2 += 1
                val //= 2
            v2_dist[v2] += 1
        print(f"    v_2(H-1): {dict(sorted(v2_dist.items(), key=lambda x: (0 if isinstance(x[0], int) else 1, x[0])))}")

# ======================================================================
# PART F: WHY +-2^k? — THE TRANSFER MATRIX EXPLANATION
# ======================================================================
print(f"\n{'='*70}")
print("PART F: WHY +-2^k? — THE TRANSFER MATRIX EXPLANATION")
print(f"{'='*70}")

print("""
  The multilinear expansion H(x) = sum_S c_S prod_{i in S} x_i
  has c_S = sum_{T subset S} (-1)^{|S|-|T|} H(T)  (Mobius inversion).

  But H(T) counts Hamiltonian paths in a tournament defined by T.
  The key insight: H is computed via a TRANSFER MATRIX of size n x n,
  where each entry is LINEAR in one arc variable x_k.

  The permanent (or path sum) of this matrix is multilinear,
  and the coefficients inherit the 2-adic structure from the
  transfer matrix entries being +-1.

  More precisely: H = sum_P prod_{i=0}^{n-2} A[P_i][P_{i+1}]
  where A[i][j] = x_k (the variable for arc (i,j)) if i<j,
  or 1-x_k if i>j.

  Each product in the sum has n-1 factors, each +-1 or 0.
  The coefficient c_S counts, with signs, the number of paths
  that USE exactly the arcs in S (in the "positive" direction).

  Since each path contributes +-1 to a specific coefficient,
  the total coefficient is an INTEGER bounded by +-(n-1)!.
  But it's always a power of 2 — WHY?

  HYPOTHESIS: The coefficients are powers of 2 because the
  transfer matrix has entries in {x, 1-x} which mod 2 are
  constant (x ≡ 1-x mod 2), making the permanent mod 2^k
  highly structured.

  Let's verify: is c_S always +-2^{v(S)} where v(S) depends
  only on the degree |S|?
""")

# Check if 2-adic valuation depends only on degree
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    print(f"\n  n={n}: v_2(c_S) by degree:")
    v2_by_deg = defaultdict(Counter)
    for S, v in ml_coeffs.items():
        deg = bin(S).count('1')
        abs_v = abs(v)
        v2 = 0
        while abs_v % 2 == 0:
            v2 += 1
            abs_v //= 2
        v2_by_deg[deg][v2] += 1

    for deg in sorted(v2_by_deg.keys()):
        print(f"    deg {deg}: v_2 distribution = {dict(sorted(v2_by_deg[deg].items()))}")

    # Is v_2(c_S) always = n - 1 - deg?
    print(f"  Check: v_2(c_S) = n-1-|S| = {n-1}-|S|?")
    for S, v in ml_coeffs.items():
        deg = bin(S).count('1')
        abs_v = abs(v)
        v2 = 0
        while abs_v % 2 == 0:
            v2 += 1
            abs_v //= 2
        expected = n - 1 - deg
        if v2 != expected:
            print(f"    FAILS: S has degree {deg}, v_2={v2}, expected {expected}")
            break
    else:
        print(f"    YES! v_2(c_S) = {n-1} - |S| for all S.")

# ======================================================================
# PART G: THE COMPLEMENT INVARIANCE H(T) = H(T^op)
# ======================================================================
print(f"\n{'='*70}")
print("PART G: COMPLEMENT INVARIANCE H(T) = H(T^op)")
print(f"{'='*70}")

print("""
  We found H(T) = H(T^op) for ALL tournaments at n=3,4,5.
  At n=6 (even n), this should FAIL.

  WHY does it hold for odd n?

  Proof: For odd n, every HP P = (v_0, v_1, ..., v_{n-1}) has
  a reversed HP P^rev = (v_{n-1}, ..., v_1, v_0) in T^op.
  The map P -> P^rev is a bijection HP(T) -> HP(T^op).
  So H(T) = H(T^op). QED.

  Wait — that's the TRIVIAL proof! It works for ALL n, not just odd.
  P is an HP of T iff P^rev is an HP of T^op.
  So H(T) = H(T^op) ALWAYS.

  Let me verify at n=6:
""")

n = 6
m = n*(n-1)//2
complement_mask = (1 << m) - 1
same_H = sum(1 for bits in all_T6 if all_T6[bits] == all_T6[bits ^ complement_mask])
print(f"  n=6: H(T) = H(T^op): {same_H}/{len(all_T6)}")

if same_H == len(all_T6):
    print(f"  YES! H(T) = H(T^op) for ALL tournaments at n=6 too!")
    print(f"  This follows from: P -> P^rev is a bijection HP(T) -> HP(T^op).")
else:
    diff = [(bits, all_T6[bits], all_T6[bits ^ complement_mask]) for bits in all_T6 if all_T6[bits] != all_T6[bits ^ complement_mask]]
    print(f"  FAILS for {len(diff)} tournaments! Example: H(T)={diff[0][1]}, H(T^op)={diff[0][2]}")

# ======================================================================
# PART H: FORMULA FOR MULTILINEAR COEFFICIENTS
# ======================================================================
print(f"\n{'='*70}")
print("PART H: FORMULA FOR MULTILINEAR COEFFICIENTS")
print(f"{'='*70}")

print("""
  From the analysis above, if v_2(c_S) = n-1-|S|, then:
  c_S = +-2^{n-1-|S|}

  For n=3: c_emptyset = 1 (v_2=2 -> 2^2=4? No, c=1)
  Hmm, let me recheck.
""")

# Recompute carefully for n=3
n = 3
m = 3
all_T3 = get_all_tournaments(3)

print(f"  n=3: Explicit multilinear coefficients:")
for S in range(1 << m):
    val = 0
    for T in range(1 << m):
        if T & S != T:
            continue
        sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
        val += sign * all_T3[T]
    if val != 0:
        deg = bin(S).count('1')
        abs_v = abs(val)
        v2 = 0
        while abs_v % 2 == 0:
            v2 += 1
            abs_v //= 2
        bits = [i for i in range(m) if S & (1 << i)]
        print(f"    S={bits}, deg={deg}, c_S={val}, |c_S|={abs(val)}, v_2={v2}")

print(f"\n  So c_emptyset = 1 (v_2 = 0), c_{{1}} = 2 (v_2 = 1), c_{{i,j}} = +-2 (v_2 = 1)")
print(f"  Expected n-1-deg: deg=0: 2, deg=1: 1, deg=2: 0")
print(f"  Actual: deg=0: 0, deg=1: 1, deg=2: 1")
print(f"  So v_2 = n-1-deg does NOT hold for the constant term!")
print(f"  But for deg >= 1: v_2(c_S) = 1 for all S at n=3.")

# Let me recheck n=4,5 more carefully
for n in [4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    # Verify pattern
    print(f"\n  n={n}: v_2(c_S) pattern")
    v2_by_deg = defaultdict(set)
    for S, v in ml_coeffs.items():
        deg = bin(S).count('1')
        abs_v = abs(v)
        v2 = 0
        while abs_v % 2 == 0:
            v2 += 1
            abs_v //= 2
        v2_by_deg[deg].add(v2)

    for deg in sorted(v2_by_deg.keys()):
        print(f"    deg {deg}: v_2 values = {sorted(v2_by_deg[deg])}")

    # The pattern is: c_S = (odd number) * 2^{v_2}
    # What is the odd part?
    print(f"  Odd parts of c_S:")
    odd_by_deg = defaultdict(Counter)
    for S, v in ml_coeffs.items():
        deg = bin(S).count('1')
        abs_v = abs(v)
        while abs_v % 2 == 0:
            abs_v //= 2
        odd_by_deg[deg][abs_v] += 1
    for deg in sorted(odd_by_deg.keys()):
        print(f"    deg {deg}: odd parts = {dict(sorted(odd_by_deg[deg].items()))}")

# ======================================================================
# PART I: WHAT SEPARATES H=11, H=13, H=15 IN (1,2,2,2,3)?
# ======================================================================
print(f"\n{'='*70}")
print("PART I: THE (1,2,2,2,3) SPLIT — FINER INVARIANTS")
print(f"{'='*70}")

n = 5
m = n*(n-1)//2

# From Part C: the stratum (1,2,2,2,3) has H in {11, 13, 15}
# and t3 distinguishes them.
# Let's find MORE invariants that might explain this.

target_score = (1, 2, 2, 2, 3)
groups = defaultdict(list)
for bits, H in all_T5.items():
    A = adj_matrix(n, bits)
    ss = score_seq(n, A)
    if ss == target_score:
        t3 = count_3cycles(n, A)

        # 5-cycle count
        t5 = 0
        for P in permutations(range(n)):
            if all(A[P[i]][P[(i+1) % n]] for i in range(n)):
                t5 += 1
        t5 //= n  # Each 5-cycle counted n times

        # Dominance graph: vertex v dominates w if v beats w AND
        # v beats all vertices beaten by w
        scores = [sum(A[i][j] for j in range(n)) for i in range(n)]

        # Subtournament H on 3-vertex subsets
        sub3_H = Counter()
        for triple in combinations(range(n), 3):
            i, j, k = triple
            sub_A = [[0]*3 for _ in range(3)]
            for a_idx, a in enumerate(triple):
                for b_idx, b in enumerate(triple):
                    if a != b:
                        sub_A[a_idx][b_idx] = A[a][b]
            sub_H = count_hp(3, sub_A)
            sub3_H[sub_H] += 1

        groups[H].append({
            'bits': bits, 't3': t3, 't5': t5,
            'sub3': dict(sub3_H), 'scores': tuple(scores)
        })

for H in sorted(groups.keys()):
    items = groups[H]
    print(f"\n  H={H}: {len(items)} tournaments")
    # Show distributions of invariants
    t3_dist = Counter(d['t3'] for d in items)
    t5_dist = Counter(d['t5'] for d in items)
    sub3_dist = Counter(tuple(sorted(d['sub3'].items())) for d in items)
    print(f"    t3: {dict(sorted(t3_dist.items()))}")
    print(f"    t5: {dict(sorted(t5_dist.items()))}")
    print(f"    3-subtournament H: {dict(sorted(sub3_dist.items()))}")

    # Show a representative
    rep = items[0]
    A = adj_matrix(n, rep['bits'])
    print(f"    Representative (bits={rep['bits']}):")
    for i in range(n):
        row = [A[i][j] for j in range(n)]
        print(f"      {row}")

# ======================================================================
# PART J: THE n=5 MASTER FORMULA
# ======================================================================
print(f"\n{'='*70}")
print("PART J: THE n=5 MASTER FORMULA H = f(score, t3)")
print(f"{'='*70}")

# We know: at n=5, score + t3 determines H (from Part C)
# Let's find the EXPLICIT FORMULA

# Build the formula
formula_data = {}
for bits, H in all_T5.items():
    A = adj_matrix(n, bits)
    ss = score_seq(n, A)
    t3 = count_3cycles(n, A)
    key = (ss, t3)
    if key not in formula_data:
        formula_data[key] = H
    else:
        assert formula_data[key] == H, f"Mismatch: {key} -> {formula_data[key]} vs {H}"

print(f"\n  Complete (score, t3) -> H table for n=5:")
for (ss, t3), H in sorted(formula_data.items()):
    print(f"    score={ss}, t3={t3} -> H={H}")

# Check: is H = a + b*t3 within each score stratum?
print(f"\n  H as linear function of t3 within score strata:")
score_groups = defaultdict(list)
for (ss, t3), H in formula_data.items():
    score_groups[ss].append((t3, H))

for ss in sorted(score_groups.keys()):
    pts = sorted(score_groups[ss])
    if len(pts) >= 2:
        # Linear fit H = a + b*t3
        t3_vals = [p[0] for p in pts]
        H_vals = [p[1] for p in pts]
        if len(set(t3_vals)) >= 2:
            # Slope
            dt3 = t3_vals[-1] - t3_vals[0]
            dH = H_vals[-1] - H_vals[0]
            slope = dH / dt3 if dt3 != 0 else 0
            intercept = H_vals[0] - slope * t3_vals[0]
            # Verify
            fits = all(abs(H - (intercept + slope * t3)) < 0.01 for t3, H in pts)
            print(f"    {ss}: H = {intercept:.0f} + {slope:.1f}*t3 {'(EXACT)' if fits else '(APPROX)'}")
        else:
            print(f"    {ss}: constant t3={pts[0][0]}, H={pts[0][1]}")
    else:
        print(f"    {ss}: single point t3={pts[0][0]}, H={pts[0][1]}")

# The most important check: is there a GLOBAL formula?
# H = f(score_seq, t3) that works across ALL strata?
print(f"\n  Testing GLOBAL linear formula H = a + b*t3 + c*something_score:")
# Try H = a + b*t3
# From above: H(1,2,2,2,3) = -1 + 4*t3 seems possible
# Check: for transitive (0,1,2,3,4), t3=0, H=1
# For (1,1,1,3,4), t3=1, H=3 -> 3 = -1 + 4*1? -1+4=3 YES
# For (0,1,3,3,3), t3=2, H=3 -> 3 = -1 + 4*2 = 7? NO

# So not purely t3. But from THM-054: tr(c_{n-3}) = 2*(n-2)!*(t3 - C(n,3)/4)
# At n=5: we know H = 1 + 2*(t3+t5) from the coefficient hierarchy

# Let's check H = 1 + 2*(t3 + t5)
print(f"\n  Testing H = 1 + 2*(t3 + t5):")
correct_count = 0
for bits, H in all_T5.items():
    A = adj_matrix(n, bits)
    t3 = count_3cycles(n, A)
    t5 = 0
    for P in permutations(range(n)):
        if all(A[P[i]][P[(i+1) % n]] for i in range(n)):
            t5 += 1
    t5 //= n
    predicted = 1 + 2*(t3 + t5)
    if predicted == H:
        correct_count += 1

print(f"  Correct: {correct_count}/{len(all_T5)}")

# From memory: H = 1 + 2*(t3+t5+t7) + 4*bc at n=7
# At n=5 there's no t7. Try H = 1 + 2*t3 + 2*t5:
# This is the same formula.

# If not exact, what's the residual?
if correct_count < len(all_T5):
    print(f"  Checking residuals:")
    residuals = Counter()
    for bits, H in list(all_T5.items())[:100]:
        A = adj_matrix(n, bits)
        t3 = count_3cycles(n, A)
        t5 = 0
        for P in permutations(range(n)):
            if all(A[P[i]][P[(i+1) % n]] for i in range(n)):
                t5 += 1
        t5 //= n
        residual = H - (1 + 2*(t3 + t5))
        residuals[residual] += 1
    print(f"    Residual distribution: {dict(sorted(residuals.items()))}")

print("\n" + "=" * 70)
print("DONE — BINARY COEFFICIENT THEOREM AND SCORE-H STRATIFICATION")
print("=" * 70)
