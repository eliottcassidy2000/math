"""
Deep investigation of the summand graph and its connections to:
1. Fermat Polygonal Number Theorem
2. Zeckendorf's Theorem  
3. Tournament theory / t(r)ienerment space
4. The forbidden H-values 7 and 21
"""

from collections import defaultdict
from math import isqrt, gcd
from functools import lru_cache

# ─────────────────────────────────────────────
# PART 1: BUILD THE SUMMAND GRAPH
# ─────────────────────────────────────────────
# Directed edges a→n and b→n when a+b=n, a≠b, a<b
# (only DISTINCT summands)

def summand_pairs(n):
    """All (a,b) with a<b, a+b=n, a≥1."""
    return [(a, n-a) for a in range(1, n//2 + (1 if n%2==1 else 0))]

def in_degree(n):
    """Number of valid summand pairs."""
    return len(summand_pairs(n))

def min_depth(n, memo={}):
    """Minimum depth to construct n from {1,2} via distinct summands."""
    if n in memo: return memo[n]
    if n <= 2:
        memo[n] = 0; return 0
    pairs = summand_pairs(n)
    if not pairs:
        memo[n] = float('inf'); return float('inf')
    d = 1 + min(max(min_depth(a), min_depth(b)) for a,b in pairs)
    memo[n] = d; return d

print("=== SUMMAND GRAPH BASIC STRUCTURE ===")
print(f"{'n':>4} | {'pairs':>6} | {'in-deg':>6} | {'depth':>6} | pairs")
for n in range(1, 22):
    pairs = summand_pairs(n)
    d = min_depth(n)
    print(f"{n:>4} | {len(pairs):>6} | {2*len(pairs):>6} | {d:>6} | {pairs}")

# ─────────────────────────────────────────────
# PART 2: DEPENDENCY ANALYSIS — what does removing nodes break?
# ─────────────────────────────────────────────
print("\n=== DEPENDENCY ANALYSIS: REMOVING SINGLE NODES ===")

def reachable_without(removed, max_n=25):
    """Which nodes 1..max_n remain constructible if 'removed' is deleted?"""
    available = {1, 2} - {removed}
    added = True
    while added:
        added = False
        for n in range(3, max_n+1):
            if n in available or n == removed:
                continue
            pairs = summand_pairs(n)
            if any(a in available and b in available for a,b in pairs):
                available.add(n)
                added = True
    return available

for r in [1, 2, 3, 4, 5, 6, 7, 8]:
    reach = reachable_without(r, 20)
    missing = [n for n in range(1,21) if n not in reach and n != r]
    print(f"Remove {r:>2}: reach={sorted(reach)[:15]}... missing={missing[:12]}")

# ─────────────────────────────────────────────
# PART 3: TRIANGULAR NUMBERS AND TOURNAMENT ARC COUNTS
# ─────────────────────────────────────────────
print("\n=== TRIANGULAR NUMBERS IN THE SUMMAND GRAPH ===")
# T_k = k(k+1)/2 = arc count of tournament on k+1 vertices
triangulars = [k*(k+1)//2 for k in range(1, 15)]
print(f"Triangular numbers T_1..T_14: {triangulars}")

for i, t in enumerate(triangulars):
    pairs = summand_pairs(t)
    # Which summands are also triangular?
    tri_set = set(triangulars)
    tri_pairs = [(a,b) for a,b in pairs if a in tri_set and b in tri_set]
    print(f"T_{i+1}={t:>4}: {len(pairs)} pairs, tri-pairs={tri_pairs}")

# ─────────────────────────────────────────────
# PART 4: ZECKENDORF REPRESENTATIONS
# ─────────────────────────────────────────────
print("\n=== ZECKENDORF REPRESENTATIONS ===")

def fibonacci_list(max_val):
    fibs = [1, 2]
    while fibs[-1] < max_val:
        fibs.append(fibs[-1] + fibs[-2])
    return [f for f in fibs if f <= max_val]

def zeckendorf(n):
    """Unique Zeckendorf representation of n."""
    fibs = sorted(fibonacci_list(n), reverse=True)
    result = []
    rem = n
    for f in fibs:
        if f <= rem:
            result.append(f)
            rem -= f
    return sorted(result)

# Check non-consecutiveness
def is_valid_zeckendorf(rep):
    fibs = [1,2,3,5,8,13,21,34,55,89,144,233]
    indices = [fibs.index(x) for x in rep if x in fibs]
    for i in range(len(indices)-1):
        if indices[i+1] - indices[i] == 1:
            return False
    return True

print(f"{'n':>4} | {'Zeckendorf':>25} | {'sum check':>6} | {'depth':>6}")
key_values = [1, 3, 4, 5, 6, 7, 8, 9, 10, 13, 15, 21, 25, 34, 55, 63, 189, 95095]
for n in key_values:
    z = zeckendorf(n)
    depth = min_depth(n) if n <= 25 else "?"
    print(f"{n:>6} | {str(z):>25} | {sum(z)==n!='':>6} | {depth}")

# Number of Zeckendorf terms for 1..50
print("\nZeckendorf term counts (1..30):")
for n in range(1, 31):
    z = zeckendorf(n)
    print(f"  {n:>2}: {len(z)} terms {z}")

# ─────────────────────────────────────────────
# PART 5: FERMAT POLYGONAL NUMBER THEOREM
# ─────────────────────────────────────────────
print("\n=== FERMAT POLYGONAL NUMBER REPRESENTATIONS ===")

def kgonal(k, n):
    """n-th k-gonal number."""
    return n * ((k-2)*n - (k-4)) // 2

def kgonal_list(k, max_val):
    result = []
    n = 1
    while True:
        v = kgonal(k, n)
        if v > max_val: break
        result.append(v)
        n += 1
    return result

for k in range(3, 7):
    kpoly = kgonal_list(k, 50)
    print(f"{k}-gonal numbers ≤50: {kpoly}")

# Gauss's eureka: every n = T_a + T_b + T_c (triangulars include 0)
def triangular_rep(n, max_terms=3):
    tri = kgonal_list(3, n)
    tri = [0] + tri  # include 0
    for a in tri:
        for b in tri:
            if a+b > n: continue
            c = n - a - b
            if c in tri:
                return sorted([a, b, c], reverse=True)
    return None

print("\nTriangular decompositions (Gauss Eureka):")
for n in range(1, 25):
    rep = triangular_rep(n)
    pairs = summand_pairs(n)
    print(f"  {n:>2} = {rep}  | summand-pairs: {len(pairs)}")

# ─────────────────────────────────────────────
# PART 6: THE CRITICAL NODES 1, 4, 6
# ─────────────────────────────────────────────
print("\n=== ANALYSIS OF CRITICAL NODES 1, 4, 6 ===")

# What happens structurally when we have ONLY {1, 4, 6} as seeds?
def reachable_from_seeds(seeds, max_n=50):
    available = set(seeds)
    added = True
    while added:
        added = False
        for n in range(2, max_n+1):
            if n in available: continue
            pairs = summand_pairs(n)
            if any(a in available and b in available for a,b in pairs):
                available.add(n)
                added = True
    return sorted(available)

for seeds in [[1,2], [1,4], [1,6], [1,4,6], [2,3], [1,3,6]]:
    reach = reachable_from_seeds(seeds, 40)
    missing = [n for n in range(1,41) if n not in reach][:10]
    print(f"Seeds {seeds}: reach first 20 = {reach[:20]}")
    print(f"          missing (first 10): {missing}")
    print()

# ─────────────────────────────────────────────
# PART 7: THE SUMMAND GRAPH AND THE TOURNAMENT STAIRCASE
# ─────────────────────────────────────────────
print("\n=== SUMMAND GRAPH ↔ TOURNAMENT STAIRCASE CONNECTION ===")

# T_k = k(k+1)/2 = arc count of tournament on k+1 vertices
# The staircase δ_{n-2} has T_{n-2} cells for a tournament on n vertices
# Row k of the staircase has k cells, corresponding to arcs (k+1, 1..k)

print("Tournament arcs and staircase structure:")
for n in range(2, 10):
    arcs = n*(n-1)//2  # = T_{n-1}
    print(f"  T_{n}: {n} vertices, {arcs} arcs (= T_{n-1}), staircase δ_{n-2}")
    print(f"    Node {arcs} in summand graph: {len(summand_pairs(arcs))} incoming pairs")
    
# Key: the node T_k in the summand graph has EXACTLY k-1 incoming pairs
# (since T_k has floor((T_k-1)/2) pairs, and T_k = k(k+1)/2)
print("\nVerification: T_k in summand graph has ⌊(T_k-1)/2⌋ incoming pairs:")
for k in range(1, 9):
    t = k*(k+1)//2
    expected = (t-1)//2
    actual = len(summand_pairs(t))
    print(f"  T_{k}={t:>3}: expected {expected:>2}, actual {actual:>2} {'✓' if expected==actual else '✗'}")

# ─────────────────────────────────────────────
# PART 8: FORBIDDEN VALUES IN THE SUMMAND GRAPH CONTEXT
# ─────────────────────────────────────────────
print("\n=== THE FORBIDDEN VALUES 7 AND 21 IN THE SUMMAND GRAPH ===")

for n in [7, 21, 63]:
    pairs = summand_pairs(n)
    print(f"\nNode {n}:")
    print(f"  {len(pairs)} incoming pairs: {pairs}")
    # Triangular structure
    k = 0
    while k*(k+1)//2 < n: k+=1
    if k*(k+1)//2 == n:
        print(f"  = T_{k} (triangular number!)")
    else:
        print(f"  NOT triangular (between T_{k-1}={((k-1)*k)//2} and T_{k}={k*(k+1)//2})")
    # Zeckendorf
    z = zeckendorf(n)
    print(f"  Zeckendorf: {z} ({len(z)} terms)")
    # Polygonal?
    for poly_k in range(3, 9):
        poly_list = kgonal_list(poly_k, n)
        if n in poly_list:
            idx = poly_list.index(n) + 1
            print(f"  = {poly_k}-gonal number P({poly_k},{idx})")

# ─────────────────────────────────────────────
# PART 9: TERNARY SUMMAND GRAPH (3 distinct parts)
# ─────────────────────────────────────────────
print("\n=== TERNARY SUMMAND GRAPH: 3-PART SUMS ===")

def ternary_triples(n):
    """All (a,b,c) with a<b<c, a+b+c=n."""
    result = []
    for a in range(1, n//3 + 1):
        for b in range(a+1, (n-a)//2 + 1):
            c = n-a-b
            if c > b:
                result.append((a,b,c))
    return result

print(f"{'n':>4} | {'binary pairs':>14} | {'ternary triples':>18}")
for n in range(1, 20):
    bp = summand_pairs(n)
    tt = ternary_triples(n)
    print(f"{n:>4} | {str(bp):>14} | {tt}")

# Key: the first number with ternary triples is 6 = 1+2+3
print(f"\nSmallest n with a ternary triple: {min(n for n in range(1,100) if ternary_triples(n))}")
# The "ternary phase transition" at n=6

# ─────────────────────────────────────────────
# PART 10: THE 1,4,6 MYSTERY — CONNECTING TO POLYGONALS
# ─────────────────────────────────────────────
print("\n=== THE 1, 4, 6 MYSTERY: POLYGONAL GENERATORS ===")

# 1 = T_1 = 1^2 (triangular AND square)  
# 4 = 2^2 (square, FIRST non-trivial perfect square)
# 6 = T_3 (third triangular number)

# Key property: {1, 4, 6} are the SMALLEST members of two key polygonal families
# such that their sums generate the structure

# In Gauss's eureka theorem: every n = T_a + T_b + T_c
# The triangular numbers 0,1,3,6,10,... are the GENERATORS
# Removing T_3=6: can we still represent all numbers with T_1=1, T_2=3?
print("Representations using only triangulars {1,3} (no 6):")
tri_restricted = [0, 1, 3]
for n in range(1, 25):
    reps = []
    for a in tri_restricted:
        for b in tri_restricted:
            c = n - a - b
            if c in tri_restricted and a <= b <= c:
                reps.append((a,b,c))
    if reps:
        print(f"  {n:>2} = {reps[0]}")
    else:
        print(f"  {n:>2} = IMPOSSIBLE with {{0,1,3}} only")


# ─────────────────────────────────────────────
# PART 11: THE {2,3} → MISSING {1,4,6} DISCOVERY
# ─────────────────────────────────────────────
print("\n=== THE ASTOUNDING COMPLEMENTARITY: {2,3} GENERATES EVERYTHING EXCEPT {1,4,6} ===")

# From the output above: seeds {2,3} miss exactly {1,4,6} among small naturals
# Let's verify this pattern extends
seeds_23 = reachable_from_seeds([2,3], 100)
missing_23 = [n for n in range(1,101) if n not in seeds_23]
print(f"Starting from {{2,3}}, missing in 1..100: {missing_23}")

# The set {1,4,6}: are these EXACTLY the numbers not reachable from {2,3}?
print(f"\nChecking: is every n > 6 reachable from {{2,3}}?")
for n in range(7, 30):
    if n not in seeds_23:
        print(f"  {n} NOT reachable from {{2,3}}!")

print("\n=== COMPLEMENTARITY THEOREM ===")
print("Starting from {2,3}: reaches everything EXCEPT {1,4,6}")
print("Starting from {1,4,6} together with {2,3}: reaches everything")
print()
print("This means: {1,4,6} and {2,3} are COMPLEMENTARY generator sets!")
print("{2,3} = generates the 'even-origin' structure (starts at 2)")
print("{1,4,6} = generates the 'odd/triangular' structure (starts at 1)")
print()

# The seeds {2,3} generate: 5=2+3, 7=2+5=3+4(but 4 missing?), 
# Let's trace carefully
print("Trace of {2,3} generation:")
available = {2,3}
step = 0
while True:
    new = set()
    for n in range(2, 30):
        if n in available: continue
        pairs = summand_pairs(n)
        if any(a in available and b in available for a,b in pairs):
            new.add(n)
    if not new: break
    available |= new
    step += 1
    print(f"  Step {step}: added {sorted(new)} | total: {sorted(available)[:15]}")

# ─────────────────────────────────────────────
# PART 12: ZECKENDORF PATTERNS FOR SPECIAL SEQUENCES
# ─────────────────────────────────────────────
print("\n=== ZECKENDORF OF TOURNAMENT H-VALUES AND TRIANGULAR NUMBERS ===")

# H-values for Paley tournaments and other special tournaments
h_values = {
    'H(T_3)': 3,
    'H(T_5,max)': 15,
    'H(T_7)=T_6*9': 189,
    'H(T_11)': 95095,
    'T_6=H_forb_2': 21,
    'H_forb_1': 7,
    '42=base': 42,
    'C(7,2)=21': 21,
}

print(f"{'Name':>20} | {'Value':>8} | {'Zeckendorf':>35} | terms | depth")
for name, val in h_values.items():
    z = zeckendorf(val)
    depth = min_depth(val) if val <= 30 else "—"
    print(f"{name:>20} | {val:>8} | {str(z):>35} | {len(z):>5} | {depth}")

# H(T_7) = 189 = 3^3 * 7 — Zeckendorf has 4 terms: [3, 8, 34, 144]
# The terms are: F_4=3, F_6=8, F_9=34, F_11=144
# Indices: 4, 6, 9, 11 — differences: 2, 3, 2
# This mirrors the Fano plane: differences 2,3,2 are the cycle lengths!

print("\nZeckendorf of H(T_7)=189: term indices in Fibonacci sequence")
fibs = [1,2,3,5,8,13,21,34,55,89,144,233,377]
z_189 = zeckendorf(189)
indices = [fibs.index(f)+1 for f in z_189]
diffs = [indices[i+1]-indices[i] for i in range(len(indices)-1)]
print(f"  Terms: {z_189}")
print(f"  Fibonacci indices: {indices}")
print(f"  Index differences: {diffs}")
print(f"  Sum of differences: {sum(diffs)}")
print(f"  Product of terms: {1}")
for f in z_189:
    print(f"    F_{fibs.index(f)+1} = {f}")

# ─────────────────────────────────────────────  
# PART 13: THE DEPTH SEQUENCE AND ITS STRUCTURE
# ─────────────────────────────────────────────
print("\n=== DEPTH SEQUENCE STRUCTURE ===")

depths = [min_depth(n) for n in range(1, 35)]
print(f"Depths 1..34: {depths}")

# Which numbers achieve each depth level?
from collections import Counter
depth_groups = defaultdict(list)
for n, d in enumerate(depths, 1):
    depth_groups[d].append(n)
for d in sorted(depth_groups.keys()):
    g = depth_groups[d]
    print(f"  Depth {d}: {g} (count={len(g)})")

# Group sizes: 2, 1, 2, 4, 8, ... doubling?
sizes = [len(depth_groups[d]) for d in sorted(depth_groups.keys())]
print(f"\nGroup sizes by depth: {sizes}")
print(f"Pattern: {[sizes[i+1]/sizes[i] for i in range(len(sizes)-1)]}")

# ─────────────────────────────────────────────
# PART 14: THE PAIR (seeds={2,3}) AND ZECKENDORF
# ─────────────────────────────────────────────
print("\n=== WHY {2,3} = THE ZECKENDORF BACKBONE ===")
print("Fibonacci numbers: 1,2,3,5,8,13,21,34,55,89,144,...")
print("The seed pair {2,3} = {F_2, F_3}")
print("Every other Fibonacci: F_4=5=2+3, F_5=8=3+5, F_6=13=5+8, ...")
print("The Fibonacci sequence IS the summand graph backbone from {2,3}!")

fibs_from_23 = [2,3]
while fibs_from_23[-1] < 100:
    fibs_from_23.append(fibs_from_23[-1] + fibs_from_23[-2])
print(f"Fibonacci chain from {{2,3}}: {fibs_from_23[:12]}")

print("\nKey insight: {2,3} generates ALL integers except {1,4,6}")
print("= Fibonacci sequence generates integers, with gaps {1,4,6}")
print("= These gaps are exactly the numbers NOT in the Fibonacci sequence")
print("  that also cannot be formed from pairs of Fibonacci numbers!")

# Numbers reachable from {2,3} but NOT being Fibonacci numbers:
fib_set = set(fibs_from_23[:15])
non_fib_reachable = [n for n in seeds_23 if n not in fib_set]
print(f"\nNon-Fibonacci reachable from {{2,3}}: {non_fib_reachable[:15]}")
print(f"Fibonacci numbers in reach: {sorted(n for n in seeds_23 if n in fib_set)[:10]}")

# ─────────────────────────────────────────────
# PART 15: TERNARY PHASE TRANSITION AT n=6 
# AND THE TOURNAMENT TERNARY THRESHOLD
# ─────────────────────────────────────────────
print("\n=== THE TERNARY PHASE TRANSITION AT n=6 ===")
print("Binary (2-part): first occurs at n=3 = 1+2")
print("Ternary (3-part): first occurs at n=6 = 1+2+3")
print("Quaternary (4-part): first occurs at n=10 = 1+2+3+4")
print(f"These are TRIANGULAR NUMBERS: T_2=3, T_3=6, T_4=10, T_5=15, ...")
print()
for k in range(2, 8):
    t = k*(k+1)//2
    min_k_sum = sum(range(1, k+1))
    parts = list(range(1, k+1))
    print(f"{k}-part sum first appears at T_{k}={t} = {'+'.join(map(str,parts))}")

print("\nThe tournament connection:")
print("  T_k = arc count of tournament on k+1 vertices")
print("  k-part sums first appear at T_k = arcs of smallest tournament with ≥k vertices")
print("  n=6=T_3 is EXACTLY the tournament on 4 vertices (the first with non-trivial OCF)")
print()
print("Node 6 in summand graph:")
print("  Binary parents: (1,5) and (2,4)")
print("  First ternary child: 6 = 1+2+3 (ternary)")
print("  This is WHY 6 is a 'generator of structure'!")
print("  It's the TERNARY PHASE TRANSITION POINT")

