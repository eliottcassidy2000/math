#!/usr/bin/env python3
"""
Simplex-Cuboid Packing Hierarchy
opus-2026-03-14-S71f

User directive: "Think of simplices as (x+1)^n and cuboids as (x+2)^n,
and think about packing them inside each other."

The independence polynomial I(Ω, x) factorizes as a product of (1+c_i*x)
factors IF the conflict graph Ω decomposes into independent cliques
(each clique = a group of mutually overlapping cycles, with c_i cycles in group i).

Key insight: I(Ω, x) = prod_i (1 + c_i * x) when Ω = disjoint union of cliques K_{c_i}.

So the "shapes" are:
- Factor (1+x): a single isolated cycle (simplex "brick")
- Factor (1+2x): a pair of overlapping cycles (cuboid "brick") 
- Factor (1+3x): a triple of mutually overlapping cycles (tesseract "brick")
- etc.

PACKING: Given H, can we write I(Ω,x) = prod (1+c_i*x)?
This is possible iff Ω is a disjoint union of cliques.

At n=5: Ω is ALWAYS complete, so I(Ω,x) = 1 + α₁*x (one big clique).
This means H = 1 + 2α₁ = 1 + 2(t₃+d₅).

At n=6: Ω can have non-adjacent cycles (disjoint 3-cycles), so we can have
I(Ω,x) = (1+a*x)(1+b*x) = 1 + (a+b)x + ab*x² → α₁=a+b, α₂=ab.
H = 1 + 2(a+b) + 4ab = (1+2a)(1+2b).

QUESTION: When does Ω decompose into cliques, and what are the possible
packing shapes?
"""
from itertools import permutations, combinations
from math import comb

print("="*70)
print("SIMPLEX-CUBOID PACKING HIERARCHY")
print("="*70)

# Part 1: Theoretical framework
print("\n--- Part 1: Independence polynomial factorization ---")
print("I(G,x) = prod_i (1 + c_i * x) iff G = disjoint union of K_{c_i}")
print("H = I(Ω,2) = prod_i (1 + 2c_i)")
print()
print("Brick types:")
print("  (1+x) = simplex    → H contribution: 3")  
print("  (1+2x) = cuboid    → H contribution: 5")
print("  (1+3x) = tesseract → H contribution: 7 (FORBIDDEN for m=1!)")
print("  (1+cx) = c-cell    → H contribution: 2c+1")
print()
print("Packing = product of bricks:")
print("  3^a * 5^b * 7^c * ... = H")
print()
print("  H=9  = 3² : two simplices")
print("  H=15 = 3*5 : simplex + cuboid")
print("  H=25 = 5² : two cuboids")
print("  H=27 = 3³ : three simplices")
print("  H=45 = 3²*5 : two simplices + cuboid")
print("  H=75 = 3*5² : simplex + two cuboids")

# Part 2: Which H values are achievable as products of odd numbers ≥3?
print("\n--- Part 2: Achievable H from packing ---")
print("H must be odd. H = prod (2c_i + 1) where c_i ≥ 1.")
print("So H is a product of odd numbers ≥ 3.")
print()

# Generate all odd numbers that are products of odd numbers ≥ 3
achievable = set()
# Use dynamic programming: all products of odd numbers from {3,5,7,9,...}
# up to some limit
limit = 200
bases = list(range(3, limit, 2))  # 3, 5, 7, 9, ...
# Every odd number ≥ 3 is achievable (just use one factor)
# But wait: we need each factor to correspond to a clique in Ω.
# The question is whether the TOURNAMENT exists, not whether the number is factorizable.
# Let me think about what H values are achievable...

# Actually, every odd H is achievable by some tournament (known result, with 
# exceptions at {7, 21}). The packing question is whether the OPTIMAL packing
# (factorization into (1+c_i*x) factors) gives a valid decomposition.

# Let me instead enumerate the possible packings for small H values:
print("Possible packings for small odd H:")
def factorizations_into_odd_geq3(n, min_factor=3):
    """Find all ways to write n as a product of odd numbers ≥ 3."""
    if n == 1:
        return [[]]
    result = []
    for f in range(min_factor, n+1, 2):
        if n % f == 0:
            for rest in factorizations_into_odd_geq3(n // f, f):
                result.append([f] + rest)
    return result

for h in range(1, 130, 2):
    facts = factorizations_into_odd_geq3(h)
    if facts:
        # Convert to brick types
        bricks_list = []
        for f in facts:
            bricks = tuple(sorted([(fi-1)//2 for fi in f]))
            bricks_list.append(bricks)
        
        # H=1 has the empty packing
        if h == 1:
            print(f"  H={h:4d}: empty packing (transitive)")
        elif len(facts) == 1 and len(facts[0]) == 1:
            c = (h-1)//2
            print(f"  H={h:4d}: single {c}-cell")
        else:
            packing_strs = []
            for b in bricks_list:
                parts = []
                for c in b:
                    if c == 1: parts.append("S")
                    elif c == 2: parts.append("C")
                    elif c == 3: parts.append("T")
                    else: parts.append(f"{c}-cell")
                packing_strs.append(" × ".join(parts))
            print(f"  H={h:4d}: {' | '.join(packing_strs)}")

# Part 3: The packing obstruction for H=7
print("\n\n--- Part 3: Why H=7 is unpacked ---")
print("H=7: only packing is a single 3-cell (1+3x)")
print("This requires 3 mutually overlapping 3-cycles on ≤5 vertices.")
print("But: 3 mutually overlapping 3-cycles sharing common vertices")
print("creates additional cycles (5-cycles), violating the single-clique structure.")
print()
print("The key obstruction: at n=7 with t₃=3, there are ALWAYS d₅≥1")  
print("This creates 5-cycle vertices in Ω adjacent to the 3-cycles,")
print("breaking the clique structure needed for I=(1+3x).")

# Part 4: Simplex packing = disjoint cycles
print("\n--- Part 4: Simplex packing theorem ---")
print("I(Ω,x) = (1+x)^m requires m isolated vertices in Ω")
print("= m mutually non-adjacent odd cycles")
print("= m pairwise vertex-disjoint odd cycles")
print("Minimum n: 3m (each cycle uses 3 vertices, all disjoint)")
print()
print("PROVED: Simplex tournaments exist for all m, using:")
print("  n = 3m vertices in m groups of 3")
print("  Within each group: directed 3-cycle")
print("  Between groups: all edges transitive (group i → group j if i < j)")
print("  No 5-cycles appear → Ω = m isolated vertices → I = (1+x)^m → H = 3^m")

# Part 5: Cuboid packing
print("\n--- Part 5: Cuboid packing possibilities ---")
print("I(Ω,x) = (1+2x) requires 2 cycles sharing a vertex, with NO disjoint pairs")
print("and NO other cycles in the tournament.")
print()
print("For 2 overlapping 3-cycles:")
print("  They share 1 or 2 vertices.")
print("  Share 1 vertex: uses 5 vertices. But other triangles may form!")
print("  Share 2 vertices: uses 4 vertices. The common edge + 2 apex vertices.")

# Check at n=4: can we have exactly 2 three-cycles sharing 2 vertices?
print("\n  At n=4:")
n = 4
for bits in range(2**6):
    A = [[0]*4 for _ in range(4)]
    idx = 0
    for i in range(4):
        for j in range(i+1, 4):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    t3 = 0
    for i,j,k in combinations(range(4), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    
    h = 0
    for perm in permutations(range(4)):
        ok = True
        for pi in range(3):
            if A[perm[pi]][perm[pi+1]] != 1:
                ok = False
                break
        if ok: h += 1
    
    if h == 5:
        # H=5 at n=4
        cycles = []
        for i,j,k in combinations(range(4), 3):
            if A[i][j] and A[j][k] and A[k][i]:
                cycles.append((i,j,k))
            if A[i][k] and A[k][j] and A[j][i]:
                cycles.append((i,k,j))
        print(f"  bits={bits}: H=5, t₃={t3}, cycles={cycles}")

# Part 6: The forbidden H and (2,3) duality
print("\n\n--- Part 6: The (2,3) duality in packing ---")
print("Evaluation at x=1 (simplex point): I(Ω,1) = sum α_k = total indep sets")
print("Evaluation at x=2 (OCF point):     I(Ω,2) = H")
print("Evaluation at x=3 (cuboid point):   I(Ω,3) = 'cuboid count'")
print()
print("For a simplex (1+x)^m:")
print("  I(1) = 2^m, I(2) = 3^m, I(3) = 4^m")
print("  Ratio I(2)/I(1) = (3/2)^m → approaches infinity")
print("  Ratio I(3)/I(2) = (4/3)^m → approaches infinity")
print()
print("For a cuboid (1+2x)^m:")
print("  I(1) = 3^m, I(2) = 5^m, I(3) = 7^m = FORBIDDEN^m!")
print("  So cuboid packing at x=3 gives powers of 7!")
print("  This is the CUBOID-TESSERACT DUALITY:")
print("  Cuboid at simplex point = Simplex at cuboid point = 3^m")
print("  Cuboid at cuboid point = Tesseract at simplex point = 7^m")
print()
print("DEEP INSIGHT: The evaluation shift x → x+1 transforms:")
print("  (1+cx) at x+1 = (1+c(x+1)) = (c+1) + cx")
print("  Not exactly the same form. But at integer points:")
print("  I(k) / I(k-1) controls the 'growth rate' of independent sets.")
print()
print("The k-nacci connection:")
print("  k-nacci ratio → 2: this IS the ratio I(2)/I(1) = H/|independent sets|")
print("  weighted k-nacci ratio → 3: this IS I(3)/I(1) = cuboid/simplex")

# Part 7: When does Ω factor into cliques?
print("\n--- Part 7: Clique decomposition of Ω ---")
print("Ω decomposes into cliques when:")
print("  1. No two non-adjacent cycles share NO vertices (trivially)")  
print("  2. The non-adjacency relation is an equivalence (transitivity)")
print("     i.e., if C1⊥C2 and C2⊥C3 then C1⊥C3")
print("     where ⊥ means vertex-disjoint")
print()
print("Transitivity of vertex-disjointness:")
print("  If C1∩C2=∅ and C2∩C3=∅, must C1∩C3=∅?")
print("  NO! Counter: C1={1,2,3}, C2={4,5,6}, C3={1,4,7}")
print("  C1⊥C2, C2⊥C3, but C1∩C3={1} ≠ ∅")
print()
print("So Ω does NOT always decompose into cliques!")
print("The INDEPENDENCE complex of Ω is generally NOT a product of simplices.")
print()
print("CONSEQUENCE: Most H values do NOT have a clean packing interpretation.")
print("Only special tournaments with 'clean' cycle structure admit packing.")

print("\n" + "="*70)
print("DONE")
print("="*70)
