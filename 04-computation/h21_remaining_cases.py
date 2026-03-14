#!/usr/bin/env python3
"""
H=21 Remaining Cases Check
opus-2026-03-14-S71h

THM-202 claimed 3 remaining cases need to be ruled out for H≠21:
  - K₆ minus 2 edges (α₁=6, α₂=2)
  - K₈ minus 1 edge (α₁=8, α₂=1)
  - K₁₀ (α₁=10, α₂=0)

But our K₃ poison analysis shows ALL graphs with I(G,2)=21
have polynomial 1+4x+3x² (α₁=4, α₂=3).

If α₁=6, α₂=2: I(G,2) = 1+12+8 = 21. This IS 21!
If α₁=8, α₂=1: I(G,2) = 1+16+4 = 21. This IS 21!
If α₁=10: I(G,2) = 1+20 = 21. This IS 21!

But the independence POLYNOMIAL is different in each case:
  - α₁=4, α₂=3: I = 1+4x+3x²
  - α₁=6, α₂=2: I = 1+6x+2x²
  - α₁=8, α₂=1: I = 1+8x+x²
  - α₁=10: I = 1+10x

So THM-202 might be wrong about the "one polynomial" claim.
Let me check: do graphs with these α vectors actually exist?
And do their I-polynomials divide by (1+3x)?
"""

from itertools import combinations

def independence_poly(n, edges):
    adj = set()
    for u, v in edges:
        adj.add((u, v))
        adj.add((v, u))
    coeffs = [0] * (n + 1)
    for mask in range(1 << n):
        verts = [i for i in range(n) if mask & (1 << i)]
        k = len(verts)
        ok = True
        for i in range(k):
            for j in range(i+1, k):
                if (verts[i], verts[j]) in adj:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            coeffs[k] += 1
    return tuple(coeffs)

def poly_at_2(coeffs):
    return sum(c * (2**k) for k, c in enumerate(coeffs))

def clean_poly(coeffs):
    c = list(coeffs)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c

print("=" * 70)
print("CHECKING THM-202 REMAINING CASES")
print("=" * 70)
print()

# Case 1: K₆ minus 2 edges
# K₆ has C(6,2)=15 edges. Removing 2 leaves 13 edges.
# α₁ = 6 means all vertices are in some independent set
# α₂ = 2 means exactly 2 independent pairs
print("Case 1: Graphs on 6 vertices with α₁=6, α₂=2")
n = 6
all_edges = list(combinations(range(n), 2))
count = 0
for missing in combinations(all_edges, 2):
    edges = [e for e in all_edges if e not in missing]
    coeffs = independence_poly(n, edges)
    if coeffs[1] == 6 and coeffs[2] == 2:
        val = poly_at_2(coeffs)
        p = clean_poly(coeffs)
        val_root = sum(c * (-1/3)**k for k, c in enumerate(p))
        poisoned = abs(val_root) < 1e-10
        count += 1
        if count <= 5:
            print(f"  edges_removed={list(missing)}, I(G,x)={p}, I(G,2)={val}, K₃-poisoned: {poisoned}")

if count == 0:
    print("  NONE FOUND")
else:
    print(f"  Total: {count}")

print()

# Case 2: K₈ minus 1 edge
# K₈ has C(8,2)=28 edges. Removing 1 leaves 27 edges.
# α₁ = 8, α₂ = 1
print("Case 2: Graphs on 8 vertices with α₁=8, α₂=1 (K₈ - 1 edge)")
n = 8
all_edges = list(combinations(range(n), 2))
count = 0
for missing in combinations(all_edges, 1):
    edges = [e for e in all_edges if e not in missing]
    coeffs = independence_poly(n, edges)
    if coeffs[1] == 8 and (len(coeffs) <= 2 or coeffs[2] == 1):
        val = poly_at_2(coeffs)
        p = clean_poly(coeffs)
        val_root = sum(c * (-1/3)**k for k, c in enumerate(p))
        poisoned = abs(val_root) < 1e-10
        count += 1
        if count <= 5:
            print(f"  edge_removed={missing[0]}, I(G,x)={p}, I(G,2)={val}, K₃-poisoned: {poisoned}")

if count == 0:
    print("  NONE FOUND")
else:
    print(f"  Total: {count}")

print()

# Case 3: K₁₀ (complete graph on 10 vertices)
# α₁ = 10, α₂ = 0
print("Case 3: K₁₀ (complete graph on 10 vertices)")
n = 10
coeffs = (1, 10) + (0,) * (n - 1)  # I(K_n, x) = 1 + nx
val = poly_at_2(coeffs)
p = clean_poly(coeffs)
val_root = sum(c * (-1/3)**k for k, c in enumerate(p))
poisoned = abs(val_root) < 1e-10
print(f"  I(K₁₀, x) = {p}, I(K₁₀, 2) = {val}, K₃-poisoned: {poisoned}")
print()

# What about K_n in general?
print("I(K_n, 2) = 1 + 2n:")
for n in range(1, 15):
    val = 1 + 2*n
    p = [1, n]
    val_root = 1 + n*(-1/3)
    poisoned = abs(val_root) < 1e-10
    is_21 = " ← H=21!" if val == 21 else ""
    print(f"  K_{n:2d}: I(K_{n},2) = {val:3d}, I(K_{n},-1/3) = {val_root:.4f}, K₃-poisoned: {poisoned}{is_21}")

print()
print("K₃ is the ONLY complete graph that is K₃-poisoned!")
print("(Because I(K₃,-1/3) = 1+3(-1/3) = 0, but I(K_n,-1/3) = 1-n/3 ≠ 0 for n≠3)")
print()

# So K₁₀ with I(K₁₀,2) = 21 is NOT K₃-poisoned.
# This means the K₃ poison alone doesn't explain why H≠21!
# We need the FULL argument:
# 1. K₃ component → blocked by THM-201
# 2. P₄ → blocked by THM-202
# 3. K₁₀ → not K₃-poisoned, but K₁₀ cannot be Ω(T) for other reasons
# 4. K₈-e → not K₃-poisoned
# 5. K₆-2e → not K₃-poisoned

print("=" * 70)
print("CORRECTION: K₃ POISON IS INCOMPLETE FOR H=21!")
print("=" * 70)
print()
print("The K₃ poison (I divisible by 1+3x) captures:")
print("  - K₃ itself (I = 1+3x)")
print("  - P₄ and K₁⊔K₃ (I = 1+4x+3x²)")
print()
print("But it MISSES graphs with different α vectors:")
print("  - K₆-2e: I = 1+6x+2x², I(-1/3) = 1-2+2/9 = -7/9 ≠ 0")
print("  - K₈-e: I = 1+8x+x², I(-1/3) = 1-8/3+1/9 = -14/9 ≠ 0")
print("  - K₁₀: I = 1+10x, I(-1/3) = 1-10/3 = -7/3 ≠ 0")
print()
print("These graphs are NOT K₃-poisoned but also have I(G,2)=21.")
print("The full H≠21 proof requires ruling out EACH of these as Ω(T).")
print()

# Wait — does my earlier analysis saying "only 1 polynomial" contradict this?
# The n≤5 exhaustive found only polynomial 1+4x+3x². But K₆-2e is on 6 vertices!
# K₈-e is on 8 vertices. K₁₀ is on 10 vertices.
# So these graphs exist at HIGHER vertex counts that I didn't check.

print("RESOLUTION: The n≤5 exhaustive search found only 1+4x+3x² because")
print("the other α vectors require n≥6 vertices. The K₃ poison is a")
print("PARTIAL explanation, covering the 4-vertex cases. The full proof")
print("of H≠21 requires additional arguments for the large-vertex cases.")
print()

# BUT: can K₆-2e actually be Ω(T)? That would mean T has exactly 6
# directed odd cycles, arranged as K₆ minus 2 edges.
# 6 directed 3-cycles sharing vertices heavily... this is very constrained.
# At n=7, there are at most C(7,3)=35 directed 3-cycles.
# Having exactly 6 and arranged as K₆-2e is very specific.

# Actually wait: the 6 VERTICES of Ω represent 6 directed odd cycles.
# Two are adjacent iff they share a vertex. K₆-2e means 6 cycles with
# exactly 2 non-adjacent pairs. Each non-adjacent pair = disjoint cycles.

print("STRUCTURAL CHECK: Can K₆-2e be Ω(T)?")
print("  K₆-2e has 6 vertices and 13 edges.")
print("  In Ω(T): vertices = directed odd cycles, edges = vertex-sharing.")
print("  K₆-2e means 6 odd cycles with exactly 2 disjoint pairs.")
print("  For 3-cycles: 2 disjoint pairs need at least 2×6=12 tournament vertices.")
print("  But the 4 non-disjoint cycles share vertices, reducing this.")
print("  This is a complex combinatorial constraint that needs case analysis.")
