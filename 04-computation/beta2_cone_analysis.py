#!/usr/bin/env python3
"""
beta2_cone_analysis.py — Cone/star structure analysis

When we add vertex v to T\v to get T, the new structure is:
- New edges: all edges incident to v
- New 2-paths through v
- New 3-paths through v

The vertex v connects to all other vertices (tournament property).
Split others into S = out-nbrs of v, S̄ = in-nbrs of v.

2-paths through v:
  (a, v, b) where a ∈ S̄, b ∈ S  — "star" paths
  (v, a, b) where a ∈ S — "cone start" paths  
  (a, b, v) where b ∈ S̄ — "cone end" paths

For each, Ω₂ membership depends on faces.

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

print("=" * 70)
print("CONE STRUCTURE: Adding vertex v to T\\v")
print("=" * 70)

# For each tournament and v=0, classify 2-paths and 3-paths
# and study how they contribute to R_p

v = 0

# Study a specific tournament first
# Transitive: 0→1→2→3→4
A = [[0, 1, 1, 1, 1],
     [0, 0, 1, 1, 1],
     [0, 0, 0, 1, 1],
     [0, 0, 0, 0, 1],
     [0, 0, 0, 0, 0]]

S = {j for j in range(n) if j != v and A[v][j]}  # out-nbrs
Sbar = {j for j in range(n) if j != v and A[j][v]}  # in-nbrs
print(f"\nTransitive T: v={v}, S(out)={S}, S̄(in)={Sbar}")

ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)

# Classify 2-paths through v
print(f"\n2-paths through v={v}:")
for path in ap2:
    a, b, c = path
    if v not in (a, b, c): continue
    if v == a:
        print(f"  ({a},{b},{c}): v-start, b={b}∈S, c={c}{'∈S' if c in S else '∈S̄'}")
    elif v == b:
        print(f"  ({a},{b},{c}): v-mid, a={a}{'∈S' if a in S else '∈S̄'}, c={c}{'∈S' if c in S else '∈S̄'}")
    elif v == c:
        print(f"  ({a},{b},{c}): v-end, a={a}, b={b}{'∈S̄' if b in Sbar else '∈S'}")

# Which are TT (transitive)?
print(f"\nTT status of v-paths:")
for path in ap2:
    a, b, c = path
    if v not in (a, b, c): continue
    tt = A[a][c]
    print(f"  ({a},{b},{c}): {'TT' if tt else 'NT'}")

# 3-paths through v
print(f"\n3-paths through v={v}:")
for path in ap3:
    a, b, c, d = path
    if v not in (a, b, c, d): continue
    pos = 'start' if a == v else ('pos1' if b == v else ('pos2' if c == v else 'end'))
    dt = A[a][c] and A[b][d]
    print(f"  ({a},{b},{c},{d}): v-{pos}, {'DT' if dt else 'nDT'}")

# Now study the RELATIVE R_p dimensions for all tournaments
print(f"\n{'='*70}")
print("RELATIVE DIMENSIONS BY OUT-DEGREE")
print("=" * 70)

# For v=0 with out-degree d: v beats d vertices, loses to n-1-d vertices.
# Paths through v are constrained by d.

for target_d in range(n):
    print(f"\n--- d_out(v=0) = {target_d} ---")
    
    # Count path types
    path_type_dist = Counter()
    
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        
        d_out = sum(A[v])
        if d_out != target_d: continue
        
        ap2 = enumerate_allowed_paths(A, n, 2)
        
        # Count v-path types
        n_star = 0  # (a, v, b): a→v, v→b, so a∈S̄, b∈S
        n_start = 0  # (v, a, b): v→a, a→b, so a∈S
        n_end = 0  # (a, b, v): a→b, b→v, so b∈S̄
        
        for path in ap2:
            a, b, c = path
            if b == v: n_star += 1
            elif a == v: n_start += 1
            elif c == v: n_end += 1
        
        path_type_dist[(n_star, n_start, n_end)] += 1
    
    for key, cnt in sorted(path_type_dist.items()):
        star, start, end = key
        total = star + start + end
        print(f"  star={star}, start={start}, end={end}, total={total}: {cnt}")

# STAR paths: (a, v, b) where a→v (a∈S̄) and v→b (b∈S).
# For this to be an allowed 2-path: a→v→b, need a≠b, a→v, v→b.
# Number: |S̄| × |S| = (n-1-d) × d where d = d_out(v).
# This is maximized at d = (n-1)/2 (n=5: d=2, star=6).

# START paths: (v, a, b) where v→a (a∈S) and a→b.
# Number: sum over a∈S of |out-nbrs of a in {j: j≠v}| 
# This depends on the internal structure of T\v.

# END paths: (a, b, v) where a→b and b→v (b∈S̄).
# Number: sum over b∈S̄ of |in-nbrs of b in {j: j≠v}|

print(f"\n{'='*70}")
print("FORMULA: star count = (n-1-d) × d")
print("=" * 70)
for d in range(n):
    star_count = (n-1-d) * d
    print(f"  d_out={d}: star count = {star_count}")

# For the RELATIVE complex R:
# R₂ = Ω₂(T) / Ω₂(T\v)
# The "extra" elements in Ω₂(T) beyond Ω₂(T\v) involve vertex v.
# These include:
# 1. All 2-paths through v that are in Ω₂(T)
# 2. Linear combinations involving v-paths and non-v-paths 
#    that are in Ω₂(T) but not decomposable into Ω₂(T\v) + v-paths

# CRITICAL INSIGHT:
# ∂₂(star path (a,v,b)) = (v,b) - (a,b) + (a,v)
# The faces (v,b) and (a,v) go through v → in R₁.
# The face (a,b) is in A₁(T\v) → in Ω₁(T\v).
# So ∂₂^R maps the star path to (v,b) + (a,v) in R₁.

# ∂₂(start path (v,a,b)) = (a,b) - (v,b) + (v,a)
# Faces: (a,b) ∈ A₁(T\v), (v,b) and (v,a) ∈ R₁.
# ∂₂^R maps to -(v,b) + (v,a) in R₁.

# ∂₂(end path (a,b,v)) = (b,v) - (a,v) + (a,b)
# Faces: (a,b) ∈ A₁(T\v), (b,v) and (a,v) ∈ R₁.
# ∂₂^R maps to (b,v) - (a,v) in R₁.

# R₁ has basis: all edges incident to v.
# These are {(v,j) : j∈S} ∪ {(j,v) : j∈S̄}. Total: n-1 edges.

# The boundary ∂₂^R maps R₂ → R₁ by taking v-incident face components.
# ker(∂₂^R) consists of 2-chains whose v-incident boundary components cancel.

# KEY OBSERVATION: for STAR paths (a,v,b):
# ∂₂^R(a,v,b) = (v,b) + (a,v)
# This is the sum of outgoing v-edge to b and incoming v-edge from a.
# For ker(∂₂^R): we need a combination of star, start, and end paths
# whose v-edge components sum to zero.

# This is like a "balanced flow" on the v-edges!

print(f"\n{'='*70}")  
print("QUOTIENT BOUNDARY ANALYSIS")
print("=" * 70)

# For star paths (a,v,b) with a∈S̄, b∈S:
# ∂₂^R = e_{(a,v)} + e_{(v,b)}  [using signed edge basis]
# Actually more carefully:
# ∂₂(a,v,b) = (v,b) - (a,b) + (a,v)
# In R₁ (modding out A₁(T\v)):
# = (v,b) + (a,v) [since (a,b) ∈ Ω₁(T\v)]
# But the signs matter: ∂₂^R(a,v,b) = +(v,b) + (a,v)

# For start paths (v,a,b) with a∈S:
# ∂₂(v,a,b) = (a,b) - (v,b) + (v,a)
# ∂₂^R = -(v,b) + (v,a)
# But wait: (v,b) might be in R₁ or not depending on whether v→b or b→v.
# If b∈S (v→b): (v,b) ∈ R₁ as an outgoing edge.
# If b∈S̄ (b→v): (v,b) ∉ A₁ (wrong direction). But (b,v) ∈ A₁.
# Hmm, this needs more care about which edges exist.

# Actually for the boundary: ∂₂(v,a,b) gives face (v,b) which is the "skip middle" face.
# (v,b) is in A₁ iff v→b. If b∈S̄ (b→v), then (v,b) ∉ A₁.
# In that case, for (v,a,b) to be in Ω₂, we need ∂₂ ∈ Ω₁ = A₁.
# (v,a,b): faces are (a,b), (v,b), (v,a). 
# (a,b) ∈ A₁ iff a→b (always in allowed path: a→b by definition).
# Wait: a→b is given since (v,a,b) is allowed: v→a, a→b. So (a,b) ∈ A₁.
# (v,a): v→a since a∈S. So (v,a) ∈ A₁.
# (v,b): v→b iff b∈S. If b∈S̄, (v,b) ∉ A₁. Then ∂₂ has a bad face.
# So (v,a,b) is in Ω₂ iff b∈S (TT: v→a→b, v→b).
# If b∈S̄: (v,a,b) is NOT individually in Ω₂, but could be in Ω₂ via cancellation.

# This means the Ω₂ structure around v is constrained by S vs S̄ structure.
# The star paths (a,v,b) with a∈S̄, b∈S: 
#   face (a,b): in A₁ iff a→b or b→a (always, since tournament).
#   But (a,b) is the "skip v" face. It's (a,b) in A₁ iff a→b.
#   For tournaments: exactly one of a→b, b→a holds.
#   So (a,b) ∈ A₁ iff a→b.
#   face (a,v): a→v (a∈S̄). ∈ A₁.
#   face (v,b): v→b (b∈S). ∈ A₁.
# ∂₂(a,v,b) = (v,b) - (a,b) + (a,v) if a→b
#            = (v,b) - ??? + (a,v) if b→a. 
# Wait: (a,b) is ALWAYS a face of ∂₂. For path (a,v,b):
# The standard face map: d₀(a,v,b) = (v,b), d₁(a,v,b) = (a,b), d₂(a,v,b) = (a,v).
# ∂₂ = d₀ - d₁ + d₂ = (v,b) - (a,b) + (a,v).
# Here (a,b) is a FORMAL generator of A₁ iff a→b.
# If b→a, then (a,b) is NOT in A₁. So ∂₂ has a term on a non-A₁ generator.
# For Ω₂ membership: ∂₂(a,v,b) must be in Ω₁ = A₁.
# If b→a: the (a,b) term is not in A₁, so (a,v,b) ∉ Ω₂ individually.
# But pairs that share the same bad face can cancel.

# Summary: star paths (a,v,b) are TT (in Ω₂) iff a→b (in addition to a→v, v→b).

# So the number of TT star paths = |{(a,b) : a∈S̄, b∈S, a→b}|.
# The number of NT star paths = |{(a,b) : a∈S̄, b∈S, b→a}|.
# Total star paths = |S̄|×|S|. TT + NT = total.

# For NT star paths (a,v,b) with b→a:
# ∂₂(a,v,b) = (v,b) - (a,b) + (a,v) where (a,b) is not in A₁.
# But (b,a) IS in A₁ (since b→a).
# Two NT star paths (a,v,b₁) and (a,v,b₂) share the bad face (a,b_i).
# Actually no: (a,v,b₁) has bad face (a,b₁) and (a,v,b₂) has bad face (a,b₂).
# These are DIFFERENT bad faces.

# But: NT star (a₁,v,b) and NT star (a₂,v,b) share bad face partner:
# (a₁,v,b) has bad face (a₁,b) and (a₂,v,b) has bad face (a₂,b).
# These are also different.

# Actually, for cancellation we need two paths with the SAME bad face.
# Same bad face (a,b) means same a and b in the d₁ position.
# So we need two paths (..., ..., ...) with the same d₁ = (a,b) face.
# For star paths (x,v,y): d₁ = (x,y). So (a₁,v,b) and (a₂,v,b) have
# d₁ = (a₁,b) and (a₂,b) — different.
# And (a,v,b₁) and (a,v,b₂) have d₁ = (a,b₁) and (a,b₂) — different.

# So NT star paths CAN'T cancel their bad faces with other star paths.
# They need START or END paths with the same bad face!

# For start path (v,a,b) with b∈S̄ (bad face (v,b)):
# d₁ = (v,b), which is also bad since b→v means (v,b) ∉ A₁.
# So these have DIFFERENT bad faces than star paths.

# Hmm, I think the cancellation mechanism is more subtle.
# Let me just verify computationally which paths are in Ω₂.

# Take the transitive tournament at n=5
A = [[0, 1, 1, 1, 1],
     [0, 0, 1, 1, 1],
     [0, 0, 0, 1, 1],
     [0, 0, 0, 0, 1],
     [0, 0, 0, 0, 0]]

ap2 = enumerate_allowed_paths(A, n, 2)
ap1 = enumerate_allowed_paths(A, n, 1)
om2 = compute_omega_basis(A, n, 2, ap2, ap1)

print("Transitive tournament: all 2-paths and Ω₂ membership")
ap2_tuples = [tuple(p) for p in ap2]
for i, path in enumerate(ap2_tuples):
    a, b, c = path
    # Check if this path appears in any Ω₂ basis vector
    in_omega = any(abs(om2[i, col]) > 1e-10 for col in range(om2.shape[1]))
    tt = "TT" if A[a][c] else "NT"
    through_v = v in path
    print(f"  ({a},{b},{c}) [{tt}] through_v={through_v} in_Ω₂={in_omega}")

# For transitive: ALL paths are TT (since all a→c for a < c).
# So Ω₂ should equal A₂, dim = 10.

# Now try a non-transitive tournament
print("\nNon-transitive example (regular n=5):")
A_reg = [[0, 1, 1, 0, 0],
         [0, 0, 1, 1, 0],
         [0, 0, 0, 1, 1],
         [1, 0, 0, 0, 1],
         [1, 1, 0, 0, 0]]

ap2 = enumerate_allowed_paths(A_reg, n, 2)
ap1 = enumerate_allowed_paths(A_reg, n, 1)
om2 = compute_omega_basis(A_reg, n, 2, ap2, ap1)
print(f"dim Ω₂ = {om2.shape[1]}")
for i, path in enumerate(ap2):
    a, b, c = tuple(path)
    tt = "TT" if A_reg[a][c] else "NT"
    in_omega = any(abs(om2[i, col]) > 1e-10 for col in range(om2.shape[1]))
    print(f"  ({a},{b},{c}) [{tt}] in_Ω₂={in_omega}")

# Show Ω₂ basis to see which NT combinations appear
print(f"\nΩ₂ basis vectors:")
for col in range(om2.shape[1]):
    v_vec = om2[:, col]
    terms = [(tuple(ap2[i]), round(v_vec[i], 3)) for i in range(len(v_vec)) if abs(v_vec[i]) > 1e-10]
    tts = [("TT" if A_reg[t[0][0]][t[0][2]] else "NT") for t in terms]
    print(f"  e{col}: {list(zip([t[0] for t in terms], [t[1] for t in terms], tts))}")

print("\nDone.")
