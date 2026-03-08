#!/usr/bin/env python3
"""
beta2_discrete_morse.py — Explore Discrete Morse Theory approach

For β₂=0, we want to show the Ω-complex is "acyclic at level 2".

Approach: Find an acyclic matching on the Ω basis elements that pairs
every Ω₂ element either with an Ω₁ or Ω₃ element.

In discrete Morse theory terms:
- A matching pairs cells of adjacent dimension
- The matching must be acyclic (no directed loops in the Hasse diagram)
- Unmatched cells are "critical" and contribute to homology
- β₂ = 0 iff there are no critical 2-cells

For the Ω complex:
- Ω₂ elements are "2-cells"
- Ω₁ elements are "1-cells"  
- Ω₃ elements are "3-cells"
- We want to match each Ω₂ element with either an Ω₁ or Ω₃ element

Note: this is not standard DMT (which uses CW complexes), but the
algebraic analogue. We need to show the chain complex is exact at Ω₂.

Actually, a better formulation: we want an explicit homotopy operator
h: Ω₂ → Ω₃ such that ∂₃ ∘ h + h' ∘ ∂₂ = id on Z₂.
This would show every 2-cycle is a 2-boundary.

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def compute_homotopy(A, n):
    """Try to find a chain homotopy s: Ω₂ → Ω₃ such that 
    ∂₃ s + s ∂₂ = id on Ω₂ (modulo im ∂₃ and ker ∂₂).
    
    Actually: we want s: Z₂ → Ω₃ such that ∂₃ ∘ s = id on Z₂.
    i.e., s is a right inverse of ∂₃ restricted to Z₂.
    """
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap3 = enumerate_allowed_paths(A, n, 3)
    
    if not ap2: return None
    
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    om1 = compute_omega_basis(A, n, 1, ap1, enumerate_allowed_paths(A, n, 0))
    d2 = om2.shape[1] if om2.ndim == 2 else 0
    if d2 == 0: return None
    
    # ∂₂ in Ω-coords
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    
    if z2_dim == 0: return {'z2': 0, 'b2': 0, 'homotopy': True}
    
    if not ap3: return {'z2': z2_dim, 'b2': 0, 'homotopy': z2_dim == 0}
    
    om3 = compute_omega_basis(A, n, 3, ap3, ap2)
    d3 = om3.shape[1] if om3.ndim == 2 else 0
    if d3 == 0: return {'z2': z2_dim, 'b2': 0, 'homotopy': z2_dim == 0}
    
    # ∂₃ in Ω-coords
    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
    rk3 = np.linalg.matrix_rank(coords3, tol=1e-8)
    
    return {
        'z2': z2_dim, 'b2': rk3,
        'homotopy': z2_dim == rk3,
        'd2': d2, 'd3': d3,
        'rk2': rk2, 'rk3': rk3,
    }


# Instead of DMT, let me try a different approach:
# EXTENDING by one vertex.
# 
# Start with a tournament T on {1,...,n-1}.
# Add vertex 0 with edges to/from all others.
# 
# The chain complex of T is: ... → Ω₂(T) → Ω₁(T) → Ω₀(T)
# The chain complex of T\0: ... → Ω₂(T\0) → Ω₁(T\0) → Ω₀(T\0)
# 
# By the SES: 0 → Ω_*(T\0) → Ω_*(T) → R_* → 0
# where R_p = Ω_p(T)/Ω_p(T\0)
# 
# The elements of R_p are p-chains in Ω_p(T) involving vertex 0.
# 
# For the LES:
# ... → H₂(T\0) → H₂(T) → H₂(T,T\0) → H₁(T\0) → H₁(T) → ...
# 
# If β₂(T\0) = 0 (induction), then H₂(T\0) = 0, so:
# 0 → H₂(T) → H₂(T,T\0) → H₁(T\0) → H₁(T) → ...
# 
# H₂(T) injects into H₂(T,T\0), and im = ker(connecting map ∂_*).
# H₂(T) = 0 iff ∂_* is injective on H₂(T,T\0).
# 
# But ∂_* maps H₂(T,T\0) → H₁(T\0) with im(∂_*) = ker(i_*).
# ∂_* is injective iff ker(∂_*) = 0 iff im(π_*) = 0.
# im(π_*) = H₂(T) (since π_* is injective). So ker(∂_*) = β₂(T).
# 
# This is circular again.
# 
# Let me try yet another approach: DIRECT proof that Z₂ = B₂.

print("=" * 70)
print("PROOF ATTEMPT: Direct Z₂ = B₂ via edge ordering")
print("=" * 70)

# Key idea: use a total order on vertices to construct explicit
# 3-chains whose boundaries are any given 2-cycle.
# 
# For a tournament T on [n], fix the vertex ordering 0 < 1 < ... < n-1.
# Define a "cone operator" C_0: Ω_p(T) → Ω_{p+1}(T) by
# C_0(a₁,...,a_{p+1}) = (0, a₁,...,a_{p+1}) if 0→a₁ and the result is allowed
#                     = 0 otherwise
# 
# Then ∂_{p+1} C_0 + C_0 ∂_p = something...
# 
# In simplicial theory: ∂(v * σ) = σ - v * ∂σ for cone vertex v.
# So ∂C + C∂ = id (the cone is a chain homotopy to a point).
# 
# But in the Ω-complex, C_0 might not preserve the Ω condition.

# Let me test: for the transitive tournament, does the cone operator work?

print("\nTransitive T: cone operator C_0")
n = 5
A = [[0, 1, 1, 1, 1],
     [0, 0, 1, 1, 1],
     [0, 0, 0, 1, 1],
     [0, 0, 0, 0, 1],
     [0, 0, 0, 0, 0]]

# For the transitive tournament, vertex 0 beats everyone.
# C_0(a,b,c) = (0,a,b,c) if 0→a (always here).
# For (0,a,b,c) to be allowed: 0→a, a→b, b→c.
# So the original (a,b,c) must be allowed. ✓
# For (0,a,b,c) to be in Ω₃: ∂₃(0,a,b,c) = (a,b,c) - (0,b,c) + (0,a,c) - (0,a,b)
# All faces must be in Ω₂. For transitive: all TT. ✓

ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)
ap1 = enumerate_allowed_paths(A, n, 1)

om2 = compute_omega_basis(A, n, 2, ap2, ap1)
om1 = compute_omega_basis(A, n, 1, ap1, enumerate_allowed_paths(A, n, 0))
om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

# ∂₂ in Ω coords
bd2 = build_full_boundary_matrix(ap2, ap1)
bd2_om = bd2 @ om2
coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
U2, S2, Vt2 = np.linalg.svd(coords2, full_matrices=True)
rk2 = int(sum(s > 1e-8 for s in S2))

z2_basis = Vt2[rk2:].T  # Ω₂ coords of Z₂
print(f"dim Z₂ = {z2_basis.shape[1]}")

# Try cone: for each Z₂ basis vector, prepend vertex 0 to get Ω₃ element
# z ∈ Z₂: z = Σ α_i (a_i, b_i, c_i)
# C_0(z) = Σ α_i (0, a_i, b_i, c_i) if 0→a_i and path is allowed

ap2_tuples = [tuple(p) for p in ap2]
ap3_tuples = [tuple(p) for p in ap3]
ap3_idx = {p: i for i, p in enumerate(ap3_tuples)}

print(f"\nApplying cone C_0 to Z₂ basis:")
for col in range(z2_basis.shape[1]):
    c_om = z2_basis[:, col]
    c_A = om2 @ c_om  # in A₂ coords
    
    # Build cone
    cone_A = np.zeros(len(ap3))
    for i in range(len(c_A)):
        if abs(c_A[i]) < 1e-10: continue
        a, b, c = ap2_tuples[i]
        cone_path = (0, a, b, c)
        if A[0][a] and cone_path in ap3_idx:
            cone_A[ap3_idx[cone_path]] = c_A[i]
        elif A[0][a]:
            print(f"    Cone path {cone_path} not in A₃!")
        else:
            print(f"    0 does not beat {a}, skip")
    
    # Check: is cone in Ω₃?
    if om3.ndim == 2 and om3.shape[1] > 0:
        coords_c, _, _, _ = np.linalg.lstsq(om3, cone_A, rcond=None)
        recon = om3 @ coords_c
        err = np.max(np.abs(cone_A - recon))
        print(f"  C_0(z{col}): Ω₃ error = {err:.2e}")
        
        # Check ∂₃(C_0(z)) = z ?
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd_cone = bd3 @ cone_A  # in A₂ coords
        z_A = om2 @ c_om
        diff = bd_cone - z_A
        err_bd = np.max(np.abs(diff))
        print(f"    ||∂₃C_0(z) - z|| = {err_bd:.2e}")
        
        if err_bd > 1e-6:
            # ∂₃(C_0(z)) ≠ z, so cone doesn't directly work
            # But ∂₃(C_0(z)) - z should be in im(∂₃) if β₂=0
            # Actually: ∂(v*σ) = σ - v*∂σ
            # For p=2: ∂₃(0,a,b,c) = (a,b,c) - (0,b,c) + (0,a,c) - (0,a,b)
            # So ∂₃(C_0(z)) = z - C_0(∂₂(z)) = z (since z is a cycle)
            # Wait — C_0(∂₂(z)) is a sum of coned 1-chains.
            # ∂₂(z) = 0 since z is a cycle. So C_0(∂₂(z)) = 0.
            # Therefore ∂₃(C_0(z)) should equal z!
            print(f"    BUG: cone formula fails!")
            # Print the discrepancy
            for i in range(len(diff)):
                if abs(diff[i]) > 1e-6:
                    print(f"      Path {ap2_tuples[i]}: diff = {diff[i]:.4f}")

# Now try non-transitive tournament where cone from sink vertex
print(f"\n{'='*70}")
print("NON-TRANSITIVE: cone from source vertex")
print("=" * 70)

# Use the regular tournament (c₃=5)
A_reg = [[0, 1, 1, 0, 0],
         [0, 0, 1, 1, 0],
         [0, 0, 0, 1, 1],
         [1, 0, 0, 0, 1],
         [1, 1, 0, 0, 0]]

# Vertex 0 beats {1,2}, loses to {3,4}
# Cone from 0: C_0(a,b,c) = (0,a,b,c) only if 0→a, i.e., a∈{1,2}

ap2 = enumerate_allowed_paths(A_reg, n, 2)
ap3 = enumerate_allowed_paths(A_reg, n, 3)
ap1 = enumerate_allowed_paths(A_reg, n, 1)

om2 = compute_omega_basis(A_reg, n, 2, ap2, ap1)
om1 = compute_omega_basis(A_reg, n, 1, ap1, enumerate_allowed_paths(A_reg, n, 0))
om3 = compute_omega_basis(A_reg, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

bd2 = build_full_boundary_matrix(ap2, ap1)
bd2_om = bd2 @ om2
coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
U2, S2, Vt2 = np.linalg.svd(coords2, full_matrices=True)
rk2 = int(sum(s > 1e-8 for s in S2))
z2_dim = om2.shape[1] - rk2

print(f"Regular T: dim Ω₂={om2.shape[1]}, Z₂={z2_dim}")
print(f"v=0 beats {[j for j in range(n) if A_reg[0][j]]}")

z2_basis = Vt2[rk2:].T
ap2_tuples = [tuple(p) for p in ap2]
ap3_tuples = [tuple(p) for p in ap3]
ap3_idx = {p: i for i, p in enumerate(ap3_tuples)}

# Cone from vertex 0
for col in range(z2_dim):
    c_A = om2 @ z2_basis[:, col]
    cone_A = np.zeros(len(ap3))
    missing = []
    for i in range(len(c_A)):
        if abs(c_A[i]) < 1e-10: continue
        a, b, c = ap2_tuples[i]
        if A_reg[0][a]:
            cone_path = (0, a, b, c)
            if cone_path in ap3_idx:
                cone_A[ap3_idx[cone_path]] = c_A[i]
            else:
                missing.append((ap2_tuples[i], c_A[i], 'not_in_A3'))
        else:
            missing.append((ap2_tuples[i], c_A[i], f'0_not_beat_{a}'))
    
    # Check ∂₃(C_0(z)) vs z
    if om3.ndim == 2 and om3.shape[1] > 0:
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd_cone = bd3 @ cone_A
        z_A = c_A
        diff = bd_cone - z_A
        err = np.max(np.abs(diff))
        print(f"\n  z{col}: ||∂₃C_0 - z|| = {err:.2e}, missing={len(missing)}")
        if missing:
            for m in missing:
                print(f"    missing: {m}")
        if err > 1e-6:
            nz_diff = [(ap2_tuples[i], round(diff[i],4)) for i in range(len(diff)) if abs(diff[i]) > 1e-8]
            print(f"    Discrepancy: {nz_diff}")

# The cone from a single vertex won't work for non-transitive tournaments
# because the vertex doesn't beat everyone.
# 
# IDEA: Use a CHAIN of cone operators from different vertices!
# Start with z ∈ Z₂. Try coning from vertex v_1 that beats the first 
# vertex in z. Then the error term can be coned from v_2, etc.
# 
# This is essentially the "divide and conquer" approach:
# Decompose z into parts, each handled by a different cone vertex.

print(f"\n{'='*70}")
print("MULTI-VERTEX CONE DECOMPOSITION")
print("=" * 70)

# For a 2-cycle z = Σ α_i (a_i, b_i, c_i):
# For each TT path (a,b,c) with a→b→c and a→c:
# Choose a vertex v that beats a (v→a). 
# Then (v,a,b,c) is a 4-path (v→a→b→c).
# For DT: need v→b. For a→c already (TT). 
# So (v,a,b,c) is DT iff v→b.

# For the transitive tournament: v=0 beats everyone, so the cone always works.
# For general tournaments: need to choose v depending on the path.

# THEOREM CANDIDATE: In a tournament T on [n], for any TT path (a,b,c),
# there exists a vertex v ≠ a,b,c with v→a and v→b (making (v,a,b,c) DT).
# 
# Proof: The set S = {u : u→a AND u→b} ∩ {u ∉ {a,b,c}}.
# Since T is a tournament on n vertices, each u ≠ a,b has a→u or u→a.
# |{u : u→a}| = n - 1 - d_out(a) (number of in-neighbors of a).
# Similarly for b.
# |S| ≥ |{u : u→a}| + |{u : u→b}| - (n - 3) by inclusion-exclusion.
# = (n-1-d_out(a)) + (n-1-d_out(b)) - (n-3) = n - 1 - d_out(a) - d_out(b).
# For n=5, d_out(a), d_out(b) ∈ {0,...,4}. Need n-1-da-db ≥ 1, i.e., da+db ≤ 3.
# This FAILS when da+db ≥ 4.

# So the simple cone approach doesn't always work. We need a more 
# sophisticated argument.

# BUT: we don't need (v,a,b,c) to be DT — just in Ω₃.
# And Ω₃ can contain non-DT paths (via cancellation).
# So maybe we can use combinations of cones.

# Key question: for any z ∈ Z₂ in a tournament, can we find x ∈ Ω₃
# with ∂₃(x) = z? This is β₂ = 0.
# We've verified this computationally. The challenge is the proof.

# Let me check: what fraction of DT paths covers Z₂?
print("\n3-chain coverage: can DT paths alone fill Z₂?")

for bits_test in [0, 2, 4, 8, 15]:
    A_t = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits_test >> idx) & 1: A_t[i][j] = 1
        else: A_t[j][i] = 1
    
    ap2 = enumerate_allowed_paths(A_t, n, 2)
    ap3 = enumerate_allowed_paths(A_t, n, 3)
    ap1 = enumerate_allowed_paths(A_t, n, 1)
    if not ap2: continue
    
    om2 = compute_omega_basis(A_t, n, 2, ap2, ap1)
    om1 = compute_omega_basis(A_t, n, 1, ap1, enumerate_allowed_paths(A_t, n, 0))
    om3 = compute_omega_basis(A_t, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))
    
    # Only DT elements of Ω₃
    d3 = om3.shape[1] if om3.ndim == 2 else 0
    if d3 == 0: continue
    
    ap3_tuples = [tuple(p) for p in ap3]
    dt_indices = [i for i, (a,b,c,d) in enumerate(ap3_tuples) if A_t[a][c] and A_t[b][d]]
    
    # DT-only Ω₃ elements: restrict om3 to rows that are DT
    dt_paths_in_om3 = om3[dt_indices, :] if dt_indices else np.zeros((0, d3))
    
    # Boundary of DT elements
    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_dt = np.zeros((len(ap2), 0))
    
    for col in range(d3):
        v = om3[:, col]
        # Is this element purely DT?
        all_dt = all(abs(v[i]) < 1e-10 or i in dt_indices for i in range(len(v)))
        if all_dt:
            bd_v = bd3 @ v
            bd3_dt = np.hstack([bd3_dt, bd_v.reshape(-1, 1)]) if bd3_dt.shape[1] > 0 else bd_v.reshape(-1, 1)
    
    # Check rank of DT-only boundaries vs Z₂
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = om2.shape[1] - rk2
    
    rk_dt = np.linalg.matrix_rank(bd3_dt, tol=1e-8) if bd3_dt.shape[1] > 0 else 0
    
    scores = sorted(sum(A_t[i]) for i in range(n))
    print(f"  T#{bits_test} (scores={scores}): Z₂={z2_dim}, DT coverage={rk_dt}, "
          f"sufficient={'YES' if rk_dt >= z2_dim else 'NO'}")

print("\nDone.")
