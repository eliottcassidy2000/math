#!/usr/bin/env python3
"""
beta2_z2_structure.py — Direct study of Z₂ = ker(∂₂) ∩ Ω₂

Goal: Understand what 2-cycles look like and why they're always boundaries.

Key structures:
- Ω₂ = transitive triples (a,b,c): a→b→c AND a→c
- Ω₃ = {(a,b,c,d): a→b→c→d, a→c, b→d} (NOT requiring a→d)
- ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
- ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)

For z ∈ Z₂: ∂₂(z) = 0 in Ω₁. We want to find w ∈ Ω₃ with ∂₃(w) = z.

Strategy: Enumerate Z₂ generators and find explicit ∂₃-preimages.

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def get_omega2_paths(A, n):
    """Return Ω₂ paths = transitive triples."""
    paths = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:  # transitive triple
                    paths.append((a,b,c))
    return paths

def get_omega3_paths(A, n):
    """Return Ω₃ paths = (a,b,c,d) with a→b→c→d, a→c, b→d."""
    paths = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b) or not A[b][c] or not A[a][c]: continue
                for d in range(n):
                    if d in (a,b,c) or not A[c][d] or not A[b][d]: continue
                    paths.append((a,b,c,d))
    return paths

def analyze_z2(A, n, verbose=False):
    """Analyze Z₂ = ker(∂₂) ∩ Ω₂ and check if Z₂ ⊆ im(∂₃|Ω₃)."""
    om2_paths = get_omega2_paths(A, n)
    om3_paths = get_omega3_paths(A, n)
    
    if not om2_paths:
        return {'z2_dim': 0, 'b2_dim': 0, 'beta2': 0}
    
    # Build ∂₂: Ω₂ → Ω₁ (restricted)
    om1_paths = [(a,b) for a in range(n) for b in range(n) if a != b and A[a][b]]
    om1_idx = {p: i for i, p in enumerate(om1_paths)}
    
    d2 = np.zeros((len(om1_paths), len(om2_paths)))
    for j, (a,b,c) in enumerate(om2_paths):
        # ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
        if (b,c) in om1_idx: d2[om1_idx[(b,c)], j] += 1
        if (a,c) in om1_idx: d2[om1_idx[(a,c)], j] -= 1
        if (a,b) in om1_idx: d2[om1_idx[(a,b)], j] += 1
    
    # But we need ∂₂ restricted to Ω₂ → Ω₁
    # Actually for path homology, ∂₂ maps Ω₂ into Ω₁ (which equals all edges for tournaments)
    # Z₂ = ker(∂₂)
    
    # Use SVD to find kernel
    U, S, Vt = np.linalg.svd(d2, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z2_dim = len(om2_paths) - rk
    
    if z2_dim == 0:
        return {'z2_dim': 0, 'b2_dim': 0, 'beta2': 0, 'n_om2': len(om2_paths), 'n_om3': len(om3_paths)}
    
    # Z₂ basis (rows of Vt beyond rank)
    z2_basis = Vt[rk:, :]  # shape (z2_dim, |Ω₂|)
    
    # Build ∂₃: Ω₃ → Ω₂
    om2_idx = {p: i for i, p in enumerate(om2_paths)}
    d3 = np.zeros((len(om2_paths), len(om3_paths)))
    for j, (a,b,c,d) in enumerate(om3_paths):
        # ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
        for sgn, face in [(1,(b,c,d)), (-1,(a,c,d)), (1,(a,b,d)), (-1,(a,b,c))]:
            if face in om2_idx:
                d3[om2_idx[face], j] += sgn
    
    # B₂ = im(∂₃)
    b2_dim = np.linalg.matrix_rank(d3, tol=1e-8)
    beta2 = z2_dim - b2_dim
    
    if verbose and z2_dim > 0:
        print(f"    |Ω₂|={len(om2_paths)}, |Ω₃|={len(om3_paths)}, z₂={z2_dim}, b₂={b2_dim}, β₂={beta2}")
        # Show Z₂ generators
        for k in range(min(z2_dim, 3)):
            z = z2_basis[k]
            nonzero = [(om2_paths[i], z[i]) for i in range(len(z)) if abs(z[i]) > 1e-8]
            nonzero.sort(key=lambda x: -abs(x[1]))
            print(f"    Z₂ gen {k}: {nonzero[:6]}...")
    
    return {
        'z2_dim': z2_dim, 'b2_dim': b2_dim, 'beta2': beta2,
        'n_om2': len(om2_paths), 'n_om3': len(om3_paths),
        'z2_basis': z2_basis if z2_dim > 0 else None,
        'om2_paths': om2_paths,
        'om3_paths': om3_paths,
    }

# ============================================================
# Study Z₂ structure at n=5
# ============================================================
print("=" * 70)
print("Z₂ STRUCTURE ANALYSIS")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

z2_dist = Counter()
type_dist = Counter()  # (z2_dim, n_om3) -> count
examples = []

for bits in range(total):
    A = build_adj(n, bits)
    info = analyze_z2(A, n)
    z2_dist[info['z2_dim']] += 1
    if 'n_om3' in info:
        type_dist[(info['z2_dim'], info.get('n_om3', 0))] += 1
    
    if info['z2_dim'] > 0 and len(examples) < 5:
        examples.append((bits, info))

print(f"\nn={n}: z₂ dimension distribution: {dict(sorted(z2_dist.items()))}")
print(f"  (z₂, |Ω₃|) distribution:")
for (z2, om3), cnt in sorted(type_dist.items()):
    if cnt >= 10:
        print(f"    z₂={z2}, |Ω₃|={om3}: {cnt}")

# Show detailed examples
print(f"\nDetailed Z₂ examples (n={n}):")
for bits, info in examples[:3]:
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
    print(f"\n  bits={bits}, scores={scores}")
    analyze_z2(A, n, verbose=True)
    
    if info['z2_basis'] is not None:
        z = info['z2_basis'][0]
        om2 = info['om2_paths']
        # What does the Z₂ generator look like?
        support = [(om2[i], round(z[i], 3)) for i in range(len(z)) if abs(z[i]) > 1e-8]
        print(f"    Full support: {support}")
        
        # Check: which faces of the support triples are in Ω₃?
        om3 = info['om3_paths']
        om3_set = set(om3)
        for (a,b,c), coeff in support:
            # What Ω₃ elements have (a,b,c) as a face?
            covers = [p for p in om3 if (a,b,c) in [(p[1],p[2],p[3]), (p[0],p[2],p[3]), (p[0],p[1],p[3]), (p[0],p[1],p[2])]]
            # More precisely: ∂₃(a,b,c,d) has face signs
            filling = []
            for p in om3:
                faces = [(1,(p[1],p[2],p[3])), (-1,(p[0],p[2],p[3])), (1,(p[0],p[1],p[3])), (-1,(p[0],p[1],p[2]))]
                for sgn, face in faces:
                    if face == (a,b,c):
                        filling.append((p, sgn))
            if filling:
                print(f"    ({a},{b},{c}) [coeff={coeff}] filled by: {filling}")


# ============================================================
# KEY QUESTION: What's the PATTERN of Z₂ generators?
# Are they always sums over 3-cycles?
# ============================================================
print(f"\n{'='*70}")
print("Z₂ GENERATOR PATTERNS")
print("=" * 70)

# For each Z₂ generator, check if its support forms a cycle pattern
pattern_dist = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    info = analyze_z2(A, n)
    if info['z2_dim'] == 0:
        continue
    
    z = info['z2_basis'][0]
    om2 = info['om2_paths']
    support = [(om2[i], z[i]) for i in range(len(z)) if abs(z[i]) > 1e-8]
    
    # Vertices involved
    verts = set()
    for (a,b,c), _ in support:
        verts.update([a,b,c])
    
    # Count of support triples
    n_support = len(support)
    
    # Check if all coefficients are ±1
    all_pm1 = all(abs(abs(c) - 1) < 0.01 for _, c in support)
    
    pattern_dist[(n_support, len(verts), all_pm1)] += 1

print(f"\nn={n}: Z₂ generator patterns (n_triples, n_vertices, all_±1):")
for (ns, nv, pm1), cnt in sorted(pattern_dist.items()):
    print(f"  {ns} triples, {nv} vertices, ±1={pm1}: {cnt}")

# ============================================================
# CRITICAL: Study the a→d condition in Ω₃
# ============================================================
print(f"\n{'='*70}")
print("Ω₃ STRUCTURE: a→d vs a←d")
print("=" * 70)

# For each Ω₃ path (a,b,c,d): a→b→c→d, a→c, b→d
# The condition a→d is NOT required! Study when it holds.
ad_stats = Counter()  # (has a→d) -> count

for bits in range(total):
    A = build_adj(n, bits)
    om3 = get_omega3_paths(A, n)
    for (a,b,c,d) in om3:
        ad_stats[A[a][d]] += 1

print(f"\nΩ₃ paths: a→d={ad_stats[1]}, a←d={ad_stats[0]} (total={sum(ad_stats.values())})")
print(f"  Fraction with a→d: {ad_stats[1]/sum(ad_stats.values()):.3f}")

# When a→d: (a,b,c,d) is FULLY transitive (a dominates b,c,d; b dominates c,d; c dominates d)
# When a←d: we have a→b→c→d→a (contains a 4-cycle!), plus a→c, b→d

# Study: are the a←d paths crucial for filling Z₂?
print(f"\n{'='*70}")
print("FILLING ANALYSIS: Which Ω₃ paths fill Z₂?")
print("=" * 70)

for bits, info in examples[:2]:
    A = build_adj(n, bits)
    if info['z2_dim'] == 0 or info.get('n_om3', 0) == 0:
        continue
    
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
    print(f"\n  bits={bits}, scores={scores}, z₂={info['z2_dim']}")
    
    om2 = info['om2_paths']
    om3 = info['om3_paths']
    z = info['z2_basis'][0]
    
    # Build ∂₃ matrix
    om2_idx = {p: i for i, p in enumerate(om2)}
    d3 = np.zeros((len(om2), len(om3)))
    for j, (a,b,c,d) in enumerate(om3):
        for sgn, face in [(1,(b,c,d)), (-1,(a,c,d)), (1,(a,b,d)), (-1,(a,b,c))]:
            if face in om2_idx:
                d3[om2_idx[face], j] += sgn
    
    # Find w with ∂₃(w) = z (least-squares)
    w = np.linalg.lstsq(d3, z, rcond=None)[0]
    resid = np.linalg.norm(d3 @ w - z)
    print(f"  Filling residual: {resid:.2e}")
    
    if resid < 1e-6:
        # Show which Ω₃ paths are used
        w_support = [(om3[i], round(w[i], 3)) for i in range(len(w)) if abs(w[i]) > 1e-8]
        print(f"  Filling uses {len(w_support)} Ω₃ paths:")
        for (a,b,c,d), coeff in w_support:
            ad_type = "a→d" if A[a][d] else "a←d"
            print(f"    ({a},{b},{c},{d}) [{ad_type}] coeff={coeff}")

print("\nDone.")
