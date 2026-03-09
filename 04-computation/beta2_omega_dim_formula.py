#!/usr/bin/env python3
"""
beta2_omega_dim_formula.py — Exact formula for dim(Ω₂) and dim(Ω₃)

Key insight: Ω₂ = TT ⊕ IT_cancelling, where:
- TT = transitive triples (individually in Ω₂)
- IT_cancelling = intransitive differences (x,b₁,y)-(x,b₂,y) where y→x

For intransitive pair (x,y) with y→x, let k(x,y) = #{b: x→b, b→y}.
Then IT contributions: max(0, k(x,y)-1) per pair.

FORMULA: dim(Ω₂) = #TT + Σ_{y→x} max(0, k(x,y)-1)

Similarly for Ω₃: extra elements come from combinations where non-Ω₂ 
face terms cancel.

Author: opus-2026-03-08-S49
"""
import sys, numpy as np
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

def formula_omega2(A, n):
    """Predicted dim(Ω₂) = #TT + Σ max(0, k(x,y)-1) for y→x."""
    tt = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b): continue
                if A[b][c] and A[a][c]:
                    tt += 1
    
    it_extra = 0
    for x in range(n):
        for y in range(n):
            if x == y: continue
            if A[x][y]: continue  # Only intransitive: y→x
            # Count k(x,y) = #{b: x→b, b→y, b≠x,y}
            k = sum(1 for b in range(n) if b != x and b != y and A[x][b] and A[b][y])
            if k >= 2:
                it_extra += k - 1
    
    return tt + it_extra

# Verify at n=5
print("=" * 70)
print("VERIFYING dim(Ω₂) FORMULA")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    ok = 0
    fail = 0
    
    for bits in range(total):
        A = build_adj(n, bits)
        
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
        actual = dim_om(om2)
        predicted = formula_omega2(A, n)
        
        if actual == predicted:
            ok += 1
        else:
            fail += 1
            if fail <= 3:
                scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
                print(f"  FAIL n={n}: bits={bits}, actual={actual}, predicted={predicted}, scores={scores}")
    
    print(f"n={n}: ok={ok}, fail={fail} {'✓' if fail==0 else '✗'}")


# Now Ω₃ formula
# A₃ path (a,b,c,d) with a→b→c→d. Non-A₂ faces:
# ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
# Face (a,c,d): allowed iff a→c (which means c→d is given, a→c→d)
# Face (a,b,d): allowed iff a→b→d (given) AND b→d (wait, this is A₂ = a→b→c paths)
#   Actually A₂ = allowed 2-paths = (x,y,z) with x→y→z (distinct vertices)
# (a,c,d): a→c→d? Need a→c. If c→a: not allowed.
# (a,b,d): a→b→d? Need a→b (given) and b→d. If d→b: not allowed.
# (b,c,d): b→c→d ✓ always (given)
# (a,b,c): a→b→c ✓ always (given)
#
# So Ω₃ individual path conditions: a→c AND b→d (doubly transitive)
# Ω₃ combinations: cancel non-allowed faces
#
# Non-allowed faces from (a,b,c,d):
# - (a,c,d) if c→a [with sign -1]
# - (a,b,d) if d→b [with sign +1]
# Both can be non-allowed simultaneously

# For Ω₃ extra elements:
# Case 1: c→a (face (a,c,d) non-allowed), b→d OK
#   Non-allowed face: -(a,c,d)
#   To cancel: need another path (a,b',c',d) with same non-allowed face (a,c,d)
#   This means: same a,c,d, different b. Path is (a,b',c,d) with a→b'→c→d
#   So we need b' with a→b', b'→c, and b'→d (so that (a,b',d) IS allowed)
#   Wait, the non-allowed face from (a,b',c,d) is also -(a,c,d) if c→a.
#   So (a,b,c,d)-(a,b',c,d) cancels the (a,c,d) face IF both have same sign.
#   Both have -(a,c,d), so their difference cancels it. ✓
#   BUT (a,b,d) and (a,b',d) need to be checked too.

# This is getting complex. Let me just try various formulas computationally.

print(f"\n{'='*70}")
print("FORMULA SEARCH FOR dim(Ω₃)")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

def count_dt4(A, n):
    """Doubly-transitive 4-paths: a→b→c→d, a→c, b→d."""
    cnt = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b) or not A[b][c] or not A[a][c]: continue
                for d in range(n):
                    if d in (a,b,c) or not A[c][d] or not A[b][d]: continue
                    cnt += 1
    return cnt

def count_a3(A, n):
    """All allowed 3-paths: a→b→c→d."""
    cnt = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b) or not A[b][c]: continue
                for d in range(n):
                    if d in (a,b,c) or not A[c][d]: continue
                    cnt += 1
    return cnt

def formula_omega3(A, n):
    """Try to predict dim(Ω₃)."""
    dt4 = count_dt4(A, n)
    
    # Extra from intransitive 3-path combinations
    # For each non-allowed face configuration, count cancelling pairs
    
    # Non-allowed face (a,c,d) from path (a,b,c,d): occurs when c→a
    # Group by (a,c,d): paths are (a,b,c,d) for different b with a→b→c, b→d
    extra = 0
    for a in range(n):
        for c in range(n):
            if c == a or A[a][c]: continue  # Need c→a (intransitive)
            for d in range(n):
                if d in (a,c) or not A[c][d]: continue
                # Count b with a→b→c→d and b→d
                bs_with_bd = []
                bs_without_bd = []
                for b in range(n):
                    if b in (a,c,d) or not A[a][b] or not A[b][c]: continue
                    if A[b][d]:
                        bs_with_bd.append(b)
                    else:
                        bs_without_bd.append(b)
                # Paths (a,b,c,d) with b→d: these have (a,c,d) as only non-allowed face
                # (since a→b→c→d, b→d OK, a→c? No, c→a)
                # Wait, also need to check (a,b,d): a→b→d AND b→d... a→b ✓, b→d ✓
                # Is (a,b,d) allowed? a→b→d, distinct vertices, YES.
                # Is (a,b,d) in A₂? Yes (a→b→d).
                # So when b→d: only non-allowed face is (a,c,d).
                k = len(bs_with_bd)
                if k >= 2:
                    extra += k - 1
    
    # Also: non-allowed face (a,b,d) from path (a,b,c,d): occurs when d→b
    # Group by (a,b,d): paths are (a,b,c,d) for different c with b→c→d, a→c
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for d in range(n):
                if d in (a,b) or A[b][d]: continue  # Need d→b (intransitive)
                # Count c with a→b→c→d and a→c
                cs_with_ac = []
                for c in range(n):
                    if c in (a,b,d) or not A[b][c] or not A[c][d]: continue
                    if A[a][c]:
                        cs_with_ac.append(c)
                k = len(cs_with_ac)
                if k >= 2:
                    extra += k - 1
    
    # But this double-counts when BOTH faces are non-allowed...
    # Need to subtract overcounting
    
    return dt4 + extra  # First attempt

for bits in range(min(total, 20)):
    A = build_adj(n, bits)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))
    actual = dim_om(om3)
    predicted = formula_omega3(A, n)
    dt4 = count_dt4(A, n)
    a3 = count_a3(A, n)
    
    match = "✓" if actual == predicted else "✗"
    print(f"  bits={bits:4d}: |A₃|={a3:2d}, dt4={dt4:2d}, predicted={predicted:2d}, actual={actual:2d} {match}")

print("\nDone.")
