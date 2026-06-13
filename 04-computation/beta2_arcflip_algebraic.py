#!/usr/bin/env python3
"""
beta2_arcflip_algebraic.py — Algebraic analysis of β₂ invariance under arc flip

PROOF STRATEGY:
1. β₂(transitive) = 0 (transitive tournament is shellable)
2. Any tournament is reachable from transitive by arc flips
3. β₂ is invariant under arc flips (THIS is what we need to prove)

Under arc flip u→v to v→u:
- V\{u,v} partitions into A,B,C,D sets
  A = {w : u→w, w→v}  (u→w→v, "transitive mediators")
  B = {w : w→u, v→w}  (v→w→u direction, "reverse mediators")
  C = {w : u→w, v→w}  (both out to w)
  D = {w : w→u, w→v}  (both in from w)

Key identity: Δ|A₃| = (n-3)·Δ|A₂|
At n=5: Δ|A₃| = 2·Δ|A₂|, verified 100%

Question: what algebraic structure makes ΔZ₂ = ΔB₂?

Author: opus-2026-03-08-S45
"""
import sys
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


def flip_arc(A, n, u, v):
    """Return copy of A with arc u→v flipped to v→u."""
    B = [row[:] for row in A]
    B[u][v] = 0
    B[v][u] = 1
    return B


def get_ABCD(A, n, u, v):
    """Partition V\{u,v} into A,B,C,D sets."""
    sets = {'A': [], 'B': [], 'C': [], 'D': []}
    for w in range(n):
        if w == u or w == v:
            continue
        uw = A[u][w]
        wv = A[w][v]
        if uw and wv:
            sets['A'].append(w)
        elif not uw and not wv:
            sets['B'].append(w)
        elif uw and not wv:
            sets['C'].append(w)
        else:
            sets['D'].append(w)
    return sets


def compute_full_chain_data(A, n):
    """Compute Z₂, B₂, dim Ω₂, dim Ω₃ and structural data."""
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    result = {
        '|A2|': len(ap2), '|A3|': len(ap3),
    }

    if not ap2:
        result.update({'dO2': 0, 'dO3': 0, 'Z2': 0, 'B2': 0, 'beta2': 0})
        return result

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 and om2.shape[0] > 0 else 0
    result['dO2'] = d2

    if d2 == 0:
        result.update({'dO3': 0, 'Z2': 0, 'B2': 0, 'beta2': 0})
        return result

    # Z₂
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2 = d2 - rk2
    result['Z2'] = z2

    # B₂
    if ap3:
        om3 = compute_omega_basis(A, n, 3, ap3, ap2)
        d3 = om3.shape[1] if om3.ndim == 2 and om3.shape[0] > 0 else 0
        result['dO3'] = d3
        if d3 > 0:
            bd3 = build_full_boundary_matrix(ap3, ap2)
            bd3_om = bd3 @ om3
            bd3_in_om2 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
            b2 = np.linalg.matrix_rank(bd3_in_om2, tol=1e-8)
        else:
            b2 = 0
    else:
        result['dO3'] = 0
        b2 = 0

    result['B2'] = b2
    result['beta2'] = z2 - b2
    return result


# ===== MAIN =====
print("=" * 70)
print("ALGEBRAIC ANALYSIS OF β₂ ARC-FLIP INVARIANCE")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Track ALL changes under all arc flips
changes = []
abcd_stats = defaultdict(list)

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    data_before = compute_full_chain_data(A, n)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue

            B = flip_arc(A, n, u, v)
            data_after = compute_full_chain_data(B, n)
            sets = get_ABCD(A, n, u, v)
            a, b, c, d = len(sets['A']), len(sets['B']), len(sets['C']), len(sets['D'])

            dZ2 = data_after['Z2'] - data_before['Z2']
            dB2 = data_after['B2'] - data_before['B2']
            dO2 = data_after['dO2'] - data_before['dO2']
            dO3 = data_after['dO3'] - data_before['dO3']
            dA2 = data_after['|A2|'] - data_before['|A2|']
            dA3 = data_after['|A3|'] - data_before['|A3|']

            changes.append({
                'bits': bits, 'u': u, 'v': v,
                'a': a, 'b': b, 'c': c, 'd': d,
                'dZ2': dZ2, 'dB2': dB2, 'dO2': dO2, 'dO3': dO3,
                'dA2': dA2, 'dA3': dA3,
                'dbeta2': dZ2 - dB2,
            })
            abcd_stats[(a,b,c,d)].append({
                'dZ2': dZ2, 'dB2': dB2, 'dO2': dO2, 'dO3': dO3,
            })

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}")

# Verify beta2 invariance
dbeta2_dist = Counter(ch['dbeta2'] for ch in changes)
print(f"\nΔβ₂ distribution: {dict(dbeta2_dist)}")
print(f"β₂ invariant: {all(ch['dbeta2'] == 0 for ch in changes)}")

# Analyze dZ2, dB2 as function of (a,b,c,d)
print(f"\n{'='*70}")
print("ΔZ₂ AND ΔB₂ AS FUNCTION OF (a,b,c,d)")
print("=" * 70)

for key in sorted(abcd_stats.keys()):
    entries = abcd_stats[key]
    dz2_vals = Counter(e['dZ2'] for e in entries)
    db2_vals = Counter(e['dB2'] for e in entries)
    do2_vals = Counter(e['dO2'] for e in entries)
    do3_vals = Counter(e['dO3'] for e in entries)
    a, b, c, d = key
    print(f"\n  (a,b,c,d)=({a},{b},{c},{d}): {len(entries)} flips")
    print(f"    ΔZ₂: {dict(dz2_vals)}")
    print(f"    ΔB₂: {dict(db2_vals)}")
    print(f"    ΔΩ₂: {dict(do2_vals)}")
    print(f"    ΔΩ₃: {dict(do3_vals)}")

# Key question: is ΔZ₂ determined by (a,b,c,d) alone?
print(f"\n{'='*70}")
print("IS ΔZ₂ DETERMINED BY (a,b,c,d)?")
print("=" * 70)

for key in sorted(abcd_stats.keys()):
    dz2_vals = set(e['dZ2'] for e in abcd_stats[key])
    if len(dz2_vals) > 1:
        print(f"  (a,b,c,d)={key}: ΔZ₂ ∈ {dz2_vals} — NOT determined")
    else:
        print(f"  (a,b,c,d)={key}: ΔZ₂ = {dz2_vals.pop()} ✓")

# Check: is ΔΩ₂ determined by (a,b,c,d)?
print(f"\nIS ΔΩ₂ DETERMINED BY (a,b,c,d)?")
for key in sorted(abcd_stats.keys()):
    do2_vals = set(e['dO2'] for e in abcd_stats[key])
    if len(do2_vals) > 1:
        print(f"  (a,b,c,d)={key}: ΔΩ₂ ∈ {do2_vals} — NOT determined")
    else:
        print(f"  (a,b,c,d)={key}: ΔΩ₂ = {do2_vals.pop()} ✓")

# Is there a FORMULA for ΔΩ₂ in terms of a,b,c,d?
print(f"\n{'='*70}")
print("FORMULA SEARCH: ΔΩ₂, ΔΩ₃, ΔZ₂ in terms of (a,b,c,d)")
print("=" * 70)

# Hypothesis: ΔΩ₂ = b - a (gained TT triples minus lost TT triples?)
# Lost 2-paths through u→v:
#   (..., u, v, ...) where u is at some position
#   At position 0: (u,v,w) for w in out(v)\{u} — these are A∪C vertices
#     ... but (u,v,w) is TT iff u→w, which is true for A,C but not D
#   Wait, this is A₂ paths, not Ω₂

# Let me just check formulas
for key in sorted(abcd_stats.keys()):
    a, b, c, d = key
    entries = abcd_stats[key]
    if not entries:
        continue
    # Average ΔΩ₂ (they might vary)
    avg_do2 = sum(e['dO2'] for e in entries) / len(entries)
    avg_do3 = sum(e['dO3'] for e in entries) / len(entries)
    avg_dz2 = sum(e['dZ2'] for e in entries) / len(entries)

    # Try formulas
    f1 = b - a  # simple antisymmetric
    f2 = b + d - a - c  # = -(out_u - out_v - 1)
    f3 = d - c  # another antisymmetric
    f4 = b*d - a*c  # product

    do2_set = set(e['dO2'] for e in entries)
    do3_set = set(e['dO3'] for e in entries)
    dz2_set = set(e['dZ2'] for e in entries)

    if len(do2_set) == 1:
        do2_val = do2_set.pop()
        print(f"  ({a},{b},{c},{d}): ΔΩ₂={do2_val}, b-a={f1}, d-c={f3}, b+d-a-c={f2}")

# Also check the simpler identity: Δ|A₂| = ?
print(f"\n{'='*70}")
print("FORMULA FOR Δ|A₂| AND Δ|A₃|")
print("=" * 70)

da2_by_abcd = defaultdict(set)
da3_by_abcd = defaultdict(set)
for ch in changes:
    key = (ch['a'], ch['b'], ch['c'], ch['d'])
    da2_by_abcd[key].add(ch['dA2'])
    da3_by_abcd[key].add(ch['dA3'])

for key in sorted(da2_by_abcd.keys()):
    a, b, c, d = key
    da2_vals = da2_by_abcd[key]
    da3_vals = da3_by_abcd[key]
    # Is Δ|A₂| = 2(b-a)?
    f_a2 = 2*(b - a)
    f_a3 = 2*(n-3)*(b - a)  # = Δ|A₂| * (n-3)

    if len(da2_vals) == 1:
        da2_val = da2_vals.pop()
        match_a2 = "✓" if da2_val == f_a2 else f"✗ (got {da2_val}, formula={f_a2})"
    else:
        match_a2 = f"varies: {da2_vals}"

    if len(da3_vals) == 1:
        da3_val = da3_vals.pop()
        match_a3 = "✓" if da3_val == f_a3 else f"✗ (got {da3_val}, formula={f_a3})"
    else:
        match_a3 = f"varies: {da3_vals}"

    print(f"  ({a},{b},{c},{d}): Δ|A₂|=2(b-a)? {match_a2}; Δ|A₃|={(n-3)}·Δ|A₂|? {match_a3}")

print("\nDone.")
