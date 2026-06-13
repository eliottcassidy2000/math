#!/usr/bin/env python3
"""
beta2_z2_formula.py — Find a formula for dim(Z₂) and dim(B₂) for tournaments

β₂=0 ⟺ dim(Z₂) = dim(B₂) = rk(∂₃)

Strategy: compute dim(Z₂), rk(∂₃), dim(Ω₂), dim(Ω₃) and various tournament
invariants, then look for an identity.

Key quantities:
- c₃ = number of 3-cycles
- TT = number of transitive triples = C(n,3) - c₃
- |A₂| = total 2-paths = Σᵢ d_i(n-1-d_i)  [where d_i = out-degree]
- J₂ = #{(a,c): c→a, A²[a,c]>0}  (from kind-pasteur)
- NT = 3*c₃ (non-transitive 2-paths)
- dim(Ω₂) = |A₂| - J₂

Author: opus-2026-03-08-S49
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


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


def compute_all(A, n):
    """Compute chain complex data for tournament A on n vertices."""
    ap = {}; om = {}
    for p in range(5):
        ap[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            om[p] = np.eye(n)
        elif ap[p]:
            om[p] = compute_omega_basis(A, n, p, ap[p], ap[p-1])
        else:
            om[p] = np.zeros((0, 0))

    dims = {p: dim_om(om[p]) for p in range(5)}
    rks = {}  # rk(∂_p) for p=1,2,3,4

    for p in range(1, 5):
        if dims[p] == 0:
            rks[p] = 0
            continue
        bd = build_full_boundary_matrix(ap[p], ap[p-1])
        bd_om = bd @ om[p]
        coords = np.linalg.lstsq(om[p-1], bd_om, rcond=None)[0]
        rks[p] = np.linalg.matrix_rank(coords, tol=1e-8)

    z = {p: dims[p] - rks.get(p, 0) for p in range(5)}
    beta = {}
    beta[0] = n - rks.get(1, 0)  # β₀ = n - rk(∂₁)
    for p in range(1, 5):
        beta[p] = z[p] - rks.get(p+1, 0)

    return dims, rks, z, beta, ap


# Tournament invariants
def c3_count(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] + A[j][i] != 1: continue
                if A[i][k] + A[k][i] != 1: continue
                if A[j][k] + A[k][j] != 1: continue
                total = A[i][j] + A[j][k] + A[k][i]
                if total == 0 or total == 3:
                    c += 1
    return c


def J2_count(A, n):
    """J₂ = #{(a,c): c→a AND A²[a,c]>0}"""
    # Compute A²
    A_np = np.array(A)
    A2 = A_np @ A_np
    j2 = 0
    for a in range(n):
        for c in range(n):
            if a == c: continue
            if A[c][a] and A2[a][c] > 0:
                j2 += 1
    return j2


def dt_count(A, n, ap3):
    """Count DT (doubly-transitive) 4-paths."""
    ct = 0
    for path in ap3:
        a, b, c, d = path
        if A[a][c] and A[b][d]:
            ct += 1
    return ct


print("=" * 70)
print("FORMULA SEARCH: dim(Z₂) and rk(∂₃)")
print("=" * 70)

for n in [4, 5, 6]:
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    total = 1 << m

    if n > 5:
        # Sample for n=6
        import random
        random.seed(42)
        samples = random.sample(range(total), min(500, total))
    else:
        samples = range(total)

    print(f"\n{'='*60}")
    print(f"n = {n} ({len(samples)} tournaments)")
    print(f"{'='*60}")

    # Collect data
    data = []
    for bits in samples:
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1

        dims, rks, z, beta, ap = compute_all(A, n)

        scores = tuple(sorted(sum(A[i][j] for j in range(n) if j != i) for i in range(n)))
        c3 = c3_count(A, n)
        j2 = J2_count(A, n)
        a2 = len(ap[2]) if ap[2] else 0
        a3 = len(ap[3]) if ap[3] else 0
        dt = dt_count(A, n, ap[3]) if ap[3] else 0

        data.append({
            'bits': bits, 'scores': scores, 'c3': c3, 'j2': j2,
            'a2': a2, 'a3': a3, 'dt': dt,
            'dims': dims, 'rks': rks, 'z': z, 'beta': beta
        })

        if len(data) % 200 == 0:
            print(f"  ... {len(data)}/{len(samples)}")

    # Test formulas for dim(Z₂) = dim(Ω₂) - rk(∂₂)
    # Since β₁ = dim(Z₁) - rk(∂₂), and dim(Z₁) = C(n-1,2):
    # rk(∂₂) = C(n-1,2) - β₁
    # dim(Z₂) = dim(Ω₂) - rk(∂₂) = (|A₂| - J₂) - (C(n-1,2) - β₁)
    #          = |A₂| - J₂ - C(n-1,2) + β₁

    print(f"\nVerify dim(Z₂) = |A₂| - J₂ - C(n-1,2) + β₁:")
    errors = 0
    for d in data:
        formula = d['a2'] - d['j2'] - (n-1)*(n-2)//2 + d['beta'][1]
        actual = d['z'][2]
        if formula != actual:
            errors += 1
            if errors <= 3:
                print(f"  ERROR: bits={d['bits']}, formula={formula}, actual={actual}")
    print(f"  Errors: {errors}/{len(data)}")

    # Test: dim(Z₂) = rk(∂₃) (i.e., β₂=0)
    print(f"\nVerify dim(Z₂) = rk(∂₃):")
    errors = 0
    for d in data:
        if d['z'][2] != d['rks'].get(3, 0):
            errors += 1
    print(f"  Errors: {errors}/{len(data)}")

    # Test: rk(∂₃) formula
    # rk(∂₃) = dim(Ω₃) - dim(Z₃)
    # Can we express rk(∂₃) = |A₂| - J₂ - C(n-1,2) + β₁?
    print(f"\nVerify rk(∂₃) = |A₂| - J₂ - C(n-1,2) + β₁:")
    errors = 0
    for d in data:
        formula = d['a2'] - d['j2'] - (n-1)*(n-2)//2 + d['beta'][1]
        actual = d['rks'].get(3, 0)
        if formula != actual:
            errors += 1
    print(f"  Errors: {errors}/{len(data)}")

    # Now look for formula for dim(Ω₃)
    print(f"\ndim(Ω₃) distribution:")
    om3_vals = Counter(d['dims'][3] for d in data)
    print(f"  {dict(sorted(om3_vals.items()))}")

    # dim(Ω₃) vs (c₃, scores)
    print(f"\ndim(Ω₃) vs c₃:")
    om3_by_c3 = defaultdict(set)
    for d in data:
        om3_by_c3[d['c3']].add(d['dims'][3])
    for c3 in sorted(om3_by_c3.keys()):
        det = "DET" if len(om3_by_c3[c3]) == 1 else "NOT DET"
        print(f"  c₃={c3}: Ω₃ ∈ {sorted(om3_by_c3[c3])} {det}")

    # Test: dim(Ω₃) = |A₃| - J₃ for some J₃?
    # First, what is |A₃|?
    a3_by_c3 = defaultdict(set)
    for d in data:
        a3_by_c3[d['c3']].add(d['a3'])
    print(f"\n|A₃| vs c₃:")
    for c3 in sorted(a3_by_c3.keys()):
        det = "DET" if len(a3_by_c3[c3]) == 1 else "NOT DET"
        print(f"  c₃={c3}: |A₃| ∈ {sorted(a3_by_c3[c3])} {det}")

    # DT count vs |A₃|
    print(f"\nDT count vs |A₃|:")
    dt_by_a3 = defaultdict(set)
    for d in data:
        dt_by_a3[d['a3']].add(d['dt'])
    for a3 in sorted(dt_by_a3.keys()):
        print(f"  |A₃|={a3}: DT ∈ {sorted(dt_by_a3[a3])}")

    # Test explicit formulas
    # For transitive: |A₂|=C(n,3)*2... no. |A₂| = Σ d_i*(n-1-d_i)
    # For n=5 transitive (scores 0,1,2,3,4): Σ d*(4-d) = 0+3+4+3+0 = 10
    # |A₃| for transitive = Σ over (a,b,c,d) with a<b<c<d = C(n,4)... let me check
    print(f"\n|A₃| for transitive (c₃=0):")
    trans = [d for d in data if d['c3'] == 0]
    if trans:
        print(f"  |A₃| = {trans[0]['a3']}, C(n,4) = {n*(n-1)*(n-2)*(n-3)//24}")
        # How many allowed 3-paths in transitive? a<b<c<d gives C(n,4).
        # But also any permutation that follows arcs: need a→b→c→d
        # In transitive: a→b iff a<b. So 3-path = increasing sequence of length 4 = C(n,4).
        # Wait, that's not right. |A₃| should be more for n=5.
        # For n=5 transitive: 5 = C(5,4). Let me check from the data.

    # |A₃| vs score sequence
    print(f"\n|A₃| by score sequence:")
    a3_by_score = defaultdict(set)
    for d in data:
        a3_by_score[d['scores']].add(d['a3'])
    for sc in sorted(a3_by_score.keys()):
        det = "DET" if len(a3_by_score[sc]) == 1 else "NOT DET"
        print(f"  {sc}: |A₃| ∈ {sorted(a3_by_score[sc])} {det}")

    # KEY: dim(Ω₃) - dim(Z₂) = dim(Z₃) ≥ 0
    # If we could show dim(Ω₃) ≥ dim(Z₂) AND the ∂₃ map is "full rank into Z₂"...
    print(f"\ndim(Ω₃) - dim(Z₂) distribution:")
    excess = Counter(d['dims'][3] - d['z'][2] for d in data)
    print(f"  {dict(sorted(excess.items()))}")
    # Min excess is dim(Z₃). If min > 0, Ω₃ always has room.

    # nullity(∂₃) = dim(Z₃)
    print(f"\ndim(Z₃) distribution:")
    z3_vals = Counter(d['z'][3] for d in data)
    print(f"  {dict(sorted(z3_vals.items()))}")

    # Test: can we express dim(Ω₃) in terms of |A₃|, c₃, scores?
    print(f"\ndim(Ω₃) = |A₃| - X₃ where X₃ is:")
    for d in data[:5]:
        x3 = d['a3'] - d['dims'][3]
        print(f"  bits={d['bits']}, |A₃|={d['a3']}, dim(Ω₃)={d['dims'][3]}, "
              f"X₃={x3}, c₃={d['c3']}")

    # Is X₃ = |A₃| - dim(Ω₃) determined by c₃?
    x3_by_c3 = defaultdict(set)
    for d in data:
        x3 = d['a3'] - d['dims'][3]
        x3_by_c3[d['c3']].add(x3)
    print(f"\nX₃ = |A₃| - dim(Ω₃) vs c₃:")
    for c3 in sorted(x3_by_c3.keys()):
        det = "DET" if len(x3_by_c3[c3]) == 1 else "NOT DET"
        print(f"  c₃={c3}: X₃ ∈ {sorted(x3_by_c3[c3])} {det}")

    # What DOES determine dim(Ω₃)?
    # Maybe dim(Ω₃) = |A₃| - J₃ where J₃ depends on non-transitive quadruples
    print(f"\nLooking for other determining factors...")
    # Compute "non-transitive 4-paths" = A₃ minus DT
    ndt_relation = Counter()
    for d in data:
        ndt = d['a3'] - d['dt']
        om3 = d['dims'][3]
        ndt_relation[(d['c3'], ndt, d['a3'], om3)] += 1

    # dim(Ω₃) by (scores, c₃)
    om3_by_score_c3 = defaultdict(set)
    for d in data:
        om3_by_score_c3[(d['scores'], d['c3'])].add(d['dims'][3])
    determined_sc = sum(1 for v in om3_by_score_c3.values() if len(v) == 1)
    total_sc = len(om3_by_score_c3)
    print(f"\ndim(Ω₃) determined by (score, c₃): {determined_sc}/{total_sc}")

print("\nDone.")
