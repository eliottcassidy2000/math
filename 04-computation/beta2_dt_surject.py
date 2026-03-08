#!/usr/bin/env python3
"""
beta2_dt_surject.py — Prove ∂₃ surjectivity onto Z₂ via DT + cancel structure

KEY STRUCTURAL FACTS:
1. DT 4-paths (a→c, b→d) are INDIVIDUALLY in Ω₃ (all faces in A₂)
2. Non-DT paths sharing a bad face give DIFFERENCES in Ω₃
3. ∂₂∘∂₃ = 0 in full complex ⟹ any A₃-chain with all faces in A₂ is in Ω₃

APPROACH: Understand exactly what determines:
  - dim(Ω₃) in terms of tournament invariants
  - The deficit dim(Z₂) - rk(∂₃|_DT)
  - Why cancel always covers the deficit

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
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


print("=" * 70)
print("DT SURJECTIVITY AND Ω₃ DECOMPOSITION")
print("=" * 70)

# ============================================================
# PART 1: Verify Ω₃ = DT + cancel (is the decomposition exact?)
# ============================================================
print("\nPART 1: Ω₃ decomposition check")

for n in [4, 5]:
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    dt_spans = 0
    dt_cancel_spans = 0
    dt_cancel_fails = 0
    total_nontrivial = 0

    for bits in range(1 << m):
        A = build_adj(n, bits)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        if not ap3: continue

        om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))
        d3 = dim_om(om3)
        if d3 == 0: continue
        total_nontrivial += 1

        # Build DT + cancel vectors in A₃ space
        dt_vecs = []
        bad_groups = defaultdict(list)

        for i, p in enumerate(ap3):
            a, b, c, d = p
            if A[a][c] and A[b][d]:
                v = np.zeros(len(ap3)); v[i] = 1
                dt_vecs.append(v)
            else:
                if not A[a][c]:
                    bad_groups[('02', a, c)].append(i)
                if not A[b][d]:
                    bad_groups[('13', b, d)].append(i)

        cancel_vecs = []
        for key, indices in bad_groups.items():
            for j in range(1, len(indices)):
                v = np.zeros(len(ap3))
                v[indices[0]] = 1; v[indices[j]] = -1
                cancel_vecs.append(v)

        all_vecs = dt_vecs + cancel_vecs
        if not all_vecs:
            dt_cancel_fails += 1
            continue

        V = np.column_stack(all_vecs)
        # Check if om3 columns are in span(V)
        coeffs = np.linalg.lstsq(V, om3, rcond=None)[0]
        residual = np.max(np.abs(om3 - V @ coeffs))

        if residual < 1e-8:
            if not cancel_vecs:
                dt_spans += 1
            else:
                dt_cancel_spans += 1
        else:
            dt_cancel_fails += 1

    print(f"\n  n={n}: {total_nontrivial} tournaments with Ω₃>0")
    print(f"    Ω₃ = DT: {dt_spans}")
    print(f"    Ω₃ = DT+cancel: {dt_cancel_spans}")
    print(f"    Ω₃ ⊋ DT+cancel: {dt_cancel_fails}")


# ============================================================
# PART 2: dim(Ω₃) formula search at n=5
# ============================================================
print(f"\n{'='*70}")
print("PART 2: dim(Ω₃) formula")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

data = []
for bits in range(1 << m):
    A = build_adj(n, bits)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    ap0 = enumerate_allowed_paths(A, n, 0)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)

    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if A[i][j]+A[j][k]+A[k][i] in (0,3))

    n_a3 = len(ap3) if ap3 else 0
    n_dt = sum(1 for p in ap3 if A[p[0]][p[2]] and A[p[1]][p[3]]) if ap3 else 0

    # Count bad face groups with cancellation
    bad_groups = defaultdict(list)
    if ap3:
        for i, p in enumerate(ap3):
            a, b, c, d = p
            if not A[a][c]: bad_groups[('02', a, c)].append(i)
            if not A[b][d]: bad_groups[('13', b, d)].append(i)
    n_cancel = sum(len(g)-1 for g in bad_groups.values() if len(g) >= 2)
    n_bad_faces = len([g for g in bad_groups.values() if len(g) >= 1])
    n_non_dt = n_a3 - n_dt

    # J₃ = dim(A₃) - dim(Ω₃) (the "junk killed" in degree 3)
    j3 = n_a3 - d3

    if d2 > 0:
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2 = d2 - rk2
    else:
        z2 = 0

    data.append({
        'bits': bits, 'scores': scores, 'c3': c3,
        'd2': d2, 'd3': d3, 'z2': z2,
        'n_a3': n_a3, 'n_dt': n_dt, 'n_non_dt': n_non_dt,
        'n_cancel': n_cancel, 'j3': j3, 'n_bad_faces': n_bad_faces,
    })

# Test: J₃ = n_non_dt - n_cancel? (i.e., j3 = singles in bad groups)
n_singles = lambda d: d['n_non_dt'] - d['n_cancel']
test1 = sum(1 for d in data if d['j3'] == n_singles(d))
print(f"\n  J₃ = n_non_dt - n_cancel: {test1}/{len(data)}")

# Test: J₃ = n_bad_faces? (each bad face contributes 1 constraint)
test2 = sum(1 for d in data if d['j3'] == d['n_bad_faces'])
print(f"  J₃ = n_bad_faces: {test2}/{len(data)}")

# Show J₃ vs various quantities
print(f"\n  J₃ distribution:")
j3_dist = Counter(d['j3'] for d in data)
print(f"    {dict(sorted(j3_dist.items()))}")

# J₃ by c₃
j3_by_c3 = defaultdict(set)
for d in data:
    j3_by_c3[d['c3']].add(d['j3'])
print(f"\n  J₃ by c₃:")
for c3 in sorted(j3_by_c3.keys()):
    det = "DET" if len(j3_by_c3[c3]) == 1 else ""
    print(f"    c₃={c3}: J₃ ∈ {sorted(j3_by_c3[c3])} {det}")

# Key relationship: dim(Ω₃) = |DT| + n_cancel?
test3 = sum(1 for d in data if d['d3'] == d['n_dt'] + d['n_cancel'])
print(f"\n  dim(Ω₃) = |DT| + n_cancel: {test3}/{len(data)}")

# What IS dim(Ω₃) - |DT|?
excess = Counter(d['d3'] - d['n_dt'] for d in data)
print(f"  dim(Ω₃) - |DT| distribution: {dict(sorted(excess.items()))}")

# Hmm, n_cancel might double-count. Let's compute the RANK of the cancel space
print(f"\n  Computing rank of cancel subspace...")
for d in data[:20]:
    if d['n_cancel'] > 0 and d['d3'] != d['n_dt']:
        bits = d['bits']
        A = build_adj(n, bits)
        ap3 = enumerate_allowed_paths(A, n, 3)

        bad_groups = defaultdict(list)
        for i, p in enumerate(ap3):
            a, b, c, dd = p
            if not A[a][c]: bad_groups[('02', a, c)].append(i)
            if not A[b][dd]: bad_groups[('13', b, dd)].append(i)

        cancel_vecs = []
        for key, indices in bad_groups.items():
            for j in range(1, len(indices)):
                v = np.zeros(len(ap3))
                v[indices[0]] = 1; v[indices[j]] = -1
                cancel_vecs.append(v)

        if cancel_vecs:
            C = np.column_stack(cancel_vecs)
            rk_cancel = np.linalg.matrix_rank(C, tol=1e-8)
        else:
            rk_cancel = 0

        om3_excess = d['d3'] - d['n_dt']
        print(f"    bits={bits}: Ω₃={d['d3']}, DT={d['n_dt']}, cancel_pairs={d['n_cancel']}, "
              f"rk_cancel={rk_cancel}, excess={om3_excess}")


# ============================================================
# PART 3: The "full rank" question — why rk(∂₃) = dim(Z₂)?
# ============================================================
print(f"\n{'='*70}")
print("PART 3: WHY rk(∂₃) = dim(Z₂)")
print("=" * 70)

# Key identity: rk(∂₃) = dim(Ω₃) - dim(Z₃)
# β₂ = 0 ⟺ dim(Ω₃) = dim(Z₂) + dim(Z₃)
#
# So the question is: WHY does this hold?
#
# From formulas:
#   dim(Z₂) = |A₂| - J₂ - C(n-1,2) + β₁
#   dim(Ω₃) = |A₃| - J₃
#   dim(Z₃) = dim(Ω₃) - rk(∂₃)
#
# If β₂ = 0: |A₃| - J₃ = (|A₂| - J₂ - C(n-1,2) + β₁) + dim(Z₃)
#
# Using |A₃| = C(n,4) + 2(n-3)c₃:
# C(n,4) + 2(n-3)c₃ - J₃ = |A₂| - J₂ - C(n-1,2) + β₁ + dim(Z₃)
#
# Rearrange: J₃ = C(n,4) + 2(n-3)c₃ - |A₂| + J₂ + C(n-1,2) - β₁ - dim(Z₃)

print("\nVerifying J₃ formula for n=5:")
for d in data:
    bits = d['bits']
    A = build_adj(n, bits)

    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    ap4 = enumerate_allowed_paths(A, n, 4)
    ap0 = enumerate_allowed_paths(A, n, 0)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3_mat = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))
    om4 = compute_omega_basis(A, n, 4, ap4, ap3) if ap4 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3_mat)
    d4 = dim_om(om4)

    # Compute all ranks
    if d3 > 0 and d2 > 0:
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd3_om = bd3 @ om3_mat
        coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
        rk3 = np.linalg.matrix_rank(coords3, tol=1e-8)
        z3 = d3 - rk3
    else:
        rk3 = 0
        z3 = d3

    if d4 > 0 and d3 > 0:
        bd4 = build_full_boundary_matrix(ap4, ap3)
        bd4_om = bd4 @ om4
        coords4 = np.linalg.lstsq(om3_mat, bd4_om, rcond=None)[0]
        rk4 = np.linalg.matrix_rank(coords4, tol=1e-8)
    else:
        rk4 = 0

    beta1 = d['z2'] - rk3 if d['d2'] > 0 else 0  # Wrong, let me recompute
    # Actually β₁ = dim(Z₁) - rk(∂₂) = C(n-1,2) - (d2 - z2) = C(n-1,2) - d2 + z2
    # Wait: rk(∂₂) = d2 - z2? No. rk(∂₂) = dim(Ω₂) - dim(Z₂) = d2 - d['z2'].
    # β₁ = dim(Z₁) - rk(∂₂) = C(n-1,2) - (d2 - d['z2'])
    beta1 = (n-1)*(n-2)//2 - (d['d2'] - d['z2'])

    a2 = len(ap2) if ap2 else 0
    j2 = a2 - d['d2']

    # J₃ predicted
    cn4 = n*(n-1)*(n-2)*(n-3)//24
    cn12 = (n-1)*(n-2)//2
    j3_predicted = cn4 + 2*(n-3)*d['c3'] - a2 + j2 + cn12 - beta1 - z3
    j3_actual = d['j3']

    d['j3_predicted'] = j3_predicted
    d['beta1'] = beta1
    d['z3'] = z3

errors = sum(1 for d in data if d['j3'] != d['j3_predicted'])
print(f"  Errors in J₃ = C(n,4)+2(n-3)c₃-|A₂|+J₂+C(n-1,2)-β₁-dim(Z₃): {errors}/{len(data)}")

if errors == 0:
    print("  ✓ Identity holds! (But it's circular — assumes β₂=0)")

# So this identity is EQUIVALENT to β₂=0, not a proof.
# Let me look for an INDEPENDENT formula for J₃.

# What actually determines J₃?
print(f"\n  J₃ by (scores, c₃, β₁):")
j3_by_triple = defaultdict(set)
for d in data:
    key = (d['scores'], d['c3'], d['beta1'])
    j3_by_triple[key].add(d['j3'])

determined = sum(1 for v in j3_by_triple.values() if len(v) == 1)
print(f"    Determined by (scores, c₃, β₁): {determined}/{len(j3_by_triple)}")

# What about by c₃ alone?
j3_by_c3_beta = defaultdict(set)
for d in data:
    j3_by_c3_beta[(d['c3'], d['beta1'])].add(d['j3'])
print(f"\n  J₃ by (c₃, β₁):")
for key in sorted(j3_by_c3_beta.keys()):
    det = "DET" if len(j3_by_c3_beta[key]) == 1 else ""
    print(f"    (c₃={key[0]}, β₁={key[1]}): J₃ ∈ {sorted(j3_by_c3_beta[key])} {det}")


# ============================================================
# PART 4: Direct formula for J₃
# ============================================================
print(f"\n{'='*70}")
print("PART 4: Direct combinatorial formula for J₃")
print("=" * 70)

# J₃ = |A₃| - dim(Ω₃) = #{3-paths not in Ω₃}
# A 3-path (a,b,c,d) is NOT in Ω₃ iff its boundary has components outside A₂
# That happens when face (a,c,d) ∉ A₂ [i.e., c→a] OR face (a,b,d) ∉ A₂ [i.e., d→b]
#
# Wait: not exactly. Being "not in Ω₃" for a SINGLE path means: the indicator vector
# e_{(a,b,c,d)} has boundary with non-A₂ components. But Ω₃ allows linear combinations
# where non-A₂ components cancel.
#
# J₃ = #{independent constraints} = #{non-A₂ faces that appear in boundaries}
#                                     - #{dependencies among those constraints}
#
# Non-A₂ face of type (a,c,d): appears when c→a (face from removing vertex 1)
# Non-A₂ face of type (a,b,d): appears when d→b (face from removing vertex 2)
#
# Wait, the four faces of (a,b,c,d) are:
# (b,c,d) with sign +1 — in A₂ iff b→c (yes, given) and c→d (yes, given)
# (a,c,d) with sign -1 — in A₂ iff a→c and c→d (c→d given; a→c iff DT cross-arc)
# (a,b,d) with sign +1 — in A₂ iff a→b (yes) and b→d (DT cross-arc)
# (a,b,c) with sign -1 — in A₂ iff a→b (yes) and b→c (yes) → ALWAYS in A₂

# So the non-A₂ faces are:
# face02 = (a,c,d) when c→a: this is the "02-bad" face
# face13 = (a,b,d) when d→b: this is the "13-bad" face

# J₃ = #{independent 02-bad constraints} + #{independent 13-bad constraints}
# Each 02-bad constraint is: for non-arc (a,c), Σ_b α_{(a,b,c,d)} = 0 for each d
# Wait, the constraint is that the coefficient of the non-A₂ face (a,c,d) in ∂₃(x)
# must be 0. This face appears from paths (a,b,c,d) where b varies.
#
# For fixed (a,c,d) with c→a and c→d:
# coeff of (a,c,d) in ∂₃(Σ α_p · p) = Σ_{b: (a,b,c,d)∈A₃} (-1) · α_{(a,b,c,d)}
# = -Σ_b α_{(a,b,c,d)}
# This must be 0: Σ_b α_{(a,b,c,d)} = 0
#
# How many such b exist? b must satisfy a→b, b→c (so b ∈ out(a)∩in(c)\{c,d,a})
# The number of valid b values = #{w: a→w, w→c, w≠a,c,d} = A²[a,c] - (1 if a→d→c else 0)
# Hmm, this is getting complicated.

# Let me just COUNT the non-A₂ faces directly.

print("\nCounting non-A₂ faces of A₃ paths for n=5:")

for d in data[:20]:
    bits = d['bits']
    A = build_adj(n, bits)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap3: continue

    # Count distinct non-A₂ faces
    bad02_faces = set()
    bad13_faces = set()
    for p in ap3:
        a, b, c, dd = p
        if not A[a][c]:
            bad02_faces.add((a, c, dd))
        if not A[b][dd]:
            bad13_faces.add((a, b, dd))

    n_bad02 = len(bad02_faces)
    n_bad13 = len(bad13_faces)
    total_bad = n_bad02 + n_bad13

    # Are there dependencies? Each bad face gives one linear constraint.
    # But constraints might be dependent if paths share multiple bad faces.
    # J₃ should be the number of INDEPENDENT constraints.

    if total_bad != d['j3']:
        print(f"  bits={bits}: J₃={d['j3']}, #bad02={n_bad02}, #bad13={n_bad13}, "
              f"total={total_bad}, diff={d['j3']-total_bad}")
    elif total_bad == d['j3'] and d['j3'] > 0:
        pass  # Match

# Test: J₃ = n_bad02_faces + n_bad13_faces?
print("\nJ₃ = #distinct_bad02 + #distinct_bad13?")
match_count = 0
mismatch_count = 0
for d in data:
    bits = d['bits']
    A = build_adj(n, bits)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap3: continue

    bad02 = set()
    bad13 = set()
    for p in ap3:
        a, b, c, dd = p
        if not A[a][c]: bad02.add((a, c, dd))
        if not A[b][dd]: bad13.add((a, b, dd))

    predicted = len(bad02) + len(bad13)
    if predicted == d['j3']:
        match_count += 1
    else:
        mismatch_count += 1
        if mismatch_count <= 5:
            print(f"  MISMATCH: bits={bits}, J₃={d['j3']}, #bad02+#bad13={predicted}")

print(f"  Match: {match_count}, Mismatch: {mismatch_count}")

if mismatch_count == 0:
    print("  ✓ J₃ = #{distinct bad02 faces} + #{distinct bad13 faces}")
    print("     Each non-A₂ face gives exactly ONE independent constraint!")

    # Now: can we count #bad02 and #bad13 combinatorially?
    # bad02 face (a,c,d): c→a, c→d, and ∃b with a→b, b→c (so A²[a,c]>0)
    # bad13 face (a,b,d): d→b, a→b, and ∃c with b→c, c→d (so A²[b,d]>0)
    #
    # Note symmetry: bad13 is the "reverse" of bad02!
    # bad02: (a,c,d) with c→a, c→d, A²[a,c]>0
    # bad13: (a,b,d) with d→b, a→b, A²[b,d]>0
    #
    # Count bad02: for each (a,c) with c→a and A²[a,c]>0 (i.e., junk pair),
    #   and for each d with c→d:
    #   this gives a bad02 face iff there exists a 3-path (a,b,c,d) in A₃
    #   i.e., ∃b with a→b, b→c, c→d → always true since A²[a,c]>0 ensures ∃b
    #   So #bad02 = Σ_{(a,c): c→a, A²[a,c]>0} |{d: c→d, d∉{a,c}}|
    #            = Σ_{(a,c)∈junk} (out-degree(c) minus 1 if c→a counted)
    #   Actually d can be anything with c→d and d∉{a,b,c}... but d∉{a,c} suffices
    #   since b will be different.
    #   Wait, d must also not equal any b. But there can be multiple b's, so d
    #   just needs d ∉ {a,c}? Actually d also needs d ∉ {b} but since b varies,
    #   the constraint is d ∉ {a,c} and c→d.
    #   No wait: for a given (a,c,d), the 3-path is (a,b,c,d) for SOME b.
    #   We need b ∉ {a,c,d}. And a→b, b→c. So we need at least one such b.
    #   A²[a,c] counts ALL b with a→b, b→c. But some of these b might equal d.
    #   So the constraint is: A²[a,c] - (1 if a→d and d→c) ≥ 1.
    #   For the bad face to exist: ∃ b∈out(a)∩in(c)\{a,c,d}.

    print(f"\n  Computing #bad02 formula...")

    # Verify: #bad02 = Σ_{(a,c): c→a, A²[a,c]>0} #{d: c→d, d∉{a,c}, ∃b∈S(a,c)\{d}}
    for d_item in data[:10]:
        bits = d_item['bits']
        A = build_adj(n, bits)
        ap3 = enumerate_allowed_paths(A, n, 3)
        if not ap3: continue

        bad02_actual = set()
        for p in ap3:
            a, b, c, dd = p
            if not A[a][c]: bad02_actual.add((a, c, dd))

        bad02_formula = set()
        A_np = np.array(A)
        A2 = (A_np @ A_np).astype(int)
        for a in range(n):
            for c in range(n):
                if a == c: continue
                if not A[c][a]: continue  # need c→a
                if A2[a][c] == 0: continue  # need A²[a,c]>0
                S_ac = [b for b in range(n) if b not in (a,c) and A[a][b] and A[b][c]]
                for dd in range(n):
                    if dd in (a,c): continue
                    if not A[c][dd]: continue  # need c→d
                    # Need ∃b ∈ S(a,c) \ {d}
                    if any(b != dd for b in S_ac):
                        bad02_formula.add((a, c, dd))

        if bad02_actual != bad02_formula:
            print(f"    MISMATCH at bits={bits}")
        else:
            print(f"    bits={bits}: #bad02={len(bad02_actual)} ✓")


print("\nDone.")
