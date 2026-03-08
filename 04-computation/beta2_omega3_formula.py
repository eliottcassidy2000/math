#!/usr/bin/env python3
"""
beta2_omega3_formula.py — Find exact formula for dim(Ω₃)

KEY RESULT from surplus_algebraic.py:
  dim(Ω₂) = |TT| + |NT| - |bad_faces₂|  (EXACT)

Now we need dim(Ω₃). From the data:
  - Individual elements of Ω₃ = ALL DT 3-paths (regardless of a↔d direction)
  - dim(Ω₃) varies even for fixed (super_DT, DT_back) when DT_back=3

The analog of the Ω₂ formula should be:
  dim(Ω₃) = |DT| + (non-DT cancellation dimension)

Non-DT 3-paths have "bad faces" — faces not in Ω₂.
The cancellation dimension = |non-DT_in_Ω₃| = number of non-DT elements
whose bad-face terms cancel.

For a 3-path (a,b,c,d) that is NOT DT:
  - Either a not→c (missing first DT condition)
  - Or b not→d (missing second DT condition)
  - Or both

The faces are:
  (b,c,d) — this face needs to be in Ω₂. It's TT if b→d.
  (a,c,d) — TT if a→d. In Ω₂ individually if TT; in Ω₂ as combo if NT.
  (a,b,d) — TT if a→d. Same.
  (a,b,c) — TT if a→c. Same.

The face (a,b,c) is in Ω₂ individually if a→c (TT triple). It's NOT in Ω₂
individually if c→a. But as part of a cancellation combination in Ω₂, it
could still work — the "bad face" (a,c) must cancel with another element.

Actually, the Ω₃ condition is: ∂₃(x) ∈ Ω₂. For an individual 3-path,
this means ALL faces must be in Ω₂. But A₂ = TT ∪ NT triples (allowed paths).
Ω₂ = TT ⊕ NT_cancel (as a subspace of span(A₂)).

A face of a 3-path is itself a 2-path. If it's an allowed path AND it's TT,
it's in Ω₂. If it's allowed but NT, it's NOT individually in Ω₂ (but could
be part of a cancellation combo).

Wait: ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c).
Each face is a 2-path. For this to be in Ω₂, each 2-path must be in A₂
(otherwise the coefficient is on a basis element not in Ω₂'s ambient space)
OR the coefficients on non-A₂ elements must cancel.

A 2-path (x,y,z) is in A₂ iff x→y AND y→z.
Face (b,c,d): b→c (from a→b→c→d path), c→d (same). So (b,c,d) ∈ A₂ always ✓
Face (a,b,c): a→b ✓, b→c ✓. So (a,b,c) ∈ A₂ always ✓
Face (a,c,d): a→c needed AND c→d ✓. So ∈ A₂ iff a→c.
Face (a,b,d): a→b ✓, b→d needed. So ∈ A₂ iff b→d.

So a 3-path (a,b,c,d) has:
  - 0 bad faces: a→c AND b→d (DT path) → ALWAYS in Ω₃
  - 1 bad face from (a,c,d): c→a AND b→d  [one-bad type 1]
  - 1 bad face from (a,b,d): a→c AND d→b  [one-bad type 2]
  - 2 bad faces: c→a AND d→b  [two-bad]

For the boundary ∂₃(p) to be in Ω₂:
  All face coefficients must project into Ω₂.
  Since (b,c,d) and (a,b,c) are always in A₂ and hence span(A₂),
  their coefficients just need to be in the Ω₂ subspace.

  But (a,c,d) when c→a: (a,c,d) is NOT in A₂. The term -(a,c,d) is a coefficient
  on a non-A₂ element. This means ∂₃(p) has a component outside span(A₂),
  which is outside Ω₂.

  UNLESS... wait. Ω₂ is a subspace of span(A₂). But ∂₃ maps into span(A₂).

  Actually no. ∂₃ maps Ω₃ → Ω₂, but the FACE MAP ∂₃ on A₃ maps into
  R^{vertices^3}, not into span(A₂). The GLMY boundary ∂₃ maps Ω₃ → Ω₂.

  But computationally (in path_homology_v2.py), the boundary matrix is built
  as a map from A₃ basis to A₂ basis. If face (a,c,d) is NOT in A₂,
  it doesn't appear in the boundary matrix.

  So the boundary matrix bd3[i,j] records the coefficient of the i-th A₂ path
  in ∂₃ of the j-th A₃ path. Non-A₂ faces are DROPPED.

  This means: the effective boundary ∂₃^A maps A₃ → span(A₂), and Ω₃ is defined
  by: p ∈ A₃ such that ∂₃^A(p) ∈ Ω₂ (where Ω₂ ⊆ span(A₂)).

  For DT paths: all faces are in A₂, so ∂₃^A includes all terms. Since all faces
  are TT (or in Ω₂), the boundary is in span(TT) ⊆ Ω₂. ✓

  For non-DT paths: some faces are NOT in A₂ and are dropped. The remaining
  faces might or might not sum to something in Ω₂.

Let me verify this computationally.

Author: opus-2026-03-08-S44
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

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

print("=" * 70)
print("Ω₃ FORMULA: DT + NON-DT CANCELLATION")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n = {n} ---")

    formula_checks = Counter()
    bad_face_analysis = defaultdict(list)

    for A in all_tournaments(n):
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)

        if not a3:
            continue

        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        # Classify 3-paths
        dt = 0
        one_bad = 0
        two_bad = 0
        for p in a3:
            a, b, c, d = p
            ac = A[a][c]  # 1 if a→c
            bd = A[b][d]  # 1 if b→d
            if ac and bd:
                dt += 1
            elif ac or bd:
                one_bad += 1
            else:
                two_bad += 1

        # Non-DT 3-paths with their "bad face" structure
        # For one-bad (c→a, b→d): bad face is (a,c,d) which is NOT in A₂
        # For one-bad (a→c, d→b): bad face is (a,b,d) which is NOT in A₂
        # For two-bad (c→a, d→b): bad faces are BOTH (a,c,d) and (a,b,d)

        # Bad face for (a,c,d): needs a→c to be in A₂. Bad iff c→a.
        # Bad face for (a,b,d): needs b→d to be in A₂. Bad iff d→b.

        # Bad face TYPE 1: pair (a,c) with c→a
        # Bad face TYPE 2: pair (b,d) with d→b

        bf_type1 = set()  # (a,c) pairs where c→a
        bf_type2 = set()  # (b,d) pairs where d→b

        for p in a3:
            a, b, c, d = p
            if A[a][c] == 0:  # c→a, face (a,c,d) is bad
                bf_type1.add((a, c))
            if A[b][d] == 0:  # d→b, face (a,b,d) is bad
                bf_type2.add((b, d))

        # Count how many non-DT paths share each bad face
        bf1_counts = defaultdict(int)
        bf2_counts = defaultdict(int)
        for p in a3:
            a, b, c, d = p
            if A[a][c] == 0:
                bf1_counts[(a, c)] += 1
            if A[b][d] == 0:
                bf2_counts[(b, d)] += 1

        # NT cancellation for Ω₃:
        # Non-DT paths that share a bad face can cancel.
        # The cancellation dimension = Σ (count - 1) over distinct bad faces
        # But it's more complex because a two-bad path has TWO bad faces.

        # Simple formula: dim(Ω₃) = |DT| + Σ_{bf} max(0, count(bf) - 1)
        # where bf ranges over all distinct bad face pairs?

        cancel_dim_1 = sum(max(0, c - 1) for c in bf1_counts.values())
        cancel_dim_2 = sum(max(0, c - 1) for c in bf2_counts.values())

        # But bad faces of type 1 and type 2 are independent...
        # A two-bad path contributes to both bf1 and bf2 counts.
        # Need to handle this carefully.

        # Actually, let's just check the formula numerically.
        # Hypothesis: dim(Ω₃) = |DT| + cancel_dim_eff
        cancel_eff = dim_om3 - dt

        bad_face_analysis[(one_bad, two_bad, len(bf_type1), len(bf_type2))].append(cancel_eff)

        # Check various formulas
        if dim_om3 == dt + cancel_dim_1:
            formula_checks['DT + cancel1'] += 1
        if dim_om3 == dt + cancel_dim_2:
            formula_checks['DT + cancel2'] += 1
        if dim_om3 == dt + cancel_dim_1 + cancel_dim_2:
            formula_checks['DT + cancel1+2'] += 1

        # Combined bad face analysis
        all_bad_faces = set()
        bf_combined_counts = defaultdict(int)
        for p in a3:
            a, b, c, d = p
            bfs = []
            if A[a][c] == 0:
                bfs.append(('type1', a, c))
            if A[b][d] == 0:
                bfs.append(('type2', b, d))
            for bf in bfs:
                bf_combined_counts[bf] += 1
                all_bad_faces.add(bf)

        cancel_combined = sum(max(0, c - 1) for c in bf_combined_counts.values())
        if dim_om3 == dt + cancel_combined:
            formula_checks['DT + cancel_combined'] += 1

        # Try: non-DT in Ω₃ = total non-DT - #distinct bad faces
        non_dt_total = one_bad + two_bad
        num_distinct_bf = len(all_bad_faces)
        if dim_om3 == dt + non_dt_total - num_distinct_bf:
            formula_checks['DT + nonDT - #bf'] += 1

        formula_checks['total'] += 1

    print(f"  Formula checks:")
    total = formula_checks['total']
    for k, v in sorted(formula_checks.items()):
        if k != 'total':
            print(f"    {k}: {v}/{total}")

    if n == 5:
        print(f"\n  Bad face structure (one_bad, two_bad, #bf1, #bf2) → cancel_eff:")
        for key in sorted(bad_face_analysis):
            vals = bad_face_analysis[key]
            uniqs = sorted(set(vals))
            print(f"    {key}: cancel_eff ∈ {uniqs} (count={len(vals)})")

# ===== VERIFY: dim(Ω₃) = |A₃| - |bad_faces_3| =====
print(f"\n{'='*70}")
print("TESTING dim(Ω₃) = |A₃| - |distinct bad faces of A₃ paths|")
print("="*70)

for n in [4, 5]:
    print(f"\n  n = {n}:")
    matches = 0
    total = 0

    for A in all_tournaments(n):
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        if not a3:
            continue

        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        # Count distinct bad faces
        all_bad_faces = set()
        for p in a3:
            a, b, c, d = p
            if A[a][c] == 0:
                all_bad_faces.add(('ac', a, c, d))  # face (a,c,d) not in A₂
            if A[b][d] == 0:
                all_bad_faces.add(('bd', a, b, d))  # face (a,b,d) not in A₂

        # Hmm, the bad face is a SPECIFIC 2-path, not just a pair.
        # (a,c,d): this is the 2-path a→c→d where the issue is a→c.
        # But if c→a, then (a,c,d) is NOT an allowed path.
        # The "bad face" as a formal element is e_{(a,c,d)} in the abstract chain.
        # But in the boundary matrix restricted to A₂, it's DROPPED.
        #
        # So the relevant object is: which non-A₂ 2-paths appear as faces?
        # Two 3-paths share a bad face if the SAME non-A₂ 2-path appears.

        # Let me redo with the ACTUAL face identification
        non_a2_faces = set()
        face_counts = defaultdict(int)
        for p in a3:
            a, b, c, d = p
            # Face (a,c,d): in A₂ iff a→c AND c→d (c→d always true from 3-path)
            # So in A₂ iff a→c
            if A[a][c] == 0:  # c→a, face (a,c,d) NOT in A₂
                non_a2_faces.add((a,c,d))
                face_counts[(a,c,d)] += 1
            # Face (a,b,d): in A₂ iff a→b AND b→d (a→b always true from 3-path)
            # So in A₂ iff b→d
            if A[b][d] == 0:  # d→b, face (a,b,d) NOT in A₂
                non_a2_faces.add((a,b,d))
                face_counts[(a,b,d)] += 1

        num_distinct_bad = len(non_a2_faces)
        predicted = len(a3) - num_distinct_bad

        if predicted == dim_om3:
            matches += 1
        total += 1

    print(f"    dim(Ω₃) = |A₃| - |distinct non-A₂ faces|: {matches}/{total}")

# ===== DEEPER: Ω₂ and Ω₃ in terms of bad face incidence =====
print(f"\n{'='*70}")
print("UNIFIED FORMULA: dim(Ω_p) = |A_p| - |distinct non-A_{p-1} faces|")
print("="*70)

for n in [4, 5]:
    print(f"\n  n = {n}:")

    for p_test in [2, 3]:
        matches = 0
        total = 0

        for A in all_tournaments(n):
            ap = enumerate_allowed_paths(A, n, p_test)
            ap_prev = enumerate_allowed_paths(A, n, p_test - 1)

            if not ap:
                continue

            omp = compute_omega_basis(A, n, p_test, ap, ap_prev)
            dim_omp = omp.shape[1] if omp.ndim == 2 else 0

            # Count distinct non-A_{p-1} faces
            ap_prev_set = set(tuple(x) for x in ap_prev)

            non_faces = set()
            for path in ap:
                for face_idx in range(len(path)):
                    face = tuple(path[:face_idx] + path[face_idx+1:])
                    if face not in ap_prev_set:
                        non_faces.add(face)

            predicted = len(ap) - len(non_faces)

            if predicted == dim_omp:
                matches += 1
            total += 1

        print(f"    p={p_test}: dim(Ω_{p_test}) = |A_{p_test}| - |non-A_{p_test-1} faces|: {matches}/{total}")

# ===== VERIFY THE UNIFIED FORMULA =====
print(f"\n{'='*70}")
print("VERIFY FOR p=2,3,4 AT n=5")
print("="*70)

n = 5
for p_test in [2, 3, 4]:
    matches = 0
    total = 0
    mismatches = []

    for A in all_tournaments(n):
        ap = enumerate_allowed_paths(A, n, p_test)
        ap_prev = enumerate_allowed_paths(A, n, p_test - 1)

        if not ap:
            continue

        omp = compute_omega_basis(A, n, p_test, ap, ap_prev)
        dim_omp = omp.shape[1] if omp.ndim == 2 else 0

        ap_prev_set = set(tuple(x) for x in ap_prev)

        non_faces = set()
        for path in ap:
            for face_idx in range(len(path)):
                face = tuple(path[:face_idx] + path[face_idx+1:])
                if face not in ap_prev_set:
                    non_faces.add(face)

        predicted = len(ap) - len(non_faces)

        if predicted == dim_omp:
            matches += 1
        else:
            if len(mismatches) < 3:
                mismatches.append((dim_omp, predicted, len(ap), len(non_faces)))
        total += 1

    print(f"  p={p_test}: {matches}/{total}")
    if mismatches:
        for mm in mismatches:
            print(f"    Mismatch: actual={mm[0]}, predicted={mm[1]}, |A_p|={mm[2]}, |non-faces|={mm[3]}")

# ===== SURPLUS IN TERMS OF BAD FACES =====
print(f"\n{'='*70}")
print("SURPLUS = dim(Ω₃) - Z₂ IN TERMS OF BAD FACES")
print("="*70)

n = 5
print(f"\n  If dim(Ω_p) = |A_p| - |BF_p| where BF_p = non-A_{p-1} faces of A_p paths,")
print(f"  then:")
print(f"    surplus = dim(Ω₃) - Z₂")
print(f"           = (|A₃| - |BF₃|) - (dim(Ω₂) - rk(∂₂))")
print(f"           = (|A₃| - |BF₃|) - (|A₂| - |BF₂|) + rk(∂₂)")
print(f"           = (|A₃| - |A₂|) - (|BF₃| - |BF₂|) + rk(∂₂)")
print()

# Verify these components
rk2_vs_bf = defaultdict(list)
surplus_formula = 0
total = 0

for A in all_tournaments(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    if not a2 or not a3:
        continue

    # Bad faces for A₂
    a1_set = set(tuple(x) for x in a1)
    bf2 = set()
    for path in a2:
        for fi in range(3):
            face = tuple(list(path[:fi]) + list(path[fi+1:]))
            if face not in a1_set:
                bf2.add(face)

    # Bad faces for A₃
    a2_set = set(tuple(x) for x in a2)
    bf3 = set()
    for path in a3:
        for fi in range(4):
            face = tuple(list(path[:fi]) + list(path[fi+1:]))
            if face not in a2_set:
                bf3.add(face)

    # Compute actual dimensions
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    om1 = compute_omega_basis(A, n, 1, a1, enumerate_allowed_paths(A, n, 0))
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    S = np.linalg.svd(coords, compute_uv=False)
    rk2 = int(sum(s > 1e-8 for s in S))
    z2 = dim_om2 - rk2

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    surplus = dim_om3 - z2

    # Formula check
    predicted = (len(a3) - len(a2)) - (len(bf3) - len(bf2)) + rk2
    if predicted == surplus:
        surplus_formula += 1
    total += 1

    rk2_vs_bf[(len(bf2), len(bf3))].append(rk2)

print(f"  surplus = (|A₃|-|A₂|) - (|BF₃|-|BF₂|) + rk(∂₂): {surplus_formula}/{total}")

# ===== WHAT ABOUT A₂ BAD FACES? =====
# For tournaments, A₁ = ALL directed edges. So ALL faces of A₂ paths are in A₁.
# That means BF₂ = 0!
print(f"\n  BF₂ for tournaments (faces of 2-paths not in A₁):")
for A in all_tournaments(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a1_set = set(tuple(x) for x in a1)
    bf2 = set()
    for path in a2:
        for fi in range(3):
            face = tuple(list(path[:fi]) + list(path[fi+1:]))
            if face not in a1_set:
                bf2.add(face)
    if bf2:
        print(f"    FOUND non-empty BF₂!")
        break
else:
    print(f"    BF₂ = ∅ for ALL tournaments (as expected: A₁ = all edges)")

print(f"\n  So for tournaments: dim(Ω₂) = |A₂| (Ω₂ = A₂)? Let's check...")

for A in all_tournaments(n):
    a2 = enumerate_allowed_paths(A, n, 2)
    a1 = enumerate_allowed_paths(A, n, 1)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 != len(a2):
        print(f"    MISMATCH: dim(Ω₂) = {dim_om2} ≠ |A₂| = {len(a2)}")
        break
else:
    print(f"    dim(Ω₂) = |A₂| for ALL tournaments at n={n}!")

# Wait — the Ω₂ formula earlier said dim(Ω₂) = TT + NT_cancel, and
# NT_cancel < NT. But we're now saying Ω₂ ⊆ span(A₂) with dim = |A₂|??
# That means Ω₂ = span(A₂) = R^{|A₂|}!

# Actually: Ω₂ = {c ∈ span(A₂) : ∂₂(c) ∈ Ω₁}.
# For tournaments, Ω₁ = span(A₁) = R^{n(n-1)/2}.
# ∂₂(a,b,c) = (b,c) - (a,c) + (a,b). ALL three terms are in A₁ for tournaments
# (since A₁ = all directed edges). So ∂₂(a,b,c) ∈ span(A₁) = Ω₁.
# Therefore EVERY element of span(A₂) has boundary in Ω₁.
# So Ω₂ = span(A₂)!

# Wait, that contradicts the earlier finding that Ω₂ = TT + NT_cancel ≠ A₂.
# Let me recheck...

print(f"\n  WAIT: Re-checking the definition of Ω₂...")
print(f"  Ω₂ = {{c ∈ span(A₂) : ∂₂(c) ∈ Ω₁}}")
print(f"  For tournaments: A₁ = all edges, Ω₁ = span(A₁)")
print(f"  ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)")
print(f"  For (a,b,c) ∈ A₂: a→b, b→c.")
print(f"  Face (a,b): a→b, so (a,b) ∈ A₁. ✓")
print(f"  Face (b,c): b→c, so (b,c) ∈ A₁. ✓")
print(f"  Face (a,c): might have a→c or c→a. Either way, the DIRECTED edge is in A₁.")
print(f"  But ∂₂ writes -(a,c), NOT -(c,a). If c→a, then (a,c) ∉ A₁.")
print(f"  Is (a,c) in A₁? NO if c→a. The edge (c,a) is in A₁, not (a,c).")
print(f"  So ∂₂(a,b,c) contains -(a,c) which is NOT in span(A₁).")
print(f"  Therefore individual NT 2-paths are NOT in Ω₂.")
print(f"  But the formula says dim(Ω₂) = |A₂|??")
print(f"\n  This is a CONTRADICTION. Let me recompute carefully...")

# Direct computation
A = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(i+1, 5):
        A[i][j] = 1  # transitive tournament

a2 = enumerate_allowed_paths(A, 5, 2)
a1 = enumerate_allowed_paths(A, 5, 1)
print(f"\n  Transitive T₅: |A₂| = {len(a2)}, |A₁| = {len(a1)}")

om2 = compute_omega_basis(A, 5, 2, a2, a1)
print(f"  dim(Ω₂) = {om2.shape[1] if om2.ndim == 2 else 0}")
print(f"  TT count = {sum(1 for p in a2 if A[p[0]][p[2]] == 1)}")
print(f"  NT count = {sum(1 for p in a2 if A[p[0]][p[2]] == 0)}")

# For transitive, ALL triples are TT (a→c for all a<c). So NT=0, Ω₂ = TT = |A₂|. ✓
# The formula dim(Ω₂) = |A₂| would only hold when there are no NT triples.
# With NT triples, dim(Ω₂) < |A₂|.

# Recheck with a cyclic tournament
A_cyc = [[0]*5 for _ in range(5)]
for i in range(5):
    A_cyc[i][(i+1)%5] = 1
    A_cyc[i][(i+2)%5] = 1

a2_cyc = enumerate_allowed_paths(A_cyc, 5, 2)
a1_cyc = enumerate_allowed_paths(A_cyc, 5, 1)
om2_cyc = compute_omega_basis(A_cyc, 5, 2, a2_cyc, a1_cyc)
tt_cyc = sum(1 for p in a2_cyc if A_cyc[p[0]][p[2]] == 1)
nt_cyc = len(a2_cyc) - tt_cyc
dim_om2_cyc = om2_cyc.shape[1] if om2_cyc.ndim == 2 else 0

print(f"\n  Cyclic C₅^{{1,2}}: |A₂| = {len(a2_cyc)}, TT = {tt_cyc}, NT = {nt_cyc}")
print(f"  dim(Ω₂) = {dim_om2_cyc}")
print(f"  |A₂| - dim(Ω₂) = {len(a2_cyc) - dim_om2_cyc}")

# So the "BF₂ = 0" assertion was WRONG. Let me re-examine.
# Faces of (a,b,c): (b,c), (a,c), (a,b) — these are 1-paths = edges.
# (a,c) is a 1-path = a directed edge from a to c.
# If c→a in the tournament, then (a,c) is NOT in A₁ (which only has edges a→c when a→c).
# So BF₂ = {(a,c) : ∃ (a,b,c) ∈ A₂ with c→a} = {(a,c) : c→a and ∃b with a→b→c}

# Wait but above the loop said BF₂ = ∅? Let me recheck the loop...

# The issue: tuple(list(path[:fi]) + list(path[fi+1:])) for a 2-path (a,b,c):
# fi=0: (b,c)  fi=1: (a,c)  fi=2: (a,b)
# face = (a,c) when fi=1. Is (a,c) in a1_set?
# a1_set contains all DIRECTED edges. (a,c) is in a1_set iff a→c.
# If c→a, then (c,a) ∈ a1_set but (a,c) ∉ a1_set.

# Hmm, let me check what enumerate_allowed_paths returns for p=1
print(f"\n  A₁ for transitive T₅: first few = {[tuple(x) for x in a1[:5]]}")
print(f"  Are they (i,j) with i→j? {all(A[p[0]][p[1]] == 1 for p in a1)}")

# So a1_set = {(i,j) : i→j}. Then (a,c) ∈ a1_set iff a→c. For NT triple with c→a:
# (a,c) ∉ a1_set. So BF₂ should be nonempty for tournaments with NT triples!

# The earlier loop must have been wrong because transitive tournaments have NO NT triples.
# Let me re-run for a non-transitive tournament:
print(f"\n  Re-testing BF₂ for NON-transitive tournament:")
A_test = [[0]*5 for _ in range(5)]
A_test[0][1] = 1; A_test[1][2] = 1; A_test[2][0] = 1  # 3-cycle
A_test[0][3] = 1; A_test[0][4] = 1
A_test[1][3] = 1; A_test[1][4] = 1
A_test[2][3] = 1; A_test[2][4] = 1
A_test[3][4] = 1

a1_t = enumerate_allowed_paths(A_test, 5, 1)
a2_t = enumerate_allowed_paths(A_test, 5, 2)
a1_set_t = set(tuple(x) for x in a1_t)
bf2_t = set()
for path in a2_t:
    for fi in range(3):
        face = tuple(list(path[:fi]) + list(path[fi+1:]))
        if face not in a1_set_t:
            bf2_t.add(face)

print(f"  |A₂| = {len(a2_t)}, |BF₂| = {len(bf2_t)}, BF₂ = {bf2_t}")
om2_t = compute_omega_basis(A_test, 5, 2, a2_t, a1_t)
dim_om2_t = om2_t.shape[1] if om2_t.ndim == 2 else 0
tt_t = sum(1 for p in a2_t if A_test[p[0]][p[2]] == 1)
print(f"  dim(Ω₂) = {dim_om2_t}, TT = {tt_t}, predicted = |A₂|-|BF₂| = {len(a2_t)-len(bf2_t)}")

print("\nDone.")
