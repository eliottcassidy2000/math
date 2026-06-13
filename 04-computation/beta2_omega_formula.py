#!/usr/bin/env python3
"""
beta2_omega_formula.py - Derive and verify explicit formula for dim(Omega_2)

KEY INSIGHT: For a tournament T on n vertices,
  dim(Omega_2) = |A_2| - J_2
where:
  |A_2| = # allowed 2-paths (a,b,c) with a->b->c, a!=c
  J_2 = |{(a,c) : c->a and exists b with a->b->c}|
      = # reverse-ordered pairs with at least one 2-step path

PROOF: Omega_2 consists of vectors in span(A_2) whose boundary lies in Omega_1 = R_1.
For (a,b,c) with a->b->c:
  d_2(a,b,c) = (b,c) - (a,c) + (a,b)
  Faces (a,b) and (b,c) are always arcs (allowed).
  Face (a,c) is allowed iff a->c (TT triple).
  If c->a (non-TT), then (a,c) is a "junk face".

For each junk pair (a,c) with c->a: let k_{a,c} = #{b : a->b, b->c, b!=a,c}.
- k_{a,c} non-TT paths all share junk face -(a,c)
- The constraint "junk cancels" is: sum of coefficients = 0
- This costs 1 dimension: k paths but only k-1 independent in Omega_2
- If k=0: no paths, no issue
- If k>=1: lose exactly 1 dimension

So dim(Omega_2) = |A_2| - J_2 where J_2 = #{(a,c) : c->a, k_{a,c}>=1}.

Then we track:
- delta(|A_2|) = 2*(d_u - d_v - 1)  [PROVED: HYP-228]
- delta(J_2) = ?  [TO FIND]
- delta(dim Om_2) = delta(|A_2|) - delta(J_2)

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def flip_arc_bits(bits, i, j, n):
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return bits ^ (1 << idx)
            idx += 1
    return bits

def compute_J2(A, n):
    """Compute J_2 = #{(a,c): c->a and exists b with a->b->c, b!=a,c}."""
    J2 = 0
    for a in range(n):
        for c in range(n):
            if a == c:
                continue
            if A[c][a] == 0:
                continue  # need c->a
            # Count intermediates b with a->b->c
            for b in range(n):
                if b == a or b == c:
                    continue
                if A[a][b] and A[b][c]:
                    J2 += 1
                    break  # only need existence
    return J2

def compute_k_values(A, n):
    """Compute k_{a,c} for all reverse pairs (a,c) with c->a."""
    k_vals = {}
    for a in range(n):
        for c in range(n):
            if a == c:
                continue
            if A[c][a]:  # c->a
                k = sum(1 for b in range(n) if b != a and b != c and A[a][b] and A[b][c])
                k_vals[(a,c)] = k
    return k_vals


print("=" * 70)
print("VERIFYING dim(Omega_2) = |A_2| - J_2 FORMULA")
print("=" * 70)

for n in [4, 5, 6]:
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs
    errors = 0

    if n <= 5:
        rng = range(total)
    else:
        import random
        random.seed(42)
        rng = [random.randint(0, total-1) for _ in range(5000)]

    for bits in rng:
        A = build_adj(n, bits)

        # Compute dim(Omega_2) via path_homology_v2
        paths2 = enumerate_allowed_paths(A, n, 2)
        paths1 = enumerate_allowed_paths(A, n, 1)
        if len(paths2) > 0 and len(paths1) > 0:
            om2 = compute_omega_basis(A, n, 2, paths2, paths1)
            dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
        else:
            dim_om2 = 0

        # Compute via formula
        A2 = len(paths2)
        J2 = compute_J2(A, n)
        formula = A2 - J2

        if dim_om2 != formula:
            errors += 1
            if errors <= 3:
                print(f"  ERROR at n={n}: bits={bits}, dim(Om2)={dim_om2}, |A2|={A2}, J2={J2}, formula={formula}")

    print(f"  n={n}: {errors} errors in {len(list(rng)) if n <= 5 else 5000} tests")

# ANALYSIS: J_2 structure
print(f"\n{'='*70}")
print("J_2 ANALYSIS")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

j2_dist = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    j2 = compute_J2(A, n)
    j2_dist[j2] += 1

print(f"n={n}: J_2 distribution: {dict(sorted(j2_dist.items()))}")

# delta(J_2) under arc flip
print(f"\n{'='*70}")
print("delta(J_2) UNDER ARC FLIP")
print("=" * 70)

dJ2_by_local = defaultdict(set)
dJ2_dist = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    j2_T = compute_J2(A, n)
    scores = [sum(row) for row in A]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc_bits(bits, i, j, n)
            A2 = build_adj(n, bits2)
            j2_T2 = compute_J2(A2, n)

            dJ2 = j2_T2 - j2_T

            if A[i][j] == 1:
                u, v = i, j
            else:
                u, v = j, i

            du, dv = scores[u], scores[v]
            puv = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[w][v])
            pvu = sum(1 for w in range(n) if w != u and w != v and A[v][w] and A[w][u])
            cout = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[v][w])
            cin = sum(1 for w in range(n) if w != u and w != v and A[w][u] and A[w][v])

            local_key = (du, dv, cout, cin, puv, pvu)
            dJ2_by_local[local_key].add(dJ2)
            dJ2_dist[dJ2] += 1

print(f"delta(J_2) distribution: {dict(sorted(dJ2_dist.items()))}")

det = sum(1 for v in dJ2_by_local.values() if len(v) == 1)
ndet = sum(1 for v in dJ2_by_local.values() if len(v) > 1)
print(f"delta(J_2) determined by (du,dv,cout,cin,puv,pvu): {det}/{det+ndet}")

if ndet > 0:
    print("Non-determined cases:")
    for key in sorted(dJ2_by_local.keys()):
        if len(dJ2_by_local[key]) > 1:
            print(f"  {key}: {sorted(dJ2_by_local[key])}")


# Track the PRECISE changes to J_2
print(f"\n{'='*70}")
print("PRECISE delta(J_2) DECOMPOSITION")
print("=" * 70)
print("Track: gain/loss of junk pairs (u,v), (v,u), and all others")

decomp_dist = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    kvals_T = compute_k_values(A, n)

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc_bits(bits, i, j, n)
            A2 = build_adj(n, bits2)
            kvals_T2 = compute_k_values(A2, n)

            if A[i][j] == 1:
                u, v = i, j
            else:
                u, v = j, i

            # Change from pair (u,v): NOT junk in T (u->v), JUNK in T' (v->u)
            # k_{u,v} in T': same as p_uv (2-step paths from u to v, unchanged by flip)
            puv = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[w][v])
            uv_change = (1 if puv >= 1 else 0)  # +1 if (u,v) enters J_2

            # Change from pair (v,u): JUNK in T (u->v means for pair (v,u): c=u, a=v, c->a=u->v=TRUE)
            # NOT junk in T' (u->v false in T')
            pvu = sum(1 for w in range(n) if w != u and w != v and A[v][w] and A[w][u])
            vu_change = -(1 if pvu >= 1 else 0)  # -1 if (v,u) leaves J_2

            # Changes in other pairs
            other_change = 0
            for a in range(n):
                for c in range(n):
                    if a == c or (a == u and c == v) or (a == v and c == u):
                        continue
                    was_junk = (kvals_T.get((a,c), 0) >= 1) if A[c][a] else False
                    is_junk = (kvals_T2.get((a,c), 0) >= 1) if A2[c][a] else False
                    if is_junk and not was_junk:
                        other_change += 1
                    elif was_junk and not is_junk:
                        other_change -= 1

            total_dJ2 = uv_change + vu_change + other_change
            decomp_dist[(uv_change, vu_change, other_change)] += 1

print("(uv_gain, vu_loss, other) distribution:")
for key in sorted(decomp_dist.keys()):
    print(f"  {key}: count={decomp_dist[key]}, total_dJ2={sum(key)}")


# Check: which "other" pairs change their J_2 status?
print(f"\n{'='*70}")
print("WHICH 'OTHER' PAIRS CHANGE J_2 STATUS?")
print("=" * 70)

change_types = Counter()
for bits in range(min(total, 200)):
    A = build_adj(n, bits)
    kvals_T = compute_k_values(A, n)

    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                u, v = i, j
            else:
                u, v = j, i

            bits2 = flip_arc_bits(bits, i, j, n)
            A2 = build_adj(n, bits2)
            kvals_T2 = compute_k_values(A2, n)

            for a in range(n):
                for c in range(n):
                    if a == c or (a == u and c == v) or (a == v and c == u):
                        continue
                    if A[c][a] != A2[c][a]:
                        continue  # junk status (c->a) didn't change for non-{u,v} pairs

                    if not A[c][a]:
                        continue  # not a junk pair at all

                    k_old = kvals_T.get((a,c), 0)
                    k_new = kvals_T2.get((a,c), 0)
                    if (k_old >= 1) != (k_new >= 1):
                        # This pair's J_2 membership changed
                        # Figure out relationship to u,v
                        if a == u:
                            rel = "a=u"
                        elif a == v:
                            rel = "a=v"
                        elif c == u:
                            rel = "c=u"
                        elif c == v:
                            rel = "c=v"
                        else:
                            rel = "disjoint"
                        direction = "entered" if k_new >= 1 else "left"
                        change_types[(rel, direction, k_old, k_new)] += 1

print("Change types for 'other' pairs:")
for key in sorted(change_types.keys()):
    rel, direction, k_old, k_new = key
    print(f"  {rel}, {direction}: k {k_old}->{k_new}, count={change_types[key]}")


# EXTENDED FORMULA: dim(Omega_3) = |A_3| - J_3,0 - J_3,1 + J_3,01 ?
print(f"\n{'='*70}")
print("FORMULA FOR dim(Omega_3)")
print("=" * 70)

# For 3-path (a,b,c,d) with a->b->c->d:
# d_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
# Face (b,c,d): always allowed (b->c->d is subpath)
# Face (a,b,c): always allowed (a->b->c is subpath)
# Face (a,c,d): JUNK iff c->a (not a->c). Junk type "02" (removing index 1, long face spans 0 to 2)
# Face (a,b,d): JUNK iff d->b (not b->d). Junk type "13" (removing index 2, long face spans 1 to 3)
#
# A 3-path can have 0, 1, or 2 junk faces.

# Count junk face types
junk_stats = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    paths3 = enumerate_allowed_paths(A, n, 3)
    for (a,b,c,d) in paths3:
        j02 = (A[c][a] == 1)  # (a,c,d) junk iff c->a
        j13 = (A[d][b] == 1)  # (a,b,d) junk iff d->b
        junk_stats[(j02, j13)] += 1

print("3-path junk face types (all n=5 tournaments):")
for key in sorted(junk_stats.keys()):
    j02, j13 = key
    label = f"junk02={'Y' if j02 else 'N'}, junk13={'Y' if j13 else 'N'}"
    print(f"  {label}: {junk_stats[key]}")


# Verify: does dim(Omega_3) follow a similar formula?
# Hypothesis: dim(Om_3) = |A_3| - J_3 where J_3 counts some constraint
# The constraint structure is more complex because of 2 possible junk faces per path

print("\nVerifying dim(Omega_3) formulas...")

def compute_J3_naive(A, n):
    """Try: J_3 = number of 'junk constraints' at degree 3."""
    paths3 = enumerate_allowed_paths(A, n, 3)

    # Group paths by their junk face types
    # For each "junk face" (i.e., non-allowed interior face),
    # we get one constraint per distinct junk face value

    # Junk type 02: face (a,c,d) with c->a
    # Constraint: for fixed (a,c,d), sum of coefficients of paths with this junk face = 0
    # But (a,c,d) also has its own junk faces...

    # Actually, the constraint is that d_3 of the combination lands in Omega_2 (not just R_2).
    # This is recursive: Omega_2 itself has constraints.

    # Simpler approach: count constraints directly.
    # For each "junk 2-face" (a non-allowed 2-path appearing as a face of some 3-path):
    # We need a constraint that the coefficient of this junk 2-face in d_3(combination) = 0.
    # BUT: this junk 2-face might also appear in Omega_2 (via NT cancellation).
    # If the junk 2-face IS in Omega_2, then it doesn't need to cancel!
    #
    # So the constraint is: d_3(combination) must be in Omega_2.
    # d_3(combination) = sum of allowed 2-paths + sum of potentially-junk 2-paths.
    # The allowed 2-paths are fine (they're in R_2, and they might or might not be in Omega_2).
    # The junk 2-paths need to be: either zero in the sum, or combinable to land in Omega_2.
    #
    # This is more complex than degree 2.
    # At degree 2: junk faces are 1-paths (a,c) with c->a. These are never in Omega_1 = R_1
    # (since R_1 only contains arcs a->c). So they MUST cancel.
    # At degree 3: junk faces are 2-paths (a,c,d) or (a,b,d) that are not in R_2.
    # These are also not in Omega_2 (since Omega_2 is a subset of R_2).
    # So they also MUST cancel.

    # Therefore: the same counting argument works.
    # For each non-R_2 face appearing in boundaries of 3-paths,
    # group the 3-paths by their junk face value, and lose 1 dimension per group.

    # But a 3-path can have TWO junk faces. So it appears in two groups.
    # This creates a coupled constraint system.

    # Let me just count the number of distinct junk 2-faces.
    junk_faces_02 = set()  # (a,c,d) with c->a
    junk_faces_13 = set()  # (a,b,d) with d->b

    for (a,b,c,d) in paths3:
        if A[c][a]:  # c->a means (a,c,d) is not in R_2
            junk_faces_02.add((a,c,d))
        if A[d][b]:  # d->b means (a,b,d) is not in R_2
            junk_faces_13.add((a,b,d))

    return len(junk_faces_02), len(junk_faces_13), len(junk_faces_02 & junk_faces_13)


for bits in range(min(total, 100)):
    A = build_adj(n, bits)

    paths3 = enumerate_allowed_paths(A, n, 3)
    paths2 = enumerate_allowed_paths(A, n, 2)
    if len(paths3) > 0 and len(paths2) > 0:
        om3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    else:
        dim_om3 = 0

    j02, j13, j_both = compute_J3_naive(A, n)
    j_total = j02 + j13 - j_both  # inclusion-exclusion?

    A3 = len(paths3)

    # Try various formulas
    f1 = A3 - j02 - j13  # subtract both
    f2 = A3 - j02 - j13 + j_both  # inclusion-exclusion
    f3 = A3 - (j02 + j13)  # just total

    if bits < 10:
        print(f"  bits={bits}: dim(Om3)={dim_om3}, |A3|={A3}, j02={j02}, j13={j13}, "
              f"j_both={j_both}, f1={f1}, f2={f2}")


# More systematic: compute dim(Om3) and test formula
print("\nSystematic test:")
errors_f1 = 0
errors_f2 = 0
for bits in range(total):
    A = build_adj(n, bits)

    paths3 = enumerate_allowed_paths(A, n, 3)
    paths2 = enumerate_allowed_paths(A, n, 2)
    if len(paths3) > 0 and len(paths2) > 0:
        om3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    else:
        dim_om3 = 0

    j02, j13, j_both = compute_J3_naive(A, n)

    if dim_om3 != len(paths3) - j02 - j13:
        errors_f1 += 1
    if dim_om3 != len(paths3) - j02 - j13 + j_both:
        errors_f2 += 1

print(f"  n=5: dim(Om3) = |A3| - j02 - j13: {errors_f1} errors")
print(f"  n=5: dim(Om3) = |A3| - j02 - j13 + j_both: {errors_f2} errors")

# If neither works, find the right formula
if errors_f1 > 0 and errors_f2 > 0:
    # Find a correction pattern
    corrections = Counter()
    for bits in range(total):
        A = build_adj(n, bits)
        paths3 = enumerate_allowed_paths(A, n, 3)
        paths2 = enumerate_allowed_paths(A, n, 2)
        if len(paths3) > 0 and len(paths2) > 0:
            om3 = compute_omega_basis(A, n, 3, paths3, paths2)
            dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
        else:
            dim_om3 = 0

        j02, j13, j_both = compute_J3_naive(A, n)
        correction = dim_om3 - (len(paths3) - j02 - j13 + j_both)
        corrections[correction] += 1

    print(f"  Correction distribution (from f2): {dict(sorted(corrections.items()))}")


print("\nDone.")
