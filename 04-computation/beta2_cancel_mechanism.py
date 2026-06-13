#!/usr/bin/env python3
"""
beta2_cancel_mechanism.py - WHY does cancellation always fill the DT deficit?

For a 3-path (a,b,c,d) in A_3:
- Face d_0 = (b,c,d)   -> always in A_2 (b->c, c->d)
- Face d_1 = (a,c,d)   -> bad if c->a ("02 bad face")
- Face d_2 = (a,b,d)   -> always in A_2 (a->b, b->d)... wait, is b->d guaranteed?
  No! b->d is NOT guaranteed. If d->b, then (a,b,d) might still be in A_2
  if a->b and b->d... which requires b->d. Hmm.
  Actually (a,b,d) is a 2-path requiring a->b and b->d. If d->b instead, this is NOT allowed.
  But wait, it IS allowed: the 2-path (a,b,d) needs a->b and b->d. If d->b, then (a,b,d) is NOT an allowed 2-path.

Wait, I need to reconsider. Let me re-derive what faces of a 3-path are in A_2.

A 3-path (a,b,c,d) requires: a->b, b->c, c->d (consecutive arcs).
Its boundary has 4 faces:
  d_0(a,b,c,d) = (b,c,d): needs b->c, c->d. YES (from original).
  d_1(a,b,c,d) = (a,c,d): needs a->c, c->d. a->c iff TT face "02".
  d_2(a,b,c,d) = (a,b,d): needs a->b, b->d. b->d iff TT face "13".
  d_3(a,b,c,d) = (a,b,c): needs a->b, b->c. YES (from original).

So the boundary of (a,b,c,d) is:
  (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)

For this to be in A_2:
- (b,c,d) is always in A_2
- (a,b,c) is always in A_2
- (a,c,d) is in A_2 iff a->c (if c->a, this is NOT an allowed 2-path... wait)
  Actually (a,c,d) as a 2-path means sequence a,c,d with a->c and c->d.
  If c->a instead of a->c, then (a,c,d) is NOT a valid 2-path.
  BUT it might still be an ELEMENT of the path space. The boundary formula
  produces coefficients on ALL sequences, not just valid paths.
  The key is: Omega_2 is spanned by VALID 2-path combinations.

Actually wait, I need to be more careful about the definition.

In GLMY path homology, the boundary maps take linear combinations of p-paths
to linear combinations of (p-1)-paths. The Omega condition ensures that the
boundary stays in the allowed path space.

For (a,b,c,d) to be in Omega_3, we need d(a,b,c,d) to be in Omega_2.
The boundary is (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c).
If (a,c,d) is not an allowed 2-path (because c->a), then the coefficient
of (a,c,d) in the boundary contributes to the "junk" part that must cancel.

So a 3-path is in Omega_3 iff its boundary has no junk component.
DT paths (a->c AND b->d) have NO junk: all 4 faces are valid.
Non-DT paths have junk that must cancel — this is what cancellation pairs do.

KEY INSIGHT: A cancellation pair is the SIMPLEST way to cancel junk.
Two paths sharing the same bad face, taken with opposite signs, cancel the junk.

QUESTION: Is there always a cancellation pair available for every junk direction
that DT misses?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
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


print("=" * 70)
print("CANCELLATION MECHANISM ANALYSIS")
print("=" * 70)

# ============================================================
# ANALYSIS 1: For deficit tournaments at n=6, what is the full
#             Omega_3 structure? Is Omega_3 = DT + cancellation?
# ============================================================
print("\nANALYSIS 1: Is Omega_3 = DT + cancellation?")
print("-" * 50)

n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

om3_eq_dtcancel = 0
om3_bigger = 0
om3_check_count = 0

for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    # Only check deficit score sequences
    if scores not in [(1,2,2,3,3,4), (2,2,2,3,3,3)]:
        continue
    om3_check_count += 1

    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1 = enumerate_allowed_paths(A, n, 1)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # Build DT + cancellation subspace
    dt_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]

    bad02 = defaultdict(list)
    bad13 = defaultdict(list)
    for i, p in enumerate(a3):
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad02[(a_, c_, d_)].append(i)
        if not A[b_][d_]:
            bad13[(a_, b_, d_)].append(i)

    all_vecs = []
    for idx in dt_idx:
        v = np.zeros(len(a3))
        v[idx] = 1
        all_vecs.append(v)
    for group_key, indices in bad02.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                all_vecs.append(v)
    for group_key, indices in bad13.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                all_vecs.append(v)

    if all_vecs:
        V = np.column_stack(all_vecs)
        rk_V = np.linalg.matrix_rank(V, tol=1e-8)
    else:
        rk_V = 0

    if rk_V == d_om3:
        om3_eq_dtcancel += 1
    else:
        om3_bigger += 1

print(f"Deficit tournaments checked: {om3_check_count}")
print(f"Omega_3 = DT + cancel: {om3_eq_dtcancel}")
print(f"Omega_3 > DT + cancel: {om3_bigger}")

if om3_bigger == 0:
    print("  ** Omega_3 IS exactly DT + cancellation for all deficit tournaments! **")


# ============================================================
# ANALYSIS 2: For ALL n=6 tournaments, is Omega_3 = DT + cancel?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: Is Omega_3 = DT + cancellation for ALL n=6?")
print("-" * 50)

om3_eq_all = 0
om3_bigger_all = 0
om3_examples = []

t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s)")

    A = build_adj(n, bits)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0
    if d_om3 == 0:
        om3_eq_all += 1
        continue

    dt_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]

    bad02 = defaultdict(list)
    bad13 = defaultdict(list)
    for i, p in enumerate(a3):
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad02[(a_, c_, d_)].append(i)
        if not A[b_][d_]:
            bad13[(a_, b_, d_)].append(i)

    all_vecs = []
    for idx in dt_idx:
        v = np.zeros(len(a3))
        v[idx] = 1
        all_vecs.append(v)
    for group_key, indices in bad02.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                all_vecs.append(v)
    for group_key, indices in bad13.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                all_vecs.append(v)

    if all_vecs:
        V = np.column_stack(all_vecs)
        rk_V = np.linalg.matrix_rank(V, tol=1e-8)
    else:
        rk_V = 0

    if rk_V == d_om3:
        om3_eq_all += 1
    else:
        om3_bigger_all += 1
        if len(om3_examples) < 3:
            scores = tuple(sorted([sum(row) for row in A]))
            om3_examples.append((bits, scores, d_om3, rk_V))

dt = time.time() - t0
print(f"\nDone in {dt:.0f}s")
print(f"Omega_3 = DT + cancel: {om3_eq_all}/{total} ({100*om3_eq_all/total:.1f}%)")
print(f"Omega_3 > DT + cancel: {om3_bigger_all}/{total} ({100*om3_bigger_all/total:.1f}%)")

if om3_bigger_all > 0:
    print(f"\nExamples where Omega_3 > DT + cancel:")
    for bits, scores, d_om3, rk_V in om3_examples:
        print(f"  bits={bits}, scores={scores}, Om3={d_om3}, DT+cancel={rk_V}")


# ============================================================
# ANALYSIS 3: What ELSE is in Omega_3 besides DT + cancel?
# ============================================================
if om3_bigger_all > 0:
    print(f"\n{'='*70}")
    print("ANALYSIS 3: What else is in Omega_3?")
    print("-" * 50)

    # For the first example, find Omega_3 elements NOT in DT+cancel
    bits = om3_examples[0][0]
    A = build_adj(n, bits)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1]

    dt_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    bad02 = defaultdict(list)
    bad13 = defaultdict(list)
    for i, p in enumerate(a3):
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad02[(a_, c_, d_)].append(i)
        if not A[b_][d_]:
            bad13[(a_, b_, d_)].append(i)

    all_vecs = []
    for idx in dt_idx:
        v = np.zeros(len(a3))
        v[idx] = 1
        all_vecs.append(v)
    for group_key, indices in bad02.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                all_vecs.append(v)
    for group_key, indices in bad13.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                all_vecs.append(v)

    V = np.column_stack(all_vecs)
    rk_V = np.linalg.matrix_rank(V, tol=1e-8)

    print(f"T#{bits}: Om3={d_om3}, DT+cancel rank={rk_V}, gap={d_om3-rk_V}")

    # Find Omega_3 vectors not in span of DT+cancel
    # Project Omega_3 basis onto complement of DT+cancel span
    U, S, Vt = np.linalg.svd(V, full_matrices=True)
    rk = int(np.sum(np.abs(S) > 1e-8))
    complement_basis = Vt[rk:].T if rk < V.shape[1] else np.zeros((V.shape[0], 0))

    # Omega_3 in A_3 coordinates
    for col in range(d_om3):
        om3_vec = om3[:, col]
        # Check if it's in DT+cancel span
        proj, _, _, _ = np.linalg.lstsq(V, om3_vec, rcond=None)
        residual = np.linalg.norm(V @ proj - om3_vec)
        if residual > 1e-6:
            # This Omega_3 element is not in DT+cancel
            nonzero = [(i, om3_vec[i]) for i in range(len(om3_vec)) if abs(om3_vec[i]) > 1e-8]
            print(f"  Extra Omega_3 element (residual={residual:.6f}):")
            for idx, coeff in sorted(nonzero, key=lambda x: abs(x[1]), reverse=True)[:6]:
                p = a3[idx]
                has_02 = A[p[0]][p[2]]
                has_13 = A[p[1]][p[3]]
                print(f"    {coeff:+.4f} * {p} [02={'Y' if has_02 else 'N'}, 13={'Y' if has_13 else 'N'}]")
            break


# ============================================================
# ANALYSIS 4: Cancellation SPAN vs Omega_3 at n=5
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: Omega_3 = DT + cancel at n=5?")
print("-" * 50)

n = 5
total5 = 1 << (n*(n-1)//2)

eq5 = 0
bigger5 = 0

for bits in range(total5):
    A = build_adj(n, bits)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0
    if d_om3 == 0:
        eq5 += 1
        continue

    dt_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    bad02 = defaultdict(list)
    bad13 = defaultdict(list)
    for i, p in enumerate(a3):
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad02[(a_, c_, d_)].append(i)
        if not A[b_][d_]:
            bad13[(a_, b_, d_)].append(i)

    all_vecs = []
    for idx in dt_idx:
        v = np.zeros(len(a3))
        v[idx] = 1
        all_vecs.append(v)
    for group_key, indices in bad02.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                all_vecs.append(v)
    for group_key, indices in bad13.items():
        if len(indices) >= 2:
            for j in range(1, len(indices)):
                v = np.zeros(len(a3))
                v[indices[0]] = 1
                v[indices[j]] = -1
                all_vecs.append(v)

    if all_vecs:
        V = np.column_stack(all_vecs)
        rk_V = np.linalg.matrix_rank(V, tol=1e-8)
    else:
        rk_V = 0

    if rk_V == d_om3:
        eq5 += 1
    else:
        bigger5 += 1

print(f"n=5: Omega_3 = DT + cancel: {eq5}/{total5}")
print(f"n=5: Omega_3 > DT + cancel: {bigger5}/{total5}")


# ============================================================
# ANALYSIS 5: Triple cancellation — three paths sharing bad face
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: Higher-order cancellations at n=6")
print("-" * 50)

# When a bad face group has >= 3 members, we get extra cancellation elements.
# For k paths sharing a bad face, we get k-1 independent cancellation elements.
# But maybe there are elements that combine THREE different bad faces?

# Count bad face group sizes at n=6
n = 6
total = 1 << (n*(n-1)//2)

group_size_dist = {'02': Counter(), '13': Counter()}
for bits in range(min(total, 5000)):
    A = build_adj(n, bits)
    a3 = enumerate_allowed_paths(A, n, 3)

    bad02 = defaultdict(list)
    bad13 = defaultdict(list)
    for i, p in enumerate(a3):
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            bad02[(a_, c_, d_)].append(i)
        if not A[b_][d_]:
            bad13[(a_, b_, d_)].append(i)

    for g in bad02.values():
        group_size_dist['02'][len(g)] += 1
    for g in bad13.values():
        group_size_dist['13'][len(g)] += 1

print(f"Bad face group sizes (first 5000 tournaments):")
print(f"  Type 02: {dict(sorted(group_size_dist['02'].items()))}")
print(f"  Type 13: {dict(sorted(group_size_dist['13'].items()))}")


print("\nDone.")
