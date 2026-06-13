#!/usr/bin/env python3
"""
beta2_deltaJ2_formula.py - Find an exact formula for delta(J_2) under arc flip

We know:
  dim(Omega_2) = |A_2| - J_2
  delta(|A_2|) = 2*(d_u - d_v - 1)  [PROVED]

If we can find delta(J_2) as a function of local data,
then delta(dim Omega_2) = delta(|A_2|) - delta(J_2) is fully determined.

Combined with delta(beta_1) and delta(rk d_3), this would prove beta_2 = 0
via arc-flip invariance from the transitive tournament.

J_2 = #{(a,c): c->a AND A^2[a,c] > 0}
     = #{reverse pairs (a,c) where there exists b with a->b->c}

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os
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


def compute_J2(A, n):
    """J_2 = #{(a,c): c->a AND exists b with a->b->c}"""
    J2 = 0
    pairs = []
    for a in range(n):
        for c in range(n):
            if a == c:
                continue
            if not A[c][a]:  # need c->a
                continue
            # Check if exists b with a->b->c
            has_path = False
            for b in range(n):
                if b != a and b != c and A[a][b] and A[b][c]:
                    has_path = True
                    break
            if has_path:
                J2 += 1
                pairs.append((a, c))
    return J2, pairs


def compute_A2_sum(A, n):
    """|A_2| = sum of off-diagonal entries of A^2"""
    total = 0
    for a in range(n):
        for c in range(n):
            if a == c:
                continue
            for b in range(n):
                if A[a][b] and A[b][c]:
                    total += 1
    return total


def compute_beta1(A, n):
    """beta_1 via full chain complex computation"""
    paths = {}
    omega = {}
    for p in range(3):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dims = {p: (omega[p].shape[1] if omega[p].ndim == 2 else 0) for p in range(3)}

    # rk(d_1)
    if dims[1] > 0:
        bd1 = build_full_boundary_matrix(paths[1], paths[0])
        rk1 = np.linalg.matrix_rank(bd1 @ omega[1], tol=1e-8)
    else:
        rk1 = 0

    # rk(d_2)
    if dims[2] > 0 and len(paths[1]) > 0 and dims[1] > 0:
        bd2 = build_full_boundary_matrix(paths[2], paths[1])
        raw2 = bd2 @ omega[2]
        coords2, _, _, _ = np.linalg.lstsq(omega[1], raw2, rcond=None)
        rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    else:
        rk2 = 0

    z1 = dims[1] - rk1
    b1 = z1 - rk2
    return b1


def flip_arc(A, n, u, v):
    """Flip arc u->v to v->u. Returns new adjacency matrix."""
    B = [row[:] for row in A]
    B[u][v] = 0
    B[v][u] = 1
    return B


def get_scores(A, n):
    return [sum(A[i]) for i in range(n)]


# ============================================================
# ANALYSIS 1: delta(J_2) under arc flips at n=5
# ============================================================
print("=" * 70)
print("ANALYSIS 1: delta(J_2) under arc flips at n=5")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

# Collect data: for each arc flip, record (d_u, d_v, delta_J2, delta_A2, delta_dimOm2)
flip_data = []

for bits in range(total):
    A = build_adj(n, bits)
    scores = get_scores(A, n)
    J2_before, _ = compute_J2(A, n)
    A2_before = compute_A2_sum(A, n)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            B = flip_arc(A, n, u, v)
            J2_after, _ = compute_J2(B, n)
            A2_after = compute_A2_sum(B, n)

            d_u = scores[u]  # out-degree of u in original
            d_v = scores[v]  # out-degree of v in original
            delta_J2 = J2_after - J2_before
            delta_A2 = A2_after - A2_before

            # Count local quantities that might determine delta_J2
            # Number of common out-neighbors of u and v
            common_out = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[v][w])
            # Number of common in-neighbors
            common_in = sum(1 for w in range(n) if w != u and w != v and A[w][u] and A[w][v])
            # Out-neighbors of u that are in-neighbors of v (u->w->v pattern through w)
            # But this is just: does u reach v in 2 steps? That's A^2[u,v]
            a2_uv = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[w][v])
            a2_vu = sum(1 for w in range(n) if w != u and w != v and A[v][w] and A[w][u])

            # The arc u->v: in original, A[u][v]=1, A[v][u]=0
            # After flip: A[v][u]=1, A[u][v]=0

            flip_data.append({
                'd_u': d_u, 'd_v': d_v,
                'delta_J2': delta_J2, 'delta_A2': delta_A2,
                'common_out': common_out, 'common_in': common_in,
                'a2_uv': a2_uv, 'a2_vu': a2_vu,
                'du_minus_dv': d_u - d_v
            })

# Check if delta_J2 is determined by local data
print(f"Total arc flips: {len(flip_data)}")

# Group by (d_u, d_v, common_out, common_in, a2_uv, a2_vu)
by_local = defaultdict(set)
for d in flip_data:
    key = (d['d_u'], d['d_v'], d['common_out'], d['common_in'], d['a2_uv'], d['a2_vu'])
    by_local[key].add(d['delta_J2'])

non_determined = 0
for key, vals in sorted(by_local.items()):
    if len(vals) > 1:
        non_determined += 1
        print(f"  NOT determined: (d_u={key[0]}, d_v={key[1]}, co={key[2]}, ci={key[3]}, "
              f"a2_uv={key[4]}, a2_vu={key[5]}): delta_J2 in {vals}")
    else:
        val = list(vals)[0]

print(f"\nNon-determined by (d_u, d_v, co, ci, a2_uv, a2_vu): {non_determined}/{len(by_local)}")

# Try simpler: is delta_J2 determined by (d_u, d_v, a2_uv, a2_vu)?
by_simple = defaultdict(Counter)
for d in flip_data:
    key = (d['d_u'], d['d_v'], d['a2_uv'], d['a2_vu'])
    by_simple[key][d['delta_J2']] += 1

print(f"\nBy (d_u, d_v, a2_uv, a2_vu):")
for key in sorted(by_simple.keys()):
    vals = by_simple[key]
    if len(vals) > 1:
        print(f"  (d_u={key[0]}, d_v={key[1]}, a2_uv={key[2]}, a2_vu={key[3]}): {dict(vals)}")

# Check a candidate formula: delta_J2 = f(d_u, d_v, a2_uv, a2_vu, common_out, common_in)
# Since J_2 counts pairs (a,c) with c->a and A^2[a,c]>0,
# flipping u->v to v->u:
#  - Changes which pairs have c->a: the pair (v,u) goes from u->v to v->u
#    So (v,u) was NOT a reverse pair (since u->v), now IS a reverse pair
#    And (u,v) was... wait, (a,c) with c->a means a is reached FROM c.
#    Original: u->v. So (a,c) = (v,u): c=u, a=v, c->a means u->v. Yes.
#    After flip: v->u. So (v,u) has c=u, a=v, c->a means u->v? No, v->u.
#    Wait: (a,c)=(v,u) means a=v, c=u, c->a means u->v.
#    Original: u->v, so (v,u) IS a reverse pair.
#    After flip: v->u, so u is NOT going to v. (v,u): c=u, a=v, c->a means u->v = 0. Not reverse.
#    And (u,v): a=u, c=v, c->a means v->u = 1 after flip. So (u,v) becomes reverse pair.
#
#    So the pair (v,u) was reverse, (u,v) becomes reverse after flip. These swap.
#    But A^2 also changes!

# Let me think about this more carefully.
# In the original, the reverse pairs are {(a,c): A[c][a]=1} = arcs of A^T.
# J_2 pairs are reverse pairs (a,c) with A^2[a,c] > 0.
#
# When we flip u->v to v->u:
# 1. Arc changes: remove u->v, add v->u
# 2. A^T changes: remove (v,u) from A^T arcs, add (u,v)
# 3. A^2[a,c] = sum_b A[a][b]*A[b][c] changes for:
#    - (a,c) where a=u or b=u or c=u (since A[u][v] changes)
#    - (a,c) where a=v or b=v or c=v (since A[v][u] changes)
# More precisely: A^2[a,c] = sum_b A[a][b]*A[b][c]
# delta(A^2[a,c]) counts paths through the changed arcs

print(f"\n{'='*70}")
print("ANALYSIS 2: Decomposing delta(J_2) into pair-level changes")
print("=" * 70)

# For each flip, track which (a,c) pairs enter/leave J_2
enter_types = Counter()  # what kind of pair enters J_2
leave_types = Counter()  # what kind of pair leaves J_2

for bits in range(min(total, 200)):  # sample
    A = build_adj(n, bits)
    J2_before, pairs_before = compute_J2(A, n)
    set_before = set(pairs_before)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            B = flip_arc(A, n, u, v)
            J2_after, pairs_after = compute_J2(B, n)
            set_after = set(pairs_after)

            entered = set_after - set_before
            left = set_before - set_after

            for (a, c) in entered:
                # Classify: is this the (u,v) pair? Adjacent to u? Adjacent to v?
                if (a, c) == (u, v):
                    enter_types['direct_uv'] += 1
                elif (a, c) == (v, u):
                    enter_types['direct_vu'] += 1
                elif a == u or c == u:
                    enter_types['adj_u'] += 1
                elif a == v or c == v:
                    enter_types['adj_v'] += 1
                else:
                    enter_types['non_adj'] += 1

            for (a, c) in left:
                if (a, c) == (u, v):
                    leave_types['direct_uv'] += 1
                elif (a, c) == (v, u):
                    leave_types['direct_vu'] += 1
                elif a == u or c == u:
                    leave_types['adj_u'] += 1
                elif a == v or c == v:
                    leave_types['adj_v'] += 1
                else:
                    leave_types['non_adj'] += 1

print(f"Pairs ENTERING J_2: {dict(enter_types)}")
print(f"Pairs LEAVING J_2: {dict(leave_types)}")
print(f"Non-adjacent changes: enter={enter_types.get('non_adj',0)}, leave={leave_types.get('non_adj',0)}")

if enter_types.get('non_adj', 0) == 0 and leave_types.get('non_adj', 0) == 0:
    print("  ** J_2 changes are LOCAL (only adjacent pairs change) **")


# ============================================================
# ANALYSIS 3: Exact formula for delta(J_2)
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: Exact formula for delta(J_2)")
print("=" * 70)

# Let's derive the formula algebraically.
# Flip u->v to v->u. Let's count which (a,c) pairs change J_2 status.
#
# A pair (a,c) is in J_2 iff:
#   (i) c->a (i.e. A[c][a] = 1)
#   (ii) exists b with a->b->c (i.e. A^2[a,c] > 0)
#
# After flip:
#   (i') changes only for (a,c) involving {u,v}:
#     - (v,u): was c->a = A[u][v]=1 (reverse), now A[u][v]=0 (not reverse)
#     - (u,v): was c->a = A[v][u]=0 (not reverse), now A[v][u]=1 (reverse)
#   (ii') A^2 changes: paths a->b->c that use arc u->v or v->u
#
# So delta(J_2) = (new J_2 pairs) - (lost J_2 pairs)
#
# Categories of change:
# Category A: (a,c) = (u,v) — was NOT reverse, becomes reverse
#   Enters J_2 iff A'^2[u,v] > 0, i.e., exists w with u->w->v in B
#   But we flipped u->v to v->u, so u->v is gone.
#   u->w->v in B: need A'[u][w]=1 and A'[w][v]=1.
#   A'[u][w] = A[u][w] for w != v. A'[w][v] = A[w][v] for w != u.
#   So: exists w != u,v with A[u][w]=1 and A[w][v]=1? That's A^2[u,v] counting only w != u,v.
#   Wait, A^2[u,v] = sum_w A[u][w]*A[w][v]. Since A[u][u]=0 and A[v][v]=0,
#   this is sum_{w != u,v} A[u][w]*A[w][v] = #{common successors of u that are predecessors of v... no}
#   = #{w: u->w AND w->v} = a2_uv (computed above).
#   After flip: u doesn't go to v, and v doesn't go to u in terms of direct arcs.
#   A'^2[u,v] = sum_w A'[u][w]*A'[w][v] = sum_{w != u,v} A[u][w]*A[w][v] = a2_uv
#   (The w=v term: A'[u][v]=0, so it contributes 0. The w=u term: A'[u][u]=0.)
#   Actually wait: in the sum, w ranges over all vertices. A'[u][w]*A'[w][v]:
#   - w=u: A'[u][u]=0 -> 0
#   - w=v: A'[u][v]=0 (flipped away) -> 0
#   - w != u,v: A'[u][w]=A[u][w], A'[w][v]=A[w][v] -> same as original
#   So A'^2[u,v] = a2_uv (same paths, since the direct arc u->v doesn't appear in 2-step paths)
#   Wait but what about A'[v][u]=1? That doesn't affect A'^2[u,v] since we need u->w->v.
#
#   So (u,v) enters J_2 iff a2_uv > 0.
#
# Category B: (a,c) = (v,u) — was reverse (A[u][v]=1), becomes NOT reverse
#   Leaves J_2 if it was in J_2, i.e., if A^2[v,u] > 0.
#   A^2[v,u] = #{w: v->w->u} = a2_vu.
#   So (v,u) leaves J_2 iff a2_vu > 0.
#
# Category C: Other pairs (a,c) where the REVERSE STATUS doesn't change,
#   but A^2[a,c] changes (so the "path exists" condition changes).
#   These are pairs where c->a status is unchanged but A^2[a,c] goes from 0 to >0 or vice versa.

# Let me verify this categorization
print("Verifying pair-level changes...")

cat_A_enters = 0  # (u,v) enters J_2
cat_A_stays_out = 0
cat_B_leaves = 0  # (v,u) leaves J_2
cat_B_stays_in = 0
cat_C_enters = 0
cat_C_leaves = 0
cat_C_net = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    J2_before, pairs_before = compute_J2(A, n)
    set_before = set(pairs_before)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            B = flip_arc(A, n, u, v)
            J2_after, pairs_after = compute_J2(B, n)
            set_after = set(pairs_after)

            entered = set_after - set_before
            left = set_before - set_after

            # Category A: (u,v) enters?
            if (u, v) in entered:
                cat_A_enters += 1
            elif (u, v) not in set_before and (u, v) not in set_after:
                cat_A_stays_out += 1

            # Category B: (v,u) leaves?
            if (v, u) in left:
                cat_B_leaves += 1
            elif (v, u) in set_before and (v, u) in set_after:
                cat_B_stays_in += 1

            # Category C: everything else
            c_entered = entered - {(u, v)}
            c_left = left - {(v, u)}
            cat_C_enters += len(c_entered)
            cat_C_leaves += len(c_left)
            cat_C_net[len(c_entered) - len(c_left)] += 1

print(f"\nCategory A (u,v) enters J_2: {cat_A_enters}")
print(f"Category A (u,v) stays out: {cat_A_stays_out}")
print(f"Category B (v,u) leaves J_2: {cat_B_leaves}")
print(f"Category B (v,u) stays in: {cat_B_stays_in}")
print(f"Category C enters: {cat_C_enters}")
print(f"Category C leaves: {cat_C_leaves}")
print(f"Category C net change distribution: {dict(sorted(cat_C_net.items()))}")


# ============================================================
# ANALYSIS 4: Category C decomposition
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: Category C — what changes A^2 threshold for non-{u,v} pairs?")
print("=" * 70)

# For a pair (a,c) not equal to (u,v) or (v,u):
# The reverse status c->a doesn't change (since only arc u->v changes).
# But A^2[a,c] changes. When does it cross the 0 threshold?
#
# A^2[a,c] = #{w: a->w->c}
# After flip: A'^2[a,c] = #{w: a->w->c in B}
#
# The only paths affected are those through u or v:
# - Path a->u->c: exists iff A[a][u]=1 and A[u][c]=1 in A.
#   In B: A'[a][u]=A[a][u], A'[u][c]=A[u][c] UNLESS c=v (A'[u][v]=0).
#   So this path is LOST iff c=v and A[a][u]=1 and A[u][v]=1.
#   But A[u][v]=1 always (we're flipping u->v). So lost iff c=v and A[a][u]=1.
#   Wait but (a,c) != (u,v) and (a,c) != (v,u). c=v is allowed if a != u.
#
# - Path a->v->c: exists iff A[a][v]=1 and A[v][c]=1 in A.
#   In B: A'[a][v]=A[a][v], A'[v][c]=A[v][c] UNLESS c=u (A'[v][u]=1 vs A[v][u]=0).
#   So this path is GAINED iff c=u and A[a][v]=1 and A[v][u]=0.
#   A[v][u]=0 always (since A[u][v]=1). So gained iff c=u and A[a][v]=1.
#   But (a,c)=(a,u): c=u, a != v (since != (v,u)).
#
# Also check: when a=v, path v->u->c uses arc u->c.
#   In A: A[v][u]=0 -> no path. In B: A'[v][u]=1. So path v->u->c gained if A[u][c]=1.
#   This means: for pair (v,c) with c != u: gained iff A[u][c]=1.
#
# When a=u, path u->v->c uses arc u->v.
#   In A: A[u][v]=1. In B: A'[u][v]=0. So path u->v->c lost if A[v][c]=1.
#   For pair (u,c) with c != v: lost iff A[v][c]=1.

# Summary of A^2 changes for non-{u,v} pairs:
# delta(A^2[a,c]):
#   If a=u, c != v: lose path u->v->c if A[v][c]=1 (-1)
#   If a != v, c=v: lose path a->u->v if A[a][u]=1 (-1)
#   If a=v, c != u: gain path v->u->c if A[u][c]=1 (+1)
#   If a != u, c=u: gain path a->v->u if A[a][v]=1 (+1)
#   Other (a,c) with a,c not in {u,v}: no change to A^2

# So Category C pairs MUST have a or c in {u,v}.
# This confirms locality!

# Let me verify by checking non-adjacent Category C changes
cat_C_nonadj = 0
for bits in range(min(total, 100)):
    A = build_adj(n, bits)
    J2_before, pairs_before = compute_J2(A, n)
    set_before = set(pairs_before)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            B = flip_arc(A, n, u, v)
            _, pairs_after = compute_J2(B, n)
            set_after = set(pairs_after)
            entered = set_after - set_before - {(u,v)}
            left = set_before - set_after - {(v,u)}
            for (a,c) in entered | left:
                if a not in (u,v) and c not in (u,v):
                    cat_C_nonadj += 1
                    print(f"  NON-ADJ C change! bits={bits}, flip ({u},{v}), pair ({a},{c})")

print(f"Non-adjacent Category C changes: {cat_C_nonadj}")
if cat_C_nonadj == 0:
    print("  ** CONFIRMED: ALL J_2 changes are local to {u,v} **")


# ============================================================
# ANALYSIS 5: Exact formula derivation
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: Exact delta(J_2) formula")
print("=" * 70)

# From the analysis above, delta(J_2) depends on:
# 1. Whether (u,v) enters J_2 (i.e., a2_uv > 0)
# 2. Whether (v,u) leaves J_2 (i.e., a2_vu > 0 AND (v,u) was in J_2)
# 3. Category C changes for pairs adjacent to {u,v}
#
# For Category C, we need to count pairs that CROSS the A^2 threshold.
# A pair enters J_2 when: it was already reverse, A^2 was 0, now A^2 > 0
# A pair leaves J_2 when: it's still reverse, A^2 was > 0, now A^2 = 0
# The A^2 ONLY changes by +1 or -1 for the affected paths.
# So threshold crossing only happens when A^2 goes from 0 to 1 (enter) or 1 to 0 (leave).

# Let me count these precisely.
# For each arc flip, compute the EXACT delta_J2 and compare to a formula.

formula_errors = 0
for bits in range(total):
    A = build_adj(n, bits)
    scores = get_scores(A, n)
    J2_before, pairs_before = compute_J2(A, n)
    set_before = set(pairs_before)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            B = flip_arc(A, n, u, v)
            J2_after, pairs_after = compute_J2(B, n)
            delta_J2 = J2_after - J2_before

            d_u = scores[u]
            d_v = scores[v]

            # Term 1: (u,v) enters?
            a2_uv = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[w][v])
            term1 = 1 if a2_uv > 0 else 0

            # Term 2: (v,u) leaves?
            # (v,u) was in J_2 iff A[u][v]=1 (yes, it's reverse) AND A^2[v,u] > 0
            a2_vu = sum(1 for w in range(n) if w != u and w != v and A[v][w] and A[w][u])
            term2 = -1 if a2_vu > 0 else 0

            # Term 3: Category C pairs
            # Pairs (u, c) with c != v, c != u: reverse means A[c][u]=1, i.e., c is in-neighbor of u
            # delta(A^2[u,c]) = -1 if A[v][c]=1 (lost path u->v->c)
            # But also: if c is an out-neighbor of u that goes through new arc v->u?
            # No: A'^2[u,c] also gains path u->... hmm, let me reconsider.
            #
            # A'^2[u,c] - A^2[u,c]:
            #   Lost: paths through old arc u->v, i.e., u->v->c where A[v][c]=1 (-1 if A[v][c]=1)
            #   Gained: paths through new arc v->u? No, that's A'[v][u] which adds path v->u->c, not u->...->c
            #   So delta(A^2[u,c]) = -(A[v][c] == 1) for c != v, c != u
            #
            # Similarly, A'^2[v,c] - A^2[v,c]:
            #   Gained: paths through new arc v->u, i.e., v->u->c where A[u][c]=1 (+1 if A[u][c]=1)
            #   Lost: nothing (no arc removed FROM v)
            #   Wait: A'[v][u]=1 (new), A'[u][c]=A[u][c]. So gained path v->u->c if A[u][c]=1.
            #   But also: path v->w->c where w uses... no, only the arcs touching {u,v} change.
            #   So delta(A^2[v,c]) = +(A[u][c] == 1) for c != u, c != v
            #
            # A'^2[a,v] - A^2[a,v]:
            #   Lost: path a->u->v where A[a][u]=1 (-1 if A[a][u]=1, since A[u][v]=1->0)
            #   So delta(A^2[a,v]) = -(A[a][u] == 1) for a != u, a != v
            #
            # A'^2[a,u] - A^2[a,u]:
            #   Gained: path a->v->u where A[a][v]=1 (+1 since A'[v][u]=1)
            #   So delta(A^2[a,u]) = +(A[a][v] == 1) for a != u, a != v

            term3 = 0
            for c in range(n):
                if c == u or c == v:
                    continue
                # Pair (u, c): reverse iff A[c][u]=1
                if A[c][u]:
                    a2_uc_before = sum(1 for w in range(n) if w != u and w != c and A[u][w] and A[w][c])
                    # delta(A^2[u,c]) = -1 if A[v][c]=1 (lost u->v->c)
                    a2_uc_after = a2_uc_before - (1 if A[v][c] else 0)
                    if a2_uc_before > 0 and a2_uc_after == 0:
                        term3 -= 1  # pair leaves J_2
                    elif a2_uc_before == 0 and a2_uc_after > 0:
                        term3 += 1  # shouldn't happen (can't increase)

                # Pair (v, c): reverse iff A[c][v]=1
                if A[c][v]:
                    a2_vc_before = sum(1 for w in range(n) if w != v and w != c and A[v][w] and A[w][c])
                    # delta(A^2[v,c]) = +1 if A[u][c]=1 (gained v->u->c)
                    a2_vc_after = a2_vc_before + (1 if A[u][c] else 0)
                    if a2_vc_before > 0 and a2_vc_after == 0:
                        term3 -= 1  # can't happen
                    elif a2_vc_before == 0 and a2_vc_after > 0:
                        term3 += 1  # pair enters J_2

            for a in range(n):
                if a == u or a == v:
                    continue
                # Pair (a, v): reverse iff A[v][a]=1
                if A[v][a]:
                    a2_av_before = sum(1 for w in range(n) if w != a and w != v and A[a][w] and A[w][v])
                    # delta(A^2[a,v]) = -1 if A[a][u]=1 (lost a->u->v)
                    a2_av_after = a2_av_before - (1 if A[a][u] else 0)
                    if a2_av_before > 0 and a2_av_after == 0:
                        term3 -= 1
                    elif a2_av_before == 0 and a2_av_after > 0:
                        term3 += 1

                # Pair (a, u): reverse iff A[u][a]=1
                if A[u][a]:
                    a2_au_before = sum(1 for w in range(n) if w != a and w != u and A[a][w] and A[w][u])
                    # delta(A^2[a,u]) = +1 if A[a][v]=1 (gained a->v->u)
                    a2_au_after = a2_au_before + (1 if A[a][v] else 0)
                    if a2_au_before > 0 and a2_au_after == 0:
                        term3 -= 1
                    elif a2_au_before == 0 and a2_au_after > 0:
                        term3 += 1

            predicted = term1 + term2 + term3
            if predicted != delta_J2:
                formula_errors += 1
                if formula_errors <= 5:
                    print(f"  ERROR: bits={bits}, flip ({u},{v}): predicted={predicted}, actual={delta_J2}")
                    print(f"    term1={term1}, term2={term2}, term3={term3}")

print(f"\nFormula errors: {formula_errors}")
if formula_errors == 0:
    print("  ** EXACT FORMULA VERIFIED for all n=5 arc flips **")


# ============================================================
# ANALYSIS 6: Simplify — is delta(J_2) determined by a small set of local invariants?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: delta(J_2) as function of local invariants")
print("=" * 70)

# The formula above is exact but complex (threshold crossings).
# Let's check if it simplifies to something nicer.

# Collect (d_u, d_v, a2_uv, a2_vu, term3, delta_J2) and see if term3 is determined
# by something simpler.

term3_by_local = defaultdict(Counter)
for bits in range(total):
    A = build_adj(n, bits)
    scores = get_scores(A, n)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue

            B = flip_arc(A, n, u, v)
            J2_a, pa = compute_J2(A, n)
            J2_b, pb = compute_J2(B, n)

            d_u = scores[u]
            d_v = scores[v]
            a2_uv = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[w][v])
            a2_vu = sum(1 for w in range(n) if w != u and w != v and A[v][w] and A[w][u])

            # Count: number of w with A[v][c]=1 AND A[c][u]=1 AND A^2[u,c]=1 (threshold pair)
            # These are the "fragile" pairs: exactly one 2-path u->?->c, and it goes through v.
            n_fragile_uc = 0
            for c in range(n):
                if c == u or c == v:
                    continue
                if A[c][u] and A[v][c]:  # reverse AND lost path through v
                    a2_uc = sum(1 for w in range(n) if w != u and w != c and A[u][w] and A[w][c])
                    if a2_uc == 1:  # exactly one path, which goes through v
                        n_fragile_uc += 1

            n_new_vc = 0
            for c in range(n):
                if c == u or c == v:
                    continue
                if A[c][v] and A[u][c]:  # reverse AND gained path through u
                    a2_vc = sum(1 for w in range(n) if w != v and w != c and A[v][w] and A[w][c])
                    if a2_vc == 0:  # was 0, now becomes 1
                        n_new_vc += 1

            n_fragile_av = 0
            for a in range(n):
                if a == u or a == v:
                    continue
                if A[v][a] and A[a][u]:
                    a2_av = sum(1 for w in range(n) if w != a and w != v and A[a][w] and A[w][v])
                    if a2_av == 1:
                        n_fragile_av += 1

            n_new_au = 0
            for a in range(n):
                if a == u or a == v:
                    continue
                if A[u][a] and A[a][v]:
                    a2_au = sum(1 for w in range(n) if w != a and w != u and A[a][w] and A[w][u])
                    if a2_au == 0:
                        n_new_au += 1

            delta_J2 = J2_b - J2_a
            term1 = 1 if a2_uv > 0 else 0
            term2 = -1 if a2_vu > 0 else 0
            expected_term3 = -n_fragile_uc + n_new_vc - n_fragile_av + n_new_au
            predicted = term1 + term2 + expected_term3

            key = (d_u, d_v)
            term3_by_local[key][delta_J2] += 1

print("delta(J_2) distribution by (d_u, d_v):")
for key in sorted(term3_by_local.keys()):
    vals = term3_by_local[key]
    print(f"  (d_u={key[0]}, d_v={key[1]}): {dict(sorted(vals.items()))}")


# ============================================================
# ANALYSIS 7: Does delta(dim Omega_2) have a clean formula?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: delta(dim Omega_2) = delta(|A_2|) - delta(J_2)")
print("=" * 70)

delta_Om2_dist = defaultdict(Counter)
delta_Om2_vals = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    scores = get_scores(A, n)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            B = flip_arc(A, n, u, v)

            A2a = compute_A2_sum(A, n)
            A2b = compute_A2_sum(B, n)
            J2a, _ = compute_J2(A, n)
            J2b, _ = compute_J2(B, n)

            d_u = scores[u]
            d_v = scores[v]

            dim_Om2_a = A2a - J2a
            dim_Om2_b = A2b - J2b
            delta_Om2 = dim_Om2_b - dim_Om2_a
            delta_A2 = A2b - A2a
            delta_J2 = J2b - J2a

            delta_Om2_dist[(d_u, d_v)][delta_Om2] += 1
            delta_Om2_vals[(d_u - d_v, delta_Om2)] += 1

print("delta(dim Omega_2) by (d_u, d_v):")
for key in sorted(delta_Om2_dist.keys()):
    vals = delta_Om2_dist[key]
    print(f"  (d_u={key[0]}, d_v={key[1]}): {dict(sorted(vals.items()))}")

print("\ndelta(dim Omega_2) by (d_u - d_v):")
by_diff = defaultdict(Counter)
for (diff, dOm), count in delta_Om2_vals.items():
    by_diff[diff][dOm] += count
for diff in sorted(by_diff.keys()):
    print(f"  d_u - d_v = {diff}: {dict(sorted(by_diff[diff].items()))}")


# ============================================================
# ANALYSIS 8: delta(beta_1) under arc flips
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS 8: delta(beta_1) under arc flips at n=5")
print("=" * 70)

delta_b1_dist = defaultdict(Counter)
delta_b1_combined = defaultdict(Counter)

for bits in range(total):
    A = build_adj(n, bits)
    scores = get_scores(A, n)
    b1_a = compute_beta1(A, n)

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            B = flip_arc(A, n, u, v)
            b1_b = compute_beta1(B, n)
            d_u = scores[u]
            d_v = scores[v]
            delta_b1 = b1_b - b1_a

            delta_b1_dist[(d_u, d_v)][delta_b1] += 1

            # Also compute delta(dim Omega_2)
            A2a = compute_A2_sum(A, n)
            A2b = compute_A2_sum(B, n)
            J2a, _ = compute_J2(A, n)
            J2b, _ = compute_J2(B, n)
            delta_Om2 = (A2b - J2b) - (A2a - J2a)

            # Key combination: delta(dim Omega_2) + delta(beta_1)
            # This SHOULD equal delta(rk d_3) if beta_2 = 0
            combo = delta_Om2 + delta_b1
            delta_b1_combined[(d_u, d_v)][combo] += 1

print("delta(beta_1) by (d_u, d_v):")
for key in sorted(delta_b1_dist.keys()):
    vals = delta_b1_dist[key]
    if len(vals) > 1:
        marker = " <-- VARIES"
    else:
        marker = ""
    print(f"  (d_u={key[0]}, d_v={key[1]}): {dict(sorted(vals.items()))}{marker}")

print("\ndelta(dim Omega_2) + delta(beta_1) by (d_u, d_v):")
for key in sorted(delta_b1_combined.keys()):
    vals = delta_b1_combined[key]
    if len(vals) > 1:
        marker = " <-- VARIES"
    else:
        marker = ""
    print(f"  (d_u={key[0]}, d_v={key[1]}): {dict(sorted(vals.items()))}{marker}")


print("\nDone.")
