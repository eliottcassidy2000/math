#!/usr/bin/env python3
"""
beta2_arcflip_algebra.py - Algebraic analysis of the arc-flip invariance

GOAL: Understand the algebraic mechanism behind delta(rk d_3) = delta(dim Z_2).

Approach: For the flip u->v to v->u, decompose Omega_2(T) and Omega_2(T')
into "stable" and "changed" parts, and track how Z_2 and im(d_3) change.

Key decomposition:
  Omega_2(T) = S + L  (stable + lost)
  Omega_2(T') = S + G  (stable + gained)

where S = span{chains in Omega_2 not involving the flipped arc},
L = elements of Omega_2(T) involving u->v,
G = elements of Omega_2(T') involving v->u.

Track:
1. dim(Z_2 cap S) -- how much of Z_2 lives in the stable part
2. How the "unstable" part of Z_2 changes

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter
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

def flip_arc(bits, i, j, n):
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return bits ^ (1 << idx)
            idx += 1
    return bits


def analyze_flip(A, n, u, v):
    """Detailed analysis of arc flip u->v to v->u.

    Returns dict with all relevant dimensions and their changes.
    """
    # Build A' (after flip)
    A2 = [row[:] for row in A]
    A2[u][v] = 0
    A2[v][u] = 1

    results = {}

    for label, adj in [('T', A), ('T\'', A2)]:
        paths = {}
        omega = {}
        for p in range(5):
            paths[p] = enumerate_allowed_paths(adj, n, p)
            if p == 0:
                omega[p] = np.eye(n)
            elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
                omega[p] = compute_omega_basis(adj, n, p, paths[p], paths[p-1])
            else:
                omega[p] = np.zeros((max(1, len(paths[p])), 0))

        dim2 = omega[2].shape[1] if omega[2].ndim == 2 else 0
        dim3 = omega[3].shape[1] if omega[3].ndim == 2 else 0

        # Identify TT triples (individually d-inv)
        tt = [(a,b,c) for (a,b,c) in paths[2] if adj[a][c]]
        # NT paths in Omega_2 (linear combos where junk cancels)
        nt_in_omega = dim2 - len(tt)

        # Boundary maps
        if dim2 > 0:
            bd2 = build_full_boundary_matrix(paths[2], paths[1])
            bd2_om = bd2 @ omega[2]
            sv2 = np.linalg.svd(bd2_om, compute_uv=False)
            rk2 = int(np.sum(np.abs(sv2) > 1e-8))
        else:
            rk2 = 0

        if dim3 > 0 and dim2 > 0:
            bd3 = build_full_boundary_matrix(paths[3], paths[2])
            bd3_om = bd3 @ omega[3]
            im3_c, _, _, _ = np.linalg.lstsq(omega[2], bd3_om, rcond=None)
            rk3 = np.linalg.matrix_rank(im3_c, tol=1e-8)
        else:
            rk3 = 0

        z2 = dim2 - rk2
        b2 = z2 - rk3

        # Count TT triples involving the flipped arc
        tt_with_uv = [t for t in tt if (t[0]==u and t[1]==v) or
                      (t[1]==u and t[2]==v) or
                      (t[0]==u and t[2]==v)]

        # Count 3-cycles involving u,v
        cycles_uv = 0
        for w in range(n):
            if w == u or w == v:
                continue
            # Check if {u,v,w} is a 3-cycle
            # In original T, arc u->v exists
            if adj[u][v] and adj[v][w] and adj[w][u]:
                cycles_uv += 1
            elif adj[u][w] and adj[w][v] and adj[v][u]:
                cycles_uv += 1

        results[label] = {
            'dim2': dim2, 'dim3': dim3, 'rk2': rk2, 'rk3': rk3,
            'z2': z2, 'b2': b2,
            'tt': len(tt), 'nt_omega': nt_in_omega,
            'tt_with_uv': len(tt_with_uv),
            'cycles_uv': cycles_uv,
        }

    return results


# Analyze all flips at n=5
n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"ALGEBRAIC DECOMPOSITION OF ARC-FLIP at n={n}")
print("=" * 70)

# Classify by (dTT, dNT, dZ2, drk3) to understand the balance
classify = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(row) for row in A]

    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                u, v = i, j
            else:
                u, v = j, i

            r = analyze_flip(A, n, u, v)
            t = r['T']
            t2 = r["T'"]

            dTT = t2['tt'] - t['tt']
            dNT = t2['nt_omega'] - t['nt_omega']
            dZ2 = t2['z2'] - t['z2']
            drk3 = t2['rk3'] - t['rk3']
            dDim2 = t2['dim2'] - t['dim2']

            # Sanity check
            assert dDim2 == dTT + dNT, f"dim2 decomp failed: {dDim2} != {dTT} + {dNT}"

            classify[(dTT, dNT, dDim2, dZ2, drk3)] += 1

print("\nClassification by (dTT, dNT, dDim2, dZ2, drk3):")
print("  dTT=change in TT triples, dNT=change in NT-cancellation dim")
print("  dDim2=dTT+dNT, dZ2=change in ker(d2), drk3=change in rk(d3)")
print()
for key in sorted(classify.keys()):
    dTT, dNT, dDim2, dZ2, drk3 = key
    ok = "OK" if dZ2 == drk3 else "FAIL"
    print(f"  dTT={dTT:+d}, dNT={dNT:+d} => dDim2={dDim2:+d}, "
          f"dZ2={dZ2:+d}, drk3={drk3:+d}: {classify[key]} flips [{ok}]")

# Key question: is dNT ALWAYS = -dTT + (p_vu - p_uv)?
# We know dDim2 = p_vu - p_uv and dTT = -(#TT lost) + (#TT gained)
# TT lost = #TT triples using arc u->v = cout + cin + p_uv
# TT gained = #TT triples using arc v->u = cout + cin + p_vu
# So dTT = p_vu - p_uv (same as dDim2!)
# Wait, that means dNT = 0 always?

print("\n--- Checking if dNT = 0 always ---")
always_zero = all(key[1] == 0 for key in classify.keys())
print(f"  dNT = 0 for all flips: {always_zero}")

if not always_zero:
    print("  SURPRISE: dNT is NOT always 0!")
    print("  Cases where dNT != 0:")
    for key, count in sorted(classify.items()):
        if key[1] != 0:
            print(f"    dTT={key[0]:+d}, dNT={key[1]:+d}, count={count}")

# Wait: from the example analysis, we saw TT went from 8 to 7 but dim(Om2) stayed at 9.
# That means dTT = -1, dNT = +1, dDim2 = 0. So dNT != 0 in general.
# But our formula says dDim2 = p_vu - p_uv. With dTT = p_vu - p_uv also?
# NO! dTT = (gained TT - lost TT). Let me recheck.

print("\n--- Detailed TT gain/loss analysis ---")
for bits in range(min(total, 50)):
    A = build_adj(n, bits)
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                u, v = i, j
            else:
                u, v = j, i

            # Count by case
            cout = cin = p_uv = p_vu = 0
            for w in range(n):
                if w == u or w == v:
                    continue
                if A[u][w] and A[v][w]:
                    cout += 1
                elif A[w][u] and A[w][v]:
                    cin += 1
                elif A[u][w] and A[w][v]:
                    p_uv += 1
                elif A[v][w] and A[w][u]:
                    p_vu += 1

            # TT triples
            tt_T = [(a,b,c) for a in range(n) for b in range(n) for c in range(n)
                    if a!=b and b!=c and a!=c and A[a][b] and A[b][c] and A[a][c]]
            A2 = [row[:] for row in A]
            A2[u][v] = 0
            A2[v][u] = 1
            tt_T2 = [(a,b,c) for a in range(n) for b in range(n) for c in range(n)
                     if a!=b and b!=c and a!=c and A2[a][b] and A2[b][c] and A2[a][c]]

            dTT = len(tt_T2) - len(tt_T)
            predicted_dTT = p_vu - p_uv  # from our earlier analysis

            if dTT != predicted_dTT and bits < 10:
                print(f"  bits={bits}, flip {u}->{v}: dTT={dTT}, predicted={predicted_dTT}, "
                      f"cout={cout}, cin={cin}, p_uv={p_uv}, p_vu={p_vu}")

# Hmm: earlier the formula was:
# Lost TT: cout + cin + p_uv (triples involving arc u->v)
# Gained TT: cout + cin + p_vu (triples involving arc v->u)
# Net dTT = (cout+cin+p_vu) - (cout+cin+p_uv) = p_vu - p_uv

# But in the example, we had dTT = 7-8 = -1 but dDim2 = 9-9 = 0.
# This means dNT = +1 compensated.
# If dTT = p_vu - p_uv = dDim2, then dNT should be 0. Contradiction!

# Let me recheck: in the example (bits=4, flip 4->3 to 3->4):
# u=4, v=3 (since A[4][3]=1 means 4->3)
# Lost TT: triples using arc 4->3
#   (a,4,3) TT: need a->4, 4->3, a->3. a->4: nobody (4->everyone). So 0.
#   (4,b,3) TT: need 4->b, b->3, 4->3. b with 4->b AND b->3:
#     b=0: 4->0 yes, 0->3 yes. (4,0,3) TT.
#     b=1: 4->1 yes, 1->3? NO (3->1). No.
#     b=2: 4->2 yes, 2->3? NO (3->2). No.
#     So just (4,0,3).
#   (4,3,c) TT: need 4->3, 3->c, 4->c. c with 3->c AND 4->c:
#     c=1: 3->1 yes, 4->1 yes. (4,3,1) TT.
#     c=2: 3->2 yes, 4->2 yes. (4,3,2) TT.
#     So (4,3,1) and (4,3,2).
# Total lost TT: (4,0,3), (4,3,1), (4,3,2) = 3.
# Lost = cout(4,3) + cin(4,3) + p_43 = ?
# For 4->3: w != 4,3:
#   w=0: 4->0, 3->0? A[3][0]=0, A[0][3]=1. So 0->3, meaning 3 does NOT go to 0.
#         Actually: A[3][0]=? Let me check: 3's out-neighbors are {1,2} (A[3][1]=1, A[3][2]=1).
#         So A[3][0]=0. And A[0][3]=1. So 0->3.
#         4->0 and 0->3: so w=0 is case p_43 (path 4->0->3). No: p_uv = paths u->w->v = 4->w->3.
#         w=0: 4->0 and 0->3? YES. So w=0 is a p_43 case.
#   w=1: 4->1 and 1->3? A[1][3]=0, A[3][1]=1. So 3->1, NOT 1->3.
#         4->1 and 3->1. So both point TO 1? No: 4->1 and 3->1 means out from {4,3} to 1.
#         Actually: u=4, v=3. A[4][1]=1 (4->1), A[3][1]=1 (3->1).
#         So both u and v go to w=1: cout(4,3) for w=1.
#   w=2: A[4][2]=1 (4->2), A[3][2]=1 (3->2). Same: cout for w=2.
#   So: cout(4,3) = 2 (w=1,2), cin(4,3) = 0, p_43 = 1 (w=0), p_34 = 0.
#   Lost = cout + cin + p_43 = 2 + 0 + 1 = 3. Matches!
#   Gained = cout + cin + p_34 = 2 + 0 + 0 = 2.
#   dTT = 2 - 3 = -1.
#   dDim2 = p_34 - p_43 = 0 - 1 = -1.
#   Wait! dDim2 should be 0 (9->9), but the formula gives -1!

# Hmm, that contradicts. Let me recheck the formula.
# Actually, I think I got the direction wrong. When we flip u->v to v->u,
# the SOURCE of the flip is u (who had the arc) and TARGET is v.
# Lost TT triples use arc u->v.
# Gained TT triples use arc v->u.
# p_vu = paths v->w->u = |{w: v->w, w->u}|
# After flip, v->u is the new arc. TT triples using v->u:
#   (a,v,u) TT: a->v, v->u, a->u
#   (v,b,u) TT: v->b, b->u, v->u
#   (v,u,c) TT: v->u, u->c, v->c

# For (v,b,u): need v->b AND b->u. This is p_vu case! v->w AND w->u.
# But wait, in T (before flip), u->v. After flip (T'), v->u.
# The vertex w must satisfy T'[v][w]=1 AND T'[w][u]=1.
# T'[v][w] = T[v][w] (unchanged for w != u)
# T'[w][u] = T[w][u] (unchanged for w != v)
# So p_vu in T = |{w: T[v][w]=1 AND T[w][u]=1}|

# In our example: v=3, u=4. T[3][w] for w=0,1,2: 3->1, 3->2 (NOT 3->0).
# T[w][4] for w=0,1,2: nobody goes to 4 (4->everyone).
# So p_34 = 0. Correct.

# p_43 = |{w: T[4][w]=1 AND T[w][3]=1}|
# T[4][w]=1 for w=0,1,2 (4->everyone except 3, well 4->3 also but w!=3)
# T[w][3]=1: A[w][3]=1 means w->3. A[0][3]=1 (0->3), A[1][3]=0, A[2][3]=0.
# So p_43 = 1 (w=0 only).

# Formula says dDim2 = p_vu - p_uv = p_34 - p_43 = 0 - 1 = -1.
# But actual dDim2 = 9 - 9 = 0!

# ERROR: my formula is wrong?! Let me recheck the verification from earlier.
# At n=5, the script verified dTT (the number of TT triples) correctly.
# But dim(Omega_2) != |TT triples|. dim(Omega_2) = |TT| + |NT cancellation|.
# So dDim2 = dTT + dNT. And dTT = p_vu - p_uv (verified).
# So dDim2 = (p_vu - p_uv) + dNT.
# In the example: dDim2 = 0, dTT = -1, so dNT = +1.

# But the first script (beta2_arcflip_exactness.py) verified that
# "delta(Omega_2)" = p_vu - p_uv. Wait, it verified delta(dim Omega_2)
# against the count of TT triples. But dim(Omega_2) is NOT |TT|!
# Let me recheck...

# Actually the first script computed delta(Omega_2) using
# compute_omega_basis (the correct dim), and separately computed
# delta(TT) using direct counting. It verified delta(TT) matches
# the formula. But it did NOT claim delta(dim Omega_2) = delta(TT)!

# SO: delta(dim Omega_2) = delta(TT) + delta(NT cancellation)
# And the question is: what is delta(NT cancellation)?

print("\n" + "=" * 70)
print("CRITICAL: delta(dim Omega_2) decomposition")
print("=" * 70)
print("dim(Omega_2) = |TT| + |NT cancellation dims|")
print("delta(dim Omega_2) = delta(TT) + delta(NT)")
print("delta(TT) = p_vu - p_uv (PROVED)")
print("Need: delta(NT cancellation) formula")

# Compute delta(NT) for all flips at n=5
dNT_by_local = {}
for bits in range(total):
    A = build_adj(n, bits)
    p1 = enumerate_allowed_paths(A, n, 2)
    o1 = compute_omega_basis(A, n, 2, p1, enumerate_allowed_paths(A, n, 1))
    dim1 = o1.shape[1] if o1.ndim == 2 else 0
    tt1 = sum(1 for (a,b,c) in p1 if A[a][c])

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            A2 = build_adj(n, bits2)
            p2 = enumerate_allowed_paths(A2, n, 2)
            o2 = compute_omega_basis(A2, n, 2, p2, enumerate_allowed_paths(A2, n, 1))
            dim2 = o2.shape[1] if o2.ndim == 2 else 0
            tt2 = sum(1 for (a,b,c) in p2 if A2[a][c])

            dTT = tt2 - tt1
            dDim = dim2 - dim1
            dNT_val = dDim - dTT

            if A[i][j] == 1:
                u, v = i, j
            else:
                u, v = j, i

            cout = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[v][w])
            cin = sum(1 for w in range(n) if w != u and w != v and A[w][u] and A[w][v])
            p_uv = sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[w][v])
            p_vu = sum(1 for w in range(n) if w != u and w != v and A[v][w] and A[w][u])

            key = (dTT, dNT_val, dDim)
            if key not in dNT_by_local:
                dNT_by_local[key] = 0
            dNT_by_local[key] += 1

print("\n(dTT, dNT, dDim2) distribution:")
for key in sorted(dNT_by_local.keys()):
    dTT, dNT_val, dDim = key
    print(f"  dTT={dTT:+d}, dNT={dNT_val:+d} => dDim2={dDim:+d}: {dNT_by_local[key]} flips")

print("\nDone.")
