#!/usr/bin/env python3
"""
beta2_cancel_structure.py — Understand WHY cancel pairs cover the DT deficit

At n=6, DT fails for 960 tournaments (deficit 1-2). Cancel pairs fix them all.
What EXACTLY does a cancel pair contribute that DT doesn't?

KEY OBSERVATION: A cancel pair (a,b₁,c,d)-(a,b₂,c,d) with shared bad face (a,c,d)
has boundary:
  [(b₁,c,d)-(b₂,c,d)] + [(a,b₁,d)-(a,b₂,d)] - [(a,b₁,c)-(a,b₂,c)]
(The (a,c,d) terms cancel.)

This is a "substitution" of b₁ for b₂ in the three good faces.
The substitution changes the TT/NT status of the faces.

QUESTION: What is the cancel pair's image in Z₂?

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


print("=" * 70)
print("CANCEL PAIR STRUCTURE ANALYSIS")
print("=" * 70)

# ============================================================
# PART 1: At n=6, analyze the deficit tournaments
# ============================================================
print("\nPART 1: n=6 deficit tournament analysis")
n = 6
m = n*(n-1)//2

deficit_cases = []
t0 = time.time()

for bits in range(1 << m):
    A = build_adj(n, bits)
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0 or not ap3: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    dt_idx = [i for i, p in enumerate(ap3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    bd3 = build_full_boundary_matrix(ap3, ap2)

    if dt_idx:
        bd3_dt = np.column_stack([bd3[:, i] for i in dt_idx])
        coords_dt = np.linalg.lstsq(om2, bd3_dt, rcond=None)[0]
        rk_dt = np.linalg.matrix_rank(coords_dt, tol=1e-8)
    else:
        rk_dt = 0

    deficit = z2_dim - rk_dt
    if deficit > 0 and len(deficit_cases) < 50:
        scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
        c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if A[i][j]+A[j][k]+A[k][i] in (0,3))

        # Find which Z₂ direction is NOT covered by DT
        z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T
        z2_proj_dt = z2_basis.T @ coords_dt if dt_idx else np.zeros((z2_dim, 0))
        # SVD to find the uncovered direction
        if z2_proj_dt.size > 0:
            U, S, Vt = np.linalg.svd(z2_proj_dt, full_matrices=True)
            # Uncovered directions: columns of U corresponding to near-zero singular values
            uncovered = U[:, sum(s > 1e-8 for s in S):]
        else:
            uncovered = np.eye(z2_dim)

        # Convert uncovered Z₂ direction to Ω₂ coords
        uncov_om2 = z2_basis @ uncovered[:, 0]
        # Convert to A₂ coords
        uncov_a2 = om2 @ uncov_om2

        # Which A₂ paths are in the uncovered cycle?
        uncov_paths = []
        ap2_list = [tuple(x) for x in ap2]
        for i, p in enumerate(ap2_list):
            if abs(uncov_a2[i]) > 1e-8:
                a, b, c = p
                is_tt = bool(A[a][c])
                uncov_paths.append((p, round(uncov_a2[i], 4), is_tt))

        # Which cancel pairs cover this direction?
        bad_groups = defaultdict(list)
        for i, p in enumerate(ap3):
            a, b, c, d = p
            if not A[a][c]: bad_groups[('02', a, c)].append(i)
            if not A[b][d]: bad_groups[('13', b, d)].append(i)

        # Find cancel pairs that contribute to the uncovered direction
        useful_cancels = []
        for key, indices in bad_groups.items():
            if len(indices) >= 2:
                for j in range(1, len(indices)):
                    cv = np.zeros(len(ap3))
                    cv[indices[0]] = 1; cv[indices[j]] = -1
                    bd_cancel = bd3 @ cv
                    # Project to Z₂
                    cancel_om2 = np.linalg.lstsq(om2, bd_cancel, rcond=None)[0]
                    cancel_z2 = z2_basis.T @ cancel_om2
                    # Does this cover the uncovered direction?
                    proj = uncovered[:, 0] @ cancel_z2
                    if abs(proj) > 1e-8:
                        p1 = tuple(ap3[indices[0]])
                        p2 = tuple(ap3[indices[j]])
                        useful_cancels.append((key, p1, p2, round(proj, 4)))

        deficit_cases.append({
            'bits': bits, 'scores': scores, 'c3': c3,
            'deficit': deficit, 'z2_dim': z2_dim, 'rk_dt': rk_dt,
            'uncov_paths': uncov_paths,
            'useful_cancels': useful_cancels,
        })

    if bits % 10000 == 0 and bits > 0:
        elapsed = time.time() - t0
        print(f"  ... {bits}/{1 << m} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"Completed in {elapsed:.0f}s")
print(f"Found {len(deficit_cases)} deficit cases")

# Analyze the deficit cases
score_dist = Counter(d['scores'] for d in deficit_cases)
print(f"\nScore distribution of deficit cases:")
for scores, count in sorted(score_dist.items()):
    print(f"  {scores}: {count}")

# Show detailed examples
for case in deficit_cases[:5]:
    print(f"\n  T#{case['bits']} scores={case['scores']}, c₃={case['c3']}")
    print(f"    deficit={case['deficit']}, Z₂={case['z2_dim']}, rk_DT={case['rk_dt']}")

    print(f"    Uncovered Z₂ cycle ({len(case['uncov_paths'])} terms):")
    for path, coeff, is_tt in sorted(case['uncov_paths'])[:8]:
        label = "TT" if is_tt else "NT"
        print(f"      {coeff:+.4f} * {path} [{label}]")
    if len(case['uncov_paths']) > 8:
        print(f"      ... ({len(case['uncov_paths'])} total)")

    print(f"    Useful cancel pairs ({len(case['useful_cancels'])}):")
    for key, p1, p2, proj in case['useful_cancels'][:5]:
        print(f"      {key}: {p1} - {p2}, projection={proj}")

# ============================================================
# PART 2: What makes the uncovered cycle special?
# ============================================================
print(f"\n{'='*70}")
print("PART 2: UNCOVERED CYCLE STRUCTURE")
print("=" * 70)

# Are the uncovered cycles predominantly NT?
tt_fraction = []
for case in deficit_cases:
    n_tt = sum(1 for _, _, is_tt in case['uncov_paths'] if is_tt)
    n_total = len(case['uncov_paths'])
    if n_total > 0:
        tt_fraction.append(n_tt / n_total)

if tt_fraction:
    avg_tt = sum(tt_fraction) / len(tt_fraction)
    print(f"  Average TT fraction in uncovered cycles: {avg_tt:.3f}")
    print(f"  Min TT fraction: {min(tt_fraction):.3f}")
    print(f"  Max TT fraction: {max(tt_fraction):.3f}")

# How many of the uncovered cycle's A₂ paths have NO DT extension?
# A path (a,b,c) has a DT extension iff ∃d with c→d, a→c, b→d.
# Since a→c is needed for DT and the path is (a,b,c) with a→b, b→c:
# For TT: a→c already. Need d with c→d, b→d (=common out-neighbor of b and c).
# For NT: c→a. Can't have a→c. So NO DT extension possible for NT paths!

print(f"\n  NT paths cannot have DT extensions (c→a blocks a→c requirement)")
print(f"  So cycles with NT terms can only be DT-filled via TT components")

# Check: in deficit cases, does the uncovered cycle ALWAYS have NT terms?
has_nt = sum(1 for case in deficit_cases if any(not is_tt for _, _, is_tt in case['uncov_paths']))
print(f"\n  Uncovered cycles with NT terms: {has_nt}/{len(deficit_cases)}")

# If the uncovered cycle has NT terms, those terms must come from the
# Ω₂ constraints. The NT paths contribute to Ω₂ when the junk cancellation
# conditions are satisfied. The cancel pairs provide "NT corrections" that
# DT paths can't.

# ============================================================
# PART 3: The cancel pair's mechanism
# ============================================================
print(f"\n{'='*70}")
print("PART 3: CANCEL PAIR MECHANISM")
print("=" * 70)

# For a cancel pair (a,b₁,c,d)-(a,b₂,c,d) sharing bad face (a,c,d) [c→a]:
# Boundary = (b₁,c,d)-(b₂,c,d) + (a,b₁,d)-(a,b₂,d) - (a,b₁,c)+(a,b₂,c)
#
# The paths (a,b₁,c) and (a,b₂,c) are both A₂ paths with a→b, b→c.
# For cancel pair of type 02 (bad face involves c→a):
#   (a,b₁,c): a→b₁, b₁→c. Is a→c? NO (c→a, that's WHY it's a bad face).
#   So (a,b₁,c) and (a,b₂,c) are BOTH NT paths!
#
# The boundary of the cancel pair includes NT DIFFERENCES: (a,b₁,c)-(a,b₂,c).
# In Ω₂, these NT differences might represent nontrivial directions that
# DT paths (whose boundaries only involve TT paths) can't reach!

print("Cancel pairs of type 02 (bad face c→a) produce NT differences in A₂.")
print("This is WHY they can fill Z₂ directions that DT paths cannot!")

# Verify: type 02 cancel pairs have NT boundary terms
n_cancel_with_nt = 0
n_cancel_total = 0
for case in deficit_cases[:10]:
    bits = case['bits']
    A = build_adj(6, bits)
    for key, p1, p2, proj in case['useful_cancels']:
        n_cancel_total += 1
        a1, b1, c1, d1 = p1
        a2, b2, c2, d2 = p2
        if key[0] == '02':
            # Shared (a,c) with c→a. Faces (a,b₁,c) and (a,b₂,c) are NT.
            # Also faces (a,b₁,d) and (a,b₂,d): check if TT or NT
            # (b₁,c,d) and (b₂,c,d): check TT/NT
            nt_terms = 0
            for face_path in [(a1,b1,c1), (a2,b2,c2)]:
                a, b, c = face_path
                if not A[a][c]:
                    nt_terms += 1
            if nt_terms > 0:
                n_cancel_with_nt += 1

if n_cancel_total > 0:
    print(f"\n  Cancel pairs with NT boundary terms: {n_cancel_with_nt}/{n_cancel_total}")


# ============================================================
# PART 4: DT boundary spans ONLY the TT part of Z₂?
# ============================================================
print(f"\n{'='*70}")
print("PART 4: DT BOUNDARY vs TT/NT DECOMPOSITION")
print("=" * 70)

# At n=5 (where DT alone works): this suggests Ω₂ = TT at n=5.
# But we know dim(Ω₂) > |TT| for many n=5 tournaments.
# How can DT fill NT-containing Z₂ cycles?
#
# Answer: The DT boundary CAN produce NT terms!
# DT(a,b,c,d) with d→a: faces (a,c,d) and (a,b,d) are NT.
# So DT paths with d→a produce NT terms in their boundary.
# At n=5, there are enough DT paths with both a→d and d→a
# to cover all Z₂ directions (both TT and NT).
# At n=6, the NT directions become harder to reach.

print("\nFor DT paths with d→a: boundary includes NT faces (a,c,d) and (a,b,d)")
print("For DT paths with a→d: ALL boundary faces are TT")
print("So DT paths with d→a provide NT contributions.")

# Count DT paths by a↔d type at n=5
n = 5
m = n*(n-1)//2
ad_dist = Counter()
for bits in range(1 << m):
    A = build_adj(n, bits)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap3: continue
    for p in ap3:
        a, b, c, d = p
        if A[a][c] and A[b][d]:  # DT
            if A[a][d]:
                ad_dist['a→d'] += 1
            else:
                ad_dist['d→a'] += 1

total_dt = sum(ad_dist.values())
print(f"\nn=5 DT paths by a↔d type:")
for key in ['a→d', 'd→a']:
    print(f"  {key}: {ad_dist[key]} ({100*ad_dist[key]/total_dt:.1f}%)")

# At n=6: same analysis
n = 6
m = n*(n-1)//2
ad_dist6 = Counter()
n_sampled = 0
import random
random.seed(42)
for bits in random.sample(range(1 << m), 2000):
    A = build_adj(n, bits)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap3: continue
    n_sampled += 1
    for p in ap3:
        a, b, c, d = p
        if A[a][c] and A[b][d]:
            if A[a][d]:
                ad_dist6['a→d'] += 1
            else:
                ad_dist6['d→a'] += 1

total_dt6 = sum(ad_dist6.values())
print(f"\nn=6 DT paths by a↔d type ({n_sampled} tournaments):")
for key in ['a→d', 'd→a']:
    print(f"  {key}: {ad_dist6[key]} ({100*ad_dist6[key]/total_dt6:.1f}%)")


# ============================================================
# PART 5: THE ALGEBRAIC KEY — decompose ∂₃ image by face type
# ============================================================
print(f"\n{'='*70}")
print("PART 5: ∂₃ IMAGE DECOMPOSITION")
print("=" * 70)

# Every DT boundary can be decomposed as:
# ∂₃(a,b,c,d) = [(b,c,d)_TT - (a,b,c)_TT] + [(a,b,d)_? - (a,c,d)_?]
#
# For a→d: both (a,b,d) and (a,c,d) are TT. Pure TT.
# For d→a: both (a,b,d) and (a,c,d) are NT.
#   So: TT_part = (b,c,d)-(a,b,c), NT_part = (a,b,d)-(a,c,d)
#
# The NT part involves "vertex substitution": in the NT pair (a,b,d) and (a,c,d),
# we substitute b for c while keeping a and d fixed.
# Since (a,b,d) is NT: a→b→d→a cycle.
# And (a,c,d) is NT: a→c→d→a cycle.
# Both share the arc d→a.
#
# This is a "b↔c swap within the d→a cycle."

# What about cancel pairs?
# Cancel pair (a,b₁,c,d)-(a,b₂,c,d) with bad face (a,c,d), c→a:
# Boundary = [(b₁,c,d)-(b₂,c,d)] + [(a,b₁,d)-(a,b₂,d)] - [(a,b₁,c)-(a,b₂,c)]
#
# The term -(a,b₁,c)+(a,b₂,c) is a "b₁↔b₂ swap within (a,*,c)" NT paths.
# Since c→a, both (a,b₁,c) and (a,b₂,c) are NT.
# So the cancel pair contributes NT DIFFERENCES.

# KEY INSIGHT: The difference between DT and cancel is:
# - DT paths with a→d give PURE TT boundaries
# - DT paths with d→a give TT + NT boundaries (NT comes in pairs)
# - Cancel pairs give NT DIFFERENCES (b₁↔b₂ substitutions)
#
# At n=5, DT-with-d→a provides enough NT to cover Z₂.
# At n=6+, some Z₂ directions need additional NT from cancel pairs.
#
# The Z₂ cycle space has a TT-part and an NT-part.
# The TT-part is always covered by DT (TT-only).
# The NT-part needs DT-with-d→a + cancel pairs.

print("DT(a→d): pure TT boundary")
print("DT(d→a): TT + NT boundary (paired NT terms)")
print("Cancel:  NT differences (b-substitution)")
print("\nZ₂ = TT_part + NT_part")
print("TT_part covered by DT(a→d)")
print("NT_part covered by DT(d→a) + cancel")


# ============================================================
# PART 6: COUNT — is |NT directions in Z₂| ≤ |DT(d→a)| + |cancel|?
# ============================================================
print(f"\n{'='*70}")
print("PART 6: DIMENSIONAL ANALYSIS")
print("=" * 70)

n = 5
m = n*(n-1)//2

for bits in range(min(50, 1 << m)):
    A = build_adj(n, bits)
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap3: continue

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    # Count DT types
    n_dt_ad = sum(1 for p in ap3 if A[p[0]][p[2]] and A[p[1]][p[3]] and A[p[0]][p[3]])
    n_dt_da = sum(1 for p in ap3 if A[p[0]][p[2]] and A[p[1]][p[3]] and A[p[3]][p[0]])

    # Count cancel pairs
    bad_groups = defaultdict(list)
    for i, p in enumerate(ap3):
        a, b, c, d = p
        if not A[a][c]: bad_groups[('02', a, c)].append(i)
        if not A[b][d]: bad_groups[('13', b, d)].append(i)
    n_cancel = sum(len(g)-1 for g in bad_groups.values() if len(g) >= 2)

    # Count NT paths in Ω₂
    ap2_list = [tuple(x) for x in ap2]
    n_nt_in_a2 = sum(1 for p in ap2_list if not A[p[0]][p[2]])

    scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if A[i][j]+A[j][k]+A[k][i] in (0,3))

    if n_dt_da + n_cancel > 0 or c3 > 0:
        print(f"  T#{bits} scores={scores}, c₃={c3}: "
              f"Z₂={z2_dim}, DT(a→d)={n_dt_ad}, DT(d→a)={n_dt_da}, "
              f"cancel={n_cancel}, NT_in_A₂={n_nt_in_a2}")


print("\nDone.")
