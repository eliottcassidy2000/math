#!/usr/bin/env python3
"""
Detailed analysis of orbit boundary maps at critical degrees d=m, m+1.

For P_7 (m=3): look at ∂_4^orb (the critical boundary) to understand why β_3 = 0.
For P_11 (m=5): look at ∂_6^orb to understand why β_5 = 1.

The orbit complex is small enough for P_7 and P_11 to look at the actual matrices.

opus-2026-03-13-S71b
"""
import numpy as np
from fractions import Fraction

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_orbit_reps(p, d):
    QR = set(pow(x, 2, p) for x in range(1, p))
    QR_list = sorted(QR)
    if d == 0:
        return [()]
    results = []
    def is_canonical(seq):
        for q in QR_list:
            if q == 1: continue
            scaled = tuple((q * s) % p for s in seq)
            if scaled < seq:
                return False
        return True
    def backtrack(seq, ps_list, ps_set):
        if len(seq) == d:
            t = tuple(seq)
            if is_canonical(t):
                results.append(t)
            return
        for s in QR_list:
            new_ps = (ps_list[-1] + s) % p
            if new_ps in ps_set:
                continue
            seq.append(s)
            ps_list.append(new_ps)
            ps_set.add(new_ps)
            backtrack(seq, ps_list, ps_set)
            seq.pop()
            ps_list.pop()
            ps_set.remove(new_ps)
    backtrack([1], [0, 1], {0, 1})
    return sorted(results)

def zm_orbit_class(seq, QR_list, p):
    best = seq
    for q in QR_list:
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best

def compute_face(seq, i, p):
    d = len(seq)
    if i == 0: return seq[1:]
    elif i == d: return seq[:-1]
    else:
        merged = (seq[i-1] + seq[i]) % p
        return seq[:i-1] + (merged,) + seq[i+1:]

def build_orbit_boundary(p, d, orbit_reps_d, orbit_reps_dm1, QR):
    reps_dm1_idx = {rep: i for i, rep in enumerate(orbit_reps_dm1)}
    B = np.zeros((len(orbit_reps_dm1), len(orbit_reps_d)), dtype=float)
    for j, sigma in enumerate(orbit_reps_d):
        for fi in range(d + 1):
            face = compute_face(sigma, fi, p)
            face_canon = zm_orbit_class(face, QR, p)
            if face_canon in reps_dm1_idx:
                B[reps_dm1_idx[face_canon], j] += (-1) ** fi
    return B

def build_orbit_constraint(p, d, reps, QR):
    QR_set = set(QR)
    junk_orbits = {}
    for sigma in reps:
        for i in range(1, d):
            merged = (sigma[i-1] + sigma[i]) % p
            if merged not in QR_set:
                face = list(sigma)
                face[i-1] = merged
                del face[i]
                ft = zm_orbit_class(tuple(face), QR, p)
                if ft not in junk_orbits:
                    junk_orbits[ft] = len(junk_orbits)
    if not junk_orbits:
        return np.zeros((0, len(reps)), dtype=float)
    C = np.zeros((len(junk_orbits), len(reps)), dtype=float)
    for j, sigma in enumerate(reps):
        for i in range(1, d):
            merged = (sigma[i-1] + sigma[i]) % p
            if merged not in QR_set:
                face = list(sigma)
                face[i-1] = merged
                del face[i]
                ft = zm_orbit_class(tuple(face), QR, p)
                C[junk_orbits[ft], j] += (-1) ** i
    return C

for p in [7, 11]:
    m = (p - 1) // 2
    QR = get_QR(p)
    print(f"{'='*60}")
    print(f"P_{p} (m={m})")
    print(f"{'='*60}\n")

    # Build orbit reps for all degrees
    orbit_reps = {}
    for d in range(2*m + 1):
        orbit_reps[d] = build_orbit_reps(p, d)

    # Compute Omega_orb
    C_orbs = {}
    Omega_orb = {}
    for d in range(2*m + 1):
        C = build_orbit_constraint(p, d, orbit_reps[d], QR)
        C_orbs[d] = C
        Omega_orb[d] = len(orbit_reps[d]) - (np.linalg.matrix_rank(C, tol=1e-8) if C.shape[0] > 0 else 0)

    print(f"Omega_orb: {[Omega_orb[d] for d in range(2*m+1)]}")

    # Build boundary maps
    B_orbs = {}
    for d in range(1, 2*m + 1):
        B_orbs[d] = build_orbit_boundary(p, d, orbit_reps[d], orbit_reps[d-1], QR)

    # Focus on d = m+1 (the critical boundary)
    d_crit = m + 1
    print(f"\n--- Critical boundary ∂_{d_crit}^orb ---")
    B = B_orbs[d_crit]
    C_src = C_orbs[d_crit]
    C_tgt = C_orbs[d_crit - 1]

    print(f"  Source: {len(orbit_reps[d_crit])} orbits, Omega={Omega_orb[d_crit]}")
    print(f"  Target: {len(orbit_reps[d_crit-1])} orbits, Omega={Omega_orb[d_crit-1]}")
    print(f"  B shape: {B.shape}")
    print(f"  rank(B) = {np.linalg.matrix_rank(B, tol=1e-8)}")

    # Stacking rank
    if C_src.shape[0] > 0:
        stacked = np.vstack([C_src, B])
    else:
        stacked = B
    R_crit = np.linalg.matrix_rank(stacked, tol=1e-8) - (np.linalg.matrix_rank(C_src, tol=1e-8) if C_src.shape[0] > 0 else 0)
    print(f"  R_{d_crit}^orb = {R_crit}")

    # Kernel of ∂_{d_crit}
    # The kernel of B restricted to ker(C_src) is what we need.
    # Use SVD of C_src to get Omega basis
    if C_src.shape[0] > 0:
        U, S, Vt = np.linalg.svd(C_src, full_matrices=True)
        null_dim = C_src.shape[1] - np.sum(S > 1e-8)
        Omega_basis = Vt[-null_dim:].T  # columns are Omega basis vectors
    else:
        Omega_basis = np.eye(len(orbit_reps[d_crit]))

    # Boundary restricted to Omega
    B_restricted = B @ Omega_basis  # each column is B applied to an Omega basis vector
    rank_B_restricted = np.linalg.matrix_rank(B_restricted, tol=1e-8)
    print(f"  rank(B|Omega) = {rank_B_restricted}")
    print(f"  ker(B|Omega) dim = {Omega_orb[d_crit] - rank_B_restricted}")

    # For P_7 (m=3): d_crit=4
    # ∂_4^orb: 13 orbits → 7 orbits
    # Omega_4 = 3, Omega_3 = 3
    # R_4 should give β_3 = Omega_3 - R_3 - R_4

    # Actually, let me compute ALL stacking ranks properly
    R_orb = {0: 0}
    for d in range(1, 2*m + 1):
        C = C_orbs[d]
        Bd = B_orbs[d]
        if C.shape[0] > 0:
            st = np.vstack([C, Bd])
        else:
            st = Bd
        R_orb[d] = np.linalg.matrix_rank(st, tol=1e-8) - (np.linalg.matrix_rank(C, tol=1e-8) if C.shape[0] > 0 else 0)
    R_orb[2*m + 1] = 0

    betti_orb = [Omega_orb[d] - R_orb[d] - R_orb[d+1] for d in range(2*m + 1)]
    print(f"\n  R_orb: {[R_orb[d] for d in range(2*m+2)]}")
    print(f"  β_orb: {betti_orb}")

    # Now look at the KERNEL of ∂_{m+1}^orb restricted to Omega
    print(f"\n--- Kernel analysis at d={m+1} ---")

    # Get Omega_{m+1} basis
    C_m1 = C_orbs[m+1]
    if C_m1.shape[0] > 0:
        U, S, Vt = np.linalg.svd(C_m1, full_matrices=True)
        null_dim = C_m1.shape[1] - np.sum(S > 1e-8)
        Om_basis = Vt[-null_dim:].T
    else:
        Om_basis = np.eye(len(orbit_reps[m+1]))

    # ∂_{m+1} restricted to Omega_{m+1}
    B_m1 = B_orbs[m+1]
    B_m1_restr = B_m1 @ Om_basis

    # Kernel of B_{m+1}|Omega
    from scipy.linalg import null_space
    ker = null_space(B_m1_restr)
    print(f"  dim(ker ∂_{m+1}|Omega) = {ker.shape[1]}")
    print(f"  = R_{m+2} + β_{m+1} = {R_orb[m+2]} + {betti_orb[m+1]}")

    # Image of ∂_{m+2} in Omega_{m+1}
    # Get Omega_{m+2} basis
    C_m2 = C_orbs[m+2]
    if C_m2.shape[0] > 0:
        U2, S2, Vt2 = np.linalg.svd(C_m2, full_matrices=True)
        null_dim2 = C_m2.shape[1] - np.sum(S2 > 1e-8)
        Om2_basis = Vt2[-null_dim2:].T
    else:
        Om2_basis = np.eye(len(orbit_reps[m+2]))

    B_m2 = B_orbs[m+2]
    im_m2 = B_m2 @ Om2_basis  # Image of ∂_{m+2} in A_{m+1}

    # Project image onto Omega_{m+1}
    # The image lands in Omega_{m+1} automatically (since ∂ maps Omega to Omega)
    # Express in Omega_{m+1} coordinates
    im_m2_in_Om = np.linalg.lstsq(Om_basis, im_m2, rcond=None)[0]
    rank_im = np.linalg.matrix_rank(im_m2_in_Om, tol=1e-8)
    print(f"  dim(im ∂_{m+2}|Omega) = {rank_im}")
    print(f"  β_{m+1} = ker - im = {ker.shape[1]} - {rank_im} = {ker.shape[1] - rank_im}")

    # ALSO: look at d=m (the other critical degree)
    print(f"\n--- Kernel analysis at d={m} ---")
    C_m = C_orbs[m]
    if C_m.shape[0] > 0:
        Um, Sm, Vtm = np.linalg.svd(C_m, full_matrices=True)
        null_dim_m = C_m.shape[1] - np.sum(Sm > 1e-8)
        Om_m_basis = Vtm[-null_dim_m:].T
    else:
        Om_m_basis = np.eye(len(orbit_reps[m]))

    B_m = B_orbs[m]
    B_m_restr = B_m @ Om_m_basis
    ker_m = null_space(B_m_restr)
    print(f"  dim(ker ∂_{m}|Omega) = {ker_m.shape[1]}")
    print(f"  = R_{m+1} + β_m = {R_orb[m+1]} + {betti_orb[m]}")

    # The kernel of ∂_m at Omega_m:
    # Elements of ker ∂_m are cycles. Those in im ∂_{m+1} are boundaries.
    # β_m = dim(ker) - dim(im)
    B_m1_full = B_orbs[m+1]
    im_m1 = B_m1_full @ Om_basis
    im_m1_in_Om_m = np.linalg.lstsq(Om_m_basis, im_m1, rcond=None)[0]
    rank_im_m1 = np.linalg.matrix_rank(im_m1_in_Om_m, tol=1e-8)
    print(f"  dim(im ∂_{m+1}|Omega) = {rank_im_m1}")
    print(f"  β_m = ker - im = {ker_m.shape[1]} - {rank_im_m1} = {ker_m.shape[1] - rank_im_m1}")

    # For P_7 (m=3): ker(∂_3) should have dim 2 (= R_4 + β_3 = 2 + 0 = 2)
    # im(∂_4) should have dim 2 (= R_4 = 2)
    # So β_3 = 0 ✓

    # For P_11 (m=5): ker(∂_5) should have dim 62 (= R_6 + β_5 = 61 + 1)
    # im(∂_6) should have dim 61 (= R_6 = 61)
    # So β_5 = 1 ✓

    # What is the CYCLE that's NOT a boundary?
    if betti_orb[m] > 0:
        print(f"\n--- The β_m cycle (not a boundary) ---")
        # Project kernel of ∂_m onto complement of im(∂_{m+1})
        # Use QR on image to get basis
        Q_im, _ = np.linalg.qr(im_m1_in_Om_m, mode='reduced')
        # Project kernel vectors
        ker_in_Om = ker_m  # already in Omega coordinates
        proj = Q_im @ (Q_im.T @ ker_in_Om)
        complement = ker_in_Om - proj
        norms = np.linalg.norm(complement, axis=0)
        sig = norms > 1e-6
        cycle_vecs = complement[:, sig]
        print(f"  Found {np.sum(sig)} independent cycle(s)")

        if np.sum(sig) > 0:
            # Express in orbit-rep coordinates
            cycle_in_A = Om_m_basis @ cycle_vecs[:, 0]
            cycle_in_A /= np.max(np.abs(cycle_in_A))
            print(f"  Cycle in orbit coordinates:")
            for j, rep in enumerate(orbit_reps[m]):
                if abs(cycle_in_A[j]) > 1e-6:
                    print(f"    [{rep}]: {cycle_in_A[j]:.6f}")

    print()
