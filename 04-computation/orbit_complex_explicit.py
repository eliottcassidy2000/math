#!/usr/bin/env python3
"""
Explicit orbit complex computation for Paley tournaments P_p.

The Z_m = QR* group acts freely on diff-seqs by component-wise multiplication.
The orbit complex has dimensions Omega_d / m and its Betti numbers satisfy
β_m^{orb} = (m-3)/2.

This script:
1. Computes diff-seqs and their Z_m orbits
2. Builds the orbit boundary maps explicitly
3. Computes orbit Betti numbers
4. Analyzes the structure of the orbit cycles

opus-2026-03-13-S71b
"""
import numpy as np
from math import gcd

def get_QR(p):
    """Quadratic residues mod p."""
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_diffseqs(p, d):
    """Build all allowed diff-sequences of length d for P_p."""
    QR = set(pow(x, 2, p) for x in range(1, p))
    QR_list = sorted(QR)
    results = []
    if d == 0:
        return [()]
    def backtrack(seq, ps_list, ps_set):
        if len(seq) == d:
            results.append(tuple(seq))
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
    backtrack([], [0], {0})
    return results

def build_omega_basis_indices(p, d, Ad):
    """Get indices of diff-seqs that form Omega_d (kernel of constraint map)."""
    QR = set(pow(x, 2, p) for x in range(1, p))
    # Build constraint matrix
    junk_faces = {}
    for seq in Ad:
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                if ft not in junk_faces:
                    junk_faces[ft] = len(junk_faces)
    if not junk_faces:
        return list(range(len(Ad)))  # All are in Omega

    n_rows = len(junk_faces)
    C = np.zeros((n_rows, len(Ad)), dtype=float)
    for j, seq in enumerate(Ad):
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                C[junk_faces[ft], j] += (-1) ** i

    # Omega_d = ker(C). For orbit computation, we need Omega-basis in A-coords.
    # Use SVD to find kernel
    if C.shape[0] == 0:
        return np.eye(len(Ad))
    U, S, Vt = np.linalg.svd(C, full_matrices=True)
    tol = 1e-8
    null_dim = C.shape[1] - np.sum(S > tol)
    if null_dim == 0:
        return np.zeros((len(Ad), 0))
    return Vt[-null_dim:].T  # shape (|A_d|, null_dim)


def zm_orbit_class(seq, QR_list, p):
    """Return canonical representative of Z_m orbit of seq."""
    m = len(QR_list)
    best = seq
    for q in QR_list:
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best


def compute_face(seq, i, p):
    """Compute face_i of diff-seq."""
    d = len(seq)
    if i == 0:
        return seq[1:]
    elif i == d:
        return seq[:-1]
    else:
        merged = (seq[i-1] + seq[i]) % p
        return seq[:i-1] + (merged,) + seq[i+1:]


def build_orbit_boundary(p, d, orbit_reps_d, orbit_reps_dm1, QR_list, Ad, Ad_prev):
    """
    Build the orbit boundary map ∂^{orb}_d.

    For each orbit representative σ in degree d, compute ∂(σ) as a chain
    in Omega_{d-1}, then express each face in terms of orbit representatives.

    The face might not be in Omega_{d-1} (it could be a "junk" face).
    But for Omega-chains, the boundary lands in Omega_{d-1}.

    Actually: we work with the FULL diff-seq complex first, then quotient.
    The boundary of a diff-seq σ is Σ (-1)^i face_i(σ).
    Some face_i are in A_{d-1}, some are not (junk faces).
    Omega_d = { chains whose boundary has no junk faces }.
    For an Omega-chain, ∂ lands in Omega_{d-1}.

    In the orbit complex, we choose orbit reps and compute:
    ∂^{orb}([σ]) = Σ_i (-1)^i [face_i(σ)]
    where [·] maps to orbit class (with possible sign/coefficient from Z_m scaling).

    Wait — for a free group action on a chain complex, the orbit boundary is:
    ∂^{orb}([σ]) = [∂(σ)]

    But ∂(σ) is a chain in Omega_{d-1}, and we need to express it in terms of
    orbit representatives. If face_i(σ) = g · τ for orbit rep τ and g ∈ Z_m,
    then [face_i(σ)] = [τ] (since the orbit class is the same).

    So ∂^{orb}([σ]) = Σ_i (-1)^i · [face_i(σ)]
    where [face_i(σ)] is the orbit class of face_i(σ).

    But WAIT: the Omega complex is a SUBCOMPLEX of the A-complex.
    The orbit complex should be computed on Omega, not on A.
    Since Omega has a basis that's not just individual diff-seqs but LINEAR COMBINATIONS,
    the orbit complex is more subtle.

    Actually, for the k=0 eigenspace, the boundary map on Omega uses the STANDARD boundary
    (no phase). And the Z_m action preserves Omega (since it permutes diff-seqs within orbits,
    and the constraint equations are Z_m-invariant).

    For the orbit complex, we need to:
    1. Choose a Z_m-equivariant basis for Omega (which exists since the action is free)
    2. The orbit boundary map is the restriction to Z_m-invariant chains

    But a simpler approach: compute the k=0 boundary ranks directly using the
    stacking trick, with the boundary matrix restricted to Z_m-invariant subspace.

    Actually simplest: just compute rank(∂^{(0)}_d | Omega_d) / m to get R_d^{orb}.
    Since the action is free, all ranks divide by m.
    """
    QR_set = set(QR_list)
    Ad_prev_idx = {seq: i for i, seq in enumerate(Ad_prev)}
    orb_dm1_idx = {}
    for i, rep in enumerate(orbit_reps_dm1):
        orb_dm1_idx[rep] = i

    n_rows = len(orbit_reps_dm1)
    n_cols = len(orbit_reps_d)
    B = np.zeros((n_rows, n_cols), dtype=float)

    for j, sigma in enumerate(orbit_reps_d):
        for i in range(d + 1):
            face = compute_face(sigma, i, p)
            if face in Ad_prev_idx:
                # Find orbit class
                face_orb = zm_orbit_class(face, QR_list, p)
                if face_orb in orb_dm1_idx:
                    B[orb_dm1_idx[face_orb], j] += (-1) ** i

    return B


def main():
    for p in [7, 11]:
        m = (p - 1) // 2
        QR = get_QR(p)
        print(f"{'='*60}")
        print(f"P_{p} (m={m}), QR = {QR}")
        print(f"{'='*60}\n")

        # Build diff-seqs
        Ad = {}
        for d in range(2*m + 1):
            Ad[d] = build_diffseqs(p, d)

        # Compute Omega dims (via constraint matrix rank)
        Omega = {}
        for d in range(2*m + 1):
            if d <= 1:
                Omega[d] = len(Ad[d])
            else:
                QR_set = set(QR)
                junk_faces = {}
                for seq in Ad[d]:
                    for i in range(1, d):
                        merged = (seq[i-1] + seq[i]) % p
                        if merged not in QR_set:
                            face = list(seq)
                            face[i-1] = merged
                            del face[i]
                            ft = tuple(face)
                            if ft not in junk_faces:
                                junk_faces[ft] = len(junk_faces)
                if not junk_faces:
                    Omega[d] = len(Ad[d])
                else:
                    C = np.zeros((len(junk_faces), len(Ad[d])), dtype=float)
                    for j, seq in enumerate(Ad[d]):
                        for i in range(1, d):
                            merged = (seq[i-1] + seq[i]) % p
                            if merged not in QR_set:
                                face = list(seq)
                                face[i-1] = merged
                                del face[i]
                                ft = tuple(face)
                                C[junk_faces[ft], j] += (-1) ** i
                    rank_C = np.linalg.matrix_rank(C, tol=1e-8)
                    Omega[d] = len(Ad[d]) - rank_C

        print(f"Omega = {[Omega[d] for d in range(2*m+1)]}")
        orb_dims = [1] + [Omega[d] // m for d in range(1, 2*m+1)]
        print(f"Orbit dims = {orb_dims}")
        print()

        # Compute Z_m orbits of diff-seqs
        orbit_reps = {}  # d -> list of canonical orbit representatives
        for d in range(2*m + 1):
            if d == 0:
                orbit_reps[d] = [()]
                continue
            seen = set()
            reps = []
            for seq in Ad[d]:
                canon = zm_orbit_class(seq, QR, p)
                if canon not in seen:
                    seen.add(canon)
                    reps.append(canon)
            orbit_reps[d] = sorted(reps)
            print(f"d={d}: {len(Ad[d])} diff-seqs, {len(reps)} orbits "
                  f"(Omega/m = {Omega[d]//m if d > 0 else 1})")
        print()

        # Now: the orbit complex approach
        # Instead of working with orbits of DIFF-SEQS, we should use the stacking trick
        # to compute the k=0 boundary ranks, then divide by m.
        # This gives R_d^{orb} = R_d^{(0)} / m.

        print("--- k=0 boundary ranks (stacking trick) ---")
        QR_set = set(QR)

        # Build constraint matrices
        C_mats = {}
        for d in range(2*m + 1):
            junk_faces = {}
            for seq in Ad[d]:
                for i in range(1, d):
                    merged = (seq[i-1] + seq[i]) % p
                    if merged not in QR_set:
                        face = list(seq)
                        face[i-1] = merged
                        del face[i]
                        ft = tuple(face)
                        if ft not in junk_faces:
                            junk_faces[ft] = len(junk_faces)
            if not junk_faces:
                C_mats[d] = np.zeros((0, len(Ad[d])), dtype=float)
            else:
                C = np.zeros((len(junk_faces), len(Ad[d])), dtype=float)
                for j, seq in enumerate(Ad[d]):
                    for i in range(1, d):
                        merged = (seq[i-1] + seq[i]) % p
                        if merged not in QR_set:
                            face = list(seq)
                            face[i-1] = merged
                            del face[i]
                            ft = tuple(face)
                            C[junk_faces[ft], j] += (-1) ** i
                C_mats[d] = C

        # Build k=0 boundary matrices and compute ranks via stacking
        # Stacking trick: rank(∂_d | Omega_d) = rank([C_d; B_d]) - rank(C_d)
        # where C_d constrains A_d (cols = A_d) and B_d maps A_d → A_{d-1} (cols = A_d)
        R = {0: 0}
        for d in range(1, 2*m + 1):
            Ad_from = Ad[d]
            Ad_to = Ad[d-1]
            to_idx = {seq: i for i, seq in enumerate(Ad_to)}
            B = np.zeros((len(Ad_to), len(Ad_from)), dtype=float)
            for j, seq in enumerate(Ad_from):
                for i in range(d + 1):
                    face = compute_face(seq, i, p)
                    if face in to_idx:
                        B[to_idx[face], j] += (-1)**i
            # Stacking rank: C_d and B both have |A_d| columns
            if C_mats[d].shape[0] > 0:
                stacked = np.vstack([C_mats[d], B])
            else:
                stacked = B
            rank_stacked = np.linalg.matrix_rank(stacked, tol=1e-8)
            rank_C = np.linalg.matrix_rank(C_mats[d], tol=1e-8) if C_mats[d].shape[0] > 0 else 0
            R[d] = rank_stacked - rank_C
        R[2*m + 1] = 0

        betti = [Omega[d] - R[d] - R[d+1] for d in range(2*m + 1)]
        print(f"  R^(0) = {[R[d] for d in range(2*m + 2)]}")
        print(f"  β^(0) = {betti}")
        print(f"  chi = {sum((-1)**d * b for d, b in enumerate(betti))}")

        # Orbit ranks and Betti
        R_orb = {d: R[d] // m if d > 0 else 0 for d in range(2*m + 2)}
        # Check divisibility
        for d in range(1, 2*m + 1):
            if R[d] % m != 0:
                print(f"  WARNING: R_{d} = {R[d]} not divisible by m={m}!")

        betti_orb = [orb_dims[d] - R_orb[d] - R_orb.get(d+1, 0) for d in range(2*m + 1)]
        print(f"\n  R^orb = {[R_orb[d] for d in range(2*m + 2)]}")
        print(f"  β^orb = {betti_orb}")
        print(f"  chi^orb = {sum((-1)**d * b for d, b in enumerate(betti_orb))}")
        print()

        # Detailed analysis: look at the orbit boundary map structure
        # For the orbit complex, use the Z_m-averaged boundary map
        print("--- Orbit boundary map (Z_m averaged) ---")
        # The Z_m-invariant boundary map acts on the space of Z_m-invariant chains.
        # A Z_m-invariant chain assigns the SAME coefficient to all elements of an orbit.
        #
        # If we choose orbit representatives σ_1, ..., σ_N, then a Z_m-invariant chain is
        # c = Σ_i a_i · (Σ_{g ∈ Z_m} g · σ_i)
        #
        # The boundary maps this to:
        # ∂c = Σ_i a_i · (Σ_{g ∈ Z_m} g · ∂(σ_i))
        #    = Σ_i a_i · (Σ_{g ∈ Z_m} ∂(g · σ_i))  [since ∂ is Z_m-equivariant]
        #
        # Now ∂(σ_i) = Σ_k (-1)^k face_k(σ_i). Each face_k(σ_i) lies in some Z_m orbit.
        # If face_k(σ_i) = h_k · τ_{j_k} for some orbit rep τ_{j_k} and h_k ∈ Z_m, then:
        #   Σ_{g ∈ Z_m} g · face_k(σ_i) = Σ_{g ∈ Z_m} g · h_k · τ_{j_k} = Σ_{g ∈ Z_m} g · τ_{j_k}
        #
        # So the orbit boundary map is:
        #   ∂^{orb}([σ_i]) = Σ_k (-1)^k [face_k(σ_i)]
        # where [·] denotes the orbit class.
        #
        # But we need to be careful: some faces might not be in A_{d-1} (junk faces).
        # For Omega-chains, the junk faces cancel. In the orbit, we need the ORBIT of
        # the junk faces to also cancel.
        #
        # Since the constraint map is also Z_m-equivariant, and Omega = ker(C),
        # the Z_m-invariant part of Omega is ker(C^{orb}) where C^{orb} is the
        # orbit constraint map.

        # Let me just directly compute the orbit boundary by:
        # 1. Taking orbit reps of A_d
        # 2. Computing boundary in A_{d-1} coordinates
        # 3. Converting to orbit coordinates (summing within each orbit)

        # This gives the orbit boundary as a matrix on A-orbits, with rows indexed by A_{d-1}-orbits.
        # Then the orbit complex Omega^{orb} = ker(C^{orb}) and ∂^{orb} acts on it.

        for d in range(1, min(2*m + 1, 7)):  # limit to manageable degrees
            reps_d = orbit_reps[d]
            reps_dm1 = orbit_reps[d-1]
            reps_dm1_idx = {rep: i for i, rep in enumerate(reps_dm1)}

            B_orb = np.zeros((len(reps_dm1), len(reps_d)), dtype=float)

            for j, sigma in enumerate(reps_d):
                for fi in range(d + 1):
                    face = compute_face(sigma, fi, p)
                    # Check if face is in A_{d-1}
                    face_canon = zm_orbit_class(face, QR, p)
                    if face_canon in reps_dm1_idx:
                        B_orb[reps_dm1_idx[face_canon], j] += (-1) ** fi

            rank_B = np.linalg.matrix_rank(B_orb, tol=1e-8)
            print(f"  d={d}: B_orb shape {B_orb.shape}, rank {rank_B}")

            if d <= 3 and len(reps_d) <= 20:
                print(f"    Orbit reps: {reps_d}")
                print(f"    B_orb:")
                for i, rep in enumerate(reps_dm1):
                    row = B_orb[i]
                    nonzero = [(j, v) for j, v in enumerate(row) if abs(v) > 1e-8]
                    if nonzero:
                        print(f"      [{rep}]: {nonzero}")

        print()

        # Now compute orbit Betti directly using orbit stacking trick
        print("--- Orbit Betti via stacking trick ---")

        # Build orbit constraint matrices
        C_orb = {}
        for d in range(2*m + 1):
            if d <= 1:
                C_orb[d] = np.zeros((0, len(orbit_reps[d])), dtype=float)
                continue

            reps = orbit_reps[d]
            reps_idx = {rep: i for i, rep in enumerate(reps)}
            QR_set_local = set(QR)

            junk_orbits = {}
            for sigma in reps:
                for i in range(1, d):
                    merged = (sigma[i-1] + sigma[i]) % p
                    if merged not in QR_set_local:
                        face = list(sigma)
                        face[i-1] = merged
                        del face[i]
                        ft = zm_orbit_class(tuple(face), QR, p)
                        if ft not in junk_orbits:
                            junk_orbits[ft] = len(junk_orbits)

            if not junk_orbits:
                C_orb[d] = np.zeros((0, len(reps)), dtype=float)
            else:
                C = np.zeros((len(junk_orbits), len(reps)), dtype=float)
                for j, sigma in enumerate(reps):
                    for i in range(1, d):
                        merged = (sigma[i-1] + sigma[i]) % p
                        if merged not in QR_set_local:
                            face = list(sigma)
                            face[i-1] = merged
                            del face[i]
                            ft = zm_orbit_class(tuple(face), QR, p)
                            C[junk_orbits[ft], j] += (-1) ** i
                C_orb[d] = C

        # Orbit Omega dims
        Omega_orb = {}
        for d in range(2*m + 1):
            if C_orb[d].shape[0] == 0:
                Omega_orb[d] = len(orbit_reps[d])
            else:
                Omega_orb[d] = len(orbit_reps[d]) - np.linalg.matrix_rank(C_orb[d], tol=1e-8)

        print(f"  Omega_orb = {[Omega_orb[d] for d in range(2*m+1)]}")
        print(f"  Expected  = {orb_dims}")

        # Build orbit boundary matrices
        B_orbs = {}
        for d in range(1, 2*m + 1):
            reps_d = orbit_reps[d]
            reps_dm1 = orbit_reps[d-1]
            reps_dm1_idx = {rep: i for i, rep in enumerate(reps_dm1)}

            B = np.zeros((len(reps_dm1), len(reps_d)), dtype=float)
            for j, sigma in enumerate(reps_d):
                for fi in range(d + 1):
                    face = compute_face(sigma, fi, p)
                    face_canon = zm_orbit_class(face, QR, p)
                    if face_canon in reps_dm1_idx:
                        B[reps_dm1_idx[face_canon], j] += (-1) ** fi
            B_orbs[d] = B

        # Stacking ranks for orbit complex
        # Same trick: rank([C_orb_d; B_orb_d]) - rank(C_orb_d), both have cols = |orbit_reps_d|
        R_orb2 = {0: 0}
        for d in range(1, 2*m + 1):
            if C_orb[d].shape[0] > 0:
                stacked = np.vstack([C_orb[d], B_orbs[d]])
            else:
                stacked = B_orbs[d]
            rank_stacked = np.linalg.matrix_rank(stacked, tol=1e-8)
            rank_C = np.linalg.matrix_rank(C_orb[d], tol=1e-8) if C_orb[d].shape[0] > 0 else 0
            R_orb2[d] = rank_stacked - rank_C
        R_orb2[2*m + 1] = 0

        betti_orb2 = [Omega_orb[d] - R_orb2[d] - R_orb2[d+1] for d in range(2*m + 1)]
        print(f"  R_orb = {[R_orb2[d] for d in range(2*m + 2)]}")
        print(f"  β_orb = {betti_orb2}")
        print(f"  chi_orb = {sum((-1)**d * b for d, b in enumerate(betti_orb2))}")

        # Verify: R_orb * m should equal R^{(0)}
        print(f"\n  Verification: R_orb * m vs R^(0):")
        match = True
        for d in range(1, 2*m + 1):
            if R_orb2[d] * m != R[d]:
                print(f"    d={d}: R_orb*m={R_orb2[d]*m} vs R^(0)={R[d]} MISMATCH")
                match = False
        if match:
            print(f"    ALL MATCH ✓")

        print()

        # KEY ANALYSIS: For P_11, look at the orbit cycle at d=m
        if p == 11:
            print("--- P_11 orbit cycle analysis ---")
            d = m  # d=5
            print(f"\nd={d}: Omega_orb={Omega_orb[d]}, R_{d}^orb={R_orb2[d]}, "
                  f"R_{d+1}^orb={R_orb2[d+1]}, β_{d}^orb={betti_orb2[d]}")

            # The cycle lives in ker(∂_d^orb) but not im(∂_{d+1}^orb)
            # ker(∂_d^orb) has dim = Omega_orb[d] - R_orb2[d]
            print(f"  ker(∂_{d}^orb) dim = {Omega_orb[d] - R_orb2[d]}")
            print(f"  im(∂_{d+1}^orb) dim = {R_orb2[d+1]}")

            # Similarly for d=m+1
            d = m + 1
            print(f"\nd={d}: Omega_orb={Omega_orb[d]}, R_{d}^orb={R_orb2[d]}, "
                  f"R_{d+1}^orb={R_orb2[d+1]}, β_{d}^orb={betti_orb2[d]}")
            print(f"  ker(∂_{d}^orb) dim = {Omega_orb[d] - R_orb2[d]}")
            print(f"  im(∂_{d+1}^orb) dim = {R_orb2[d+1]}")

        print(f"\n{'='*60}\n")


if __name__ == "__main__":
    main()
