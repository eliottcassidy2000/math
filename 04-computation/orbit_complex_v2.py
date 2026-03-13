#!/usr/bin/env python3
"""
Orbit complex computation for Paley tournaments P_p.
Computes the Z_m orbit complex directly (without full k=0 verification).

opus-2026-03-13-S71b
"""
import numpy as np
import time

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_diffseqs(p, d):
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

def zm_orbit_class(seq, QR_list, p):
    best = seq
    for q in QR_list:
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best

def compute_face(seq, i, p):
    d = len(seq)
    if i == 0:
        return seq[1:]
    elif i == d:
        return seq[:-1]
    else:
        merged = (seq[i-1] + seq[i]) % p
        return seq[:i-1] + (merged,) + seq[i+1:]

def main():
    for p in [7, 11]:
        m = (p - 1) // 2
        QR = get_QR(p)
        print(f"{'='*60}")
        print(f"P_{p} (m={m}), QR = {QR}")
        print(f"{'='*60}\n")

        t0 = time.time()

        # Build diff-seqs and find orbits
        Ad = {}
        orbit_reps = {}
        for d in range(2*m + 1):
            Ad[d] = build_diffseqs(p, d)
            if d == 0:
                orbit_reps[d] = [()]
            else:
                seen = set()
                reps = []
                for seq in Ad[d]:
                    canon = zm_orbit_class(seq, QR, p)
                    if canon not in seen:
                        seen.add(canon)
                        reps.append(canon)
                orbit_reps[d] = sorted(reps)
            print(f"d={d}: |A|={len(Ad[d])}, orbits={len(orbit_reps[d])}", flush=True)

        print(f"\nDiff-seqs built in {time.time()-t0:.1f}s\n")

        # Build orbit constraint matrices (junk faces mapped to orbit classes)
        QR_set = set(QR)
        C_orb = {}
        Omega_orb = {}
        for d in range(2*m + 1):
            reps = orbit_reps[d]
            if d <= 1:
                C_orb[d] = np.zeros((0, len(reps)), dtype=float)
                Omega_orb[d] = len(reps)
                continue

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
                C_orb[d] = np.zeros((0, len(reps)), dtype=float)
                Omega_orb[d] = len(reps)
            else:
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
                C_orb[d] = C
                rank_C = np.linalg.matrix_rank(C, tol=1e-8)
                Omega_orb[d] = len(reps) - rank_C

        print(f"Omega_orb = {[Omega_orb[d] for d in range(2*m+1)]}")
        chi_orb = sum((-1)**d * Omega_orb[d] for d in range(2*m+1))
        print(f"chi_orb = {chi_orb}\n")

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

        # Stacking ranks
        R_orb = {0: 0}
        for d in range(1, 2*m + 1):
            t1 = time.time()
            if C_orb[d].shape[0] > 0:
                stacked = np.vstack([C_orb[d], B_orbs[d]])
            else:
                stacked = B_orbs[d]
            rank_stacked = np.linalg.matrix_rank(stacked, tol=1e-8)
            rank_C = np.linalg.matrix_rank(C_orb[d], tol=1e-8) if C_orb[d].shape[0] > 0 else 0
            R_orb[d] = rank_stacked - rank_C
            dt = time.time() - t1
            if dt > 0.5:
                print(f"  d={d}: R_orb={R_orb[d]} ({dt:.1f}s)", flush=True)
        R_orb[2*m + 1] = 0

        betti_orb = [Omega_orb[d] - R_orb[d] - R_orb[d+1] for d in range(2*m + 1)]
        print(f"\nR_orb = {[R_orb[d] for d in range(2*m + 2)]}")
        print(f"β_orb = {betti_orb}")
        print(f"chi_orb = {sum((-1)**d * b for d, b in enumerate(betti_orb))}")
        print()

        # Analysis of the cycle structure
        print("--- Orbit cycle analysis ---")
        for d in range(2*m + 1):
            if betti_orb[d] > 0:
                ker_dim = Omega_orb[d] - R_orb[d]
                im_dim = R_orb[d+1]
                print(f"  d={d}: β_orb={betti_orb[d]}, ker={ker_dim}, im={im_dim}")

        # For small orbit complexes, print the boundary matrices
        if p == 7:
            print("\n--- P_7 orbit boundary matrices ---")
            for d in range(1, 2*m + 1):
                reps_d = orbit_reps[d]
                reps_dm1 = orbit_reps[d-1]
                print(f"\n  ∂_{d}^orb: {len(reps_d)} → {len(reps_dm1)}")
                if len(reps_d) <= 15 and len(reps_dm1) <= 15:
                    print(f"    cols (d={d} reps): {reps_d}")
                    print(f"    rows (d={d-1} reps): {reps_dm1}")
                    B = B_orbs[d]
                    for i in range(B.shape[0]):
                        row = [int(B[i,j]) if B[i,j] == int(B[i,j]) else B[i,j]
                               for j in range(B.shape[1])]
                        print(f"    {reps_dm1[i]}: {row}")

        # For P_11, look at boundary near d=m
        if p == 11:
            print("\n--- P_11: orbit boundary near d=m ---")
            for d in [m-1, m, m+1, m+2]:
                reps_d = orbit_reps[d]
                reps_dm1 = orbit_reps[d-1]
                B = B_orbs[d]
                # Report shape and rank
                rank_B = np.linalg.matrix_rank(B, tol=1e-8)
                print(f"  ∂_{d}^orb: ({len(reps_dm1)} × {len(reps_d)}), "
                      f"rank(B)={rank_B}")

                # Rank of stacked [C; B]
                if C_orb[d].shape[0] > 0:
                    stacked = np.vstack([C_orb[d], B])
                    rank_s = np.linalg.matrix_rank(stacked, tol=1e-8)
                    rank_c = np.linalg.matrix_rank(C_orb[d], tol=1e-8)
                    print(f"    rank([C;B])={rank_s}, rank(C)={rank_c}, "
                          f"R_orb={rank_s - rank_c}")
                else:
                    print(f"    R_orb = {rank_B} (no constraints)")

        print(f"\n{'='*60}\n")


if __name__ == "__main__":
    main()
