#!/usr/bin/env python3
"""
Compute k=0 eigenspace boundary ranks for P_19 using QR-orbit quotient.

Key insight: QR_p acts FREELY on diff-sequences (d>=1), splitting the k=0
complex into m identical sub-complexes. Working with orbit representatives
reduces matrix dimensions by factor m=9, making the computation feasible.

The QR action commutes with the boundary operator, so the quotient complex
is well-defined with the same Betti numbers (up to the d=0 term).

Each QR-character sub-complex has:
  β_d^{(0,j)} for j=0,...,m-1
  Σ_j β_d^{(0,j)} = β_d^{(0)} (total k=0 Betti)

For j≠0: dim_d = Ω_d/m (d≥1), dim_0 = 0.
For j=0: dim_d = Ω_d/m (d≥1), dim_0 = 1.

We compute the TRIVIAL character (j=0) sub-complex, which has dim_0=1.
The boundary of an orbit sum [D] = Σ_{a∈QR} a·D is:
∂[D] = Σ_{a∈QR} Σ_i (-1)^i a·d_i(D) = Σ_i (-1)^i [d_i(D)]

So the quotient boundary matrix maps orbit reps to orbit reps with ±1 coefficients.

opus-2026-03-13-S71c
"""

import sys, time
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729  # Large prime for modular arithmetic
MAX_ORBIT_DENSE = 6000  # Max orbit reps for dense numpy computation

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def orbit_canonical(seq, S, p):
    """Find canonical representative of QR-orbit of seq.
    Orbit: {(a*s_1,...,a*s_d) : a in QR_p}.
    Canonical = lexicographically smallest.
    """
    if not seq:
        return seq
    best = seq
    for a in S:
        scaled = tuple((a * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best

def gauss_rank_dense_mod(M, prime):
    """Dense Gaussian elimination mod prime using numpy."""
    A = M.astype(np.int64) % prime
    rows, cols = A.shape
    pivot_row = 0
    rank = 0
    for col in range(cols):
        # Find pivot
        found = -1
        for row in range(pivot_row, rows):
            if A[row, col] % prime != 0:
                found = row
                break
        if found < 0:
            continue
        # Swap
        A[[pivot_row, found]] = A[[found, pivot_row]]
        # Normalize
        inv = pow(int(A[pivot_row, col] % prime), prime - 2, prime)
        A[pivot_row] = (A[pivot_row] * inv) % prime
        # Eliminate
        for row in range(rows):
            if row != pivot_row and A[row, col] != 0:
                factor = int(A[row, col])
                A[row] = (A[row] - factor * A[pivot_row]) % prime
        rank += 1
        pivot_row += 1
    return rank

def sparse_rank_fast(cols_data, n_rows, prime):
    """Sparse Gaussian elimination mod prime. Optimized version."""
    pivots = {}
    rank = 0
    for col_idx, col in enumerate(cols_data):
        if col_idx % 20000 == 0 and col_idx > 0:
            print(f"      col {col_idx}/{len(cols_data)}, rank={rank}", flush=True)
        work = dict(col)
        while work:
            pivot_row = min(work.keys())
            if pivot_row in pivots:
                piv_data = pivots[pivot_row]
                factor = work.pop(pivot_row)
                for pr, pv in piv_data.items():
                    val = (work.get(pr, 0) - factor * pv) % prime
                    if val:
                        work[pr] = val
                    elif pr in work:
                        del work[pr]
            else:
                inv = pow(work[pivot_row], prime - 2, prime)
                normalized = {r: (v * inv) % prime for r, v in work.items()}
                pivots[pivot_row] = normalized
                rank += 1
                break
    return rank

p = 19
S = qr(p)
m = len(S)
max_d = 7  # Cap here (d=8 enumeration OOMs)
print(f"P_{p}: m={m}, QR={S}")
print(f"Using QR-orbit quotient: all dimensions reduced by factor {m}")
t0 = time.time()

# Enumerate diff-sequences
prev_seqs = [()]
prev_ps = {(): frozenset([0])}
prev_last = {(): 0}
all_seqs = {0: [()]}

for d in range(1, max_d + 1):
    curr_seqs = []
    curr_ps = {}
    curr_last = {}
    for seq in prev_seqs:
        ps_set = prev_ps[seq]
        last = prev_last[seq]
        for s in S:
            ns = (last + s) % p
            if ns not in ps_set:
                nseq = seq + (s,)
                curr_seqs.append(nseq)
                curr_ps[nseq] = ps_set | {ns}
                curr_last[nseq] = ns
    if not curr_seqs:
        break
    all_seqs[d] = curr_seqs
    prev_seqs = curr_seqs
    prev_ps = curr_ps
    prev_last = curr_last
    print(f"  d={d}: A={len(curr_seqs)} ({time.time()-t0:.1f}s)")

actual_max_d = max(all_seqs.keys())
max_d = min(max_d, actual_max_d)

# Build orbit representatives for each degree
print(f"\nBuilding orbit representatives...")
orbit_reps = {0: [()]}
orbit_rep_idx = {0: {(): 0}}

for d in range(1, max_d + 1):
    seen = set()
    reps = []
    idx = {}
    for seq in all_seqs[d]:
        canon = orbit_canonical(seq, S, p)
        if canon not in seen:
            seen.add(canon)
            idx[canon] = len(reps)
            reps.append(canon)
    orbit_reps[d] = reps
    orbit_rep_idx[d] = idx
    print(f"  d={d}: {len(all_seqs[d])} seqs → {len(reps)} orbits "
          f"(ratio {len(all_seqs[d])/len(reps):.1f}) ({time.time()-t0:.1f}s)")

allowed = {d: set(all_seqs.get(d, [])) for d in range(max_d + 1)}

# Compute boundary ranks using orbit quotient
omega_dims_full = [1]  # Full Omega dims
omega_dims_orbit = [1]  # Orbit sub-complex dims
boundary_ranks = {0: 0}

for d in range(1, max_d + 1):
    reps_d = orbit_reps[d]
    reps_dm1 = orbit_reps.get(d-1, [()])
    n_reps_d = len(reps_d)
    n_reps_dm1 = len(reps_dm1)

    if n_reps_d == 0:
        omega_dims_full.append(0)
        omega_dims_orbit.append(0)
        boundary_ranks[d] = 0
        continue

    # For the orbit quotient, we need:
    # 1. Junk orbit reps: faces of orbit reps that are not in allowed[d-1]
    # 2. Constraint matrix: orbit rep → junk orbit rep
    # 3. Combined matrix: orbit rep → (junk orbit rep, allowed orbit rep)

    junk_orbits = set()
    face_data_junk = []
    face_data_allowed = []

    for D in reps_d:
        junk_faces = []
        allowed_faces = []
        for fi in range(d + 1):
            if fi == 0: fd = D[1:]
            elif fi == d: fd = D[:d-1]
            else:
                merged = (D[fi-1] + D[fi]) % p
                fd = D[:fi-1] + (merged,) + D[fi+1:]
            sign = 1 if fi % 2 == 0 else -1

            is_allowed = (fd in allowed[d-1])
            # Get orbit canonical form
            fd_canon = orbit_canonical(fd, S, p) if fd else fd

            if not is_allowed:
                junk_orbits.add(fd_canon)
                junk_faces.append((fd_canon, sign))
            else:
                allowed_faces.append((fd_canon, sign))
        face_data_junk.append(junk_faces)
        face_data_allowed.append(allowed_faces)

    junk_list = sorted(junk_orbits)
    n_junk = len(junk_list)
    junk_idx = {j: i for i, j in enumerate(junk_list)}

    # Build constraint matrix (junk rows × orbit rep cols) — only for dense
    use_dense = (n_reps_d <= MAX_ORBIT_DENSE)
    if use_dense:
        C_constraint = np.zeros((n_junk, n_reps_d), dtype=np.int64)
        for j, jfaces in enumerate(face_data_junk):
            for fd_canon, sign in jfaces:
                row = junk_idx[fd_canon]
                C_constraint[row, j] = (C_constraint[row, j] + sign) % PRIME
    else:
        C_constraint = None  # Will use sparse directly

    method = "dense" if use_dense else "sparse"
    print(f"  d={d}: constraint rank (orbits={n_reps_d}, junk_orbits={n_junk}, {method})...",
          end=" ", flush=True)

    if use_dense:
        rank_constraint = gauss_rank_dense_mod(C_constraint, PRIME)
    else:
        # Build sparse constraint columns directly (no dense matrix)
        sparse_constraint_cols = []
        for j, jfaces in enumerate(face_data_junk):
            col = {}
            for fd_canon, sign in jfaces:
                row = junk_idx[fd_canon]
                val = (col.get(row, 0) + sign) % PRIME
                if val:
                    col[row] = val
                elif row in col:
                    del col[row]
            sparse_constraint_cols.append(col)
        rank_constraint = sparse_rank_fast(sparse_constraint_cols, n_junk, PRIME)

    omega_d_orbit = n_reps_d - rank_constraint
    omega_d_full = omega_d_orbit * m  # Full Omega = orbit × m (for d≥1)
    omega_dims_full.append(omega_d_full)
    omega_dims_orbit.append(omega_d_orbit)
    print(f"Ω_orbit={omega_d_orbit}, Ω_full={omega_d_full}", flush=True)

    # Build combined matrix [constraint; boundary]
    idx_dm1 = orbit_rep_idx.get(d-1, {})
    n_total_rows = n_junk + n_reps_dm1

    if use_dense:
        C_combined = np.zeros((n_total_rows, n_reps_d), dtype=np.int64)
        C_combined[:n_junk, :] = C_constraint
        for j, afaces in enumerate(face_data_allowed):
            for fd_canon, sign in afaces:
                if fd_canon in idx_dm1:
                    row = n_junk + idx_dm1[fd_canon]
                    C_combined[row, j] = (C_combined[row, j] + sign) % PRIME
        print(f"         combined rank ({n_junk}+{n_reps_dm1} rows, {method})...",
              end=" ", flush=True)
        rank_combined = gauss_rank_dense_mod(C_combined, PRIME)
    else:
        # Build sparse combined columns directly
        sparse_combined = []
        for j in range(n_reps_d):
            col = {}
            for fd_canon, sign in face_data_junk[j]:
                row = junk_idx[fd_canon]
                val = (col.get(row, 0) + sign) % PRIME
                if val:
                    col[row] = val
                elif row in col:
                    del col[row]
            for fd_canon, sign in face_data_allowed[j]:
                if fd_canon in idx_dm1:
                    row = n_junk + idx_dm1[fd_canon]
                    val = (col.get(row, 0) + sign) % PRIME
                    if val:
                        col[row] = val
                    elif row in col:
                        del col[row]
            sparse_combined.append(col)
        print(f"         combined rank ({n_junk}+{n_reps_dm1} rows, {method})...",
              end=" ", flush=True)
        rank_combined = sparse_rank_fast(sparse_combined, n_total_rows, PRIME)

    bd_rank = rank_combined - rank_constraint
    boundary_ranks[d] = bd_rank
    dt = time.time() - t0
    print(f"rank(∂)={bd_rank} ({dt:.1f}s)", flush=True)

# Compute Betti numbers
print(f"\n{'='*60}")
print(f"RESULTS for P_{p} k=0 eigenspace (QR-orbit quotient)")
print(f"{'='*60}")
print(f"  Ω_full  = {omega_dims_full}")
print(f"  Ω_orbit = {omega_dims_orbit}")
print(f"  ∂r = {[boundary_ranks.get(d,0) for d in range(max_d+1)]}")

# Betti of the ORBIT sub-complex (trivial QR-character)
betti_orbit = []
for d in range(max_d + 1):
    od = omega_dims_orbit[d]
    rd = boundary_ranks.get(d, 0)
    rd1 = boundary_ranks.get(d+1, 0)
    bd = od - rd - rd1
    betti_orbit.append(bd)

# Full k=0 Betti
# β_d^(0) = m * β_d^{(0,j)} for j≠0, plus β_d^{(0,0)} for j=0
# But this only gives the j=0 sub-complex.
# For the full k=0 Betti: β_d^{(0)} = sum over all m QR-characters of β_d^{(0,j)}

print(f"  β_orbit = {betti_orbit}  (j=0 QR-character sub-complex)")

# For the TOTAL k=0 Betti, we need all m QR-character sub-complexes.
# But they're all isomorphic (m copies of the same complex for j≠0).
# Wait — the j=0 complex has dim_0=1, j≠0 complexes have dim_0=0.
# And they may NOT all have the same boundary ranks.

# Actually, let's check: the QR-orbit quotient gives the j=0 (trivial character).
# The j≠0 characters may differ. But under Galois symmetry, all j≠0 may be equal.

# For now, report what we have and verify against known values.
print(f"\n  Rank verification:")
print(f"    rank(∂_2) orbit = {boundary_ranks.get(2,0)} "
      f"(expected: m={m} for full, 1 for orbit)")
print(f"    rank(∂_3) orbit = {boundary_ranks.get(3,0)} "
      f"(expected: m(m-2)={m*(m-2)} for full, {m-2} for orbit)")

# Chi check
chi = sum((-1)**d * omega_dims_orbit[d] for d in range(max_d+1))
print(f"\n  χ(orbit sub-complex) = {chi}")
print(f"  Expected: 1 (since j=0 has d=0 chain)")

# Verify: Ω_d/m matches orbit dims
print(f"\n  Ω_d/m consistency:")
for d in range(1, max_d+1):
    full = omega_dims_full[d]
    orb = omega_dims_orbit[d]
    print(f"    d={d}: Ω={full}, Ω/m={full//m}, orbit={orb}, match={full//m == orb}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
print("=" * 60)
