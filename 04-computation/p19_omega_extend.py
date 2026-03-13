#!/usr/bin/env python3
"""
Extend P_19 orbit Omega computation to d=8,9 using sparse modular rank.

Known: Ω_orb through d=7: [1, 1, 8, 60, 417, 2648, 15140, 76474]
Need: Ω_8, Ω_9 (and ideally more) to compute β_9.

Strategy: enumerate orbit reps at d=8, build sparse constraint matrix,
compute rank mod q to get Ω_8.

opus-2026-03-13-S71b
"""

import time

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_orbit_reps_fast(p, d, QR_list):
    """Build orbit reps starting with s_1=1 (canonical under Z_m)."""
    QR_set = set(QR_list)
    m = len(QR_list)
    if d == 0:
        return [()]
    if d == 1:
        return [(1,)]

    results = []
    def backtrack(seq, ps_set, last_ps):
        if len(seq) == d:
            results.append(tuple(seq))
            return
        for s in QR_list:
            new_ps = (last_ps + s) % p
            if new_ps in ps_set or new_ps == 0:
                continue
            seq.append(s)
            ps_set.add(new_ps)
            backtrack(seq, ps_set, new_ps)
            seq.pop()
            ps_set.discard(new_ps)

    backtrack([1], {0, 1}, 1)
    return results  # Don't sort — too slow for large lists

def zm_orbit_class(seq, QR_list, p):
    best = seq
    for q in QR_list:
        if q == 1:
            continue
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best

def build_constraint_cols_sparse(p, d, reps, QR_list, q_mod):
    """Build sparse constraint columns for orbit reps.
    Returns (cols, n_junk) where cols[j] = dict of {row: value mod q}."""
    QR_set = set(QR_list)
    junk_map = {}
    cols = []

    for j, sigma in enumerate(reps):
        col = {}
        for i in range(1, d):  # internal faces
            merged = (sigma[i-1] + sigma[i]) % p
            if merged not in QR_set and merged != 0:
                # Junk face
                face = list(sigma)
                face[i-1] = merged
                del face[i]
                ft = zm_orbit_class(tuple(face), QR_list, p)
                if ft not in junk_map:
                    junk_map[ft] = len(junk_map)
                row = junk_map[ft]
                col[row] = (col.get(row, 0) + (-1)**i) % q_mod
        # Clean zeros
        col = {k: v for k, v in col.items() if v != 0}
        cols.append(col)

        if (j+1) % 50000 == 0:
            print(f"    Built {j+1}/{len(reps)} constraint columns, {len(junk_map)} junk orbits")

    return cols, len(junk_map)

def sparse_rank_mod(cols, q, n_rows):
    """Sparse Gaussian elimination mod q."""
    pivot_col = {}
    reduced = []

    for j, col_orig in enumerate(cols):
        c = dict(col_orig)
        while c:
            r = min(c.keys())
            v = c[r] % q
            if v == 0:
                del c[r]
                continue
            if r in pivot_col:
                pv = reduced[pivot_col[r]]
                pv_val = pv[r]
                inv_pv = pow(pv_val, q - 2, q)
                factor = (v * inv_pv) % q
                for pr, pval in pv.items():
                    c[pr] = (c.get(pr, 0) - factor * pval) % q
                c = {k: v % q for k, v in c.items() if v % q != 0}
            else:
                pivot_col[r] = len(reduced)
                reduced.append(c)
                break

        if (j+1) % 50000 == 0:
            print(f"    Rank computation: {j+1}/{len(cols)} columns processed, rank={len(reduced)}")

    return len(reduced)

# ============================================================
p = 19
m = (p - 1) // 2  # m = 9
QR = get_QR(p)
q_mod = 997

print(f"P_{p} (m={m}), QR = {QR}")
print(f"Working mod {q_mod}\n")

# Known Omega values
omega_known = [1, 1, 8, 60, 417, 2648, 15140, 76474]
print(f"Known Ω_orb: {omega_known}")
print(f"Need: Ω_8, Ω_9 to compute β_9\n")

# Bottom recursion with known values
R = [0, 0]
for d in range(2, 8):
    R.append(omega_known[d-1] - R[-1])
print(f"R_orb (from bottom): {R}")
print(f"R_8 = Ω_7 - R_7 = {omega_known[7]} - {R[7]} = {omega_known[7] - R[7]}")
R.append(omega_known[7] - R[7])
print(f"R = {R}")
print()

# Need Ω_8 to get R_9 = Ω_8 - R_8
# And Ω_9 to get Budget = Ω_9 - R_9

# Compute d=8 orbit reps
print(f"=== Computing d=8 orbit reps ===")
t0 = time.time()
reps_8 = build_orbit_reps_fast(p, 8, QR)
t1 = time.time()
n_reps_8 = len(reps_8)
print(f"  {n_reps_8} orbit reps in {t1-t0:.1f}s")

# Build constraint
print(f"\n  Building constraint matrix...")
t0 = time.time()
C_cols_8, n_junk_8 = build_constraint_cols_sparse(p, 8, reps_8, QR, q_mod)
t1 = time.time()
print(f"  {n_junk_8} junk orbits in {t1-t0:.1f}s")

# Compute rank
print(f"\n  Computing rank...")
t0 = time.time()
rank_C_8 = sparse_rank_mod(C_cols_8, q_mod, n_junk_8)
t1 = time.time()
omega_8 = n_reps_8 - rank_C_8
print(f"  rank(C_8) = {rank_C_8} in {t1-t0:.1f}s")
print(f"  Ω_8^orb = {n_reps_8} - {rank_C_8} = {omega_8}")

R_9 = omega_8 - R[8]
print(f"  R_9 = Ω_8 - R_8 = {omega_8} - {R[8]} = {R_9}")

# Free memory
del reps_8, C_cols_8

# Now d=9
print(f"\n=== Computing d=9 orbit reps ===")
t0 = time.time()
reps_9 = build_orbit_reps_fast(p, 9, QR)
t1 = time.time()
n_reps_9 = len(reps_9)
print(f"  {n_reps_9} orbit reps in {t1-t0:.1f}s")

print(f"\n  Building constraint matrix...")
t0 = time.time()
C_cols_9, n_junk_9 = build_constraint_cols_sparse(p, 9, reps_9, QR, q_mod)
t1 = time.time()
print(f"  {n_junk_9} junk orbits in {t1-t0:.1f}s")

print(f"\n  Computing rank...")
t0 = time.time()
rank_C_9 = sparse_rank_mod(C_cols_9, q_mod, n_junk_9)
t1 = time.time()
omega_9 = n_reps_9 - rank_C_9
print(f"  rank(C_9) = {rank_C_9} in {t1-t0:.1f}s")
print(f"  Ω_9^orb = {n_reps_9} - {rank_C_9} = {omega_9}")

# Now compute Budget and β_9
Budget_9 = omega_9 - R_9
print(f"\n=== RESULT ===")
print(f"Budget_9 = Ω_9 - R_9 = {omega_9} - {R_9} = {Budget_9}")
print(f"β_9^orb = Budget_9 - R_10 (need R_10 from top recursion)")
print(f"Predicted β_9^orb = (m-3)/2 = {(m-3)//2}")
print(f"Predicted β_9 = m(m-3)/2 = {m*(m-3)//2}")

# Complete Omega sequence so far
omega_all = omega_known + [omega_8, omega_9]
print(f"\nΩ_orb = {omega_all}")
print(f"χ partial = {sum((-1)**d * omega_all[d] for d in range(len(omega_all)))}")
