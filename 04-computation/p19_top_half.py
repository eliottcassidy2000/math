#!/usr/bin/env python3
"""
Compute TOP HALF of P_19 orbit Omega sequence: d=10,...,18.

Key insight: for d > m = 9, the orbit rep counts SHRINK rapidly
(constraints become very tight), so these are MUCH easier than Ω_9.

From the top recursion, R_10 gives β_9 = Ω_9 - R_9 - R_10.
But we don't need Ω_9! We can use the Euler characteristic:
  χ = p = Σ (-1)^d Ω_d
So Ω_9 = p - Σ_{d≠9} (-1)^d Ω_d, and then β_9 = Ω_9 - R_9 - R_10.

opus-2026-03-13-S71b
"""

import time

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_orbit_reps_fast(p, d, QR_list):
    """Build orbit reps starting with s_1=1 (canonical under Z_m)."""
    QR_set = set(QR_list)
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
    return results

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
    """Build sparse constraint columns for orbit reps."""
    QR_set = set(QR_list)
    junk_map = {}
    cols = []

    for j, sigma in enumerate(reps):
        col = {}
        for i in range(1, d):  # internal faces
            merged = (sigma[i-1] + sigma[i]) % p
            if merged not in QR_set and merged != 0:
                face = list(sigma)
                face[i-1] = merged
                del face[i]
                ft = zm_orbit_class(tuple(face), QR_list, p)
                if ft not in junk_map:
                    junk_map[ft] = len(junk_map)
                row = junk_map[ft]
                col[row] = (col.get(row, 0) + (-1)**i) % q_mod
        col = {k: v for k, v in col.items() if v != 0}
        cols.append(col)

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

    return len(reduced)

# ============================================================
p = 19
m = (p - 1) // 2  # m = 9
QR = get_QR(p)
q_mod = 997

print(f"P_{p} (m={m}), QR = {QR}")
print(f"Working mod {q_mod}\n")

# Known bottom half: Ω_orb[0..8]
omega_bottom = [1, 1, 8, 60, 417, 2648, 15140, 76474, 331958]
print(f"Known bottom: Ω_orb[0..8] = {omega_bottom}")

# Bottom recursion
R_bottom = [0, 0]
for d in range(2, 10):
    R_bottom.append(omega_bottom[d-1] - R_bottom[-1])
print(f"R_orb (bottom): {R_bottom}")
print(f"R_9 = {R_bottom[9]}")

# Compute top half: d=10,11,...,18
omega_top = {}
t_total = time.time()

for d in range(18, 9, -1):  # 18, 17, 16, ..., 10
    print(f"\n=== d={d} ===")
    t0 = time.time()
    reps = build_orbit_reps_fast(p, d, QR)
    t1 = time.time()
    n_reps = len(reps)
    print(f"  {n_reps} orbit reps in {t1-t0:.1f}s")

    if n_reps == 0:
        omega_top[d] = 0
        print(f"  Ω_{d}^orb = 0")
        continue

    t0 = time.time()
    cols, n_junk = build_constraint_cols_sparse(p, d, reps, QR, q_mod)
    t1 = time.time()
    print(f"  {n_junk} junk orbits, constraints built in {t1-t0:.1f}s")

    t0 = time.time()
    rank_C = sparse_rank_mod(cols, q_mod, n_junk)
    t1 = time.time()
    omega_d = n_reps - rank_C
    omega_top[d] = omega_d
    print(f"  rank(C_{d}) = {rank_C}, Ω_{d}^orb = {omega_d} in {t1-t0:.1f}s")

    del reps, cols  # free memory

# Now compute Ω_9 from Euler characteristic
print(f"\n{'='*60}")
print(f"ASSEMBLING FULL SEQUENCE")
print(f"{'='*60}")

# Build full omega array
omega_full = list(omega_bottom)  # d=0..8
# d=9 unknown, d=10..18 from top
for d in range(10, 19):
    omega_full.append(omega_top.get(d, 0))

# chi = p = sum (-1)^d Omega_d
# Omega_9 = (p - sum_{d!=9} (-1)^d Omega_d) / (-1)^9
# = (p - sum_{d!=9} (-1)^d Omega_d) * (-1)
chi_partial = sum((-1)**d * omega_full[d] for d in range(9))
chi_partial += sum((-1)**d * omega_top.get(d, 0) for d in range(10, 19))
# chi = chi_partial + (-1)^9 * Omega_9 = chi_partial - Omega_9
# p = chi_partial - Omega_9
# Omega_9 = chi_partial - p
omega_9 = chi_partial - p
print(f"\nchi partial (d≠9) = {chi_partial}")
print(f"Ω_9^orb = chi_partial - p = {chi_partial} - {p} = {omega_9}")

# Insert into full sequence
omega_all = list(omega_bottom) + [omega_9] + [omega_top.get(d, 0) for d in range(10, 19)]
print(f"\nFull Ω_orb = {omega_all}")
chi = sum((-1)**d * omega_all[d] for d in range(len(omega_all)))
print(f"χ = {chi} (should be {p})")

# Top recursion: R_18 = Ω_18, R_17 = Ω_17 - R_18, ...
R_top = [0] * 19
R_top[18] = omega_all[18] if 18 < len(omega_all) else 0
# Actually: ∂_d: Ω_d → Ω_{d-1}, so R_d = rank(∂_d).
# For d ≥ m+2: β_d = 0, so Ω_d = R_d + R_{d+1}
# Top recursion: R_{2m} = Ω_{2m}, R_{2m-1} = Ω_{2m-1} - R_{2m}, etc.
# (since β_d = 0 for d ≥ m+2)

for d in range(2*m, m+1, -1):  # d = 18, 17, ..., 11
    if d == 2*m:
        R_top[d] = omega_all[d]
    else:
        R_top[d] = omega_all[d] - R_top[d+1]
    print(f"R_{d} = {R_top[d]}")

# Now: β_{m+1} = Ω_{m+1} - R_{m+1} - R_{m+2}
# But β_{m+1} might be nonzero! We need R_{m+1} = R_10.
# From top recursion: R_{m+2} = Ω_{m+2} - R_{m+3} (since β_{m+2} = 0)
# R_{m+1} = Ω_{m+1} - R_{m+2} - β_{m+1}
# This is circular. But from the Euler characteristic formula:
# R_10 = alternating sum from top: R_10 = Ω_10 - Ω_11 + Ω_12 - ... + Ω_18
# (only if β_d = 0 for d ≥ 11)

# Actually from top: for d > m+1 (d ≥ 11): β_d = 0
# So R_11 = Ω_11 - R_12, R_12 = Ω_12 - R_13, etc.
# R_11 = alternating sum of Ω_11,...,Ω_18
R_11 = sum((-1)**(d-11) * omega_all[d] for d in range(11, 19))
print(f"\nR_11 = {R_11} (from top recursion)")

# β_10 = Ω_10 - R_10 - R_11
# But β_10 = β_m+1. And R_10 from bottom recursion = Ω_9 - R_9 = ...wait
# From bottom: R_10 = Ω_9 - R_9 - β_9 (since β_9 might be nonzero)
# From top: R_10 = Ω_10 - R_11 - β_10 (β_10 might be nonzero)

# Key relation from THM-130: β_9 = β_m, β_10 = β_{m+1}
# And χ = 1 - β_m + β_{m+1} = p
# So β_{m+1} = p - 1 + β_m

# From bottom at d=m: Budget_m = Ω_m - R_m
# β_m + R_{m+1} = Budget_m
# From top at d=m+1: Budget_{m+1} = Ω_{m+1} - R_{m+1}  (via top recursion to get R_{m+2})
# β_{m+1} + R_{m+2}...wait this gets complicated.

# Simpler: use the FULL recursion.
# From top: R_{m+2} through R_{2m} are determined (β=0 for d≥m+2).
# R_{m+1} from top: R_{m+1} = Ω_{m+1} - R_{m+2} - β_{m+1}
# But we don't know β_{m+1} yet.

# Use the system of equations:
# β_m = Ω_m - R_m - R_{m+1}
# β_{m+1} = Ω_{m+1} - R_{m+1} - R_{m+2}  [where R_{m+2} is known from top]
# β_{m+1} - β_m = p - 1  [from χ = p]

# From the third equation: β_{m+1} = β_m + p - 1
# Substitute into second: β_m + p - 1 = Ω_{m+1} - R_{m+1} - R_{m+2}
# From first: R_{m+1} = Ω_m - R_m - β_m
# Substitute: β_m + p - 1 = Ω_{m+1} - (Ω_m - R_m - β_m) - R_{m+2}
# β_m + p - 1 = Ω_{m+1} - Ω_m + R_m + β_m - R_{m+2}
# p - 1 = Ω_{m+1} - Ω_m + R_m - R_{m+2}
# This is an identity that should hold regardless of β_m!

# Check: p - 1 = Ω_10 - Ω_9 + R_9 - R_11
check_lhs = p - 1
check_rhs = omega_all[10] - omega_9 + R_bottom[9] - R_11
print(f"\nConsistency check: p-1 = Ω_10 - Ω_9 + R_9 - R_11")
print(f"  LHS = {check_lhs}")
print(f"  RHS = {omega_all[10]} - {omega_9} + {R_bottom[9]} - {R_11} = {check_rhs}")
print(f"  Match: {'✓' if check_lhs == check_rhs else '✗'}")

# To find β_m, use the first equation: β_m = Ω_m - R_m - R_{m+1}
# And R_{m+1} = Ω_m - R_m - β_m → circular again.
# We need one more equation. Use:
# R_{m+1} from the top: R_{m+1} = Ω_{m+1} - R_{m+2} - β_{m+1} = Ω_{m+1} - R_{m+2} - (β_m + p - 1)
# Set equal to bottom: Ω_m - R_m - β_m = Ω_{m+1} - R_{m+2} - β_m - (p - 1)
# Ω_m - R_m = Ω_{m+1} - R_{m+2} - (p - 1)
# Already shown this holds. So we have ONE free parameter: β_m.

# Budget_m = Ω_m - R_m = omega_9 - R_bottom[9]
Budget_m = omega_9 - R_bottom[9]
print(f"\nBudget_m = Ω_9 - R_9 = {omega_9} - {R_bottom[9]} = {Budget_m}")

# Budget_{m+1} = Ω_{m+1} - R_{m+2}
R_m_plus_2 = R_11  # Wait: R_{m+2} = R_11 from top recursion
Budget_m_plus_1 = omega_all[10] - R_11
print(f"Budget_{m+1} = Ω_10 - R_11 = {omega_all[10]} - {R_11} = {Budget_m_plus_1}")

# Budget_m = β_m + R_{m+1}
# Budget_{m+1} = β_{m+1} + R_{m+1}  [NO! R_{m+1} is INCOMING to d=m, not d=m+1]
# Wait, let me be careful about notation.
# R_d = rank(∂_d: Ω_d → Ω_{d-1}) = image of d landing in degree d-1
# β_d = Ω_d - R_d - R_{d+1}
# where R_d is the rank of boundary FROM degree d, R_{d+1} is rank FROM degree d+1

# So: β_m = Ω_m - R_m - R_{m+1}
# β_{m+1} = Ω_{m+1} - R_{m+1} - R_{m+2}

# Budget_m = Ω_m - R_m = β_m + R_{m+1}
# Budget_{m+1} = Ω_{m+1} - R_{m+2} = β_{m+1} + R_{m+1}

# So: Budget_{m+1} - Budget_m = β_{m+1} - β_m = p - 1
print(f"\nBudget_{m+1} - Budget_m = {Budget_m_plus_1} - {Budget_m} = {Budget_m_plus_1 - Budget_m}")
print(f"Should equal p - 1 = {p - 1}: {'✓' if Budget_m_plus_1 - Budget_m == p - 1 else '✗'}")

# Budget_m = β_m + R_{m+1}
# Budget_{m+1} = β_{m+1} + R_{m+1} = (β_m + p - 1) + R_{m+1} = Budget_m + p - 1 ✓

# We STILL can't determine β_m without more info. But THM-130 predicts:
# β_m^orb = (m-3)/2 = 3
# β_m = m·(m-3)/2 = 27
# β_{m+1} = β_m + p - 1 = 27 + 18 = 45 = m(m+1)/2 = C(10,2)

beta_m_pred = m * (m - 3) // 2  # 27
R_m_plus_1_pred = Budget_m - beta_m_pred

print(f"\n{'='*60}")
print(f"PREDICTED (THM-130):")
print(f"{'='*60}")
print(f"β_m^orb = (m-3)/2 = {(m-3)//2}")
print(f"β_m = m(m-3)/2 = {beta_m_pred}")
print(f"β_{{m+1}} = β_m + p-1 = {beta_m_pred + p - 1} = C(m+1,2) = {m*(m+1)//2}")
print(f"R_{{m+1}} = Budget_m - β_m = {Budget_m} - {beta_m_pred} = {R_m_plus_1_pred}")
print(f"R_{{m+1}}^orb = R_{{m+1}}/m = {R_m_plus_1_pred/m}")

# Check: β_{m+1} = Budget_{m+1} - R_{m+1}
beta_mp1_check = Budget_m_plus_1 - R_m_plus_1_pred
print(f"\nβ_{{m+1}} check = Budget_{{m+1}} - R_{{m+1}} = {Budget_m_plus_1} - {R_m_plus_1_pred} = {beta_mp1_check}")
print(f"Should be {beta_m_pred + p - 1}: {'✓' if beta_mp1_check == beta_m_pred + p - 1 else '✗'}")

# Summary
t_end = time.time()
print(f"\n{'='*60}")
print(f"SUMMARY")
print(f"{'='*60}")
print(f"Full Ω_orb = {omega_all}")
print(f"χ = {chi}")
print(f"Total time: {t_end - t_total:.1f}s")
