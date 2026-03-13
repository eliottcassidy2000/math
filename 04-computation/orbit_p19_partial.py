#!/usr/bin/env python3
"""
P_19 orbit complex: partial computation (bottom half only).
m=9, so we need d=0..9 for bottom recursion, d=10..18 for top.
d=0..7 is feasible, d=8+ gets too large.

Compute Omega_orb for d=0..7, then use the partial data to verify patterns.
Also compute the FULL Omega values (Omega_d) since we need those.

opus-2026-03-13-S71b
"""
import time

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

def sparse_rank_mod(cols_sparse, q):
    pivot_col = {}
    reduced = []
    for col in cols_sparse:
        c = dict(col)
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
                c = {k: v for k, v in c.items() if v % q != 0}
            else:
                pivot_col[r] = len(reduced)
                reduced.append(c)
                break
    return len(reduced)

p = 19
m = 9
QR = get_QR(p)
q = 997

print(f"P_{p} (m={m}), QR = {QR}\n")

# Build orbit reps up to d=7
max_d = 7
orbit_reps = {}
for d in range(max_d + 1):
    t1 = time.time()
    orbit_reps[d] = build_orbit_reps(p, d)
    dt = time.time() - t1
    print(f"d={d}: {len(orbit_reps[d])} reps ({dt:.1f}s)", flush=True)

QR_set = set(QR)

# Orbit constraint
C_orb_cols = {}
C_orb_n_rows = {}
Omega_orb = {}
for d in range(max_d + 1):
    reps = orbit_reps[d]
    if d <= 1:
        C_orb_cols[d] = []
        C_orb_n_rows[d] = 0
        Omega_orb[d] = len(reps)
        continue
    junk_orbits = {}
    cols = []
    for sigma in reps:
        col = {}
        for i in range(1, d):
            merged = (sigma[i-1] + sigma[i]) % p
            if merged not in QR_set:
                face = list(sigma)
                face[i-1] = merged
                del face[i]
                ft = zm_orbit_class(tuple(face), QR, p)
                if ft not in junk_orbits:
                    junk_orbits[ft] = len(junk_orbits)
                r = junk_orbits[ft]
                sign = ((-1) ** i) % q
                col[r] = (col.get(r, 0) + sign) % q
        col = {k: v for k, v in col.items() if v % q != 0}
        cols.append(col)
    C_orb_cols[d] = cols
    C_orb_n_rows[d] = len(junk_orbits)
    rank_C = sparse_rank_mod(cols, q)
    Omega_orb[d] = len(reps) - rank_C

print(f"\nOmega_orb (d=0..{max_d}): {[Omega_orb[d] for d in range(max_d+1)]}")
print(f"Omega (full, multiply by m): {[Omega_orb[d]*m if d > 0 else 1 for d in range(max_d+1)]}")

# Stacking ranks
R_orb = {0: 0}
for d in range(1, max_d + 1):
    t1 = time.time()
    reps_d = orbit_reps[d]
    reps_dm1 = orbit_reps[d-1]
    reps_dm1_idx = {rep: i for i, rep in enumerate(reps_dm1)}
    n_junk = C_orb_n_rows[d]
    stacked_cols = []
    for j, sigma in enumerate(reps_d):
        col = {}
        if j < len(C_orb_cols[d]):
            for r, v in C_orb_cols[d][j].items():
                col[r] = v
        for fi in range(d + 1):
            face = compute_face(sigma, fi, p)
            face_canon = zm_orbit_class(face, QR, p)
            if face_canon in reps_dm1_idx:
                r = n_junk + reps_dm1_idx[face_canon]
                sign = ((-1) ** fi) % q
                col[r] = (col.get(r, 0) + sign) % q
        col = {k: v for k, v in col.items() if v % q != 0}
        stacked_cols.append(col)
    rank_stacked = sparse_rank_mod(stacked_cols, q)
    rank_C = sparse_rank_mod(C_orb_cols[d], q) if C_orb_cols[d] else 0
    R_orb[d] = rank_stacked - rank_C
    dt = time.time() - t1
    print(f"  d={d}: R_orb={R_orb[d]} ({dt:.1f}s)", flush=True)

print(f"\nR_orb (d=0..{max_d}): {[R_orb[d] for d in range(max_d+1)]}")
print(f"R (full, multiply by m): {[R_orb[d]*m if d > 0 else 0 for d in range(max_d+1)]}")

# Bottom recursion check
print(f"\nBottom recursion check:")
R_rec = {0: 0}
R_rec[1] = 0  # ∂_1 = 0
for d in range(2, max_d + 1):
    R_rec[d] = Omega_orb[d-1] - R_rec[d-1]
    match = "✓" if R_rec[d] == R_orb[d] else "✗"
    print(f"  R_{d} = Omega_{d-1} - R_{d-1} = {Omega_orb[d-1]} - {R_rec[d-1]} = {R_rec[d]} {match}")

# Extrapolate R_m (d=9) using recursion
print(f"\nExtrapolation to d=m={m}:")
print(f"  Need Omega_orb for d={max_d+1}..{m-1}")
print(f"  From d=0..{max_d} recursion, next would be:")
print(f"  R_{max_d+1} = Omega_{max_d} - R_{max_d} = ... (need Omega_{max_d})")

# Known Omega from p19_omega scripts: Omega = [1, 9, 72, 540, 3753, ?, ?, ...]
# Wait, we only have through d=4.
# Our orbit gives: Omega_full = m * Omega_orb for d >= 1
# d=0: 1, d=1: 9, d=2: 72, d=3: 540, d=4: 3753
print(f"\n  Full Omega from orbit: {[Omega_orb[d]*m if d > 0 else 1 for d in range(max_d+1)]}")
print(f"  Known from scripts: [1, 9, 72, 540, 3753, ?, ...]")
print(f"  Match through d=4: "
      f"{'✓' if [1, 9, 72, 540, 3753] == [Omega_orb[d]*m if d > 0 else 1 for d in range(5)] else '✗'}")

# New Omega values from orbit computation!
print(f"\n  NEW Omega values for P_19:")
for d in range(max_d + 1):
    val = Omega_orb[d] * m if d > 0 else 1
    print(f"    Omega_{d} = {val}")
