#!/usr/bin/env python3
"""
Extract the explicit non-boundary orbit m-cycle for P_11 (m=5).

β_orb = 1 means there's exactly ONE independent orbit 5-cycle that
is NOT in the image of ∂_6. We extract it using modular arithmetic
and examine its structure to understand the diagonal-type correspondence.

For P_11: m=5, QR = {1,3,4,5,9}, the single diagonal type is gap 2.

opus-2026-03-13-S71b
"""

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_orbit_reps(p, d, QR_list):
    """Build canonical orbit representatives for A_d.
    Fix first entry = 1 (since Z_m acts freely, every orbit has
    exactly one rep starting with 1)."""
    QR_set = set(QR_list)
    m = len(QR_list)
    if d == 0:
        return [()]
    if d == 1:
        return [(1,)]  # single orbit

    results = []
    def backtrack(seq, ps_set):
        if len(seq) == d:
            results.append(tuple(seq))
            return
        last_ps = sum(seq) % p
        for s in QR_list:
            new_ps = (last_ps + s) % p
            if new_ps in ps_set or new_ps == 0:
                continue
            seq.append(s)
            ps_set.add(new_ps)
            backtrack(seq, ps_set)
            seq.pop()
            ps_set.discard(new_ps)

    # Start with s_1 = 1
    backtrack([1], {0, 1})
    return sorted(results)

def compute_face(seq, i, p):
    d = len(seq)
    if i == 0:
        return seq[1:]
    elif i == d:
        return seq[:-1]
    else:
        merged = (seq[i-1] + seq[i]) % p
        return seq[:i-1] + (merged,) + seq[i+1:]

def zm_orbit_class(seq, QR_list, p):
    """Canonical orbit representative = lexicographic min over Z_m scaling."""
    best = seq
    for q in QR_list:
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best

def orbit_index(seq, QR_list, p):
    """Map a sequence to its canonical orbit rep, returning the scaling factor."""
    best = seq
    best_q = 1
    for q in QR_list:
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
            best_q = q
    return best, best_q

def build_sparse_boundary_and_constraint(p, d, reps_d, reps_dm1, QR_list):
    """Build sparse boundary and constraint columns for orbit complex."""
    QR_set = set(QR_list)
    m = len(QR_list)

    # Index d-1 reps
    rep_idx = {rep: i for i, rep in enumerate(reps_dm1)}

    # Also need to track junk faces for constraint
    junk_faces = {}

    boundary_cols = []
    constraint_cols = []

    for j, sigma in enumerate(reps_d):
        b_col = {}
        c_col = {}

        for fi in range(d + 1):
            face = compute_face(sigma, fi, p)
            # Check if face is in A_{d-1}
            face_canon = zm_orbit_class(face, QR_list, p)

            if face_canon in rep_idx:
                # Valid face
                row = rep_idx[face_canon]
                b_col[row] = b_col.get(row, 0) + (-1)**fi
            else:
                # Junk face — contributes to constraint
                if face_canon not in junk_faces:
                    junk_faces[face_canon] = len(junk_faces)
                row = junk_faces[face_canon]
                c_col[row] = c_col.get(row, 0) + (-1)**fi

        boundary_cols.append(b_col)
        constraint_cols.append(c_col)

    return boundary_cols, constraint_cols, len(junk_faces)

def gauss_elim_mod(cols, q, n_rows):
    """Sparse Gaussian elimination mod q. Returns (rank, pivot_info, reduced_matrix)."""
    pivot_col = {}
    reduced = []
    col_order = []  # which original column each reduced col came from

    for j, col in enumerate(cols):
        c = dict(col)
        # Reduce
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
                col_order.append(j)
                break

    return len(reduced), pivot_col, reduced, col_order

def find_kernel_mod(cols, q, n_rows, n_cols):
    """Find kernel vectors of the column matrix mod q.
    Returns list of kernel vectors as dicts {col_idx: value}."""
    # Track reduction: each column becomes a linear combination
    # We track: reduced_col = original_col - sum(factor * pivot_cols)
    pivot_col = {}
    reduced = []
    col_order = []
    history = []  # For each reduced col, the combination coefficients

    non_pivot_cols = []
    col_combos = {}  # col j -> dict of {pivot_col_idx: factor} used to reduce it

    for j in range(n_cols):
        c = dict(cols[j])
        combo = {}  # tracks: current = col_j - sum(combo[k] * col_{col_order[k]})

        while c:
            r = min(c.keys())
            v = c[r] % q
            if v == 0:
                del c[r]
                continue
            if r in pivot_col:
                pidx = pivot_col[r]
                pv = reduced[pidx]
                pv_val = pv[r]
                inv_pv = pow(pv_val, q - 2, q)
                factor = (v * inv_pv) % q
                combo[pidx] = (combo.get(pidx, 0) + factor) % q
                for pr, pval in pv.items():
                    c[pr] = (c.get(pr, 0) - factor * pval) % q
                c = {k: v % q for k, v in c.items() if v % q != 0}
            else:
                pivot_col[r] = len(reduced)
                reduced.append(c)
                col_order.append(j)
                history.append(combo)
                break
        else:
            # Column reduced to zero — it's in the span
            # col_j = sum(combo[k] * col_{col_order[k]})
            # So col_j - sum(...) = 0 is a kernel relation
            non_pivot_cols.append((j, combo))

    # Build kernel vectors
    kernel = []
    for j, combo in non_pivot_cols:
        kvec = {j: 1}
        for pidx, factor in combo.items():
            orig_col = col_order[pidx]
            kvec[orig_col] = (-factor) % q
        kernel.append(kvec)

    return kernel

# ============================================================
# MAIN: Extract the cycle for P_11
# ============================================================

p = 11
m = (p - 1) // 2
QR = get_QR(p)
q = 997  # working modulus

print(f"P_{p} (m={m}), QR = {QR}")
print(f"Working mod {q}\n")

# Build orbit reps for degrees m-1 through m+2
print("Building orbit reps...")
reps = {}
for d in range(m + 3):
    reps[d] = build_orbit_reps(p, d, QR)
    print(f"  d={d}: {len(reps[d])} orbit reps")

# Build boundary and constraint for degree m
print(f"\nBuilding boundary ∂_{m} and constraint C_{m}...")
B_cols_m, C_cols_m, n_junk_m = build_sparse_boundary_and_constraint(
    p, m, reps[m], reps[m-1], QR)
print(f"  {len(reps[m])} source orbits, {len(reps[m-1])} target orbits, {n_junk_m} junk orbits")

# Build boundary and constraint for degree m+1
print(f"Building boundary ∂_{m+1} and constraint C_{m+1}...")
B_cols_m1, C_cols_m1, n_junk_m1 = build_sparse_boundary_and_constraint(
    p, m+1, reps[m+1], reps[m], QR)
print(f"  {len(reps[m+1])} source orbits, {len(reps[m])} target orbits, {n_junk_m1} junk orbits")

# Compute ranks using stacking trick
print(f"\n=== Rank computations ===")

# R_m = rank([C_m; B_m]) - rank(C_m)
n_src_m = len(reps[m])
n_tgt_m = len(reps[m-1])

# Stack constraint on top, boundary below
stacked_m = []
for j in range(n_src_m):
    col = dict(C_cols_m[j])
    for r, v in B_cols_m[j].items():
        col[r + n_junk_m] = col.get(r + n_junk_m, 0) + v
    stacked_m.append(col)

rank_C_m, _, _, _ = gauss_elim_mod(C_cols_m, q, n_junk_m)
rank_stacked_m, _, _, _ = gauss_elim_mod(stacked_m, q, n_junk_m + n_tgt_m)
R_m = rank_stacked_m - rank_C_m
Omega_m = n_src_m - rank_C_m

print(f"rank(C_{m}) = {rank_C_m}")
print(f"rank([C_{m}; B_{m}]) = {rank_stacked_m}")
print(f"R_{m} = {R_m}")
print(f"Ω_{m}^orb = {Omega_m}")

# R_{m+1} = rank([C_{m+1}; B_{m+1}]) - rank(C_{m+1})
n_src_m1 = len(reps[m+1])
n_tgt_m1 = len(reps[m])

stacked_m1 = []
for j in range(n_src_m1):
    col = dict(C_cols_m1[j])
    for r, v in B_cols_m1[j].items():
        col[r + n_junk_m1] = col.get(r + n_junk_m1, 0) + v
    stacked_m1.append(col)

rank_C_m1, _, _, _ = gauss_elim_mod(C_cols_m1, q, n_junk_m1)
rank_stacked_m1, _, _, _ = gauss_elim_mod(stacked_m1, q, n_junk_m1 + n_tgt_m1)
R_m1 = rank_stacked_m1 - rank_C_m1
Omega_m1 = n_src_m1 - rank_C_m1

print(f"\nrank(C_{m+1}) = {rank_C_m1}")
print(f"rank([C_{m+1}; B_{m+1}]) = {rank_stacked_m1}")
print(f"R_{m+1} = {R_m1}")
print(f"Ω_{m+1}^orb = {Omega_m1}")

beta_m = Omega_m - R_m - R_m1
print(f"\nβ_{m}^orb = Ω_{m} - R_{m} - R_{m+1} = {Omega_m} - {R_m} - {R_m1} = {beta_m}")

# Now: find the kernel of ∂_m restricted to Ω_m
# We need ker(∂_m) ∩ Ω_m = ker([C_m; B_m]) ... no.
# Actually: ker(∂_m|Ω_m) = {v ∈ ker(C_m) : B_m·v = 0}
# = ker([C_m; B_m]) (joint kernel)

print(f"\n=== Extracting kernel of [C_{m}; B_{m}] ===")

# Build the combined matrix
combined_m = []
for j in range(n_src_m):
    col = dict(C_cols_m[j])
    for r, v in B_cols_m[j].items():
        col[r + n_junk_m] = col.get(r + n_junk_m, 0) + v
    combined_m.append(col)

kernel_vecs = find_kernel_mod(combined_m, q, n_junk_m + n_tgt_m, n_src_m)
print(f"dim(ker(∂_{m}|Ω_{m})) = {len(kernel_vecs)}")
print(f"Expected: R_{m+1} + β_{m} = {R_m1} + {beta_m} = {R_m1 + beta_m}")

# Now find image of ∂_{m+1} in Ω_m
# im(∂_{m+1}|Ω_{m+1}) lands in Ω_m
# We need to find which kernel vectors are NOT in the image

# The image vectors are B_{m+1} applied to each Ω_{m+1} basis vector
# But we need Ω_{m+1} basis first = ker(C_{m+1})

print(f"\n=== Finding Ω_{m+1} basis (ker C_{m+1}) ===")
omega_m1_kernel = find_kernel_mod(C_cols_m1, q, n_junk_m1, n_src_m1)
print(f"dim(Ω_{m+1}) = {len(omega_m1_kernel)}, expected {Omega_m1}")

# Apply ∂_{m+1} to each Ω_{m+1} basis vector
# Result is in the space of reps[m]
print(f"\n=== Computing image of ∂_{m+1} in A_{m} ===")
image_cols = []
for kvec in omega_m1_kernel:
    # kvec is a linear combination of reps[m+1] columns
    # Apply B_{m+1}: each column of B_{m+1} corresponds to a rep[m+1] column
    img = {}
    for j, coeff in kvec.items():
        for r, v in B_cols_m1[j].items():
            img[r] = (img.get(r, 0) + coeff * v) % q
    img = {k: v % q for k, v in img.items() if v % q != 0}
    image_cols.append(img)

# Now we need to check: which kernel vectors of [C_m; B_m] are NOT
# in the span of image_cols (when projected to A_m)?
# Actually, the kernel vectors live in A_m (they're coefficients of reps[m]).
# The image vectors also live in A_m.
# We need to find kernel vectors that are NOT in the column span of image_cols.

# Strategy: stack image_cols, then try adding each kernel vector.
# If the rank increases, that kernel vector is NOT in the image.

print(f"Number of image columns: {len(image_cols)}")
rank_image, _, _, _ = gauss_elim_mod(image_cols, q, n_tgt_m)
print(f"rank(image) = {rank_image}, expected R_{m+1} = {R_m1}")

# Find kernel vectors not in image
print(f"\n=== Identifying non-boundary cycles ===")
non_boundary = []
for i, kvec in enumerate(kernel_vecs):
    # Check if kvec is in span of image_cols
    test_cols = image_cols + [{k: v % q for k, v in kvec.items() if v % q != 0}]
    test_rank, _, _, _ = gauss_elim_mod(test_cols, q, n_tgt_m)
    if test_rank > rank_image:
        non_boundary.append((i, kvec))
        print(f"  Kernel vector {i} is NOT a boundary! (rank increased to {test_rank})")

print(f"\nFound {len(non_boundary)} non-boundary cycle(s) (expected β_{m} = {beta_m})")

# Print the non-boundary cycles
for idx, (i, kvec) in enumerate(non_boundary):
    print(f"\n=== Non-boundary cycle #{idx+1} ===")
    # kvec maps orbit rep indices to coefficients
    print(f"  Support size: {len(kvec)} orbit reps")
    for j, coeff in sorted(kvec.items()):
        rep = reps[m][j]
        # Compute partial sums
        ps = [0]
        for s in rep:
            ps.append((ps[-1] + s) % p)
        print(f"    coeff={coeff:4d}  rep={rep}  partial_sums={ps}")

    # Analyze the structure
    print(f"\n  Structure analysis:")
    total_reps = len(kvec)
    entries = set()
    for j in kvec:
        for s in reps[m][j]:
            entries.add(s)
    print(f"    QR elements appearing: {sorted(entries)}")

    # Check: do the reps have any pattern related to "gap 2" diagonal?
    # In the m-gon interpretation, gap 2 means connecting vertex i to i+2
    # In QR terms, consecutive QR elements (by index) differ by gap g
    print(f"\n  Diff-seq patterns:")
    for j, coeff in sorted(kvec.items()):
        rep = reps[m][j]
        # What are the consecutive ratios?
        ratios = [(rep[k+1] * pow(rep[k], p-2, p)) % p for k in range(len(rep)-1)]
        print(f"    rep={rep}, ratios={ratios}, coeff={coeff}")
