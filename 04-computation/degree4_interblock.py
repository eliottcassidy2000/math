"""
Degree-4 inter-block dependency analysis.
Each block individually has full row rank, but the stacked matrix loses m² rank.
The lost rank comes from columns that hit multiple blocks.

Key insight: a column (s,t,u,v) hits blocks 1,2,3 when faces 1,2,3 are junk.
A column hitting blocks i and j creates a "coupling" between those blocks.

opus-2026-03-13-S71b
"""
import numpy as np

def analyze_interblock(p):
    m = (p - 1) // 2
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))

    # Build A_4 and constraint matrix
    A4 = []
    QR_list = sorted(QR)

    for s in QR_list:
        for t in QR_list:
            st = (s + t) % p
            if st == 0: continue
            for u in QR_list:
                stu = (st + u) % p
                tu = (t + u) % p
                if stu == 0 or stu == s or tu == 0: continue
                for v in QR_list:
                    stuv = (stu + v) % p
                    tuv = (tu + v) % p
                    uv = (u + v) % p
                    if stuv == 0 or tuv == 0 or uv == 0: continue
                    A4.append((s, t, u, v))

    # Classify columns by which blocks they hit
    multi_block_cols = {
        '13': [], '12': [], '23': [], '123': []
    }

    for j, (s, t, u, v) in enumerate(A4):
        st = (s + t) % p
        tu = (t + u) % p
        uv = (u + v) % p
        b1 = st not in QR
        b2 = tu not in QR
        b3 = uv not in QR
        key = ''
        if b1: key += '1'
        if b2: key += '2'
        if b3: key += '3'
        if key in multi_block_cols:
            multi_block_cols[key].append(j)

    print(f"P_{p} (m={m}):")
    print(f"  Multi-block columns: 12={len(multi_block_cols['12'])}, "
          f"13={len(multi_block_cols['13'])}, 23={len(multi_block_cols['23'])}, "
          f"123={len(multi_block_cols['123'])}")

    # The nullity comes from: when we stack the blocks, columns hitting
    # multiple blocks create equations like:
    #   entry_in_block_i + entry_in_block_j = ... (coupling)
    # With individual blocks having full rank, the coupling columns
    # create a "gluing" structure.

    # Let's look at this differently. Consider the column space:
    # Block 1 rows have entries only from columns where s+t ∉ QR
    # Block 2 rows have entries only from columns where t+u ∉ QR
    # Block 3 rows have entries only from columns where u+v ∉ QR

    # The rank loss happens because some linear combinations of
    # multi-block columns are zero across ALL blocks simultaneously.

    # Let's look at the COUPLING matrix: for each multi-block column,
    # what row does it hit in each block?

    # For a (1,3) column: s+t ∈ NQR, u+v ∈ NQR, t+u ∈ QR
    # In block 1: hits row (s+t, u, v) with coeff -1
    # In block 3: hits row (s, t, u+v) with coeff -1
    # (Same sign!)

    # For a (1,2) column: s+t ∈ NQR, t+u ∈ NQR, u+v ∈ QR
    # In block 1: hits row (s+t, u, v) with coeff -1
    # In block 2: hits row (s, t+u, v) with coeff +1
    # (Opposite signs!)

    # For a (2,3) column: t+u ∈ NQR, u+v ∈ NQR, s+t ∈ QR
    # In block 2: hits row (s, t+u, v) with coeff +1
    # In block 3: hits row (s, t, u+v) with coeff -1
    # (Opposite signs!)

    # For a (1,2,3) column: all three NQR
    # In block 1: coeff -1, block 2: coeff +1, block 3: coeff -1

    # The coupling pattern is key. Let me look at the block-1 and block-3
    # columns specifically, since the "13" category dominates.

    # Focus: what's the rank of the COUPLING MATRIX between blocks 1 and 3?
    # For each (1,3) column j, define the pair:
    #   (row_in_block1, row_in_block3) = ((s+t, u, v), (s, t, u+v))

    coupling_13 = []
    for j in multi_block_cols['13']:
        s, t, u, v = A4[j]
        st = (s + t) % p
        uv = (u + v) % p
        coupling_13.append(((st, u, v), (s, t, uv)))

    # How many distinct row-1 values and row-3 values?
    rows1_in_13 = set(c[0] for c in coupling_13)
    rows3_in_13 = set(c[1] for c in coupling_13)
    print(f"\n  (1,3) coupling: {len(coupling_13)} columns, "
          f"{len(rows1_in_13)} distinct B1 rows, {len(rows3_in_13)} distinct B3 rows")

    # Build the coupling bipartite graph / matrix
    b1_list = sorted(rows1_in_13)
    b3_list = sorted(rows3_in_13)
    b1_idx = {r: i for i, r in enumerate(b1_list)}
    b3_idx = {r: i for i, r in enumerate(b3_list)}

    # Each coupling column maps B1_row -> B3_row with both coeff -1
    # So the coupling is: for each column, row_b1 and row_b3 get the same coefficient
    # In the combined system, this means row_b1 - row_b3 = 0 (cancellation)
    # forming a m²-dimensional null space?

    # Actually, let's think about it more carefully.
    # In the stacked matrix [C1; C2; C3], a column that hits block 1 at row r1
    # with coeff a1 and block 3 at row r3 with coeff a3 has nonzero entries at
    # positions (r1_in_block1, a1) and (r3_in_block3, a3).
    # The row rank loss occurs because there exist linear combinations of rows
    # that are zero.

    # Since each block individually has full row rank, the null space of [C1; C2; C3]^T
    # (left null space) consists of vectors [v1; v2; v3] where v1*C1 + v2*C2 + v3*C3 = 0
    # but v1 ≠ 0 or v2 ≠ 0 or v3 ≠ 0.

    # For a column that only hits block i, the constraint is v_i[r_i] * coeff = 0.
    # Since block i individually has full rank, this means v_i is determined by
    # the multi-block column constraints.

    # Let me compute the left null space directly for small p.
    if p <= 11:
        # Build the full stacked constraint matrix
        all_junk_f1 = set()
        all_junk_f2 = set()
        all_junk_f3 = set()
        for s, t, u, v in A4:
            st = (s + t) % p
            tu = (t + u) % p
            uv = (u + v) % p
            if st not in QR: all_junk_f1.add((st, u, v))
            if tu not in QR: all_junk_f2.add((s, tu, v))
            if uv not in QR: all_junk_f3.add((s, t, uv))

        f1_list = sorted(all_junk_f1)
        f2_list = sorted(all_junk_f2)
        f3_list = sorted(all_junk_f3)
        f1_idx = {r: i for i, r in enumerate(f1_list)}
        f2_idx = {r: i for i, r in enumerate(f2_list)}
        f3_idx = {r: i for i, r in enumerate(f3_list)}

        n1, n2, n3 = len(f1_list), len(f2_list), len(f3_list)
        C_full = np.zeros((n1 + n2 + n3, len(A4)), dtype=np.float64)

        for j, (s, t, u, v) in enumerate(A4):
            st = (s + t) % p
            tu = (t + u) % p
            uv = (u + v) % p
            if st not in QR:
                C_full[f1_idx[(st, u, v)], j] = -1
            if tu not in QR:
                C_full[n1 + f2_idx[(s, tu, v)], j] = 1
            if uv not in QR:
                C_full[n1 + n2 + f3_idx[(s, t, uv)], j] = -1

        # SVD to get left null space
        U, S_vals, Vt = np.linalg.svd(C_full, full_matrices=True)
        tol = 1e-8
        null_dim = np.sum(S_vals < tol) if len(S_vals) < n1+n2+n3 else (n1+n2+n3 - np.sum(S_vals > tol))
        actual_rank = np.sum(S_vals > tol)
        print(f"\n  SVD: rank = {actual_rank}, left null dim = {n1+n2+n3 - actual_rank}")

        # Look at the null vectors
        null_start = actual_rank
        null_vecs = U[:, null_start:]
        print(f"  Null space dimension: {null_vecs.shape[1]}")

        # Check which blocks the null vectors span
        if null_vecs.shape[1] > 0 and null_vecs.shape[1] <= 30:
            # For each null vector, check support in each block
            for i in range(min(null_vecs.shape[1], 5)):
                v = null_vecs[:, i]
                s1 = np.max(np.abs(v[:n1]))
                s2 = np.max(np.abs(v[n1:n1+n2]))
                s3 = np.max(np.abs(v[n1+n2:]))
                print(f"  Null vec {i}: block1_max={s1:.4f}, block2_max={s2:.4f}, block3_max={s3:.4f}")

    return m

for p in [7, 11, 19]:
    print(f"\n{'='*60}")
    analyze_interblock(p)
    print(f"{'='*60}")
