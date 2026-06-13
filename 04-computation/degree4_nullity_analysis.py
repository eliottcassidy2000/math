"""
Degree-4 constraint matrix structure analysis for Paley tournaments.
Goal: Understand WHY nullity = m² at degree 4.

A 4-diff-seq (s,t,u,v) has faces:
  face 0 = (t,u,v) — always in A_3 (all entries QR, partial sums OK)
  face 1 = (s+t, u, v) — junk iff s+t ∉ QR
  face 2 = (s, t+u, v) — junk iff t+u ∉ QR
  face 3 = (s, t, u+v) — junk iff u+v ∉ QR
  face 4 = (s,t,u) — always in A_3

The constraint equation for (s,t,u,v):
  C(s,t,u,v) = -junk(face1) + junk(face2) - junk(face3)

where junk(face_i) means the row corresponding to that face in the junk space.

opus-2026-03-13-S71b
"""
import numpy as np
from itertools import product

def analyze_degree4(p):
    m = (p - 1) // 2
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))

    # Build A_4: all (s,t,u,v) with s,t,u,v ∈ QR and all partial sums distinct
    A4 = []
    QR_list = sorted(QR)

    for s in QR_list:
        for t in QR_list:
            st = (s + t) % p
            if st == 0:
                continue
            for u in QR_list:
                stu = (st + u) % p
                if stu == 0 or stu == s:
                    continue
                # Check: partial sums 0, s, s+t, s+t+u must be distinct
                # Already: s≠0, s+t≠0. Need s+t≠s (t≠0 ✓), s+t+u≠0, s+t+u≠s (t+u≠0)
                tu = (t + u) % p
                if tu == 0:
                    continue
                for v in QR_list:
                    stuv = (stu + v) % p
                    if stuv == 0:
                        continue
                    # s+t+u+v ≠ s means t+u+v ≠ 0
                    tuv = (tu + v) % p
                    if tuv == 0:
                        continue
                    # s+t+u+v ≠ s+t means u+v ≠ 0
                    uv = (u + v) % p
                    if uv == 0:
                        continue
                    A4.append((s, t, u, v))

    print(f"P_{p} (m={m}): |A_4| = {len(A4)}, |A_4|/p = {len(A4)//p}")

    # Classify junk faces
    junk_profiles = {'none': 0, 'f1': 0, 'f2': 0, 'f3': 0,
                     'f12': 0, 'f13': 0, 'f23': 0, 'f123': 0}

    # Track junk face tuples
    junk_f1_set = set()  # (r, u, v) where r = s+t ∈ NQR
    junk_f2_set = set()  # (s, r, v) where r = t+u ∈ NQR
    junk_f3_set = set()  # (s, t, r) where r = u+v ∈ NQR

    for s, t, u, v in A4:
        st = (s + t) % p
        tu = (t + u) % p
        uv = (u + v) % p

        f1 = st not in QR  # face 1 junk
        f2 = tu not in QR  # face 2 junk
        f3 = uv not in QR  # face 3 junk

        if f1: junk_f1_set.add((st, u, v))
        if f2: junk_f2_set.add((s, tu, v))
        if f3: junk_f3_set.add((s, t, uv))

        key = ''
        if not f1 and not f2 and not f3:
            key = 'none'
        elif f1 and not f2 and not f3:
            key = 'f1'
        elif not f1 and f2 and not f3:
            key = 'f2'
        elif not f1 and not f2 and f3:
            key = 'f3'
        elif f1 and f2 and not f3:
            key = 'f12'
        elif f1 and not f2 and f3:
            key = 'f13'
        elif not f1 and f2 and f3:
            key = 'f23'
        else:
            key = 'f123'
        junk_profiles[key] += 1

    print(f"\nJunk face profiles (counts / p):")
    for k, v in junk_profiles.items():
        print(f"  {k:6s}: {v//p:8d}  (total {v})")

    # Check overlaps between junk face sets
    f12_overlap = junk_f1_set & junk_f2_set
    f13_overlap = junk_f1_set & junk_f3_set
    f23_overlap = junk_f2_set & junk_f3_set
    f123_overlap = junk_f1_set & junk_f2_set & junk_f3_set

    print(f"\nJunk face tuple counts:")
    print(f"  |junk_f1| = {len(junk_f1_set)}")
    print(f"  |junk_f2| = {len(junk_f2_set)}")
    print(f"  |junk_f3| = {len(junk_f3_set)}")
    print(f"  |f1 ∩ f2| = {len(f12_overlap)}")
    print(f"  |f1 ∩ f3| = {len(f13_overlap)}")
    print(f"  |f2 ∩ f3| = {len(f23_overlap)}")
    print(f"  |f1 ∩ f2 ∩ f3| = {len(f123_overlap)}")
    total_junk = len(junk_f1_set | junk_f2_set | junk_f3_set)
    print(f"  |f1 ∪ f2 ∪ f3| = {total_junk}")

    # Now build the actual constraint matrix
    # Rows: all junk 3-tuples (from the union)
    # For each column (s,t,u,v) in A_4, the constraint is:
    #   -δ_{face1_junk}(row for face1) + δ_{face2_junk}(row for face2) - δ_{face3_junk}(row for face3)

    all_junk = sorted(junk_f1_set | junk_f2_set | junk_f3_set)
    junk_idx = {j: i for i, j in enumerate(all_junk)}

    n_rows = len(all_junk)
    n_cols = len(A4)

    print(f"\nConstraint matrix: {n_rows} rows × {n_cols} cols")

    # Build sparse constraint matrix
    C = np.zeros((n_rows, n_cols), dtype=np.int8)
    for j, (s, t, u, v) in enumerate(A4):
        st = (s + t) % p
        tu = (t + u) % p
        uv = (u + v) % p

        if st not in QR:
            C[junk_idx[(st, u, v)], j] -= 1  # face 1 sign = (-1)^1
        if tu not in QR:
            C[junk_idx[(s, tu, v)], j] += 1   # face 2 sign = (-1)^2
        if uv not in QR:
            C[junk_idx[(s, t, uv)], j] -= 1   # face 3 sign = (-1)^3

    # Compute rank
    rank = np.linalg.matrix_rank(C.astype(np.float64))
    nullity = n_rows - rank

    print(f"  rank(C_4) = {rank}")
    print(f"  nullity = {n_rows} - {rank} = {nullity}")
    print(f"  nullity / m² = {nullity / m**2:.4f}")
    print(f"  Omega_4/p = {len(A4)//p} - {rank} = {len(A4)//p - rank}")

    # Analyze overlap structure more carefully
    # Key question: what are the junk face tuples in the overlap?
    print(f"\n--- Overlap Analysis ---")

    # f1 faces have type (NQR, QR, QR)
    # f2 faces have type (QR, NQR, QR)
    # f3 faces have type (QR, QR, NQR)
    # So overlaps require: position has to be BOTH QR and NQR — impossible!
    # Wait, that's the KEY: let me check the types.

    # A face in junk_f1 is (r, u, v) where r ∈ {0} ∪ NQR, u ∈ QR, v ∈ QR
    # A face in junk_f2 is (s, r', v) where s ∈ QR, r' ∈ {0} ∪ NQR, v ∈ QR
    # A face in junk_f3 is (s, t, r'') where s ∈ QR, t ∈ QR, r'' ∈ {0} ∪ NQR

    # For a tuple to be in BOTH f1 and f2: first entry must be both NQR and QR
    # => first entry is 0? No, s+t = 0 requires t = -s, but -1 ∉ QR so -s ∉ QR
    # Wait: r can be 0 too (if s+t = 0 mod p). But for QR entries, s+t ≠ 0.
    # Actually for Paley: s,t ∈ QR and -1 ∉ QR means s+t ≠ 0 always.
    # So r = s+t is NEVER 0 for (s,t) ∈ QR×QR.

    # Similarly for t+u and u+v.

    # So: junk_f1 tuples have first entry in NQR
    #     junk_f2 tuples have first entry in QR, second in NQR
    #     junk_f3 tuples have first entry in QR, second in QR, third in NQR

    # f1 ∩ f2: need first entry in NQR ∩ QR = ∅. So overlap is EMPTY.
    # f1 ∩ f3: need first entry in NQR and first entry in QR? NO!
    # Wait: f3 tuple is (s, t, r'') with s ∈ QR. f1 tuple is (r, u, v) with r ∈ NQR.
    # For overlap: r = s (needs to be both QR and NQR) => empty.
    # f2 ∩ f3: f2 = (s, r', v) with second ∈ NQR. f3 = (s, t, r'') with second ∈ QR.
    # For overlap: need second entry both NQR and QR => empty.

    # ALL pairwise overlaps should be EMPTY!

    print(f"f1 ∩ f2 = {len(f12_overlap)} (should be 0)")
    print(f"f1 ∩ f3 = {len(f13_overlap)} (should be 0)")
    print(f"f2 ∩ f3 = {len(f23_overlap)} (should be 0)")

    # So the three blocks are completely disjoint!
    # Then why doesn't the constraint matrix have full row rank?
    # The issue must be WITHIN the blocks (intra-block dependencies).

    # Let's compute rank of each block separately
    f1_rows = [junk_idx[j] for j in sorted(junk_f1_set)]
    f2_rows = [junk_idx[j] for j in sorted(junk_f2_set)]
    f3_rows = [junk_idx[j] for j in sorted(junk_f3_set)]

    # For each column, which block(s) does it contribute to?
    # A column contributes to block 1 if face 1 is junk (s+t ∉ QR)
    # A column contributes to block 2 if face 2 is junk (t+u ∉ QR)
    # A column contributes to block 3 if face 3 is junk (u+v ∉ QR)
    # A column can contribute to multiple blocks!

    # Extract sub-matrices
    C1 = C[f1_rows, :]
    C2 = C[f2_rows, :]
    C3 = C[f3_rows, :]

    rank1 = np.linalg.matrix_rank(C1.astype(np.float64))
    rank2 = np.linalg.matrix_rank(C2.astype(np.float64))
    rank3 = np.linalg.matrix_rank(C3.astype(np.float64))

    print(f"\nPer-block analysis:")
    print(f"  Block 1 (face 1): {len(f1_rows)} rows, rank = {rank1}, nullity = {len(f1_rows)-rank1}")
    print(f"  Block 2 (face 2): {len(f2_rows)} rows, rank = {rank2}, nullity = {len(f2_rows)-rank2}")
    print(f"  Block 3 (face 3): {len(f3_rows)} rows, rank = {rank3}, nullity = {len(f3_rows)-rank3}")
    print(f"  Sum of block ranks = {rank1+rank2+rank3}")
    print(f"  Full rank = {rank}")
    print(f"  Discrepancy = {rank1+rank2+rank3 - rank}")

    # Since the blocks are disjoint in rows, rank(C) = rank([C1; C2; C3])
    # But the columns can overlap between blocks!
    # rank([C1; C2; C3]) ≤ rank(C1) + rank(C2) + rank(C3)
    # with equality iff the column spaces are "compatible"

    print(f"\n  Total nullity = {nullity}")
    print(f"  Sum of block nullities = {len(f1_rows)-rank1 + len(f2_rows)-rank2 + len(f3_rows)-rank3}")

    # Count how many columns hit each combination of blocks
    col_block_counts = {'none': 0, '1': 0, '2': 0, '3': 0,
                        '12': 0, '13': 0, '23': 0, '123': 0}
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
        if not key: key = 'none'
        col_block_counts[key] += 1

    print(f"\nColumn-to-block distribution (counts / p):")
    for k, v in col_block_counts.items():
        print(f"  blocks {k:4s}: {v//p:8d}")

    return rank, nullity, m

for p in [7, 11, 19]:
    print(f"\n{'='*60}")
    r, n, m = analyze_degree4(p)
    print(f"\nSUMMARY: nullity = {n}, m² = {m**2}, match = {n == m**2}")
    print(f"{'='*60}")
