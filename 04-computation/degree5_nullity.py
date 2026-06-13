"""
Degree-5 nullity analysis for Paley tournaments.
At degree 5, there are 4 interior faces with signs -1,+1,-1,+1.
Same-sign pairs: (1,3) and (2,4).

From the degree-4 analysis, each same-sign pair should contribute m²
coupling nullity via the tensor product mechanism.

PREDICTION: nullity_5 involves both (1,3) and (2,4) couplings.
If independent: nullity = 2*m². But may be more complex.

opus-2026-03-13-S71b
"""
import numpy as np

def analyze_degree5(p):
    m = (p - 1) // 2
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    QR_list = sorted(QR)

    # Build A_5: all (s,t,u,v,w) with s,t,u,v,w ∈ QR and all partial sums distinct
    A5 = []
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
                    # Also: stuv ≠ s (t+u+v ≠ 0), stuv ≠ st (u+v ≠ 0 ✓)
                    if (t + u + v) % p == 0: continue
                    for w in QR_list:
                        stuvw = (stuv + w) % p
                        tuvw = (tuv + w) % p
                        uvw = (uv + w) % p
                        vw = (v + w) % p
                        if stuvw == 0 or tuvw == 0 or uvw == 0 or vw == 0: continue
                        # Check all partial sum distinctness
                        # PS: 0, s, s+t, s+t+u, s+t+u+v, s+t+u+v+w
                        # Need stuvw ≠ s,st,stu: t+u+v+w≠0, u+v+w≠0, v+w≠0 (checked)
                        if (t + u + v + w) % p == 0: continue
                        if (u + v + w) % p == 0: continue
                        A5.append((s, t, u, v, w))

    print(f"P_{p} (m={m}): |A_5| = {len(A5)}, |A_5|/p = {len(A5)//p}")

    # Interior faces: face 1 through face 4
    # face 1 = (s+t, u, v, w) — junk iff s+t ∉ QR
    # face 2 = (s, t+u, v, w) — junk iff t+u ∉ QR
    # face 3 = (s, t, u+v, w) — junk iff u+v ∉ QR
    # face 4 = (s, t, u, v+w) — junk iff v+w ∉ QR

    junk_sets = [set(), set(), set(), set()]  # for faces 1-4
    signs = [-1, +1, -1, +1]  # (-1)^i for i=1,2,3,4

    for s, t, u, v, w in A5:
        st = (s + t) % p
        tu = (t + u) % p
        uv = (u + v) % p
        vw = (v + w) % p

        if st not in QR: junk_sets[0].add((st, u, v, w))
        if tu not in QR: junk_sets[1].add((s, tu, v, w))
        if uv not in QR: junk_sets[2].add((s, t, uv, w))
        if vw not in QR: junk_sets[3].add((s, t, u, vw))

    for i in range(4):
        print(f"  |junk_f{i+1}| = {len(junk_sets[i])}")

    # Check pairwise overlaps (should all be 0 by QR/NQR type argument)
    # f1: (NQR, QR, QR, QR), f2: (QR, NQR, QR, QR),
    # f3: (QR, QR, NQR, QR), f4: (QR, QR, QR, NQR)
    for i in range(4):
        for j in range(i+1, 4):
            overlap = junk_sets[i] & junk_sets[j]
            if overlap:
                print(f"  WARNING: f{i+1} ∩ f{j+1} = {len(overlap)}")

    all_junk = set()
    for s in junk_sets:
        all_junk |= s
    total_junk = len(all_junk)
    print(f"  Total junk rows = {total_junk}")

    # Build constraint matrix
    all_junk_list = sorted(all_junk)
    junk_idx = {j: i for i, j in enumerate(all_junk_list)}

    n_rows = len(all_junk_list)
    n_cols = len(A5)
    C = np.zeros((n_rows, n_cols), dtype=np.int8)

    for j, (s, t, u, v, w) in enumerate(A5):
        st = (s + t) % p
        tu = (t + u) % p
        uv = (u + v) % p
        vw = (v + w) % p

        faces = [(st, u, v, w), (s, tu, v, w), (s, t, uv, w), (s, t, u, vw)]
        for i, face in enumerate(faces):
            if face in junk_idx:
                C[junk_idx[face], j] += signs[i]

    rank = np.linalg.matrix_rank(C.astype(np.float64))
    nullity = n_rows - rank

    print(f"\n  Constraint matrix: {n_rows} × {n_cols}")
    print(f"  rank(C_5) = {rank}")
    print(f"  nullity = {nullity}")
    print(f"  m² = {m**2}")
    print(f"  2*m² = {2*m**2}")
    print(f"  nullity / m² = {nullity / m**2:.4f}")
    print(f"  Omega_5/p = {len(A5)//p} - {rank} = {len(A5)//p - rank}")

    # Per-block rank analysis
    for i in range(4):
        rows_i = [junk_idx[j] for j in sorted(junk_sets[i])]
        Ci = C[rows_i, :]
        ri = np.linalg.matrix_rank(Ci.astype(np.float64))
        print(f"  Block {i+1}: {len(rows_i)} rows, rank = {ri}, nullity = {len(rows_i)-ri}")

    return rank, nullity, m

for p in [7, 11]:
    print(f"\n{'='*60}")
    r, n, m = analyze_degree5(p)
    print(f"\nSUMMARY: nullity={n}, m²={m**2}, 2m²={2*m**2}")
    print(f"  ratio nullity/m² = {n/m**2:.4f}")
    print(f"{'='*60}")
