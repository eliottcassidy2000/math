"""
Degree-4: the (1,3) coupling creates exactly m² dependencies.
All null vectors live in blocks 1 and 3 (same sign: both -1).

A (1,3) column has s+t ∈ NQR, u+v ∈ NQR, t+u ∈ QR.
In block 1: hits (s+t, u, v) with coeff -1
In block 3: hits (s, t, u+v) with coeff -1

So the coupling maps block-1 row (r, u, v) to block-3 row (s, t, r')
where r = s+t ∈ NQR and r' = u+v ∈ NQR.

This is a bipartite graph. The null space consists of vectors v1, v3 such that
v1[r,u,v] = v3[s,t,r'] for each coupling edge.

Let's understand this bipartite graph structure.

opus-2026-03-13-S71b
"""
import numpy as np
from collections import defaultdict

def coupling_analysis(p):
    m = (p - 1) // 2
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    NQR = set(range(1, p)) - QR
    QR_list = sorted(QR)

    # Find all (1,3) columns: s+t ∈ NQR, t+u ∈ QR, u+v ∈ NQR
    coupling_edges = []  # (block1_row, block3_row)
    for s in QR_list:
        for t in QR_list:
            st = (s + t) % p
            if st == 0 or st in QR:
                continue  # need s+t ∈ NQR
            for u in QR_list:
                tu = (t + u) % p
                if tu == 0 or tu not in QR:
                    continue  # need t+u ∈ QR
                stu = (st + u) % p
                if stu == 0 or stu == s:
                    continue
                for v in QR_list:
                    uv = (u + v) % p
                    if uv == 0 or uv in QR:
                        continue  # need u+v ∈ NQR
                    stuv = (stu + v) % p
                    tuv = (tu + v) % p
                    if stuv == 0 or tuv == 0:
                        continue
                    # This is a (1,3) column
                    b1_row = (st, u, v)   # in NQR × QR × QR
                    b3_row = (s, t, uv)   # in QR × QR × NQR
                    coupling_edges.append((b1_row, b3_row))

    print(f"P_{p} (m={m}): {len(coupling_edges)} coupling edges")

    # Also need (1,2,3) columns: s+t, t+u, u+v all ∈ NQR
    # Block 1: (s+t, u, v) with -1, Block 2: (s, t+u, v) with +1, Block 3: (s, t, u+v) with -1
    # For the (1,3) null space, these contribute with same-sign entries in B1 and B3
    coupling_edges_123 = []
    for s in QR_list:
        for t in QR_list:
            st = (s + t) % p
            if st == 0 or st in QR: continue
            for u in QR_list:
                tu = (t + u) % p
                if tu == 0 or tu in QR: continue  # need t+u ∈ NQR too
                stu = (st + u) % p
                if stu == 0 or stu == s: continue
                for v in QR_list:
                    uv = (u + v) % p
                    if uv == 0 or uv in QR: continue
                    stuv = (stu + v) % p
                    tuv = (tu + v) % p
                    if stuv == 0 or tuv == 0: continue
                    b1_row = (st, u, v)
                    b3_row = (s, t, uv)
                    coupling_edges_123.append((b1_row, b3_row))

    print(f"  Plus {len(coupling_edges_123)} coupling edges from (1,2,3) columns")
    all_coupling = coupling_edges + coupling_edges_123
    print(f"  Total coupling edges between B1 and B3: {len(all_coupling)}")

    # Build bipartite incidence
    b1_rows = sorted(set(e[0] for e in all_coupling))
    b3_rows = sorted(set(e[1] for e in all_coupling))
    print(f"  B1 rows involved: {len(b1_rows)}, B3 rows involved: {len(b3_rows)}")

    # Check: are these ALL the rows in blocks 1 and 3?
    # Block 1 has all (r, u, v) with r ∈ NQR, u ∈ QR, v ∈ QR
    # But not all such tuples are valid junk faces (need the face to actually appear)
    # For a face (r, u, v) to be in block 1, need ∃ (s,t) with s+t=r, s,t ∈ QR
    # and the full partial sums valid.

    # The bipartite graph structure:
    # B1 node (r, u, v) connects to B3 node (s, t, r') when
    #   s + t = r, u + v = r', and (s,t,u,v) is a valid 4-diff-seq

    # Key: the transformation (r, u, v) -> (s, t, r') is:
    #   Given r ∈ NQR, u,v ∈ QR with r' = u+v ∈ NQR:
    #   Choose any (s,t) ∈ QR×QR with s+t = r
    #   Then b3_row = (s, t, r')

    # So the coupling factorizes!
    # (r, u, v) connects to (s, t, u+v) for each decomposition s+t = r

    # The bipartite graph has a PRODUCT structure:
    # B1 side: parametrized by (r, u, v) — but we can write as (r, (u,v))
    # B3 side: parametrized by (s, t, r') — but we can write as ((s,t), r')
    # Edge exists when r = s+t and r' = u+v

    # So the bipartite adjacency is:
    #   A[(r,(u,v)), ((s,t),r')] = 1 iff r = s+t AND r' = u+v

    # This is the TENSOR PRODUCT of two bipartite graphs:
    #   G1: r ↔ (s,t) where s+t = r (sum decomposition graph)
    #   G2: (u,v) ↔ r' where u+v = r' (sum decomposition graph)

    # G1 connects NQR value r to QR pairs (s,t) with s+t=r.
    #   N(r) = (m+1)/2 for each r ∈ NQR.
    # G2 connects QR pair (u,v) to NQR value r' = u+v.
    #   Each QR pair (u,v) maps to exactly one value r' (which is in NQR if u+v ∈ NQR).
    #   But wait: not every (u,v) ∈ QR² has u+v ∈ NQR.
    #   #{(u,v) : u+v ∈ NQR} = m(m+1)/2.
    #   Each NQR value r' receives N(r') = (m+1)/2 pairs.

    # The TENSOR product has the property that
    # rank(A ⊗ B) = rank(A) × rank(B)
    # and similarly for the null space.

    # Wait, but the coupling matrix is not exactly A ⊗ B.
    # The rows of the coupling matrix are indexed by B1 = (r, u, v)
    # and columns by "coupling columns" = (s,t,u,v).
    # Let me think of it differently.

    # The key insight: there's a FACTORIZATION structure.
    # Group B1 rows by r-value: for fixed r, the rows are {(r,u,v) : ...}
    # Group B3 rows by r'-value: for fixed r', the rows are {(s,t,r') : ...}

    # For fixed (r, r'), the coupling sub-block connects
    #   {(r, u, v) : u+v = r'} to {(s, t, r') : s+t = r}
    # through columns (s,t,u,v) with both conditions met.

    # Within this sub-block, each column (s,t,u,v) puts -1 at row (r,u,v) in B1
    # and -1 at row (s,t,r') in B3.
    # So the sub-block is: column (s,t,u,v) = -e_{(r,u,v)} in B1 and -e_{(s,t,r')} in B3.
    # The sub-block in the STACKED [B1; B3] matrix has rank = max(#B1_rows, #B3_rows)
    # where #B1_rows = #{(u,v) : u+v = r'} = N(r') = (m+1)/2
    #       #B3_rows = #{(s,t) : s+t = r} = N(r) = (m+1)/2
    # and #columns = #B1_rows × #B3_rows = ((m+1)/2)²

    # Actually in the stacked matrix, each column contributes exactly 2 nonzero entries.
    # This is an m×m bipartite graph (B1_rows ↔ B3_rows) where every pair is connected.
    # The rank of such a complete bipartite matrix is... let me think.

    # The sub-block for (r, r') has:
    # - Rows in B1: (u,v) pairs with u+v = r', indexed 1..N(r')
    # - Rows in B3: (s,t) pairs with s+t = r, indexed 1..N(r)
    # - Columns: all (s,t,u,v) = all pairs, so N(r)×N(r') columns
    # - Column (s_i, t_j, u_k, v_l) has -1 at B1 row k and -1 at B3 row i

    # Stacked: matrix is [I_B1 ⊗ 1_{B3}^T; 1_{B1}^T ⊗ I_B3] (up to signs)
    # Actually each column is -e_{k} (in B1 space) stacked with -e_{i} (in B3 space)
    # The column for (i,k) is: [-e_k; -e_i]
    # Rank of this matrix = rank of [-I_a ⊗ 1_b^T; -1_a^T ⊗ I_b]
    # = a + b - 1 (standard complete bipartite rank)
    # where a = N(r') and b = N(r)

    # So the rank of the (r,r') sub-block is N(r)+N(r')-1 = (m+1)/2 + (m+1)/2 - 1 = m
    # The number of rows is N(r)+N(r') = m+1
    # Nullity per sub-block = 1

    # Number of (r, r') sub-blocks = |NQR|² = m²
    # Total nullity = m² × 1 = m²  !!!

    print(f"\n  TENSOR PRODUCT ANALYSIS:")
    print(f"  Each (r, r') ∈ NQR × NQR gives a sub-block:")
    print(f"    B1 rows: N(r') = (m+1)/2 = {(m+1)//2}")
    print(f"    B3 rows: N(r) = (m+1)/2 = {(m+1)//2}")
    print(f"    Columns: N(r)×N(r') = {((m+1)//2)**2}")
    print(f"    Rank of sub-block = N(r)+N(r')-1 = {m}")
    print(f"    Nullity per sub-block = 1")
    print(f"    Number of sub-blocks: m² = {m**2}")
    print(f"    TOTAL NULLITY = m² = {m**2}  ✓")

    # Verify: check that the sub-blocks are independent
    # (i.e., different (r,r') sub-blocks don't share rows)
    # B1 rows for (r,r') are {(r,u,v) : u+v=r'} — different r' give different rows (fixed r).
    # B3 rows for (r,r') are {(s,t,r') : s+t=r} — different r give different rows (fixed r').
    # So the sub-blocks are DISJOINT in rows. ✓

    # Verify this computationally for small p
    if p <= 11:
        # Build the coupling matrix restricted to B1 and B3
        b1_idx = {r: i for i, r in enumerate(b1_rows)}
        b3_idx = {r: i for i, r in enumerate(b3_rows)}
        n1 = len(b1_rows)
        n3 = len(b3_rows)

        # Coupling matrix: each edge maps to [-e_{b1}; -e_{b3}]
        C_coupling = np.zeros((n1 + n3, len(all_coupling)), dtype=np.float64)
        for j, (b1, b3) in enumerate(all_coupling):
            C_coupling[b1_idx[b1], j] = -1
            C_coupling[n1 + b3_idx[b3], j] = -1

        rank_coupling = np.linalg.matrix_rank(C_coupling)
        nullity_coupling = (n1 + n3) - rank_coupling
        print(f"\n  Coupling matrix verification:")
        print(f"    Rows: {n1} (B1) + {n3} (B3) = {n1+n3}")
        print(f"    Cols: {len(all_coupling)}")
        print(f"    Rank: {rank_coupling}")
        print(f"    Nullity: {nullity_coupling}")
        print(f"    Expected nullity: m² = {m**2}")
        print(f"    Match: {nullity_coupling == m**2}")

        # Also check per-(r,r') sub-block
        NQR_list = sorted(NQR)
        total_nullity = 0
        for r in NQR_list[:3]:  # just check first few
            for rp in NQR_list[:3]:
                # B1 rows with first entry r and u+v = rp
                b1_sub = [(r, u, v) for u in QR_list for v in QR_list
                          if (u+v)%p == rp and (r,u,v) in set(b1_rows)]
                # B3 rows with third entry rp and s+t = r
                b3_sub = [(s, t, rp) for s in QR_list for t in QR_list
                          if (s+t)%p == r and (s,t,rp) in set(b3_rows)]

                if b1_sub and b3_sub:
                    print(f"    Sub-block (r={r}, r'={rp}): "
                          f"|B1|={len(b1_sub)}, |B3|={len(b3_sub)}")

        # Count total sub-block sizes
        total_b1 = 0
        total_b3 = 0
        for r in NQR_list:
            for rp in NQR_list:
                n_b1 = sum(1 for u in QR_list for v in QR_list
                           if (u+v)%p == rp and (r,u,v) in set(b1_rows))
                n_b3 = sum(1 for s in QR_list for t in QR_list
                           if (s+t)%p == r and (s,t,rp) in set(b3_rows))
                total_b1 += n_b1
                total_b3 += n_b3
                total_nullity += (1 if n_b1 > 0 and n_b3 > 0 else 0)

        print(f"    Total non-empty sub-blocks: {total_nullity}")

    return m

for p in [7, 11]:
    print(f"\n{'='*60}")
    coupling_analysis(p)
    print(f"{'='*60}")
