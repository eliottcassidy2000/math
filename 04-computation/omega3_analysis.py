#!/usr/bin/env python3
"""
omega3_analysis.py — opus-2026-03-13-S71b

Analyze the constraint structure at degree 3 for Paley tournaments.
A 3-diff-seq (s,t,u) with s,t,u ∈ QR has 4 faces:
  face 0 = (t,u)           — always in A_2 (both entries QR)
  face 1 = ((s+t)%p, u)    — in A_2 iff s+t ∈ QR
  face 2 = (s, (t+u)%p)    — in A_2 iff t+u ∈ QR
  face 3 = (s,t)            — always in A_2

The constraint for (s,t,u) involves:
  +face0 - face1 + face2 - face3 = 0 in A_2
  
Only faces 1 and 2 can fail to be in A_2:
  face 1 junk iff s+t ∉ QR (i.e., s+t ∈ NQR ∪ {0})
  face 2 junk iff t+u ∉ QR (i.e., t+u ∈ NQR ∪ {0})

Cases:
  (a) Both in A_2: no constraint (column is zero)
  (b) Only face 1 junk: constraint has entry -1 at junk row for (s+t, u)
  (c) Only face 2 junk: constraint has entry +1 at junk row for (s, t+u)
  (d) Both junk: constraint has entries -1 at (s+t, u) and +1 at (s, t+u)

Goal: compute rank(C_3) and verify rank(C_3)/|A_3| = 4/p.
"""

import sys
sys.path.insert(0, '04-computation')

def analyze_degree3(p):
    """Analyze degree-3 constraint structure for P_p."""
    QR = set()
    for x in range(1, p):
        QR.add((x*x) % p)
    m = len(QR)
    NQR = set(range(1, p)) - QR
    
    S = sorted(QR)
    
    # Count A_3 (3-diff-seqs with distinct partial sums)
    # Partial sums: 0, s, s+t, s+t+u all distinct mod p
    A3 = []
    for s in S:
        for t in S:
            st = (s + t) % p
            if st == 0:  # s+t = 0 mod p
                continue
            for u in S:
                stu = (st + u) % p
                if stu == 0 or stu == s:  # s+t+u = 0 or s+t+u = s (i.e. t+u=0)
                    continue
                A3.append((s, t, u))
    
    print(f"P_{p}: m={m}, |A_3|/p = {len(A3)}")
    
    # Classify constraints
    case_a = 0  # both faces in A_2
    case_b = 0  # only face 1 junk
    case_c = 0  # only face 2 junk
    case_d = 0  # both faces junk
    
    # Build constraint matrix
    # Junk faces are indexed by (merged_value, remaining_value, which_face)
    # Actually, junk faces are tuples NOT in A_2
    # A tuple (a,b) is in A_2 iff a ∈ QR, b ∈ QR, and a+b ≠ 0 mod p
    # (and partial sums 0, a, a+b distinct, which requires a≠0, b≠0, a+b≠0)
    # Since a,b ∈ QR implies a≠0,b≠0, only need a+b≠0
    # But for p≡3 mod 4, -1 ∉ QR, so a+b=0 iff b=-a ∉ QR. So all QR×QR are in A_2.
    # Wait... but face 1 = ((s+t)%p, u). If s+t ∉ QR, then face 1 ∉ A_2.
    # If s+t ∈ QR, face 1 = (s+t, u) with both in QR → in A_2. ✓
    
    junk_faces = {}  # maps junk face tuple → index
    columns = []
    
    for (s, t, u) in A3:
        st = (s + t) % p
        tu = (t + u) % p
        
        f1_junk = st not in QR  # s+t ∉ QR
        f2_junk = tu not in QR  # t+u ∉ QR
        
        if not f1_junk and not f2_junk:
            case_a += 1
            columns.append({})
            continue
        
        col = {}
        
        if f1_junk:
            # face 1 = (st, u), sign = (-1)^1 = -1
            jf = ('f1', st, u)
            if jf not in junk_faces:
                junk_faces[jf] = len(junk_faces)
            col[junk_faces[jf]] = -1
        
        if f2_junk:
            # face 2 = (s, tu), sign = (-1)^2 = +1
            jf = ('f2', s, tu)
            if jf not in junk_faces:
                junk_faces[jf] = len(junk_faces)
            col[junk_faces[jf]] = 1
        
        if f1_junk and not f2_junk:
            case_b += 1
        elif not f1_junk and f2_junk:
            case_c += 1
        else:
            case_d += 1
        
        columns.append(col)
    
    print(f"  Case (a) both OK:       {case_a}")
    print(f"  Case (b) only f1 junk:  {case_b}")
    print(f"  Case (c) only f2 junk:  {case_c}")
    print(f"  Case (d) both junk:     {case_d}")
    print(f"  Total junk rows:        {len(junk_faces)}")
    
    # Now compute rank
    # Use Gaussian elimination mod a prime
    PRIME = 191
    n_rows = len(junk_faces)
    n_cols = len(columns)
    
    # Build dense matrix for small cases
    if n_rows < 5000:
        import numpy as np
        M = np.zeros((n_rows, n_cols), dtype=int)
        for j, col in enumerate(columns):
            for r, v in col.items():
                M[r, j] = v % PRIME
        
        # Gaussian elimination mod PRIME
        rank = 0
        used_rows = set()
        for j in range(n_cols):
            # Find pivot
            pivot = -1
            for r in range(n_rows):
                if r not in used_rows and M[r, j] % PRIME != 0:
                    pivot = r
                    break
            if pivot == -1:
                continue
            rank += 1
            used_rows.add(pivot)
            # Eliminate
            inv = pow(int(M[pivot, j]), PRIME - 2, PRIME)
            for r in range(n_rows):
                if r != pivot and M[r, j] % PRIME != 0:
                    factor = (M[r, j] * inv) % PRIME
                    for c in range(n_cols):
                        M[r, c] = (M[r, c] - factor * M[pivot, c]) % PRIME
        
        print(f"  rank(C_3) = {rank}")
        print(f"  rank/|A_3| = {rank}/{len(A3)} = {rank/len(A3):.6f}")
        print(f"  4/p = {4/p:.6f}")
        print(f"  Match: {abs(rank/len(A3) - 4/p) < 1e-10}")
        
        # Analyze the junk face structure more carefully
        # How many junk rows come from face 1 vs face 2?
        f1_rows = sum(1 for k in junk_faces if k[0] == 'f1')
        f2_rows = sum(1 for k in junk_faces if k[0] == 'f2')
        print(f"\n  Junk rows from face 1: {f1_rows}")
        print(f"  Junk rows from face 2: {f2_rows}")
        
        # Analyze N(r) for face 1: #{(s,t,u) : s+t = r, r ∉ QR}
        # For fixed r ∉ QR: #{s ∈ QR : r-s ∈ QR} × m (choices for u)
        # The #{s ∈ QR : r-s ∈ QR} = N(r) from degree-2 analysis
        print(f"\n  Face 1 junk detail:")
        for r in sorted(NQR):
            count = sum(1 for (s,t,u) in A3 if (s+t)%p == r)
            n_r = sum(1 for s in QR if (r - s) % p in QR)
            # But we also need (s,t,u) to be a valid 3-diff-seq
            # i.e., partial sums 0, s, s+t=r, s+t+u distinct
            # s ≠ 0 ✓ (QR), r ≠ 0 ✓ (NQR), s+t+u ≠ 0 and ≠ s
            print(f"    r={r}: count={count}, N(r)={n_r}, N(r)×m_adj=?")
        
        # Count case (d) more carefully
        print(f"\n  Case (d) analysis: both s+t ∉ QR AND t+u ∉ QR")
        # For (s,t,u): s+t ∈ NQR∪{0} AND t+u ∈ NQR∪{0}
        # s+t = 0 impossible (since -s ∉ QR for p ≡ 3 mod 4)
        # t+u = 0 impossible (same reason)
        # So: s+t ∈ NQR AND t+u ∈ NQR
        
        return rank, len(A3)


for p in [7, 11, 19, 23]:
    r, a = analyze_degree3(p)
    print(f"\n  Omega_3/p = {a - r}")
    print(f"  Expected: see known values")
    print()
    print("=" * 60)
    print()
