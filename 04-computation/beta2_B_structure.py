#!/usr/bin/env python3
"""beta2_B_structure.py - Deep structural analysis of the unfiltered B matrix

The B matrix maps T' 2-paths to swap cycle coordinates.
For each T' path (a,b,c), the column of B is:
  d_3(v,a,b,c) + d_3(a,b,c,v) restricted to swap-relevant coordinates.

Key structural questions:
1. What is the EXACT rank of B in terms of graph parameters?
2. Why is rank(B) >= swap_dim always?
3. Is there a natural decomposition of B into blocks?

The swap cycle space lives in span{e_{(a,b,v)} - e_{(v,a,b)} : (a,b) in arcs(P->Q)}.
The constraint is zero row sums (for a in P) and zero col sums (for b in Q).

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def all_tournaments_gen(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def analyze_B_detailed(A, n, v):
    """Return detailed B matrix structure for vertex v."""
    d_out = sum(A[v])
    d_in = n - 1 - d_out
    P = [a for a in range(n) if a != v and A[v][a] == 1]  # out-neighbors
    Q = [b for b in range(n) if b != v and A[b][v] == 1]  # in-neighbors

    if not P or not Q:
        return None

    arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
    if len(arcs_PQ) < 2:
        return None

    # Build bipartite constraint matrix C
    m = len(arcs_PQ)
    rows = []
    row_labels = []
    for a in P:
        row = [0] * m
        has_arc = False
        for j, (a2, b2) in enumerate(arcs_PQ):
            if a2 == a:
                row[j] = 1
                has_arc = True
        if has_arc:
            rows.append(row)
            row_labels.append(('P', a))
    for b in Q:
        row = [0] * m
        has_arc = False
        for j, (a2, b2) in enumerate(arcs_PQ):
            if b2 == b:
                row[j] = 1
                has_arc = True
        if has_arc:
            rows.append(row)
            row_labels.append(('Q', b))

    if not rows:
        return None

    C_mat = np.array(rows, dtype=float)
    Sc = np.linalg.svd(C_mat, compute_uv=False)
    rank_C = int(sum(s > 1e-8 for s in Sc))
    ker_dim = m - rank_C
    if ker_dim == 0:
        return None

    # Get T' paths
    others = [x for x in range(n) if x != v]
    B_sub, vlist = get_induced(A, n, others)
    paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
    Tp_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

    # Build B matrix in the "swap arc" coordinates
    # Columns indexed by T' paths, rows indexed by arcs_PQ
    # For T' path (a,b,c):
    #   d_3(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)
    #   d_3(a,b,c,v) = (b,c,v) - (a,c,v) + (a,b,v) - (a,b,c)
    #   Sum = -(v,b,c)+(v,a,c)-(v,a,b)+(b,c,v)-(a,c,v)+(a,b,v)
    # The swap-relevant parts are (x,y,v) and (v,x,y) for (x,y) in arcs_PQ

    arc_idx = {arc: i for i, arc in enumerate(arcs_PQ)}

    # B_swap: rows = arcs_PQ (via coefficient of e_{(a,b,v)} - e_{(v,a,b)})
    # But actually swap cycles use BOTH (a,b,v) and (v,a,b) coordinates
    # Let's project onto the swap coordinates: coefficient of [(a,b,v) - (v,a,b)]

    B_swap = np.zeros((m, len(Tp_orig)))

    for j, (a, b, c) in enumerate(Tp_orig):
        # From d_3(v,a,b,c) + d_3(a,b,c,v):
        # Terms with (x,y,v): +(b,c,v), -(a,c,v), +(a,b,v)
        # Terms with (v,x,y): -(v,b,c), +(v,a,c), -(v,a,b)

        # Coefficient in swap coordinate (x,y) means:
        # coeff of (x,y,v) contributes +coeff to swap coord
        # coeff of (v,x,y) contributes -coeff to swap coord (because swap = (x,y,v)-(v,x,y))

        # (a,b,v): +1 in (a,b,v), so +1 to swap coord (a,b) if (a,b) in arcs_PQ
        if (a, b) in arc_idx:
            B_swap[arc_idx[(a, b)], j] += 1

        # (a,c,v): -1 in (a,c,v), so -1 to swap coord (a,c) if (a,c) in arcs_PQ
        if (a, c) in arc_idx:
            B_swap[arc_idx[(a, c)], j] -= 1

        # (b,c,v): +1 in (b,c,v), so +1 to swap coord (b,c) if (b,c) in arcs_PQ
        if (b, c) in arc_idx:
            B_swap[arc_idx[(b, c)], j] += 1

        # (v,a,b): -1 in (v,a,b), so -(-1)=+1 to swap coord (a,b) if (a,b) in arcs_PQ
        if (a, b) in arc_idx:
            B_swap[arc_idx[(a, b)], j] += 1

        # (v,a,c): +1 in (v,a,c), so -(+1)=-1 to swap coord (a,c) if (a,c) in arcs_PQ
        if (a, c) in arc_idx:
            B_swap[arc_idx[(a, c)], j] -= 1

        # (v,b,c): -1 in (v,b,c), so -(-1)=+1 to swap coord (b,c) if (b,c) in arcs_PQ
        if (b, c) in arc_idx:
            B_swap[arc_idx[(b, c)], j] += 1

    # Wait - let me reconsider. The (a,b) terms get doubled:
    # (a,b,v) gives +1 and (v,a,b) gives -(-1)=+1, total +2
    # (a,c) terms: (a,c,v) gives -1 and (v,a,c) gives -(+1)=-1, total -2
    # (b,c) terms: (b,c,v) gives +1 and (v,b,c) gives -(-1)=+1, total +2
    #
    # So B_swap[arc, T'path(a,b,c)] = 2*delta_{arc=(a,b)} - 2*delta_{arc=(a,c)} + 2*delta_{arc=(b,c)}
    #
    # That's just 2 * [incidence of arc in {(a,b),(b,c)} minus (a,c)]
    # = 2 * (boundary-like structure)

    # Actually let me recompute more carefully.
    B_swap2 = np.zeros((m, len(Tp_orig)))
    for j, (a, b, c) in enumerate(Tp_orig):
        if (a, b) in arc_idx:
            B_swap2[arc_idx[(a, b)], j] += 2   # from (a,b,v) and (v,a,b)
        if (a, c) in arc_idx:
            B_swap2[arc_idx[(a, c)], j] -= 2   # from (a,c,v) and (v,a,c)
        if (b, c) in arc_idx:
            B_swap2[arc_idx[(b, c)], j] += 2   # from (b,c,v) and (v,b,c)

    rank_B_swap = np.linalg.matrix_rank(B_swap2, tol=1e-8)

    # Now compute the image of B_swap restricted to ker(C)
    # swap_dim = dim(ker(C)) since ker_dim = swap_dim (from S42 findings)

    # Also compute B_swap * (something) = constraint satisfaction
    # C * B_swap should be related to something...
    CB = C_mat @ B_swap2

    # Classify T' paths by their relationship to v
    path_types = defaultdict(int)
    for (a, b, c) in Tp_orig:
        a_in_P = a in P
        b_in_P = b in P
        c_in_P = c in P
        a_in_Q = a in Q
        b_in_Q = b in Q
        c_in_Q = c in Q
        ptype = (
            'P' if a_in_P else 'Q',
            'P' if b_in_P else 'Q',
            'P' if c_in_P else 'Q'
        )
        path_types[ptype] += 1

    # Count arcs between P and Q
    PtoQ = sum(1 for a in P for b in Q if A[a][b] == 1)
    QtoP = sum(1 for b in Q for a in P if A[b][a] == 1)
    PtoP = sum(1 for a1 in P for a2 in P if a1 != a2 and A[a1][a2] == 1)
    QtoQ = sum(1 for b1 in Q for b2 in Q if b1 != b2 and A[b1][b2] == 1)

    return {
        'd_out': d_out, 'd_in': d_in,
        'P_size': len(P), 'Q_size': len(Q),
        'arcs_PQ_count': len(arcs_PQ),
        'PtoQ': PtoQ, 'QtoP': QtoP,
        'PtoP': PtoP, 'QtoQ': QtoQ,
        'rank_C': rank_C,
        'ker_dim': ker_dim,
        'num_Tp': len(Tp_orig),
        'rank_B_swap': rank_B_swap,
        'surplus': rank_B_swap - ker_dim,
        'path_types': dict(path_types),
        'CB_rank': np.linalg.matrix_rank(CB, tol=1e-8),
        'CB_norm': np.max(np.abs(CB)),
    }


# ============================================================
# Part 1: n=5 exhaustive detailed analysis
# ============================================================
print("=" * 70)
print("PART 1: B MATRIX STRUCTURE at n=5 (exhaustive)")
print("=" * 70)

n = 5
data5 = []
for A in all_tournaments_gen(n):
    for v in range(n):
        result = analyze_B_detailed(A, n, v)
        if result:
            data5.append(result)

print(f"n=5: {len(data5)} (T,v) cases")

# Summary statistics
surpluses = [d['surplus'] for d in data5]
print(f"  surplus: min={min(surpluses)}, max={max(surpluses)}")
print(f"  rank_B_swap: {sorted(set(d['rank_B_swap'] for d in data5))}")
print(f"  ker_dim: {sorted(set(d['ker_dim'] for d in data5))}")
print(f"  num_Tp: {sorted(set(d['num_Tp'] for d in data5))}")

# Is CB always zero? (Does B map into ker(C)?)
print(f"\n  C*B analysis:")
cb_norms = [d['CB_norm'] for d in data5]
print(f"  max |C*B|: {max(cb_norms):.6f}")
print(f"  Is C*B always zero? {max(cb_norms) < 1e-8}")

# Path type distribution
print(f"\n  T' path types (P/Q classification):")
all_types = defaultdict(int)
for d in data5:
    for k, v2 in d['path_types'].items():
        all_types[k] += v2
for k in sorted(all_types.keys()):
    print(f"    {k}: {all_types[k]}")

# Degree distribution
print(f"\n  By (d_out, d_in):")
for dout in sorted(set(d['d_out'] for d in data5)):
    subset = [d for d in data5 if d['d_out'] == dout]
    din = n - 1 - dout
    surp = [d['surplus'] for d in subset]
    print(f"    ({dout},{din}): {len(subset)} cases, surplus={min(surp)}-{max(surp)}, "
          f"arcs_PQ={sorted(set(d['arcs_PQ_count'] for d in subset))}")


# ============================================================
# Part 2: Test key hypothesis: rank(B_swap) = rank(B_swap restricted to ker(C)) + rank_C?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: RANK DECOMPOSITION ANALYSIS")
print("=" * 70)

n = 5
count_decomp = 0
count_total = 0

for A in all_tournaments_gen(n):
    for v in range(n):
        d_out = sum(A[v])
        d_in = n - 1 - d_out
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        arc_idx = {arc: i for i, arc in enumerate(arcs_PQ)}

        # Build C
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = int(sum(s > 1e-8 for s in Sc))
        ker_dim = m - rank_C
        if ker_dim == 0:
            continue

        # Build B_swap
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

        B_swap = np.zeros((m, len(Tp_orig)))
        for j, (a, b, c) in enumerate(Tp_orig):
            if (a, b) in arc_idx:
                B_swap[arc_idx[(a, b)], j] += 2
            if (a, c) in arc_idx:
                B_swap[arc_idx[(a, c)], j] -= 2
            if (b, c) in arc_idx:
                B_swap[arc_idx[(b, c)], j] += 2

        # B_swap / 2 has integer entries in {-1, 0, +1, +2}
        B_int = B_swap / 2

        rank_B = np.linalg.matrix_rank(B_swap, tol=1e-8)

        # Project B_swap onto ker(C) and im(C^T)
        _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
        ker_basis = Vt[rank_C:]  # rows are ker(C) basis vectors
        im_basis = Vt[:rank_C]   # rows are im(C^T) basis vectors

        # B in ker(C) coordinates
        B_ker = ker_basis @ B_swap  # ker_dim x #Tp
        rank_B_ker = np.linalg.matrix_rank(B_ker, tol=1e-8)

        # B in im(C^T) coordinates
        B_im = im_basis @ B_swap   # rank_C x #Tp
        rank_B_im = np.linalg.matrix_rank(B_im, tol=1e-8)

        count_total += 1
        if rank_B == rank_B_ker + rank_B_im:
            count_decomp += 1
        else:
            if count_total - count_decomp <= 3:
                print(f"  RANK DECOMP FAIL: v={v}, rank_B={rank_B}, "
                      f"rank_B_ker={rank_B_ker}, rank_B_im={rank_B_im}, "
                      f"sum={rank_B_ker + rank_B_im}")

print(f"\nn=5: rank decomposition holds for {count_decomp}/{count_total}")
print(f"  (rank(B) = rank(proj_ker B) + rank(proj_im B))")


# ============================================================
# Part 3: Is C*B = 0? (Does image of B lie in ker(C)?)
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: DOES im(B) SUBSET ker(C)?")
print("=" * 70)

for n in [5, 6]:
    if n <= 6:
        gen = list(all_tournaments_gen(n))

    cb_zero = 0
    cb_nonzero = 0
    total = 0

    for A in gen:
        for v in range(n):
            d_out = sum(A[v])
            P = [a for a in range(n) if a != v and A[v][a] == 1]
            Q = [b for b in range(n) if b != v and A[b][v] == 1]
            arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

            if len(arcs_PQ) < 2:
                continue

            m = len(arcs_PQ)
            arc_idx = {arc: i for i, arc in enumerate(arcs_PQ)}

            rows = []
            for a in P:
                row = [0] * m
                for j, (a2, b2) in enumerate(arcs_PQ):
                    if a2 == a:
                        row[j] = 1
                if any(r != 0 for r in row):
                    rows.append(row)
            for b in Q:
                row = [0] * m
                for j, (a2, b2) in enumerate(arcs_PQ):
                    if b2 == b:
                        row[j] = 1
                if any(r != 0 for r in row):
                    rows.append(row)

            C_mat = np.array(rows, dtype=float)
            Sc = np.linalg.svd(C_mat, compute_uv=False)
            rank_C = int(sum(s > 1e-8 for s in Sc))
            ker_dim = m - rank_C
            if ker_dim == 0:
                continue

            others = [x for x in range(n) if x != v]
            B_sub, vlist = get_induced(A, n, others)
            paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
            Tp_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

            B_swap = np.zeros((m, len(Tp_orig)))
            for j, (a, b, c) in enumerate(Tp_orig):
                if (a, b) in arc_idx:
                    B_swap[arc_idx[(a, b)], j] += 2
                if (a, c) in arc_idx:
                    B_swap[arc_idx[(a, c)], j] -= 2
                if (b, c) in arc_idx:
                    B_swap[arc_idx[(b, c)], j] += 2

            CB = C_mat @ B_swap
            total += 1
            if np.max(np.abs(CB)) < 1e-8:
                cb_zero += 1
            else:
                cb_nonzero += 1
                if cb_nonzero <= 3:
                    print(f"  CB nonzero: n={n}, v={v}, max|CB|={np.max(np.abs(CB)):.4f}")

    print(f"\nn={n}: C*B=0 for {cb_zero}/{total}, nonzero for {cb_nonzero}/{total}")


# ============================================================
# Part 4: Analyze B_swap column structure
# For T' path (a,b,c), the column is:
#   +2*e_{(a,b)} - 2*e_{(a,c)} + 2*e_{(b,c)}  (in arc_PQ coords)
# This is like a "boundary" of the path in the arc_PQ space!
# ============================================================
print(f"\n{'=' * 70}")
print("PART 4: COLUMN STRUCTURE — B as 'boundary' in arc space")
print("=" * 70)

print("""
Each column of B_swap (for T' path (a,b,c)) is:
  +2*delta_{(a,b) in arcs_PQ} - 2*delta_{(a,c) in arcs_PQ} + 2*delta_{(b,c) in arcs_PQ}

This is the BOUNDARY MAP of a 1-chain in the bipartite graph G(P,Q)!

Specifically: if we think of arcs_PQ as edges of a bipartite graph,
then (a,b,c) defines a "2-simplex" with edges (a,b), (a,c), (b,c).
The column is the (signed) boundary of this simplex restricted to arcs_PQ.

KEY INSIGHT: B = 2 * (partial boundary of T' 2-paths in the bipartite arc graph)
""")

# Let's verify: how many of (a,b), (a,c), (b,c) are in arcs_PQ for each T' path?
n = 5
for A in [list(all_tournaments_gen(n))[42]]:  # one example
    for v in [2]:
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
        arc_set = set(arcs_PQ)

        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

        print(f"Example: n={n}, v={v}, |P|={len(P)}, |Q|={len(Q)}, "
              f"|arcs_PQ|={len(arcs_PQ)}")
        print(f"P={P}, Q={Q}")
        print(f"arcs_PQ = {arcs_PQ}")
        print(f"\nT' paths and their arc_PQ memberships:")
        for (a, b, c) in Tp_orig:
            ab = (a, b) in arc_set
            ac = (a, c) in arc_set
            bc = (b, c) in arc_set
            n_hits = sum([ab, ac, bc])
            a_type = 'P' if a in P else 'Q'
            b_type = 'P' if b in P else 'Q'
            c_type = 'P' if c in P else 'Q'
            col = []
            if ab: col.append(f"+e_{(a,b)}")
            if ac: col.append(f"-e_{(a,c)}")
            if bc: col.append(f"+e_{(b,c)}")
            print(f"  ({a},{b},{c}) [{a_type}{b_type}{c_type}]: "
                  f"hits={n_hits}, col={'  '.join(col) if col else '0'}")


# ============================================================
# Part 5: Key test — is B_swap restricted to ker(C) always surjective onto ker(C)?
# i.e., does im(B_swap) contain ker(C)?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 5: Does im(B_swap) contain ker(C)?")
print("=" * 70)

for n in [5, 6]:
    if n <= 6:
        gen = list(all_tournaments_gen(n))

    contains = 0
    not_contains = 0
    total = 0

    for A in gen:
        for v in range(n):
            P = [a for a in range(n) if a != v and A[v][a] == 1]
            Q = [b for b in range(n) if b != v and A[b][v] == 1]
            arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

            if len(arcs_PQ) < 2:
                continue

            m = len(arcs_PQ)
            arc_idx = {arc: i for i, arc in enumerate(arcs_PQ)}

            rows = []
            for a in P:
                row = [0] * m
                for j, (a2, b2) in enumerate(arcs_PQ):
                    if a2 == a:
                        row[j] = 1
                if any(r != 0 for r in row):
                    rows.append(row)
            for b in Q:
                row = [0] * m
                for j, (a2, b2) in enumerate(arcs_PQ):
                    if b2 == b:
                        row[j] = 1
                if any(r != 0 for r in row):
                    rows.append(row)

            C_mat = np.array(rows, dtype=float)
            Sc = np.linalg.svd(C_mat, compute_uv=False)
            rank_C = int(sum(s > 1e-8 for s in Sc))
            ker_dim = m - rank_C
            if ker_dim == 0:
                continue

            others = [x for x in range(n) if x != v]
            B_sub, vlist = get_induced(A, n, others)
            paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
            Tp_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

            B_swap = np.zeros((m, len(Tp_orig)))
            for j, (a, b, c) in enumerate(Tp_orig):
                if (a, b) in arc_idx:
                    B_swap[arc_idx[(a, b)], j] += 2
                if (a, c) in arc_idx:
                    B_swap[arc_idx[(a, c)], j] -= 2
                if (b, c) in arc_idx:
                    B_swap[arc_idx[(b, c)], j] += 2

            # Check if im(B_swap) >= ker(C)
            # This means: for each basis vector of ker(C), it should be in im(B_swap)
            _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
            ker_basis = Vt[rank_C:]  # rows are ker(C) basis

            rank_B = np.linalg.matrix_rank(B_swap, tol=1e-8)

            # Check containment: add ker_basis rows to B_swap rows, see if rank increases
            combined = np.vstack([B_swap, ker_basis])
            rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)

            total += 1
            # im(B_swap) contains ker(C) iff rank_combined = rank(B_swap) (adding ker_basis doesn't increase rank)
            # Wait, this checks if rows of ker_basis are in the row span of B_swap.
            # But we want columns of ker_basis^T in column span of B_swap.
            # Let me redo: we want ker(C) ⊆ im(B_swap) as subsets of R^m (column space)

            rank_B_col = np.linalg.matrix_rank(B_swap, tol=1e-8)
            combined_col = np.hstack([B_swap, ker_basis.T])  # m x (#Tp + ker_dim)
            rank_combined_col = np.linalg.matrix_rank(combined_col, tol=1e-8)

            if rank_combined_col == rank_B_col:
                contains += 1
            else:
                not_contains += 1

    print(f"n={n}: ker(C) subset im(B_swap) for {contains}/{total}")
    if not_contains > 0:
        print(f"  FAILS {not_contains} times!")


print("\n\nDone.")
