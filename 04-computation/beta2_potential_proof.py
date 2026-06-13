#!/usr/bin/env python3
r"""beta2_potential_proof.py - Vertex potential approach to b1 <= 1

KEY IDEA: The TT cocycle condition z(a,b) + z(b,c) = z(a,c) for TT triples
means that for each vertex a, z(a, .) restricted to N+(a) has "potential"
structure: z(b,c) = z(a,c) - z(a,b) for all b,c in N+(a) with b->c.

This means: the out-edges from a DETERMINE all edges within N+(a).
Similarly, the in-edges to a determine edges within N-(a).

Combined with the 1-cycle condition (kirchhoff at each vertex), this
is highly constrained. Let's count degrees of freedom.

For a tournament on n vertices:
- n(n-1)/2 directed edges
- dim(Z_1) = (n-1)(n-2)/2 (number of edges minus rank of d_1 = n-1)
- dim(cocycle space) = dim(Z_1) - rk(TT boundaries in Z_1) = b1

POTENTIAL APPROACH:
Define phi_a(b) = z(a,b) for b in N+(a).
TT condition: z(b,c) = phi_a(c) - phi_a(b) for b,c in N+(a) with b->c.
So z on edges WITHIN N+(a) is determined by phi_a.

But z(b,c) is also constrained by phi_b (for c in N+(b)),
and by phi_c (for b in N-(c)), etc.

CONSISTENCY: phi_a(c) - phi_a(b) = phi_b(c) - phi_b(b)... wait,
phi_b(b) doesn't exist. Let me be more careful.

For vertex a with N+(a) = {b1,...,bk}:
  z(a, bi) = phi_a(bi)
  z(bi, bj) = phi_a(bj) - phi_a(bi)  if bi->bj (both in N+(a))

For vertex bi with N+(bi) = {...}:
  z(bi, bj) = phi_{bi}(bj)  if bj in N+(bi)
  z(bj, w) = phi_{bi}(w) - phi_{bi}(bj)  if bj->w, bj,w in N+(bi)

From both: phi_{bi}(bj) = phi_a(bj) - phi_a(bi)  (for bj in N+(a) cap N+(bi))

This gives: phi_{bi}(bj) = phi_a(bj) - phi_a(bi)

Similarly if c in N+(a) cap N+(bi), then:
  phi_{bi}(c) = phi_a(c) - phi_a(bi)

So phi_bi is DETERMINED on N+(a) cap N+(bi) by phi_a.

This propagation could make phi globally determined up to 1 DOF.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def compute_cocycle_via_potential(A, n):
    """Express the TT-cocycle condition as a potential consistency system."""
    # Edges
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    num_edges = len(edges)

    # Out-neighborhoods
    Nplus = [[] for _ in range(n)]
    Nminus = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                Nplus[i].append(j)
                Nminus[j].append(i)

    # The potential system:
    # Variables: z(a,b) for each edge a->b (num_edges variables)
    # Constraints:
    #   1. TT: z(a,b) + z(b,c) = z(a,c) for each TT triple
    #   2. 1-cycle (Kirchhoff): sum_{u: u->v} z(u,v) = sum_{w: v->w} z(v,w) for each v

    # TT constraints
    tt_constraints = []
    for a in range(n):
        for b in Nplus[a]:
            for c in Nplus[b]:
                if A[a][c]:  # TT: a->b->c, a->c
                    # z(a,b) + z(b,c) - z(a,c) = 0
                    row = np.zeros(num_edges)
                    row[edge_idx[(a,b)]] = 1
                    row[edge_idx[(b,c)]] = 1
                    row[edge_idx[(a,c)]] = -1
                    tt_constraints.append(row)

    # Kirchhoff constraints
    kirchhoff = []
    for v in range(n):
        row = np.zeros(num_edges)
        for u in Nminus[v]:
            row[edge_idx[(u,v)]] = 1
        for w in Nplus[v]:
            row[edge_idx[(v,w)]] = -1
        kirchhoff.append(row)

    # Combined constraint matrix
    if tt_constraints:
        C_tt = np.array(tt_constraints)
    else:
        C_tt = np.zeros((0, num_edges))

    C_kirch = np.array(kirchhoff)
    C_all = np.vstack([C_tt, C_kirch])

    # Rank of combined system
    sv = np.linalg.svd(C_all, compute_uv=False)
    rk_all = int(sum(s > 1e-8 for s in sv))

    # Null space dimension = b1 (+ dim of B_1 space)
    # Actually the null space of [TT + Kirchhoff] is ker(TT) cap ker(Kirchhoff)
    # = (cycles that are TT-cocycles)
    # = b1 + dim(boundaries that are TT-cocycles)
    # Hmm, boundaries ARE TT boundaries so they're automatically TT.
    # Actually no: B_1 = im(d_0) which is trivial (d_0: R^{n} -> R^{edges} by d(v) = sum of out-edges - sum of in-edges)
    # Wait, that's d_1. B_1 is not relevant here since we're looking at H_1.

    # Let me reconsider. We want dim of:
    # {z in R^{edges} : z in Z_1 AND z perp TT boundaries}
    # = {z : Kirchhoff(z)=0 AND TT-cocycle(z)=0}
    # = null(C_all)

    dim_null = num_edges - rk_all

    # But this includes the boundary space B_1.
    # B_1 = im(d_2), and we showed TT boundaries span im(d_2).
    # Wait no, B_1 in path homology is im(d_2: Omega_2 -> Omega_1).
    # And Z_1 = ker(d_1: Omega_1 -> Omega_0).
    # b_1 = dim(Z_1/B_1) = dim(Z_1) - dim(B_1).

    # The null space of C_all gives Z_1 cap TT-cocycle-space.
    # Since TT boundaries span all of im(d_2) = B_1,
    # the TT-cocycle-space in Z_1 is exactly H_1 = Z_1/B_1.
    # Wait, that's not right either. Let me think again.

    # Z_1 = {z in R^{edges} : Kirchhoff = 0}. dim = (n-1)(n-2)/2.
    # TT-cocycle condition: z perp im(M_TT).
    # Since im(M_TT) = im(d_2) = B_1, this gives:
    # z in Z_1 AND z perp B_1 = complement of B_1 in Z_1.
    # This is the "harmonic" representative of H_1 (using the standard inner product).
    # Its dimension = b_1 = dim(H_1).

    # So null(C_all) should have dimension = b_1.
    # But C_all includes Kirchhoff AND TT constraints.

    # Actually: C_tt * z = 0 means z is in the "TT kernel".
    # C_kirch * z = 0 means z is a 1-cycle.
    # Together: z is a 1-cycle that's also a TT-cocycle.
    # Since TT boundaries span B_1, and the cocycle condition is
    # z perp TT-boundaries, we have: z perp B_1 AND z in Z_1.
    # This is exactly the harmonic space, dim = b_1.

    # So dim_null should equal b_1!

    # But wait, I need to also account for the Omega_1 structure.
    # In GLMY path homology, the chain space is Omega_1 (not all of R^{edges}).
    # Omega_1 is a SUBSPACE of R^{edges}.

    # Hmm, the Omega_1 constraint is that for each non-allowed 1-path,
    # the coefficient is 0. But in a tournament, ALL edges are allowed 1-paths.
    # So Omega_1 = R^{edges} for a tournament!

    # OK good. So dim_null = b_1.

    return dim_null, rk_all, num_edges, len(tt_constraints), n


# ============================================================
# Part 1: Verify dim_null = b_1
# ============================================================
print("=" * 70)
print("VERIFY: null(TT + Kirchhoff) = b_1")
print("=" * 70)

for n in [3, 4, 5]:
    total = 2 ** (n*(n-1)//2)
    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        dim_null, _, num_e, num_tt, _ = compute_cocycle_via_potential(A, n)

        # Compare with direct b1 computation
        from beta2_good_vertex_final import compute_b1
        b1 = compute_b1(A, n)

        if dim_null != b1:
            print(f"  MISMATCH at n={n}, bits={bits}: null={dim_null}, b1={b1}")
            break
    else:
        b1_dist = Counter()
        for bits in range(total):
            A = [[0] * n for _ in range(n)]
            idx = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (bits >> idx) & 1:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    idx += 1
            dim_null, _, _, _, _ = compute_cocycle_via_potential(A, n)
            b1_dist[dim_null] += 1
        print(f"  n={n}: VERIFIED. b1 distribution: {dict(sorted(b1_dist.items()))}")


# ============================================================
# Part 2: How many TT constraints are independent?
# ============================================================
print(f"\n{'=' * 70}")
print("TT CONSTRAINT RANK vs KIRCHHOFF RANK")
print("=" * 70)

for n in [4, 5]:
    total = 2 ** (n*(n-1)//2)
    tt_rk_dist = Counter()
    kirch_rk = n - 1  # always

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        # Compute TT constraint matrix
        edges = []
        for i in range(n):
            for j in range(n):
                if A[i][j]:
                    edges.append((i, j))
        edge_idx = {e: k for k, e in enumerate(edges)}
        num_edges = len(edges)

        tt_rows = []
        for a in range(n):
            for b in range(n):
                if b == a or not A[a][b]: continue
                for c in range(n):
                    if c == a or c == b or not A[b][c]: continue
                    if A[a][c]:
                        row = np.zeros(num_edges)
                        row[edge_idx[(a,b)]] = 1
                        row[edge_idx[(b,c)]] = 1
                        row[edge_idx[(a,c)]] = -1
                        tt_rows.append(row)

        if tt_rows:
            C_tt = np.array(tt_rows)
            sv = np.linalg.svd(C_tt, compute_uv=False)
            rk_tt = int(sum(s > 1e-8 for s in sv))
        else:
            rk_tt = 0

        # Combined rank
        C_kirch = np.zeros((n, num_edges))
        for v in range(n):
            for j in range(n):
                if A[j][v]:
                    C_kirch[v, edge_idx[(j,v)]] = 1
                if A[v][j]:
                    C_kirch[v, edge_idx[(v,j)]] = -1

        C_all = np.vstack([C_tt, C_kirch]) if tt_rows else C_kirch
        sv_all = np.linalg.svd(C_all, compute_uv=False)
        rk_all = int(sum(s > 1e-8 for s in sv_all))

        sv_k = np.linalg.svd(C_kirch, compute_uv=False)
        rk_k = int(sum(s > 1e-8 for s in sv_k))

        c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]))
        num_tt = len(tt_rows)
        dim_null = num_edges - rk_all

        tt_rk_dist[(c3, rk_tt, rk_k, rk_all, dim_null)] += 1

    print(f"\nn={n}: num_edges={n*(n-1)//2}")
    print(f"  (c3, rk_TT, rk_Kirch, rk_combined, null_dim) => count")
    for key, cnt in sorted(tt_rk_dist.items()):
        print(f"  c3={key[0]}, rk_TT={key[1]}, rk_K={key[2]}, rk_all={key[3]}, b1={key[4]}: {cnt}")


# ============================================================
# Part 3: Potential propagation analysis
# ============================================================
print(f"\n{'=' * 70}")
print("POTENTIAL PROPAGATION")
print("=" * 70)

# For each vertex a, out-edges determine edges in N+(a).
# How many edge variables are DETERMINED by choosing out-edges from each vertex?
#
# Out-edges from v: (n-1)/2 edges (for regular tournament)
# Total out-edges: n * (n-1)/2 = n(n-1)/2 = total edges (each edge counted once as out-edge of source)
# So out-edges from ALL vertices = ALL edges. No new info from potential.
#
# BUT: edges WITHIN N+(a) are determined by potential:
# z(b,c) = z(a,c) - z(a,b) for b,c in N+(a) with b->c
# So the WITHIN-N+(a) edges are NOT free variables.

# Let me count:
# Variables: z(a,b) for each edge. Total: n(n-1)/2.
# Constraints from vertex a's out-potential:
#   For each edge (b,c) with b,c in N+(a): z(b,c) = z(a,c) - z(a,b).
#   Number of such constraints: #edges within N+(a).

# The "free" variables after TT constraints are:
# out-edges from each vertex minus overlaps
# This is tricky because an edge (b,c) might be constrained by MULTIPLE vertices a.

# Let me just count: how many edges in N+(a) for each a?
for n in [5, 7]:
    random.seed(42)
    if n == 5:
        A = random_tournament(n)
    else:
        A = random_tournament(n)

    print(f"\nn={n}:")
    total_within = 0
    for a in range(n):
        Np = [j for j in range(n) if A[a][j]]
        within = sum(1 for b in Np for c in Np if b != c and A[b][c])
        total_within += within
        print(f"  v={a}: |N+(v)|={len(Np)}, #edges within N+(v)={within}")

    print(f"  Total within-constraints: {total_within}")
    print(f"  Total edges: {n*(n-1)//2}")
    print(f"  These constraints have overlaps (edge in multiple N+(a))")


# ============================================================
# Part 4: Direct DOF analysis with potential variables
# ============================================================
print(f"\n{'=' * 70}")
print("DEGREES OF FREEDOM: POTENTIAL PARAMETERIZATION")
print("=" * 70)

# Alternative parameterization: instead of z(a,b) for each edge,
# define phi(a) for each vertex (n variables), plus a "cycle correction" sigma.
# Then z(a,b) = phi(b) - phi(a) + sigma * delta(a,b)
# where delta(a,b) = 1 if (a,b) is a 3-cycle edge, 0 otherwise.
#
# This doesn't quite work because the cycle-edge indicator is not well-defined.
# But the TT-cocycle IS a coboundary + cycle-correction.
#
# Let me test: can the cocycle be written as phi(b) - phi(a) + sigma * f(a,b)
# where f is some explicit function on edges?

for n in [4, 5]:
    total = 2 ** (n*(n-1)//2)
    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        r = compute_cocycle_via_potential(A, n)
        if r[0] != 1:  # b1 != 1
            continue

        # Compute the cocycle explicitly
        edges = []
        for i in range(n):
            for j in range(n):
                if A[i][j]:
                    edges.append((i, j))
        edge_idx = {e: k for k, e in enumerate(edges)}
        num_edges = len(edges)

        # Build full constraint matrix
        tt_rows = []
        for a in range(n):
            for b in range(n):
                if b == a or not A[a][b]: continue
                for c in range(n):
                    if c == a or c == b or not A[b][c]: continue
                    if A[a][c]:
                        row = np.zeros(num_edges)
                        row[edge_idx[(a,b)]] = 1
                        row[edge_idx[(b,c)]] = 1
                        row[edge_idx[(a,c)]] = -1
                        tt_rows.append(row)

        C_kirch = np.zeros((n, num_edges))
        for v in range(n):
            for j in range(n):
                if A[j][v]:
                    C_kirch[v, edge_idx[(j,v)]] = 1
                if A[v][j]:
                    C_kirch[v, edge_idx[(v,j)]] = -1

        if tt_rows:
            C_all = np.vstack([np.array(tt_rows), C_kirch])
        else:
            C_all = C_kirch

        U, S, Vt = np.linalg.svd(C_all, full_matrices=True)
        rk = int(sum(s > 1e-8 for s in S))
        z = Vt[rk]  # null vector = cocycle

        # Try to decompose: z(a,b) = phi(b) - phi(a) + sigma * g(a,b)
        # where g(a,b) is related to 3-cycle membership
        #
        # For TT edge (a,b) with a->b, a->c, b->c: z(a,b) = phi(b)-phi(a)
        # For 3-cycle edge (a,b) with a->b->c->a: z(a,b) = phi(b)-phi(a) + sigma*g(a,b)
        #
        # From TT constraint: z(a,b)+z(b,c)-z(a,c) = 0 on TT triples.
        # If z = coboundary + correction:
        #   [phi(b)-phi(a)+sigma*g(a,b)] + [phi(c)-phi(b)+sigma*g(b,c)] - [phi(c)-phi(a)+sigma*g(a,c)] = 0
        #   sigma * [g(a,b) + g(b,c) - g(a,c)] = 0
        #   So g must be a cocycle on TT triples too: g(a,b)+g(b,c)=g(a,c) for TT.

        # The simplest g: count the number of 3-cycles containing (a,b).
        # Or: g(a,b) = 1/(n-2) * (number of 3-cycles containing edge a->b)

        # Let me compute the 3-cycle participation for each edge
        cyc_count = np.zeros(num_edges)
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if A[i][j] and A[j][k] and A[k][i]:
                        cyc_count[edge_idx[(i,j)]] += 1
                        cyc_count[edge_idx[(j,k)]] += 1
                        cyc_count[edge_idx[(k,i)]] += 1
                    elif A[i][k] and A[k][j] and A[j][i]:
                        cyc_count[edge_idx[(i,k)]] += 1
                        cyc_count[edge_idx[(k,j)]] += 1
                        cyc_count[edge_idx[(j,i)]] += 1

        # Check if z is proportional to some simple function of cyc_count
        # z vs cyc_count correlation?
        z_norm = z / np.max(np.abs(z))
        cc_norm = cyc_count / max(np.max(cyc_count), 1)

        # Check if z = alpha * cyc_count + beta * ones
        # (Least squares: z = a * cyc_count + b * 1 + residual)
        X = np.column_stack([cyc_count, np.ones(num_edges)])
        coeff, res, _, _ = np.linalg.lstsq(X, z, rcond=None)
        pred = X @ coeff
        residual = np.max(np.abs(z - pred))

        if bits <= 30 or residual < 1e-6:
            print(f"\n  n={n}, bits={bits}:")
            print(f"    z = {coeff[0]:.6f} * cyc_count + {coeff[1]:.6f} * 1")
            print(f"    Max residual: {residual:.2e}")
            if residual < 1e-6:
                print(f"    EXACT FIT!")

        # Check z vs coboundary + cycle_count
        # Try: z(a,b) = phi(b) - phi(a)
        # Least squares for phi
        D1 = np.zeros((num_edges, n))
        for k, (a, b) in enumerate(edges):
            D1[k, b] = 1
            D1[k, a] = -1
        phi, _, _, _ = np.linalg.lstsq(D1, z, rcond=None)
        coboundary = D1 @ phi
        correction = z - coboundary

        # Is correction proportional to cyc_count?
        if np.max(np.abs(correction)) > 1e-10 and np.max(cyc_count) > 0:
            X2 = cyc_count.reshape(-1, 1)
            coeff2, _, _, _ = np.linalg.lstsq(X2, correction, rcond=None)
            pred2 = X2.flatten() * coeff2[0]
            res2 = np.max(np.abs(correction - pred2))
            if bits <= 30 or res2 < 1e-6:
                print(f"    Correction = {coeff2[0]:.6f} * cyc_count, residual={res2:.2e}")

        if bits > 30 and residual > 1e-6:
            continue

    break  # Only n=4 for detailed output


# ============================================================
# Part 5: Systematic check: z = alpha * cyc_count + beta
# ============================================================
print(f"\n{'=' * 70}")
print("z = alpha * cyc_count + beta ?")
print("=" * 70)

for n in [4, 5]:
    total = 2 ** (n*(n-1)//2)
    exact_fit = 0
    total_b1_1 = 0

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        r = compute_cocycle_via_potential(A, n)
        if r[0] != 1:
            continue
        total_b1_1 += 1

        edges = []
        for i in range(n):
            for j in range(n):
                if A[i][j]:
                    edges.append((i, j))
        edge_idx = {e: k for k, e in enumerate(edges)}
        num_edges = len(edges)

        tt_rows = []
        for a in range(n):
            for b in range(n):
                if b == a or not A[a][b]: continue
                for c in range(n):
                    if c == a or c == b or not A[b][c]: continue
                    if A[a][c]:
                        row = np.zeros(num_edges)
                        row[edge_idx[(a,b)]] = 1
                        row[edge_idx[(b,c)]] = 1
                        row[edge_idx[(a,c)]] = -1
                        tt_rows.append(row)

        C_kirch = np.zeros((n, num_edges))
        for v in range(n):
            for j in range(n):
                if A[j][v]:
                    C_kirch[v, edge_idx[(j,v)]] = 1
                if A[v][j]:
                    C_kirch[v, edge_idx[(v,j)]] = -1

        if tt_rows:
            C_all = np.vstack([np.array(tt_rows), C_kirch])
        else:
            C_all = C_kirch

        U, S, Vt = np.linalg.svd(C_all, full_matrices=True)
        rk = int(sum(s > 1e-8 for s in S))
        z = Vt[rk]

        # Cycle count per edge
        cyc_count = np.zeros(num_edges)
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if A[i][j] and A[j][k] and A[k][i]:
                        cyc_count[edge_idx[(i,j)]] += 1
                        cyc_count[edge_idx[(j,k)]] += 1
                        cyc_count[edge_idx[(k,i)]] += 1
                    elif A[i][k] and A[k][j] and A[j][i]:
                        cyc_count[edge_idx[(i,k)]] += 1
                        cyc_count[edge_idx[(k,j)]] += 1
                        cyc_count[edge_idx[(j,i)]] += 1

        X = np.column_stack([cyc_count, np.ones(num_edges)])
        coeff, _, _, _ = np.linalg.lstsq(X, z, rcond=None)
        pred = X @ coeff
        residual = np.max(np.abs(z - pred))
        if residual < 1e-6:
            exact_fit += 1

    print(f"  n={n}: {exact_fit}/{total_b1_1} b1=1 tournaments have z = a*cyc_count + b*1")


# ============================================================
# Part 6: What IS the cocycle? Try more basis functions
# ============================================================
print(f"\n{'=' * 70}")
print("COCYCLE DECOMPOSITION: MORE BASIS FUNCTIONS")
print("=" * 70)

# Besides cyc_count, try:
# - out_degree of source: deg+(a)
# - in_degree of target: deg-(b)
# - number of common out-neighbors of a,b
# - score difference: score(b) - score(a)

for n in [5]:
    total = 2 ** (n*(n-1)//2)
    fits = Counter()

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        r = compute_cocycle_via_potential(A, n)
        if r[0] != 1:
            continue

        edges = []
        for i in range(n):
            for j in range(n):
                if A[i][j]:
                    edges.append((i, j))
        edge_idx = {e: k for k, e in enumerate(edges)}
        num_edges = len(edges)

        # Build constraint matrix
        tt_rows = []
        for a in range(n):
            for b in range(n):
                if b == a or not A[a][b]: continue
                for c in range(n):
                    if c == a or c == b or not A[b][c]: continue
                    if A[a][c]:
                        row = np.zeros(num_edges)
                        row[edge_idx[(a,b)]] = 1
                        row[edge_idx[(b,c)]] = 1
                        row[edge_idx[(a,c)]] = -1
                        tt_rows.append(row)

        C_kirch = np.zeros((n, num_edges))
        for v in range(n):
            for j in range(n):
                if A[j][v]:
                    C_kirch[v, edge_idx[(j,v)]] = 1
                if A[v][j]:
                    C_kirch[v, edge_idx[(v,j)]] = -1

        C_all = np.vstack([np.array(tt_rows), C_kirch]) if tt_rows else C_kirch
        U, S, Vt = np.linalg.svd(C_all, full_matrices=True)
        rk = int(sum(s > 1e-8 for s in S))
        z = Vt[rk]

        scores = [sum(A[i]) for i in range(n)]

        # Build feature matrix for edges
        features = []
        for (a, b) in edges:
            cyc = 0
            for k in range(n):
                if k == a or k == b:
                    continue
                if A[a][b] and A[b][k] and A[k][a]:
                    cyc += 1
                if A[a][k] and A[k][b] and A[b][a]:  # can't happen (A[a][b]=1)
                    pass
                # Actually: directed 3-cycles containing edge a->b
                # Need a->b, then b->c->a for some c: A[b][c] and A[c][a]
                # OR: d->a->b for some d: A[d][a] and then A[b][d]
                # Hmm, let me just count properly
            cyc_ab = 0
            for k in range(n):
                if k == a or k == b:
                    continue
                if A[b][k] and A[k][a]:  # a->b->k->a
                    cyc_ab += 1
            features.append([
                1,  # constant
                cyc_ab,  # 3-cycles through this edge
                scores[a],  # out-degree of source
                scores[b],  # out-degree of target
                scores[b] - scores[a],  # score difference
            ])

        F = np.array(features)
        coeff, res, rank_F, _ = np.linalg.lstsq(F, z, rcond=None)
        pred = F @ coeff
        residual = np.max(np.abs(z - pred))

        if residual < 1e-6:
            fits['exact'] += 1
        else:
            fits['inexact'] += 1

    print(f"  n=5: {fits.get('exact', 0)}/{fits.get('exact', 0) + fits.get('inexact', 0)} exact fits with [1, cyc, deg+, deg-, diff]")


# ============================================================
# Part 7: Understand cocycle at n=5, c3=4 disconnected case
# ============================================================
print(f"\n{'=' * 70}")
print("DISCONNECTED CASE ANALYSIS (n=5, c3=4)")
print("=" * 70)

n = 5
for bits in [40]:  # known disconnected case
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    print(f"\nbits={bits}:")
    print(f"Adjacency:")
    for i in range(n):
        row = []
        for j in range(n):
            row.append(str(A[i][j]))
        print(f"  {' '.join(row)}")

    scores = [sum(A[i]) for i in range(n)]
    print(f"Scores: {scores}")

    # Find 3-cycles
    from beta2_cocycle_algebra import get_3cycles
    cycles = get_3cycles(A, n)
    print(f"3-cycles: {cycles}")

    # Find edges of each cycle
    for (i, j, k) in cycles:
        es = []
        if A[i][j]: es.append(f"{i}->{j}")
        if A[j][k]: es.append(f"{j}->{k}")
        if A[k][i]: es.append(f"{k}->{i}")
        if A[j][i]: es.append(f"{j}->{i}")
        if A[k][j]: es.append(f"{k}->{j}")
        if A[i][k]: es.append(f"{i}->{k}")
        print(f"  Cycle {(i,j,k)}: {', '.join(es)}")

    # Compute cocycle
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edges.append((i, j))
    edge_idx = {e: k for k, e in enumerate(edges)}
    num_edges = len(edges)

    tt_rows = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                if A[a][c]:
                    row = np.zeros(num_edges)
                    row[edge_idx[(a,b)]] = 1
                    row[edge_idx[(b,c)]] = 1
                    row[edge_idx[(a,c)]] = -1
                    tt_rows.append(row)

    C_kirch = np.zeros((n, num_edges))
    for v in range(n):
        for j in range(n):
            if A[j][v]:
                C_kirch[v, edge_idx[(j,v)]] = 1
            if A[v][j]:
                C_kirch[v, edge_idx[(v,j)]] = -1

    C_all = np.vstack([np.array(tt_rows), C_kirch])
    U, S, Vt = np.linalg.svd(C_all, full_matrices=True)
    rk = int(sum(s > 1e-8 for s in S))
    z = Vt[rk]
    z_norm = z / np.max(np.abs(z))

    print(f"\nCocycle (normalized):")
    for k, e in enumerate(edges):
        if abs(z_norm[k]) > 1e-8:
            print(f"  {e[0]}->{e[1]}: {z_norm[k]:.6f}")

    # Kirchhoff check
    print(f"\nKirchhoff check:")
    for v in range(n):
        in_sum = sum(z_norm[edge_idx[(u,v)]] for u in range(n) if A[u][v])
        out_sum = sum(z_norm[edge_idx[(v,w)]] for w in range(n) if A[v][w])
        print(f"  v={v}: in={in_sum:.6f}, out={out_sum:.6f}, diff={in_sum-out_sum:.2e}")


print("\n\nDone.")
