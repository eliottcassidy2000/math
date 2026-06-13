"""
omega2_exact_formula.py — Verify exact formula for dim(Omega_2) and filling ratio.

PROVEN (THM-119): dim(Omega_2) = |A_2| - #NA_faces
where #NA_faces = #{(a,c) : c->a in T, and exists b with a->b->c in T}

Key formulas to verify:
  |A_2| = C(n,3) + 2*c3    (HYP-268, already confirmed)
  #NA = #{backward pairs in support of A^2}
  isolated = #{(a,c): c->a, A^2[a,c]=0} = backward pairs NOT in any 3-cycle
  #NA = C(n,2) - n*(n-1)/2 + #{forward pairs in A^2>0}  ... need to work this out

From THM-119:
  dim(Omega_2) = |A_2| - #NA
  = C(n,3) + 2*c3 - #NA

We want: #NA as a function of simple tournament invariants.

Let B = #{backward pairs (a,c) with c->a AND A^2[a,c] > 0}
    = #{(a,c): c->a, exists b: a->b, b->c}
    = sum over (a,c) with c->a of 1_{A^2[a,c]>0}

Note: A^2[a,c] = #{b: a->b->c} = out-neighborhood of a intersected with in-neighborhood of c.

For a backward pair (a,c) with c->a:
A^2[a,c] > 0 iff there exists a path a->b->c, i.e., (a,c) sits in some 3-cycle a->b->c->a.

So: #NA = #{backward pairs in some 3-cycle}
    = #{edges of T participating in 3-cycles, counted as backward relative to some ordering}

Actually: #NA = #{(a,c) ordered: c->a AND a,c share a 3-cycle}
         = total backward pairs - isolated backward pairs
         = C(n,2) - #{forward pairs} - #{isolated backward}
Wait: C(n,2) total unordered pairs. Each pair has one direction.
Let's define:
  E = n(n-1)/2 = number of edges (unordered)
  Forward pairs = |A_1| = E (all edges are allowed 1-paths in both directions)

No. Let me be more careful. For ordered pair (a,c) with a<c (say):
  Either a->c (forward) or c->a (backward).

For the NA count, we need ordered pairs (a,c) (not necessarily a<c) where:
  c->a (so the edge goes c to a) AND exists b with a->b->c.

But since each unordered pair {a,c} has exactly one direction, and the NA condition
treats (a,c) as the ordered pair where c->a, we can count:

#NA = #{unordered pairs {a,c} where the LOSER (say c->a direction) has A^2[a,c]>0}

This equals the number of EDGES that participate in at least one 3-cycle,
where we count each edge once (the backward direction of some ordered pair).

Actually every 3-cycle {a,b,c} with a->b->c->a contributes 3 such pairs:
  (a,c) with c->a, A^2[a,c]>=1 (via b)
  (b,a) with a->b... wait no. Let me think again.

The NA face condition: face (a,c) is non-allowed in boundary of 2-path (a,b,c)
iff c->a. The 2-path (a,b,c) means a->b->c.

So #NA = #{distinct ordered (a,c) : c->a AND exists b with a->b->c}
       = #{edges e=(c->a) such that e is in some 3-cycle}

Each 3-cycle a->b->c->a contributes 3 backward edges: c->a, a->b...
no, let me count from the 3-cycle perspective.

3-cycle: a->b->c->a. The backward edges (relative to 2-paths) are:
  2-path (a,b,c): backward edge c->a ✓  => NA pair (a,c)
  2-path (b,c,a): backward edge a->b ✓  => NA pair (b,a)
  2-path (c,a,b): backward edge b->c ✓  => NA pair (c,b)

So each 3-cycle contributes 3 NA pairs. But pairs can appear from multiple 3-cycles.

An edge (c->a) is an NA pair iff it's in at LEAST one 3-cycle.
Let e_cyc = #{edges in at least one 3-cycle}.
Then #NA = e_cyc.

And: dim(Omega_2) = C(n,3) + 2*c3 - e_cyc.

Is e_cyc a simple function of c3 and n? Let's check.

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            c3 += 1
    return c3

def edges_in_3cycles(A, n):
    """Count edges that participate in at least one directed 3-cycle."""
    in_cycle = set()
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            in_cycle.add((i,j))
            in_cycle.add((j,k))
            in_cycle.add((k,i))
        elif A[i][k] and A[k][j] and A[j][i]:
            in_cycle.add((i,k))
            in_cycle.add((k,j))
            in_cycle.add((j,i))
    return len(in_cycle)

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def compute_omega_dim(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return 0
    if p == 0: return dim_Ap
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return dim_Ap
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    sv = np.linalg.svd(P, compute_uv=False)
    rank = int(sum(s > 1e-10 for s in sv))
    return dim_Ap - rank

def count_na_faces(A, n, allowed_p, allowed_pm1):
    """Count distinct non-allowed (p-1)-faces appearing in boundaries of p-paths."""
    allowed_pm1_set = set(allowed_pm1)
    na_faces = set()
    for path in allowed_p:
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                na_faces.add(face)
    return len(na_faces)


def main():
    print("=" * 70)
    print("EXACT FORMULA FOR dim(Omega_2)")
    print("=" * 70)

    # Part 1: Verify dim(Omega_2) = |A_2| - #NA = C(n,3) + 2*c3 - e_cyc
    print("\n--- Part 1: Exhaustive check dim(Omega_2) = C(n,3) + 2*c3 - e_cyc ---")

    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        violations = 0
        total = 0

        for bits in range(N):
            A = bits_to_adj(bits, n)
            c3 = count_3cycles(A, n)
            e_cyc = edges_in_3cycles(A, n)

            a1 = enumerate_allowed_paths(A, n, 1)
            a2 = enumerate_allowed_paths(A, n, 2)
            na = count_na_faces(A, n, a2, a1)
            dim_om2 = compute_omega_dim(A, n, 2, a2, a1)

            # Check: |A_2| = C(n,3) + 2*c3
            assert len(a2) == comb(n, 3) + 2*c3, f"A2 formula fail at bits={bits}"

            # Check: #NA = e_cyc
            if na != e_cyc:
                violations += 1
                if violations <= 3:
                    print(f"  MISMATCH at bits={bits}: #NA={na}, e_cyc={e_cyc}")

            # Check: dim(Omega_2) = |A_2| - #NA (from THM-119)
            if dim_om2 != len(a2) - na:
                print(f"  *** THM-119 FAIL at bits={bits}: dim={dim_om2}, |A2|-NA={len(a2)-na}")

            total += 1

        if violations == 0:
            print(f"  n={n}: CONFIRMED #NA = e_cyc for all {total} tournaments")
        else:
            print(f"  n={n}: {violations}/{total} mismatches")

    # Part 2: What is e_cyc as a function of c3 and n?
    print("\n--- Part 2: e_cyc vs c3 relationship ---")

    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        data = defaultdict(list)

        for bits in range(N):
            A = bits_to_adj(bits, n)
            c3 = count_3cycles(A, n)
            e_cyc = edges_in_3cycles(A, n)
            data[c3].append(e_cyc)

        print(f"\n  n={n}: e_cyc grouped by c3:")
        print(f"  {'c3':>4} | {'min':>4} | {'max':>4} | {'mean':>6} | {'values':>20}")
        for c3 in sorted(data.keys()):
            vals = data[c3]
            unique = sorted(set(vals))
            v_str = str(unique) if len(unique) <= 6 else f"{unique[:3]}...{unique[-2:]}"
            print(f"  {c3:4d} | {min(vals):4d} | {max(vals):4d} | {np.mean(vals):6.1f} | {v_str}")

    # Part 3: Is e_cyc = 3*c3 when ALL edges are in cycles?
    # For a regular tournament, every edge is in exactly (n-3)/2 3-cycles,
    # so all edges are in cycles => e_cyc = n(n-1)/2... no, e_cyc counts DIRECTED edges.
    # Actually e_cyc = #{directed edges in at least one 3-cycle}.
    # Total directed edges = n(n-1)/2 (one per unordered pair in a tournament... no!
    # n(n-1) directed edges total, n(n-1)/2 actual tournament arcs).
    #
    # Wait: we count ORDERED pairs (a,c) with A[c][a]=1 AND A^2[a,c]>0.
    # That's the number of tournament ARCS (directed edges) in at least one 3-cycle.
    # Total arcs = n(n-1)/2.

    print("\n--- Part 3: e_cyc / (total arcs) ratio ---")
    for n in [5, 6, 7]:
        rng = np.random.RandomState(42 + n)
        total_arcs = n * (n - 1) // 2
        N = min(2**(n*(n-1)//2), 500)

        ecyc_vals = []
        for trial in range(N):
            if N == 2**(n*(n-1)//2):
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)
            e_cyc = edges_in_3cycles(A, n)
            ecyc_vals.append(e_cyc)

        print(f"  n={n}: total arcs={total_arcs}, "
              f"e_cyc range=[{min(ecyc_vals)},{max(ecyc_vals)}], "
              f"mean={np.mean(ecyc_vals):.1f}, "
              f"transitive e_cyc={0}")

    # Part 4: Formula for dim(Omega_2) at n=5,6 exhaustive
    print("\n--- Part 4: dim(Omega_2) = C(n,3) + 2*c3 - e_cyc ---")
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        dim_vals = defaultdict(list)

        for bits in range(N):
            A = bits_to_adj(bits, n)
            c3 = count_3cycles(A, n)
            e_cyc = edges_in_3cycles(A, n)
            dim_formula = comb(n, 3) + 2*c3 - e_cyc
            dim_vals[(c3, e_cyc)].append(dim_formula)

        print(f"\n  n={n}: dim(Omega_2) = {comb(n,3)} + 2*c3 - e_cyc")
        for key in sorted(dim_vals.keys()):
            c3, e_cyc = key
            vals = dim_vals[key]
            dim_om2 = comb(n, 3) + 2*c3 - e_cyc
            print(f"    c3={c3}, e_cyc={e_cyc}: dim(Omega_2)={dim_om2}, count={len(vals)}")

    # Part 5: Defect rate U-shape investigation
    # The "defect rate" = probability that beta_p > 0 for random tournaments.
    # It shows a U-shape: high at low p and high p, low in the middle.
    # This is related to the filling ratio: when dim(Omega_p) is close to C(n,p+1),
    # there's little room for homology, so defect rate is low.
    print("\n--- Part 5: Defect rate U-shape at n=7,8 ---")
    for n in [7, 8]:
        rng = np.random.RandomState(42 + n)
        N = 200

        nonzero_counts = defaultdict(int)

        for trial in range(N):
            A = random_tournament(n, rng)
            paths_cache = {}
            for p in range(n):
                paths_cache[p] = enumerate_allowed_paths(A, n, p)

            for p in range(1, n-1):  # beta_0 always 1, beta_{n-1} special
                ob = compute_omega_dim(A, n, p, paths_cache[p], paths_cache[p-1])
                # Need ker(d_p) and im(d_{p+1}) to get beta_p
                # For now just track dim(Omega_p)
                pass

            # Actually compute full betti vector for smaller n
            if n <= 7:
                from numpy.linalg import svd
                om_dims = []
                for p in range(n):
                    om_dim = compute_omega_dim(A, n, p, paths_cache[p],
                                              paths_cache[p-1] if p > 0 else [])
                    om_dims.append(om_dim)

                # Build boundary matrices in Omega basis
                # This is expensive, so let's use a simpler approach
                # Just check if dim(Omega_p) > 0 and track

        # Actually let me just compute betti numbers properly for n=7
        if n == 7:
            print(f"\n  n={n}: Computing full Betti vectors for {N} tournaments...")
            betti_nonzero = defaultdict(int)
            dim_omega_sum = defaultdict(float)

            for trial in range(N):
                A = random_tournament(n, rng)

                # Cache all allowed paths
                ap = {}
                for p in range(n):
                    ap[p] = enumerate_allowed_paths(A, n, p)

                # Compute omega bases
                omega_bases = {}
                for p in range(n):
                    dim_Ap = len(ap[p])
                    if dim_Ap == 0:
                        omega_bases[p] = np.zeros((0, 0))
                        continue
                    if p == 0:
                        omega_bases[p] = np.eye(dim_Ap)
                        continue
                    apm1_set = set(ap[p-1])
                    non_allowed = {}
                    na_count = 0
                    for j, path in enumerate(ap[p]):
                        for sign, face in boundary_coeffs(path):
                            if len(set(face)) == len(face) and face not in apm1_set:
                                if face not in non_allowed:
                                    non_allowed[face] = na_count
                                    na_count += 1
                    if na_count == 0:
                        omega_bases[p] = np.eye(dim_Ap)
                    else:
                        P = np.zeros((na_count, dim_Ap))
                        for j, path in enumerate(ap[p]):
                            for sign, face in boundary_coeffs(path):
                                if face in non_allowed:
                                    P[non_allowed[face], j] += sign
                        U, S, Vt = np.linalg.svd(P, full_matrices=True)
                        rank = int(sum(s > 1e-10 for s in S))
                        ns = Vt[rank:].T
                        omega_bases[p] = ns if ns.shape[1] > 0 else np.zeros((dim_Ap, 0))

                    dim_omega_sum[p] += omega_bases[p].shape[1]

                # Compute boundary maps in Omega basis
                betti = []
                for p in range(n):
                    dim_p = omega_bases[p].shape[1] if omega_bases[p].ndim == 2 else 0

                    # ker(d_p)
                    if p == 0 or dim_p == 0:
                        ker_dp = dim_p
                    else:
                        bd = np.zeros((len(ap[p-1]), len(ap[p])))
                        idx_pm1 = {path: i for i, path in enumerate(ap[p-1])}
                        for j, path in enumerate(ap[p]):
                            for sign, face in boundary_coeffs(path):
                                if face in idx_pm1:
                                    bd[idx_pm1[face], j] += sign
                        dp_om = bd @ omega_bases[p]
                        sv = np.linalg.svd(dp_om, compute_uv=False)
                        rank_dp = int(sum(s > 1e-8 for s in sv))
                        ker_dp = dim_p - rank_dp

                    # im(d_{p+1})
                    if p >= n - 1:
                        im_dp1 = 0
                    else:
                        dim_p1 = omega_bases[p+1].shape[1] if omega_bases[p+1].ndim == 2 else 0
                        if dim_p1 == 0:
                            im_dp1 = 0
                        else:
                            bd1 = np.zeros((len(ap[p]), len(ap[p+1])))
                            idx_p = {path: i for i, path in enumerate(ap[p])}
                            for j, path in enumerate(ap[p+1]):
                                for sign, face in boundary_coeffs(path):
                                    if face in idx_p:
                                        bd1[idx_p[face], j] += sign
                            dp1_om = bd1 @ omega_bases[p+1]
                            sv1 = np.linalg.svd(dp1_om, compute_uv=False)
                            im_dp1 = int(sum(s > 1e-8 for s in sv1))

                    bp = ker_dp - im_dp1
                    betti.append(bp)
                    if bp > 0:
                        betti_nonzero[p] += 1

                if trial % 50 == 0:
                    print(f"    trial {trial}: betti={betti}", flush=True)

            print(f"\n  n={n}: Defect rates (P(beta_p > 0)):")
            for p in range(n):
                rate = betti_nonzero[p] / N
                avg_dim = dim_omega_sum[p] / N if p > 0 else n
                print(f"    p={p}: rate={rate:.3f} ({betti_nonzero[p]}/{N}), avg dim(Omega_p)={avg_dim:.1f}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
