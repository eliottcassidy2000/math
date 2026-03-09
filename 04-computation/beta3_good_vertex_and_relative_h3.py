"""
beta3_good_vertex_and_relative_h3.py — Two key ingredients for beta_3 <= 1 proof

Ingredient 1: Good vertex existence for beta_3
  For every tournament T with beta_3(T) > 0, does there exist v with beta_3(T\\v) = 0?

Ingredient 2: Relative homology bound
  When beta_3(T\\v) = 0, H_3(T) ≅ H_3(T, T\\v). Is dim H_3(T, T\\v) <= 1?

  H_3(T, T\\v) = ker(d_3^rel) / im(d_4^rel) where d^rel are the relative boundary maps
  on the quotient complex C_*(T) / C_*(T\\v).

  C_p^rel = span of allowed p-paths in T that pass through v.

Author: kind-pasteur-S46 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter, defaultdict
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

def compute_omega_basis(ap, p):
    if not ap.get(p, []):
        return np.zeros((0, 0)), 0
    if p == 0:
        return np.eye(len(ap[p])), len(ap[p])

    apm1_set = set(ap.get(p-1, []))
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in apm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(len(ap[p])), len(ap[p])

    P = np.zeros((na_count, len(ap[p])))
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = int(sum(s > 1e-10 for s in S))
    ns = Vt[rank:].T
    if ns.shape[1] > 0:
        return ns, ns.shape[1]
    else:
        return np.zeros((len(ap[p]), 0)), 0


def compute_beta3(A, n):
    """Compute beta_3 for tournament A."""
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p)

    omega_bases = {}
    omega_dims = {}
    for p in range(min(6, n)):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

    dim3 = omega_dims.get(3, 0)
    if dim3 == 0:
        return 0

    # Build d_3
    bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_coeffs(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3
    sv3 = np.linalg.svd(d3_om, compute_uv=False)
    rank_d3 = int(sum(s > 1e-8 for s in sv3))
    ker_d3 = dim3 - rank_d3

    if ker_d3 == 0:
        return 0

    # Build d_4
    dim4 = omega_dims.get(4, 0)
    if dim4 == 0:
        return ker_d3

    bd4 = np.zeros((len(ap[3]), len(ap.get(4, []))))
    idx3 = {path: i for i, path in enumerate(ap[3])}
    for j, path in enumerate(ap.get(4, [])):
        for sign, face in boundary_coeffs(path):
            if face in idx3:
                bd4[idx3[face], j] += sign

    O4 = omega_bases[4]
    d4_om = bd4 @ O4
    O3pinv = np.linalg.pinv(O3)
    d4_omega3 = O3pinv @ d4_om
    sv4 = np.linalg.svd(d4_omega3, compute_uv=False)
    rank_d4 = int(sum(s > 1e-8 for s in sv4))

    return ker_d3 - rank_d4


def compute_full_homology(A, n, max_p=5):
    """Compute all Betti numbers and return chain complex data."""
    ap = {}
    for p in range(min(max_p+2, n)):
        ap[p] = enumerate_allowed_paths(A, n, p)

    omega_bases = {}
    omega_dims = {}
    for p in range(min(max_p+2, n)):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

    return ap, omega_bases, omega_dims


def deletion_adj(A, n, v):
    vertices = [i for i in range(n) if i != v]
    n1 = len(vertices)
    A1 = np.zeros((n1, n1), dtype=int)
    for i, vi in enumerate(vertices):
        for j, vj in enumerate(vertices):
            A1[i][j] = A[vi][vj]
    return A1, n1


def compute_relative_h3(A, n, v):
    """
    Compute dim H_3(T, T\\v) = dim of relative homology.

    Relative chain complex: C_p^rel = Omega_p(T) / Omega_p(T\\v)
    Concretely: Omega_p^rel = span of Omega_p basis elements that use vertex v.

    We compute Omega_p(T) and Omega_p(T\\v), find the quotient,
    then compute kernel/image of induced boundary maps.
    """
    # Full complex
    ap_full = {}
    for p in range(min(6, n)):
        ap_full[p] = enumerate_allowed_paths(A, n, p)

    # Deletion complex (T\\v)
    A1, n1 = deletion_adj(A, n, v)
    ap_del = {}
    for p in range(min(6, n1)):
        ap_del[p] = enumerate_allowed_paths(A1, n1, p)

    # Map deletion paths to full paths (relabel)
    verts_del = [i for i in range(n) if i != v]
    def relabel_path(path):
        return tuple(verts_del[i] for i in path)

    # Omega bases for full and deletion
    ob_full = {}
    od_full = {}
    ob_del = {}
    od_del = {}
    for p in range(min(6, n)):
        if not ap_full.get(p, []):
            ob_full[p] = np.zeros((0, 0))
            od_full[p] = 0
        else:
            ob_full[p], od_full[p] = compute_omega_basis(ap_full, p)
    for p in range(min(6, n1)):
        if not ap_del.get(p, []):
            ob_del[p] = np.zeros((0, 0))
            od_del[p] = 0
        else:
            ob_del[p], od_del[p] = compute_omega_basis(ap_del, p)

    # For each p, identify "through-v" paths = those in Omega_p(T) but NOT in Omega_p(T\\v)
    # We need relative Omega, which is quotient Omega(T) / Omega(T\\v).

    # Method: Express Omega(T\\v) inside Omega(T) via relabeling, then take quotient.
    # Paths through v: allowed p-paths in T that include vertex v in their vertex set.

    # Simple approach: among Omega_p(T) basis vectors, find those that are NOT in span of Omega_p(T\\v).
    # The relative chain group is the quotient.

    # Actually, let me use a direct computation:
    # C_p^rel has basis = allowed p-paths through v in Omega_p(T), modulo Omega relations.

    # Simpler: just compute boundary maps on the quotient directly.
    # A path "involves v" if v appears in it.

    # For the boundary d_3^rel: C_3^rel -> C_2^rel
    # and d_4^rel: C_4^rel -> C_3^rel
    # Then H_3^rel = ker(d_3^rel) / im(d_4^rel)

    # Identify "through-v" allowed paths for each degree
    through_v = {}
    for p in range(min(6, n)):
        through_v[p] = [path for path in ap_full[p] if v in path]

    # Build Omega_p^rel: the through-v part of Omega_p
    # First, get Omega_p(T) basis. Among these basis vectors, some project nontrivially onto through-v paths.
    # The relative Omega_p = image of Omega_p(T) in the quotient by non-v paths.

    # Let's index through-v paths
    tv_idx = {}
    for p in range(min(6, n)):
        tv_idx[p] = {path: i for i, path in enumerate(through_v[p])}

    # Project Omega_p(T) basis onto through-v coordinates
    rel_bases = {}
    rel_dims = {}
    for p in range(min(6, n)):
        if od_full.get(p, 0) == 0 or not through_v[p]:
            rel_bases[p] = np.zeros((0, 0))
            rel_dims[p] = 0
            continue

        # Full Omega basis: columns of ob_full[p], indexed by ap_full[p]
        num_full = len(ap_full[p])
        num_tv = len(through_v[p])

        # Projection matrix: extract through-v rows
        proj = np.zeros((num_tv, num_full))
        full_idx = {path: i for i, path in enumerate(ap_full[p])}
        for path, ti in tv_idx[p].items():
            proj[ti, full_idx[path]] = 1.0

        # Project Omega basis
        proj_omega = proj @ ob_full[p]  # (num_tv x omega_dim)

        # Find rank = dim of relative Omega
        if proj_omega.size == 0:
            rel_bases[p] = np.zeros((0, 0))
            rel_dims[p] = 0
            continue
        sv = np.linalg.svd(proj_omega, compute_uv=False)
        rank = int(sum(s > 1e-10 for s in sv))
        rel_dims[p] = rank

        if rank == 0:
            rel_bases[p] = np.zeros((num_tv, 0))
            continue

        # Get basis for the image
        U, S, Vt = np.linalg.svd(proj_omega, full_matrices=False)
        rel_bases[p] = U[:, :rank]  # (num_tv x rank) orthonormal basis

    # Now build relative boundary maps
    # d_3^rel: C_3^rel -> C_2^rel
    dim3_rel = rel_dims.get(3, 0)
    dim2_rel = rel_dims.get(2, 0)

    if dim3_rel == 0:
        return 0, rel_dims

    # Boundary matrix on through-v paths: bd3_tv[i2, j3] where i2 indexes through_v[2], j3 indexes through_v[3]
    # But we need to include ALL faces, not just through-v ones.
    # Wait — the relative boundary goes from rel to rel. A through-v 3-path has faces that may or may not go through v.
    # In the relative complex, the non-v faces are zero (they're in C_*(T\\v)).
    # So d^rel only keeps the through-v faces!

    # Build raw relative boundary: through_v[3] -> through_v[2]
    bd3_tv = np.zeros((len(through_v.get(2, [])), len(through_v[3])))
    for j, path in enumerate(through_v[3]):
        for sign, face in boundary_coeffs(path):
            if face in tv_idx.get(2, {}):
                bd3_tv[tv_idx[2][face], j] += sign

    # Restrict to relative Omega bases
    if dim2_rel > 0 and dim3_rel > 0:
        d3_rel = rel_bases[2].T @ bd3_tv @ rel_bases[3]
    elif dim3_rel > 0:
        # d3_rel maps to 0
        d3_rel = np.zeros((0, dim3_rel))

    sv3 = np.linalg.svd(d3_rel, compute_uv=False) if d3_rel.size > 0 else np.array([])
    rank_d3_rel = int(sum(s > 1e-8 for s in sv3))
    ker_d3_rel = dim3_rel - rank_d3_rel

    if ker_d3_rel == 0:
        return 0, rel_dims

    # d_4^rel: C_4^rel -> C_3^rel
    dim4_rel = rel_dims.get(4, 0)
    if dim4_rel == 0:
        return ker_d3_rel, rel_dims

    bd4_tv = np.zeros((len(through_v.get(3, [])), len(through_v.get(4, []))))
    for j, path in enumerate(through_v.get(4, [])):
        for sign, face in boundary_coeffs(path):
            if face in tv_idx.get(3, {}):
                bd4_tv[tv_idx[3][face], j] += sign

    if dim3_rel > 0 and dim4_rel > 0:
        d4_rel = rel_bases[3].T @ bd4_tv @ rel_bases[4]
    else:
        d4_rel = np.zeros((dim3_rel, 0))

    sv4 = np.linalg.svd(d4_rel, compute_uv=False) if d4_rel.size > 0 else np.array([])
    rank_d4_rel = int(sum(s > 1e-8 for s in sv4))

    h3_rel = ker_d3_rel - rank_d4_rel
    return h3_rel, rel_dims


def main():
    print("=" * 70)
    print("BETA_3 GOOD VERTEX & RELATIVE H_3 — Key ingredients for beta_3 <= 1")
    print("=" * 70)

    # Part 1: Good vertex existence at n=6 (exhaustive)
    print("\n--- Part 1: Good vertex existence at n=6 (exhaustive) ---")
    n = 6
    total_b3 = 0
    has_good = 0
    no_good = 0

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue
        total_b3 += 1

        found_good = False
        for v_del in range(n):
            A1, n1 = deletion_adj(A, n, v_del)
            b3v = compute_beta3(A1, n1)
            if b3v == 0:
                found_good = True
                break

        if found_good:
            has_good += 1
        else:
            no_good += 1

    print(f"  beta_3>0 tournaments: {total_b3}")
    print(f"  Has good vertex (beta_3(T\\v)=0): {has_good}")
    print(f"  NO good vertex: {no_good}")

    # Part 2: Good vertex at n=7 (sampled)
    print("\n--- Part 2: Good vertex existence at n=7 (sampled) ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 500
    b3_count = 0
    good_count = 0
    no_good_count = 0
    num_good_dist = Counter()

    for trial in range(N):
        A = random_tournament(n, rng)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue
        b3_count += 1

        num_good = 0
        for v_del in range(n):
            A1, n1 = deletion_adj(A, n, v_del)
            b3v = compute_beta3(A1, n1)
            if b3v == 0:
                num_good += 1

        num_good_dist[num_good] += 1
        if num_good > 0:
            good_count += 1
        else:
            no_good_count += 1

        if (trial + 1) % 100 == 0:
            print(f"  n=7: {trial+1}/{N} done, {b3_count} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count}/{N}")
    print(f"  Has good vertex: {good_count}")
    print(f"  NO good vertex: {no_good_count}")
    print(f"  #good vertices distribution: {dict(sorted(num_good_dist.items()))}")

    # Part 3: Good vertex at n=8 (sampled)
    print("\n--- Part 3: Good vertex existence at n=8 (sampled) ---")
    n = 8
    rng8 = np.random.RandomState(99)
    N8 = 200
    b3_count8 = 0
    good_count8 = 0
    no_good8 = 0
    num_good_dist8 = Counter()

    for trial in range(N8):
        A = random_tournament(n, rng8)
        try:
            b3 = compute_beta3(A, n)
        except:
            continue
        if b3 == 0:
            continue
        b3_count8 += 1

        num_good = 0
        for v_del in range(n):
            A1, n1 = deletion_adj(A, n, v_del)
            b3v = compute_beta3(A1, n1)
            if b3v == 0:
                num_good += 1

        num_good_dist8[num_good] += 1
        if num_good > 0:
            good_count8 += 1
        else:
            no_good8 += 1

        if (trial + 1) % 50 == 0:
            print(f"  n=8: {trial+1}/{N8} done, {b3_count8} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count8}/{N8}")
    print(f"  Has good vertex: {good_count8}")
    print(f"  NO good vertex: {no_good8}")
    print(f"  #good vertices distribution: {dict(sorted(num_good_dist8.items()))}")

    # Part 4: Relative H_3(T, T\\v) computation at n=6
    print("\n--- Part 4: dim H_3(T, T\\v) at n=6 (exhaustive for beta_3>0) ---")
    n = 6
    rel_h3_max = 0
    rel_h3_dist = Counter()
    rel_h3_when_good = Counter()

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue

        for v_del in range(n):
            A1, n1 = deletion_adj(A, n, v_del)
            b3v = compute_beta3(A1, n1)

            h3_rel, rel_dims = compute_relative_h3(A, n, v_del)
            rel_h3_dist[h3_rel] += 1
            rel_h3_max = max(rel_h3_max, h3_rel)

            if b3v == 0:
                rel_h3_when_good[h3_rel] += 1

    print(f"  H_3(T,T\\v) distribution (all v, all T with beta_3>0):")
    for h, cnt in sorted(rel_h3_dist.items()):
        print(f"    dim={h}: {cnt}")
    print(f"  H_3(T,T\\v) when beta_3(T\\v)=0:")
    for h, cnt in sorted(rel_h3_when_good.items()):
        print(f"    dim={h}: {cnt}")
    print(f"  Max dim H_3(T,T\\v): {rel_h3_max}")

    # Part 5: Relative H_3 at n=7 (sampled)
    print("\n--- Part 5: dim H_3(T, T\\v) at n=7 (sampled) ---")
    n = 7
    rng5 = np.random.RandomState(321)
    N5 = 200
    rel_h3_dist7 = Counter()
    rel_h3_good7 = Counter()
    b3_count5 = 0

    for trial in range(N5):
        A = random_tournament(n, rng5)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue
        b3_count5 += 1

        for v_del in range(n):
            A1, n1 = deletion_adj(A, n, v_del)
            b3v = compute_beta3(A1, n1)

            h3_rel, rel_dims = compute_relative_h3(A, n, v_del)
            rel_h3_dist7[h3_rel] += 1

            if b3v == 0:
                rel_h3_good7[h3_rel] += 1

        if (trial + 1) % 50 == 0:
            print(f"  n=7: {trial+1}/{N5} done, {b3_count5} with beta_3>0", flush=True)

    print(f"\n  beta_3>0: {b3_count5}/{N5}")
    print(f"  H_3(T,T\\v) distribution:")
    for h, cnt in sorted(rel_h3_dist7.items()):
        print(f"    dim={h}: {cnt}")
    print(f"  H_3(T,T\\v) when beta_3(T\\v)=0:")
    for h, cnt in sorted(rel_h3_good7.items()):
        print(f"    dim={h}: {cnt}")

    # Part 6: LES verification — does H_3(T) = H_3(T,T\\v) when beta_3(T\\v)=0?
    print("\n--- Part 6: LES isomorphism check at n=6 ---")
    n = 6
    iso_ok = 0
    iso_fail = 0

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        b3 = compute_beta3(A, n)
        if b3 == 0:
            continue

        for v_del in range(n):
            A1, n1 = deletion_adj(A, n, v_del)
            b3v = compute_beta3(A1, n1)

            if b3v == 0:
                h3_rel, _ = compute_relative_h3(A, n, v_del)
                if h3_rel == b3:
                    iso_ok += 1
                else:
                    iso_fail += 1
                    if iso_fail <= 5:
                        print(f"  FAIL: bits={bits}, v={v_del}, beta_3(T)={b3}, H_3^rel={h3_rel}")

    print(f"  LES isomorphism (beta_3(T\\v)=0 => beta_3=H_3^rel): OK={iso_ok}, FAIL={iso_fail}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
