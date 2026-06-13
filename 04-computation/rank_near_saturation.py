"""
rank_near_saturation.py -- Why is rank(d_4) always ker(d_3) or ker(d_3)-1?

If we can understand this, we prove beta_3 in {0,1} (THM-098).

Key question: what is the cokernel of the map d_4: Omega_4 -> ker(d_3)?
When dim(coker) = 1, there's exactly one independent element in ker(d_3)
that cannot be hit by d_4. What makes this element special?

Approach:
1. Collect detailed rank data for d_3, d_4 at n=6 (exhaustive), n=7 (sampled)
2. When beta_3=1: find the unique cokernel element and analyze its structure
3. Look for patterns: does cokernel relate to score sequence? cycle structure?

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
    """Compute basis of Omega_p as null space of non-allowed face constraints."""
    if not ap[p]:
        return np.zeros((0, 0)), 0
    if p == 0:
        return np.eye(len(ap[p])), len(ap[p])

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

def compute_chain_data(A, n, max_p=5):
    """Compute detailed chain complex data up to level max_p."""
    ap = {}
    for p in range(max_p + 2):
        if p <= n - 1:
            ap[p] = enumerate_allowed_paths(A, n, p)
        else:
            ap[p] = []

    omega_bases = {}
    omega_dims = {}
    for p in range(max_p + 2):
        if p > n - 1 or not ap[p]:
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

    # Build boundary maps d_p: Omega_p -> (A_{p-1} space, then restrict to Omega_{p-1})
    results = {}
    for p in range(1, max_p + 1):
        dim_p = omega_dims.get(p, 0)
        dim_pm1 = omega_dims.get(p-1, 0)
        if dim_p == 0 or dim_pm1 == 0:
            results[p] = {'rank': 0, 'ker': dim_p, 'dim': dim_p}
            continue

        bd = np.zeros((len(ap[p-1]), len(ap[p])))
        idx_pm1 = {path: i for i, path in enumerate(ap[p-1])}
        for j, path in enumerate(ap[p]):
            for sign, face in boundary_coeffs(path):
                if face in idx_pm1:
                    bd[idx_pm1[face], j] += sign

        dp_om = bd @ omega_bases[p]
        sv = np.linalg.svd(dp_om, compute_uv=False)
        rank = int(sum(s > 1e-8 for s in sv))
        results[p] = {'rank': rank, 'ker': dim_p - rank, 'dim': dim_p}

    # Compute Betti numbers
    betti = []
    for p in range(max_p + 1):
        ker_p = results.get(p, {}).get('ker', omega_dims.get(p, 0) if p == 0 else 0)
        if p == 0:
            ker_p = omega_dims[0]
            im_1 = results.get(1, {}).get('rank', 0)
            betti.append(ker_p - im_1)
        else:
            im_next = results.get(p+1, {}).get('rank', 0)
            betti.append(ker_p - im_next)

    return {
        'ap': ap,
        'omega_bases': omega_bases,
        'omega_dims': omega_dims,
        'boundary': results,
        'betti': betti
    }


def main():
    print("=" * 70)
    print("RANK NEAR-SATURATION ANALYSIS")
    print("Why is rank(d_4) always = ker(d_3) or ker(d_3)-1?")
    print("=" * 70)

    # Part 1: n=6 exhaustive — collect (ker_d3, rank_d4, beta_3) triples
    print("\n--- Part 1: n=6 exhaustive ---")
    n = 6
    N = 2**(n*(n-1)//2)
    gap_counter = Counter()  # gap = ker(d_3) - rank(d_4)
    detail_by_gap = defaultdict(list)  # gap -> list of (bits, scores, dims)
    total_b3 = 0

    for bits in range(N):
        A = bits_to_adj(bits, n)
        data = compute_chain_data(A, n, max_p=4)
        betti = data['betti']
        b3 = betti[3] if len(betti) > 3 else 0

        # Extract key data
        ker_d3 = data['boundary'].get(3, {}).get('ker', 0)
        rank_d4 = data['boundary'].get(4, {}).get('rank', 0)
        gap = ker_d3 - rank_d4

        gap_counter[gap] += 1
        if gap > 0:
            total_b3 += 1
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            dims = {p: data['omega_dims'].get(p, 0) for p in range(6)}
            ranks = {p: data['boundary'].get(p, {}).get('rank', 0) for p in range(1, 5)}
            detail_by_gap[gap].append({
                'bits': bits,
                'scores': scores,
                'dims': dims,
                'ranks': ranks,
                'ker_d3': ker_d3,
                'rank_d4': rank_d4,
                'betti': betti[:5]
            })

        if (bits + 1) % 5000 == 0:
            print(f"  n=6: {bits+1}/{N} done, {total_b3} with beta_3>0", flush=True)

    print(f"\n  RESULTS n=6:")
    print(f"  Total tournaments: {N}")
    print(f"  Gap distribution (ker(d_3) - rank(d_4)):")
    for gap in sorted(gap_counter.keys()):
        cnt = gap_counter[gap]
        pct = 100*cnt/N
        print(f"    gap={gap}: {cnt} ({pct:.1f}%)")
    print(f"  beta_3>0 (gap>0): {total_b3}")

    # Analyze gap=1 cases
    if 1 in detail_by_gap:
        g1 = detail_by_gap[1]
        print(f"\n  Gap=1 cases ({len(g1)} tournaments):")
        score_dist = Counter(d['scores'] for d in g1)
        for sc, cnt in score_dist.most_common():
            print(f"    scores {sc}: {cnt}")

        # Show dimension profiles
        dim_profiles = Counter(
            tuple(d['dims'][p] for p in range(6))
            for d in g1
        )
        print(f"\n  Omega dimension profiles for gap=1:")
        for prof, cnt in dim_profiles.most_common():
            print(f"    dims[0..5]={prof}: {cnt}")

        # Show ker/rank profiles
        print(f"\n  Boundary rank profiles for gap=1:")
        rank_profiles = Counter(
            (d['ranks'][1], d['ranks'][2], d['ranks'][3], d['ranks'].get(4, 0))
            for d in g1
        )
        for prof, cnt in rank_profiles.most_common():
            print(f"    (rk_d1, rk_d2, rk_d3, rk_d4)={prof}: {cnt}")

    # Part 2: n=7 sampled — same analysis at larger n
    print("\n--- Part 2: n=7 sampled ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 1000
    gap_counter7 = Counter()
    detail_by_gap7 = defaultdict(list)

    for trial in range(N):
        A = random_tournament(n, rng)
        data = compute_chain_data(A, n, max_p=5)
        betti = data['betti']
        b3 = betti[3] if len(betti) > 3 else 0

        ker_d3 = data['boundary'].get(3, {}).get('ker', 0)
        rank_d4 = data['boundary'].get(4, {}).get('rank', 0)
        gap = ker_d3 - rank_d4

        gap_counter7[gap] += 1
        if gap > 0:
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            dims = {p: data['omega_dims'].get(p, 0) for p in range(7)}
            ranks = {p: data['boundary'].get(p, {}).get('rank', 0) for p in range(1, 6)}
            detail_by_gap7[gap].append({
                'trial': trial,
                'scores': scores,
                'dims': dims,
                'ranks': ranks,
                'ker_d3': ker_d3,
                'rank_d4': rank_d4,
                'betti': betti[:6]
            })

        if (trial + 1) % 200 == 0:
            print(f"  n=7: {trial+1}/{N} done", flush=True)

    print(f"\n  RESULTS n=7:")
    print(f"  Gap distribution:")
    for gap in sorted(gap_counter7.keys()):
        cnt = gap_counter7[gap]
        pct = 100*cnt/N
        print(f"    gap={gap}: {cnt} ({pct:.1f}%)")

    if 1 in detail_by_gap7:
        g1_7 = detail_by_gap7[1]
        print(f"\n  Gap=1 cases ({len(g1_7)} tournaments):")
        score_dist7 = Counter(d['scores'] for d in g1_7)
        for sc, cnt in score_dist7.most_common():
            print(f"    scores {sc}: {cnt}")

        dim_profiles7 = Counter(
            tuple(d['dims'][p] for p in range(7))
            for d in g1_7
        )
        print(f"\n  Omega dimension profiles for gap=1:")
        for prof, cnt in dim_profiles7.most_common(10):
            print(f"    dims[0..6]={prof}: {cnt}")

    # Part 3: Relationship between ker(d_3) and dim(Omega_4)
    # If dim(Omega_4) = 0 then rank(d_4) = 0, so gap = ker(d_3)
    # When does dim(Omega_4) = 0?
    print("\n--- Part 3: When is Omega_4 = 0 at n=6? ---")
    n = 6
    om4_zero = 0
    om4_zero_b3 = 0
    N6 = 2**(n*(n-1)//2)

    for bits in range(N6):
        A = bits_to_adj(bits, n)
        ap4 = enumerate_allowed_paths(A, n, 4)
        if not ap4:
            om4_zero += 1
            # Quick beta_3 check
            data = compute_chain_data(A, n, max_p=4)
            b3 = data['betti'][3] if len(data['betti']) > 3 else 0
            if b3 > 0:
                om4_zero_b3 += 1

    print(f"  Omega_4 = 0: {om4_zero}/{N6} ({100*om4_zero/N6:.1f}%)")
    print(f"  Of those with Omega_4=0, beta_3>0: {om4_zero_b3}")

    # Part 4: For beta_3=1 cases at n=6, what is dim(Omega_4)?
    print("\n--- Part 4: Omega_4 dimension when beta_3=1 at n=6 ---")
    if 1 in detail_by_gap:
        om4_dims = Counter(d['dims'][4] for d in detail_by_gap[1])
        print(f"  dim(Omega_4) distribution for beta_3=1:")
        for dim_val in sorted(om4_dims.keys()):
            print(f"    dim(Omega_4)={dim_val}: {om4_dims[dim_val]}")

    # Part 5: n=8 sampled — does gap ever exceed 1?
    print("\n--- Part 5: n=8 sampled ---")
    n = 8
    rng8 = np.random.RandomState(99)
    N8 = 300
    gap_counter8 = Counter()

    for trial in range(N8):
        A = random_tournament(n, rng8)
        try:
            data = compute_chain_data(A, n, max_p=5)
            betti = data['betti']
            ker_d3 = data['boundary'].get(3, {}).get('ker', 0)
            rank_d4 = data['boundary'].get(4, {}).get('rank', 0)
            gap = ker_d3 - rank_d4
            gap_counter8[gap] += 1
        except Exception as e:
            pass

        if (trial + 1) % 100 == 0:
            print(f"  n=8: {trial+1}/{N8} done", flush=True)

    print(f"\n  RESULTS n=8:")
    print(f"  Gap distribution:")
    for gap in sorted(gap_counter8.keys()):
        cnt = gap_counter8[gap]
        pct = 100*cnt/N8
        print(f"    gap={gap}: {cnt} ({pct:.1f}%)")

    # Part 6: Focus on the cokernel — what element of ker(d_3) is NOT in im(d_4)?
    print("\n--- Part 6: Cokernel element structure at n=6 ---")
    n = 6
    # Pick a few beta_3=1 tournaments and examine the cokernel
    count = 0
    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        data = compute_chain_data(A, n, max_p=4)
        betti = data['betti']
        if len(betti) > 3 and betti[3] == 1:
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))

            # Find the cokernel element
            ap = data['ap']
            O3 = data['omega_bases'][3]
            dim3 = data['omega_dims'][3]
            dim4 = data['omega_dims'][4]

            if dim3 == 0 or dim4 == 0:
                continue

            # Build d_3 in Omega coordinates
            bd3 = np.zeros((len(ap[2]), len(ap[3])))
            idx2 = {path: i for i, path in enumerate(ap[2])}
            for j, path in enumerate(ap[3]):
                for sign, face in boundary_coeffs(path):
                    if face in idx2:
                        bd3[idx2[face], j] += sign
            d3_om = bd3 @ O3
            U3, S3, Vt3 = np.linalg.svd(d3_om, full_matrices=True)
            rank_d3 = int(sum(s > 1e-8 for s in S3))
            ker_basis = Vt3[rank_d3:].T

            if ker_basis.shape[1] == 0:
                continue

            # Build d_4 in Omega coordinates
            O4 = data['omega_bases'][4]
            bd4 = np.zeros((len(ap[3]), len(ap[4])))
            idx3 = {path: i for i, path in enumerate(ap[3])}
            for j, path in enumerate(ap[4]):
                for sign, face in boundary_coeffs(path):
                    if face in idx3:
                        bd4[idx3[face], j] += sign
            d4_om = bd4 @ O4

            # Map d4 image to Omega_3 coordinates
            O3pinv = np.linalg.pinv(O3)
            d4_omega3 = O3pinv @ d4_om

            # Find cokernel: project ker(d_3) onto complement of im(d_4)
            U4, S4, Vt4 = np.linalg.svd(d4_omega3, full_matrices=True)
            rank_d4 = int(sum(s > 1e-8 for s in S4))
            im_basis = U4[:, :rank_d4]

            # The cokernel element: ker(d_3) projected orthogonally to im(d_4)
            for i in range(ker_basis.shape[1]):
                v = ker_basis[:, i]
                if rank_d4 > 0:
                    proj = im_basis @ (im_basis.T @ v)
                    residual = v - proj
                    if np.linalg.norm(residual) > 1e-6:
                        coker = residual / np.linalg.norm(residual)
                        # Express in A_3 basis
                        gen_A3 = O3 @ coker

                        # Analyze: which 3-paths have nonzero coefficient?
                        nonzero = [(ap[3][j], gen_A3[j])
                                   for j in range(len(ap[3]))
                                   if abs(gen_A3[j]) > 1e-6]

                        # Vertex sets used
                        vsets = set(frozenset(p) for p, c in nonzero)
                        # Signs pattern
                        signs = [np.sign(c) for p, c in nonzero]

                        print(f"\n  bits={bits}, scores={scores}")
                        print(f"  dims: Om3={dim3}, Om4={dim4}")
                        print(f"  ker(d3)={ker_basis.shape[1]}, rank(d4)={rank_d4}")
                        print(f"  Cokernel element: {len(nonzero)} nonzero 3-paths")
                        print(f"  Vertex sets: {len(vsets)} (sizes: {Counter(len(vs) for vs in vsets)})")
                        print(f"  Signs: {Counter(int(s) for s in signs)}")

                        # Check: does the cokernel element use ALL vertices?
                        all_verts = set()
                        for vs in vsets:
                            all_verts.update(vs)
                        print(f"  Vertices used: {sorted(all_verts)} ({len(all_verts)}/{n})")

                        count += 1
                        break

        if count >= 10:
            break

    print("\nDONE.")


if __name__ == '__main__':
    main()
