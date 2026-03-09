"""
Prove beta_1 * beta_3 = 0 for ALL tournaments (HYP-299).

Confirmed computationally: n=6 exhaustive (0/32768), n=7 (0/500), n=8 (0/300).

Strategy: Show that beta_1 > 0 and beta_3 > 0 cannot coexist.

Key facts:
- beta_1 <= 1 for all tournaments (THM-103)
- beta_3 <= 1 for all tournaments through n=8 (HYP-300)
- beta_2 = 0 for all tournaments (HYP-249, essentially proved)
- chi = 1 - beta_1 + beta_2 - beta_3 + ... = 1 - beta_1 - beta_3 + higher

If chi >= 0 always and beta_1, beta_3 <= 1, then beta_1 + beta_3 <= 1 + higher.
But we need chi >= 0 AND higher terms to be small enough.

Actually from the data:
- chi in {0, 1} for n <= 7
- This means 1 - beta_1 - beta_3 + (beta_4 - beta_5 + ...) in {0, 1}

If beta_4 = beta_5 = ... = 0 (generic case): chi = 1 - beta_1 - beta_3
  chi = 0: beta_1 + beta_3 = 1, so exactly one is 1
  chi = 1: beta_1 = beta_3 = 0

So for generic tournaments: beta_1 * beta_3 = 0 iff chi >= 0.

For non-generic (beta_4 > 0): need to check that beta_1 = beta_3 = 0 always.
From data: beta_4 > 0 implies beta_1 = beta_3 = 0 (Paley T_7, etc.)

Approach: Show that the chain complex structure forces mutual exclusivity.
Specifically: if there's a 1-cycle (beta_1=1), the boundary map d_4
saturates im(d_4) enough to kill beta_3.

Let's investigate the chain complex structure more carefully.
"""
import numpy as np
from math import comb
from collections import defaultdict
from itertools import combinations

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

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    ns = Vt[rank:].T
    return ns if ns.shape[1] > 0 else np.zeros((dim_Ap, 0))

def build_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx:
                M[idx[face], j] += sign
    return M

def full_chain_complex(A, n, max_dim=None):
    """Compute the full chain complex: all Omega dims, boundary maps, ranks."""
    if max_dim is None: max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

    omega = {}
    omega_dims = []
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        omega_dims.append(omega[p].shape[1] if omega[p].ndim == 2 else 0)

    # Boundary maps restricted to Omega
    bd_ranks = []
    kernels = []
    images = []
    betti = []

    for p in range(max_dim + 1):
        dim_om = omega_dims[p]
        if dim_om == 0:
            bd_ranks.append(0)
            kernels.append(0)
            images.append(0)
            betti.append(0)
            continue

        # d_p: Omega_p -> Omega_{p-1}
        bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_om = bd_p @ omega[p]
        if bd_p_om.size > 0:
            sv = np.linalg.svd(bd_p_om, compute_uv=False)
            rk = sum(s > 1e-8 for s in sv)
        else:
            rk = 0
        bd_ranks.append(rk)
        ker = dim_om - rk
        kernels.append(ker)

        # im(d_{p+1}): Omega_{p+1} -> Omega_p
        dim_om1 = omega_dims[p+1]
        if dim_om1 > 0:
            bd1 = build_boundary_matrix(allowed[p+1], allowed[p])
            bd1_om = bd1 @ omega[p+1]
            sv1 = np.linalg.svd(bd1_om, compute_uv=False)
            im = sum(s > 1e-8 for s in sv1)
        else:
            im = 0
        images.append(im)
        betti.append(ker - im)

    return {
        'omega_dims': omega_dims[:max_dim+1],
        'bd_ranks': bd_ranks,
        'kernels': kernels,
        'images': images,
        'betti': betti,
    }

def is_strongly_connected(A, n):
    visited = set()
    stack = [0]
    while stack:
        v = stack.pop()
        if v in visited: continue
        visited.add(v)
        for u in range(n):
            if A[v][u] and u not in visited:
                stack.append(u)
    if len(visited) != n: return False
    visited = set()
    stack = [0]
    while stack:
        v = stack.pop()
        if v in visited: continue
        visited.add(v)
        for u in range(n):
            if A[u][v] and u not in visited:
                stack.append(u)
    return len(visited) == n

def main():
    print("=" * 70)
    print("BETA_1 * BETA_3 = 0: Proof investigation")
    print("=" * 70)

    # Part 1: At n=6, study the chain complex structure for beta_1>0 vs beta_3>0
    n = 6
    total = 2**(n*(n-1)//2)

    b1_tours = []
    b3_tours = []
    trivial_tours = []

    for bits in range(total):
        A = bits_to_adj(bits, n)
        cc = full_chain_complex(A, n)
        betti = cc['betti']

        b1 = betti[1] if len(betti) > 1 else 0
        b3 = betti[3] if len(betti) > 3 else 0

        if b1 > 0:
            b1_tours.append((bits, cc))
        elif b3 > 0:
            b3_tours.append((bits, cc))
        else:
            trivial_tours.append((bits, cc))

    print(f"\n  n=6: {len(b1_tours)} beta_1>0, {len(b3_tours)} beta_3>0, {len(trivial_tours)} trivial")

    # Part 2: Compare chain complex structure
    print("\n--- Chain complex structure comparison ---")

    def summarize_cc(tours, label):
        ods = defaultdict(int)
        kers = defaultdict(int)
        ims = defaultdict(int)
        for _, cc in tours[:50]:
            od_key = tuple(cc['omega_dims'])
            ods[od_key] += 1
            ker_key = tuple(cc['kernels'])
            kers[ker_key] += 1
            im_key = tuple(cc['images'])
            ims[im_key] += 1

        print(f"\n  {label} ({len(tours)} tours):")
        print(f"    Omega dims distribution (top 5):")
        for od, cnt in sorted(ods.items(), key=lambda x: -x[1])[:5]:
            print(f"      {list(od)}: {cnt}")
        print(f"    ker(d_p) distribution (top 5):")
        for k, cnt in sorted(kers.items(), key=lambda x: -x[1])[:5]:
            print(f"      {list(k)}: {cnt}")
        print(f"    im(d_{'{p+1}'}) distribution (top 5):")
        for im, cnt in sorted(ims.items(), key=lambda x: -x[1])[:5]:
            print(f"      {list(im)}: {cnt}")

    summarize_cc(b1_tours, "beta_1 > 0")
    summarize_cc(b3_tours, "beta_3 > 0")
    summarize_cc(trivial_tours[:100], "trivial (first 100)")

    # Part 3: KEY QUESTION — what differs between beta_1>0 and beta_3>0?
    # Focus on Omega_3 dimension and ker(d_3) vs im(d_4)
    print("\n--- Part 3: Structural difference ---")

    print("\n  beta_1>0 tours:")
    for bits, cc in b1_tours[:10]:
        print(f"    bits={bits}: Omega={cc['omega_dims']}, ker={cc['kernels']}, im={cc['images']}, betti={cc['betti']}")

    print("\n  beta_3>0 tours:")
    for bits, cc in b3_tours[:10]:
        print(f"    bits={bits}: Omega={cc['omega_dims']}, ker={cc['kernels']}, im={cc['images']}, betti={cc['betti']}")

    # Part 4: Look at the DEFICIENCY that creates beta_1
    # beta_1 = ker(d_1) - im(d_2) = C(n-1,2) - rank(d_2|Omega_2)
    # So beta_1 = 1 means rank(d_2|Omega_2) = C(n-1,2) - 1
    # And beta_3 = ker(d_3) - im(d_4)
    # When beta_1=1: is ker(d_3) = im(d_4) forced? Why?
    print("\n--- Part 4: Why beta_1=1 forces beta_3=0 ---")

    # Check: for beta_1>0 tournaments, what is ker(d_3) - im(d_4)?
    gap3_b1 = []
    for bits, cc in b1_tours:
        gap = cc['kernels'][3] - cc['images'][3]
        gap3_b1.append(gap)

    gap3_b3 = []
    for bits, cc in b3_tours:
        gap = cc['kernels'][3] - cc['images'][3]
        gap3_b3.append(gap)

    print(f"  beta_1>0: ker(d_3)-im(d_4) values: {sorted(set(gap3_b1))}")
    print(f"  beta_3>0: ker(d_3)-im(d_4) values: {sorted(set(gap3_b3))}")

    # Part 5: Does beta_1>0 imply ker(d_3) = im(d_4)?
    print(f"\n  beta_1>0: ALWAYS ker(d_3) = im(d_4)? {all(g == 0 for g in gap3_b1)}")
    print(f"  beta_3>0: ALWAYS ker(d_1) = im(d_2)? ", end="")
    gap1_b3 = [cc['kernels'][1] - cc['images'][1] for _, cc in b3_tours]
    print(f"{all(g == 0 for g in gap1_b3)}")

    # Part 6: Track im(d_2) for beta_1>0 vs beta_3>0
    print("\n--- Part 5: im(d_2) comparison ---")
    im2_b1 = [cc['images'][1] for _, cc in b1_tours]
    im2_b3 = [cc['images'][1] for _, cc in b3_tours]
    im2_triv = [cc['images'][1] for _, cc in trivial_tours[:200]]

    print(f"  beta_1>0: im(d_2) values: {sorted(set(im2_b1))}, always = {comb(n-1,2)-1}? {all(v == comb(n-1,2)-1 for v in im2_b1)}")
    print(f"  beta_3>0: im(d_2) values: {sorted(set(im2_b3))}")
    print(f"  trivial: im(d_2) values: {sorted(set(im2_triv))}")

    # Part 7: Is there a structural constraint linking im(d_2) and ker(d_3)?
    # If im(d_2) is lower (fewer boundaries killing 1-cycles), does that
    # free up ker(d_3)?
    print("\n--- Part 6: Correlation between im(d_2) and ker(d_3) ---")
    all_data = []
    for bits, cc in b1_tours + b3_tours + trivial_tours[:500]:
        im2 = cc['images'][1]
        ker3 = cc['kernels'][3] if len(cc['kernels']) > 3 else 0
        im4 = cc['images'][3] if len(cc['images']) > 3 else 0
        b1 = cc['betti'][1]
        b3 = cc['betti'][3] if len(cc['betti']) > 3 else 0
        all_data.append((im2, ker3, im4, b1, b3))

    by_im2 = defaultdict(list)
    for im2, ker3, im4, b1, b3 in all_data:
        by_im2[im2].append((ker3, im4, b1, b3))

    for im2_val in sorted(by_im2.keys()):
        entries = by_im2[im2_val]
        b1_rate = sum(1 for _, _, b1, _ in entries if b1 > 0) / len(entries)
        b3_rate = sum(1 for _, _, _, b3 in entries if b3 > 0) / len(entries)
        avg_ker3 = np.mean([k for k, _, _, _ in entries])
        avg_im4 = np.mean([i for _, i, _, _ in entries])
        if len(entries) >= 3:
            print(f"  im(d_2)={im2_val}: {len(entries)} tours, P(b1>0)={100*b1_rate:.1f}%, P(b3>0)={100*b3_rate:.1f}%, avg ker(d_3)={avg_ker3:.1f}, avg im(d_4)={avg_im4:.1f}")

    # Part 8: Strong connectivity constraint
    print("\n--- Part 7: Strong connectivity role ---")
    sc_b1 = sum(1 for bits, _ in b1_tours if is_strongly_connected(bits_to_adj(bits, n), n))
    sc_b3 = sum(1 for bits, _ in b3_tours if is_strongly_connected(bits_to_adj(bits, n), n))
    sc_triv = sum(1 for bits, _ in trivial_tours if is_strongly_connected(bits_to_adj(bits, n), n))

    print(f"  beta_1>0: {100*sc_b1/len(b1_tours):.1f}% SC")
    print(f"  beta_3>0: {100*sc_b3/len(b3_tours):.1f}% SC")
    print(f"  trivial:  {100*sc_triv/len(trivial_tours):.1f}% SC")

if __name__ == '__main__':
    main()
