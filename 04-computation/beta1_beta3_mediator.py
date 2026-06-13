"""
Beta_1 * Beta_3 = 0: The im(d_2) mediator hypothesis

At n=6 we discovered:
  beta_1 > 0 <=> im(d_2) = 9 (drops by 1 from max of 10)
  beta_3 > 0 => im(d_2) = 10 (saturated)

Since im(d_2) can't be both 9 and 10, mutual exclusion follows.

Question: Does this generalize to n=7, 8?
Specifically: is there a single rank quantity that separates
beta_1>0 and beta_3>0 tournaments at every n?

Also: what is im(d_2) in terms of the tournament?
im(d_2) = rank of d_2 restricted to Omega_2.
d_2: Omega_2 -> Omega_1.
The rank deficit ker(d_2)/im(d_3) gives beta_2.

Since beta_2=0 always (at n<=8), we have ker(d_2)=im(d_3).
So rank(d_2|Omega_2) = dim(Omega_2) - ker(d_2) = dim(Omega_2) - im(d_3).

The chain complex is: ... -> Omega_3 -d3-> Omega_2 -d2-> Omega_1 -d1-> Omega_0 -> 0

beta_1 = ker(d_1) - im(d_2)
beta_3 = ker(d_3) - im(d_4)

If beta_2=0: ker(d_2)=im(d_3), so
  im(d_2) = dim(Omega_2) - im(d_3)

The critical insight: im(d_2) = dim(Omega_2) - im(d_3).
beta_1>0 iff im(d_2) < ker(d_1)
  iff dim(Omega_2) - im(d_3) < ker(d_1)

And ker(d_1) = dim(Omega_1) always for p=1 (since Omega_0 = vertices, d_1 is
surjective iff... actually d_1: Omega_1 -> Omega_0 has im(d_1) = n-1 for
connected tournaments). Actually:
beta_0 = dim(Omega_0) - im(d_1) = n - im(d_1)
For strongly connected: beta_0 = 1, so im(d_1) = n-1.
For non-SC with k components: beta_0 = k.

So ker(d_1) = dim(Omega_1) - im(d_1) = dim(Omega_1) - (n - beta_0).

And beta_1 = ker(d_1) - im(d_2)
          = dim(Omega_1) - n + beta_0 - im(d_2)
          = dim(Omega_1) - n + beta_0 - dim(Omega_2) + im(d_3)

This is getting algebraic. Let me just check computationally at n=7.
"""

import numpy as np
from math import comb

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

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

def full_chain_complex(A, n, max_p=None):
    """Compute full chain complex data: omega dims, ker, im, betti."""
    if max_p is None:
        max_p = n - 1

    allowed = {}
    for p in range(-1, max_p + 2):
        if p < 0:
            allowed[p] = []
        elif p < n:
            allowed[p] = enumerate_allowed_paths(A, n, p)
        else:
            allowed[p] = []

    omega_bases = {}
    omega_dims = []
    for p in range(max_p + 1):
        ob = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        omega_bases[p] = ob
        omega_dims.append(ob.shape[1] if ob.ndim == 2 else 0)

    # Compute boundary maps restricted to Omega
    # d_p: Omega_p -> Omega_{p-1}
    kers = []
    ims = []
    bettis = []

    for p in range(max_p + 1):
        # Compute im(d_{p+1}) and ker(d_p)
        # ker(d_p): rank of d_p restricted to Omega_p
        if p == 0:
            # d_0 = 0, ker(d_0) = Omega_0
            ker_p = omega_dims[0]
        else:
            bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
            if omega_dims[p] > 0:
                bd_p_om = bd_p @ omega_bases[p]
                sv = np.linalg.svd(bd_p_om, compute_uv=False)
                rk = sum(s > 1e-8 for s in sv)
            else:
                rk = 0
            ker_p = omega_dims[p] - rk
        kers.append(ker_p)

        # im(d_{p+1}): rank of d_{p+1} restricted to Omega_{p+1}
        if p + 1 <= max_p:
            bd_p1 = build_boundary_matrix(allowed[p+1], allowed[p])
            ob_p1 = omega_bases.get(p+1, np.zeros((0, 0)))
            dim_p1 = ob_p1.shape[1] if ob_p1.ndim == 2 else 0
            if dim_p1 > 0:
                bd_p1_om = bd_p1 @ ob_p1
                sv = np.linalg.svd(bd_p1_om, compute_uv=False)
                im_p = sum(s > 1e-8 for s in sv)
            else:
                im_p = 0
        else:
            im_p = 0
        ims.append(im_p)

        bettis.append(ker_p - im_p)

    return {
        'omega_dims': omega_dims,
        'kers': kers,
        'ims': ims,
        'bettis': bettis,
        'allowed': allowed
    }

def main():
    print("=" * 70)
    print("BETA_1 * BETA_3 = 0: im(d_2) MEDIATOR HYPOTHESIS")
    print("=" * 70)

    # Part 1: Verify at n=6 (exhaustive) and extend to n=7 (sampled)
    for n in [6, 7, 8]:
        print(f"\n{'='*60}")
        print(f"n={n}")
        print(f"{'='*60}")

        if n <= 6:
            N = 2**(n*(n-1)//2)
            exhaustive = True
        else:
            N = 2000
            exhaustive = False

        rng = np.random.RandomState(42 + n)

        # Track im(d_2) for each category
        from collections import defaultdict
        im2_b1 = []  # im(d_2) for beta_1 > 0
        im2_b3 = []  # im(d_2) for beta_3 > 0
        im2_triv = []  # im(d_2) for trivial
        im2_both = []  # im(d_2) for both > 0 (should be empty)

        # Also track: ker(d_1), ker(d_3)
        ker1_b1 = []
        ker1_b3 = []
        ker3_b1 = []
        ker3_b3 = []

        # Also: dim(Omega_2), im(d_3)
        om2_b1 = []
        om2_b3 = []
        imd3_b1 = []
        imd3_b3 = []

        cnt = 0
        b1_cnt = 0
        b3_cnt = 0

        for trial in range(N):
            if exhaustive:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            data = full_chain_complex(A, n, max_p=min(n-1, 5))
            b1 = data['bettis'][1] if len(data['bettis']) > 1 else 0
            b3 = data['bettis'][3] if len(data['bettis']) > 3 else 0

            im2 = data['ims'][1] if len(data['ims']) > 1 else 0
            ker1 = data['kers'][1] if len(data['kers']) > 1 else 0
            ker3 = data['kers'][3] if len(data['kers']) > 3 else 0
            im3 = data['ims'][2] if len(data['ims']) > 2 else 0
            om2 = data['omega_dims'][2] if len(data['omega_dims']) > 2 else 0

            if b1 > 0 and b3 > 0:
                im2_both.append(im2)
            elif b1 > 0:
                im2_b1.append(im2)
                ker1_b1.append(ker1)
                ker3_b1.append(ker3)
                om2_b1.append(om2)
                imd3_b1.append(im3)
                b1_cnt += 1
            elif b3 > 0:
                im2_b3.append(im2)
                ker1_b3.append(ker1)
                ker3_b3.append(ker3)
                om2_b3.append(om2)
                imd3_b3.append(im3)
                b3_cnt += 1
            else:
                im2_triv.append(im2)

            cnt += 1
            if not exhaustive and cnt % 500 == 0:
                print(f"  {cnt}/{N}: b1={b1_cnt}, b3={b3_cnt}")

        print(f"\n  Total: {cnt}, beta_1>0: {b1_cnt}, beta_3>0: {b3_cnt}, both>0: {len(im2_both)}")

        if im2_b1:
            print(f"\n  im(d_2) for beta_1>0: {sorted(set(im2_b1))}")
        if im2_b3:
            print(f"  im(d_2) for beta_3>0: {sorted(set(im2_b3))}")
        if im2_triv:
            print(f"  im(d_2) for trivial: {sorted(set(im2_triv))}")
        if im2_both:
            print(f"  *** BOTH > 0: im(d_2) = {sorted(set(im2_both))}")

        # Check: do the sets overlap?
        set_b1 = set(im2_b1)
        set_b3 = set(im2_b3)
        set_triv = set(im2_triv)

        overlap = set_b1 & set_b3
        print(f"\n  im(d_2) overlap between b1>0 and b3>0: {overlap if overlap else 'NONE (disjoint!)'}")
        print(f"  im(d_2) overlap between b1>0 and triv: {set_b1 & set_triv if set_b1 & set_triv else 'NONE'}")
        print(f"  im(d_2) overlap between b3>0 and triv: {set_b3 & set_triv if set_b3 & set_triv else 'NONE'}")

        # Maximum possible im(d_2)
        trans_ker1 = comb(n, 2) - (n - 1)  # for transitive: ker(d_1) = C(n,2) - (n-1)
        print(f"\n  Transitive im(d_2) would be: ker(d_1) = C(n,2) - (n-1) = {trans_ker1}")

        if im2_b1:
            print(f"\n  ker(d_1) for beta_1>0: {sorted(set(ker1_b1))}")
            print(f"  ker(d_3) for beta_1>0: {sorted(set(ker3_b1))[:10]}")
            print(f"  dim(Omega_2) for beta_1>0: {sorted(set(om2_b1))[:10]}")
        if im2_b3:
            print(f"\n  ker(d_1) for beta_3>0: {sorted(set(ker1_b3))}")
            print(f"  ker(d_3) for beta_3>0: {sorted(set(ker3_b3))[:10]}")
            print(f"  dim(Omega_2) for beta_3>0: {sorted(set(om2_b3))[:10]}")

        # Part 2: Deeper look — what determines im(d_2)?
        # im(d_2) = rank(d_2: Omega_2 -> Omega_1)
        # = dim(Omega_2) - ker(d_2)
        # Since beta_2=0 always (checked): ker(d_2)=im(d_3)
        # So im(d_2) = dim(Omega_2) - im(d_3)
        if im2_b1 and imd3_b1:
            # Check: im(d_2) = dim(Omega_2) - im(d_3)?
            check = all(im2_b1[i] == om2_b1[i] - imd3_b1[i] for i in range(len(im2_b1)))
            print(f"\n  Verification: im(d_2) = dim(Omega_2) - im(d_3) for beta_1>0: {check}")

        # Part 3: Does ker(d_1) vary or is it fixed?
        print(f"\n  ker(d_1) values across ALL tournaments:")
        all_ker1 = ker1_b1 + ker1_b3
        if im2_triv:
            # Need to also track ker1 for trivial
            pass
        if all_ker1:
            print(f"    beta_1>0: {sorted(set(ker1_b1))}")
            print(f"    beta_3>0: {sorted(set(ker1_b3))}")

    # Part 4: The algebraic constraint
    print("\n" + "=" * 70)
    print("ALGEBRAIC INTERPRETATION")
    print("=" * 70)
    print("""
At n=6:
  ker(d_1) = 10 always (= C(6,2) - 5 = 10, since all n=6 are SC? No...)
  Actually: ker(d_1) = dim(Omega_1) - im(d_1) = 15 - im(d_1)

  For beta_1 > 0: im(d_2) = 9, ker(d_1) = 10, so beta_1 = 10 - 9 = 1
  For beta_3 > 0: im(d_2) = 10, ker(d_1) = 10, so beta_1 = 10 - 10 = 0

  The KEY is: why does beta_3>0 force im(d_2) = 10 (maximal)?

  Hypothesis: beta_3>0 means ker(d_3) > im(d_4), which means d_3 has
  a larger kernel, which means im(d_3) is smaller, which forces
  im(d_2) = dim(Omega_2) - im(d_3) to be larger (since beta_2=0).

  In the chain Omega_3 -d3-> Omega_2 -d2-> Omega_1:
  - beta_2=0 means ker(d_2) = im(d_3) exactly
  - So rank(d_2|Omega_2) = dim(Omega_2) - im(d_3)
  - If im(d_3) drops (due to beta_3>0 stealing kernel from d_3),
    then im(d_2) increases!
  - And im(d_2) increasing makes beta_1 = ker(d_1) - im(d_2) smaller.

  This is a SEESAW EFFECT through the exact sequence at Omega_2!
  beta_3 up => ker(d_3) up => im(d_3) down => im(d_2) up => beta_1 down

  The constraint is: beta_2 = 0.
  If beta_2 were nonzero, the seesaw would have slack and both could
  be nonzero simultaneously.

  THEOREM CANDIDATE: If beta_2(T) = 0 for all tournaments T on n vertices,
  then beta_1 * beta_3 = 0 for all such T.
""")

if __name__ == '__main__':
    main()
