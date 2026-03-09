"""
Why is beta_2 = 0 for all tournaments?

Beta_2 = ker(d_2) - im(d_3) in the Omega chain complex.

The chain complex Omega_* for tournaments has a remarkable property:
ALL even Betti numbers appear to vanish. This needs explanation.

Key structural observation: for tournaments, Omega_0 = Z^n (vertices),
Omega_1 = Z^{C(n,2)} (all directed edges are allowed since tournaments
are complete), and Omega_2 depends on the orientation.

For the TRANSITIVE tournament: Omega_p = A_p (all paths are allowed),
dim(Omega_p) = C(n,p+1), and all Betti numbers are 0 (contractible).

For a general tournament T, the defects from transitivity create homology.
The question is: why only in odd dimensions?

Possible explanations:
1. Poincare duality type argument (but path homology isn't simplicial)
2. Transfer/symmetry argument
3. The chi(T) = 1 constraint forces even Betti to vanish
   (chi = sum (-1)^p beta_p, and if chi=1 and only odd bettis nonzero,
    then sum (-1)^{2k+1} beta_{2k+1} = 1 - 1 = 0,
    but chi=1 doesn't force this by itself)

Let me investigate the STRUCTURE of ker(d_2) and im(d_3) directly.
"""

import numpy as np
from math import comb
from collections import defaultdict

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

def main():
    print("=" * 70)
    print("BETA_2 VANISHING: WHY?")
    print("=" * 70)

    # Part 1: Detailed chain complex at n=5,6
    for n in [5, 6]:
        print(f"\n{'='*60}")
        print(f"n={n}: Detailed chain complex analysis for beta_2")
        print(f"{'='*60}")

        if n <= 5:
            N = 2**(n*(n-1)//2)
            exhaustive = True
        else:
            N = 2**(n*(n-1)//2)
            exhaustive = True

        # Track ker(d_2), im(d_3), and their relation
        ker2_vals = []
        im3_vals = []
        om2_vals = []
        om3_vals = []
        rank_d2_vals = []

        for trial in range(N):
            A = bits_to_adj(trial, n)

            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            # Compute Omega bases
            ob = {}
            for p in range(n):
                ob[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

            dim_om = [ob[p].shape[1] if ob[p].ndim == 2 else 0 for p in range(n)]

            # d_2: Omega_2 -> Omega_1
            bd2 = build_boundary_matrix(allowed[2], allowed[1])
            if dim_om[2] > 0:
                bd2_om = bd2 @ ob[2]
                sv2 = np.linalg.svd(bd2_om, compute_uv=False)
                rank_d2 = sum(s > 1e-8 for s in sv2)
            else:
                rank_d2 = 0
            ker_d2 = dim_om[2] - rank_d2

            # d_3: Omega_3 -> Omega_2
            if dim_om[3] > 0 and 3 < n:
                bd3 = build_boundary_matrix(allowed[3], allowed[2])
                bd3_om = bd3 @ ob[3]
                sv3 = np.linalg.svd(bd3_om, compute_uv=False)
                im_d3 = sum(s > 1e-8 for s in sv3)
            else:
                im_d3 = 0

            beta2 = ker_d2 - im_d3
            if beta2 != 0:
                print(f"  *** BETA_2 NONZERO: trial {trial}, beta_2 = {beta2}")

            ker2_vals.append(ker_d2)
            im3_vals.append(im_d3)
            om2_vals.append(dim_om[2])
            om3_vals.append(dim_om[3] if 3 < n else 0)
            rank_d2_vals.append(rank_d2)

        print(f"\n  {N} tournaments checked, ALL beta_2 = 0: {all(k == i for k, i in zip(ker2_vals, im3_vals))}")
        print(f"\n  ker(d_2) distribution: {sorted(set(ker2_vals))}")
        print(f"  im(d_3) distribution: {sorted(set(im3_vals))}")
        print(f"  dim(Omega_2) distribution: {sorted(set(om2_vals))}")
        print(f"  rank(d_2) distribution: {sorted(set(rank_d2_vals))}")

        # Key question: is ker(d_2) = im(d_3) by a STRUCTURAL reason?
        # Compute the actual relationship
        pairs = defaultdict(int)
        for k, i in zip(ker2_vals, im3_vals):
            pairs[(k, i)] += 1
        print(f"\n  (ker(d_2), im(d_3)) pairs:")
        for (k, i), cnt in sorted(pairs.items()):
            print(f"    ({k}, {i}): {cnt} times, beta_2 = {k-i}")

    # Part 2: Is there a dimensional reason?
    # In the transitive tournament: ker(d_p) = im(d_{p+1}) for all p (acyclicity).
    # For a general tournament, defects from transitivity create non-trivial homology.
    # But why only at ODD p?
    #
    # Consider the Euler characteristic:
    # chi = sum (-1)^p beta_p = sum (-1)^p dim(Omega_p) - 2 * sum (-1)^p im(d_{p+1})
    # No, that's not right. chi = sum (-1)^p dim(Omega_p) always.
    # And chi = sum (-1)^p beta_p.
    #
    # If chi = 1 (as verified for n<=7):
    # 1 = beta_0 - beta_1 + beta_2 - beta_3 + ...
    # = 1 - beta_1 + 0 - beta_3 + 0 - beta_5 + ...
    # = 1 - (beta_1 + beta_3 + beta_5 + ...)
    # => beta_1 + beta_3 + beta_5 + ... = 0
    #
    # But wait, Betti numbers are non-negative! So chi = 1 AND all even betti = 0
    # would force ALL odd betti = 0 too. That contradicts the fact that beta_1 and
    # beta_3 CAN be nonzero.
    #
    # So chi is NOT always 1 when beta_p > 0 for some odd p.
    # From S44: chi in {0, 1} at n<=7, with chi=0 iff some odd beta > 0.
    # So: chi = 1 - (beta_1 + beta_3 + ...) if even betti vanish.
    # Since we observe chi in {0,1}: beta_1 + beta_3 + ... in {0, 1}.
    # Combined with beta_1*beta_3 = 0: at most ONE odd betti is nonzero, and it equals 1.

    print("\n" + "=" * 70)
    print("Part 3: THE BIG PICTURE")
    print("=" * 70)
    print("""
If all even Betti vanish (beta_0 = 1 excepted, beta_2 = beta_4 = ... = 0):

  chi = 1 - beta_1 - beta_3 - beta_5 - ...

Observed: chi in {0, 1} at n <= 7.
  chi = 0 <=> exactly ONE odd beta = 1
  chi = 1 <=> all beta = 0

Combined with seesaw (adjacent odd bettis exclusive):
  AT MOST ONE odd beta is nonzero, and when nonzero it equals 1.

This means the path homology of a tournament is EXTRAORDINARILY SIMPLE:
  H_p(T) = Z  for exactly one odd p (or none)
  H_0(T) = Z  always (connected)
  All other H_p = 0

The tournament path homology is either CONTRACTIBLE or has a
SINGLE odd-dimensional hole of size 1.
""")

    # Part 3: Verify this at n=7
    print("Verification at n=7 (sampled):")
    n = 7
    rng = np.random.RandomState(42)
    N = 2000
    betti_profiles = defaultdict(int)

    for trial in range(N):
        A = random_tournament(n, rng)

        allowed = {}
        for p in range(-1, n):
            allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

        ob = {}
        for p in range(n):
            ob[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        dim_om = [ob[p].shape[1] if ob[p].ndim == 2 else 0 for p in range(n)]

        bettis = []
        for target_p in range(n):
            if dim_om[target_p] == 0:
                bettis.append(0)
                continue
            bd_p = build_boundary_matrix(allowed[target_p], allowed[target_p-1])
            bd_p_om = bd_p @ ob[target_p]
            if bd_p_om.size > 0:
                sv = np.linalg.svd(bd_p_om, compute_uv=False)
                rk = sum(s > 1e-8 for s in sv)
            else:
                rk = 0
            ker = dim_om[target_p] - rk

            if target_p + 1 < n and dim_om[target_p+1] > 0:
                bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
                bd1_om = bd1 @ ob[target_p+1]
                sv1 = np.linalg.svd(bd1_om, compute_uv=False)
                im = sum(s > 1e-8 for s in sv1)
            else:
                im = 0
            bettis.append(ker - im)

        profile = tuple(bettis)
        betti_profiles[profile] += 1

    print(f"  {N} random n=7 tournaments:")
    for profile, cnt in sorted(betti_profiles.items(), key=lambda x: -x[1]):
        nonzero = [(p, b) for p, b in enumerate(profile) if b != 0 and p > 0]
        desc = "contractible" if not nonzero else f"H_{nonzero[0][0]}=Z^{nonzero[0][1]}" if len(nonzero) == 1 else str(nonzero)
        chi = sum((-1)**p * b for p, b in enumerate(profile))
        print(f"    {list(profile)}: {cnt} ({100*cnt/N:.1f}%), chi={chi}, type={desc}")

if __name__ == '__main__':
    main()
