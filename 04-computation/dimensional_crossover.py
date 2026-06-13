"""
Dimensional crossover: Why does beta_3 overtake beta_1?

The crossover between P(beta_1>0) and P(beta_3>0) happens at d=6 (n=7).
This is the SAME dimension where Paley T_7 first exhibits beta_4 > 0.

Hypothesis: As dimension increases, the path complex has "room" for
higher-dimensional homological features but fewer 1-dimensional holes.

Structural investigation:
1. What fraction of 3-element subsets form 3-cycles? (= P(cyclic triple))
   This determines |A_2| - C(n,3) = excess paths at dimension 2.
2. How does the "cyclic fraction" change with n?
3. Does beta_3 correlate with the cyclic fraction?
4. What is the relationship between beta_1 and strong connectivity?

Key observation from MEMORY:
- beta_1 requires SC (100% at n=6)
- beta_3 does NOT require SC
- At larger n, almost all tournaments are SC => beta_1 should be MORE
  prevalent, not less. So the decrease in beta_1 is NOT about SC.

The decrease in beta_1 might be about the Omega_1 constraint becoming
MORE restrictive as n grows, killing 1-cycles.
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

def betti_single(A, n, target_p):
    allowed = {}
    for p in [target_p - 1, target_p, target_p + 1]:
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega_p = compute_omega_basis(A, n, target_p, allowed[target_p], allowed[target_p-1])
    dim_om = omega_p.shape[1] if omega_p.ndim == 2 else 0
    if dim_om == 0:
        return 0, 0, 0, 0  # beta, dim_omega, ker, im
    bd_p = build_boundary_matrix(allowed[target_p], allowed[target_p-1])
    bd_p_om = bd_p @ omega_p
    if bd_p_om.size > 0:
        sv = np.linalg.svd(bd_p_om, compute_uv=False)
        rk = sum(s > 1e-8 for s in sv)
    else:
        rk = 0
    ker = dim_om - rk
    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed[target_p+1], allowed[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = sum(s > 1e-8 for s in sv1)
    else:
        im = 0
    return ker - im, dim_om, ker, im

def count_cycles(A, n, k):
    count = 0
    for subset in combinations(range(n), k):
        # Check all cyclic permutations
        from itertools import permutations
        for perm in permutations(subset):
            if all(A[perm[i]][perm[(i+1) % k]] for i in range(k)):
                count += 1
    return count // k

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
    print("DIMENSIONAL CROSSOVER: beta_1 vs beta_3 mechanism")
    print("=" * 70)

    # Part 1: Track ker, im decomposition of beta_1 and beta_3
    # beta_p = ker(d_p) - im(d_{p+1}) within Omega_p
    print("\n--- Part 1: ker/im decomposition across dimensions ---")

    for n in range(5, 9):
        rng = np.random.RandomState(42 + n)
        N = min(300, 2**(n*(n-1)//2))

        b1_stats = {'pos': 0, 'dim_omega': [], 'ker': [], 'im': []}
        b3_stats = {'pos': 0, 'dim_omega': [], 'ker': [], 'im': []}
        sc_count = 0
        c3_vals = []
        total = 0

        for _ in range(N):
            if n <= 5 and _ < 2**(n*(n-1)//2):
                A = bits_to_adj(_, n)
            else:
                A = random_tournament(n, rng)

            sc = is_strongly_connected(A, n)
            if sc: sc_count += 1

            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1
            c3_vals.append(c3)

            b1, dim_o1, ker1, im1 = betti_single(A, n, 1)
            b1_stats['dim_omega'].append(dim_o1)
            b1_stats['ker'].append(ker1)
            b1_stats['im'].append(im1)
            if b1 > 0: b1_stats['pos'] += 1

            if n >= 6:
                b3, dim_o3, ker3, im3 = betti_single(A, n, 3)
                b3_stats['dim_omega'].append(dim_o3)
                b3_stats['ker'].append(ker3)
                b3_stats['im'].append(im3)
                if b3 > 0: b3_stats['pos'] += 1

            total += 1

        print(f"\n  n={n} (d={n-1}), {total} samples:")
        print(f"    SC fraction: {100*sc_count/total:.1f}%")
        print(f"    avg c3: {np.mean(c3_vals):.1f}, c3/C(n,3) = {np.mean(c3_vals)/comb(n,3):.3f}")

        print(f"    beta_1: P(>0) = {100*b1_stats['pos']/total:.2f}%")
        print(f"      avg dim(Omega_1) = {np.mean(b1_stats['dim_omega']):.1f}")
        print(f"      avg ker(d_1) = {np.mean(b1_stats['ker']):.2f}")
        print(f"      avg im(d_2) = {np.mean(b1_stats['im']):.2f}")
        print(f"      avg ker-im = {np.mean([k-i for k,i in zip(b1_stats['ker'], b1_stats['im'])]):.3f}")

        if n >= 6:
            print(f"    beta_3: P(>0) = {100*b3_stats['pos']/total:.2f}%")
            print(f"      avg dim(Omega_3) = {np.mean(b3_stats['dim_omega']):.1f}")
            print(f"      avg ker(d_3) = {np.mean(b3_stats['ker']):.2f}")
            print(f"      avg im(d_4) = {np.mean(b3_stats['im']):.2f}")
            print(f"      avg ker-im = {np.mean([k-i for k,i in zip(b3_stats['ker'], b3_stats['im'])]):.3f}")

    # Part 2: The "dimensional surplus" — how many extra constraints
    # does Omega_1 have vs the simplex?
    print("\n" + "=" * 70)
    print("Part 2: Omega_1 constraints vs simplex dimension")
    print("=" * 70)
    print("For transitive: dim(Omega_1) = C(n,2). For general: dim(Omega_1) = C(n,2).")
    print("(All edges ARE allowed 1-paths in a tournament.)")
    print("So Omega_1 dimension is ALWAYS C(n,2) regardless of tournament.")
    print("beta_1 = ker(d_1) - im(d_2).")
    print("d_1: Omega_1 -> Omega_0 = R^n has matrix n x C(n,2).")
    print("rank(d_1) = n-1 always (d maps edge (a,b) to b-a).")
    print("ker(d_1) = C(n,2) - (n-1) = C(n-1,2).")
    print("im(d_2) = rank(d_2|Omega_2) = dim(Omega_2) - ker(d_2|Omega_2).")
    print()

    for n in range(3, 9):
        cn2 = comb(n, 2)
        ker_d1 = cn2 - (n - 1)  # = C(n-1, 2)
        print(f"  n={n}: dim(Omega_1)={cn2}, ker(d_1)={ker_d1} = C({n-1},2)")

    print("\nbeta_1 = ker(d_1) - im(d_2) = C(n-1,2) - rank(d_2|Omega_2)")
    print("beta_1 > 0 iff rank(d_2|Omega_2) < C(n-1,2)")
    print("As n grows, C(n-1,2) grows as n^2/2.")
    print("If rank(d_2|Omega_2) also grows as n^2/2 (saturating),")
    print("then beta_1 = 0 becomes more likely.")
    print()

    # Part 3: Why does beta_3 increase?
    print("=" * 70)
    print("Part 3: Why beta_3 increases with dimension")
    print("=" * 70)
    print("beta_3 = ker(d_3|Omega_3) - im(d_4|Omega_4)")
    print("As n grows:")
    print("  - dim(Omega_3) grows (more non-transitive 4-tuples)")
    print("  - ker(d_3) could grow faster than im(d_4)")
    print("  - More 'room' for 3-cycles in the chain complex")
    print()

    # Verify: compute dim(Omega_p), rank(d_p) explicitly at each n
    for n in range(5, 9):
        rng = np.random.RandomState(42 + n)
        N = min(100, 2**(n*(n-1)//2))

        omega_avg = defaultdict(float)
        cnt = 0
        for _ in range(N):
            if n <= 5 and _ < 2**(n*(n-1)//2):
                A = bits_to_adj(_, n)
            else:
                A = random_tournament(n, rng)

            # Compute all Omega dims
            allowed = {}
            for p in range(-1, n):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
            for p in range(n):
                apm1 = allowed[p-1]
                ap = allowed[p]
                omega_p = compute_omega_basis(A, n, p, ap, apm1)
                od = omega_p.shape[1] if omega_p.ndim == 2 else 0
                omega_avg[p] += od
            cnt += 1

        print(f"\n  n={n}: average Omega dimensions (simplex C(n,p+1) in parens)")
        for p in range(n):
            c = comb(n, p+1)
            avg_o = omega_avg[p] / cnt
            ratio = avg_o / c
            surplus = avg_o - c
            print(f"    p={p}: <Omega>={avg_o:.1f}, C={c}, ratio={ratio:.3f}, surplus={surplus:.1f}")

if __name__ == '__main__':
    main()
