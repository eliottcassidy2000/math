"""
wild_ideas_S71m.py -- opus-2026-03-14-S71m
WILD IDEAS SESSION: p-adic H, quantum tournaments, tropical H,
tournament zeta functions, and more

Investigating radically new perspectives on tournament structure.
"""

import sys
import numpy as np
from math import factorial, comb, gcd, log, log2, sqrt, pi
from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def compute_H_dp(adj, n):
    """Count Hamiltonian paths via DP."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def all_tournaments(n):
    """Generate all tournaments on n vertices with H values."""
    m = n*(n-1)//2
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    results = []
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        H = compute_H_dp(adj, n)
        results.append((bits, H, adj))
    return results

def v2(n):
    """2-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v

def vp(n, p):
    """p-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def main():
    print("=" * 70)
    print("WILD IDEAS SESSION")
    print("opus-2026-03-14-S71m")
    print("=" * 70)

    # Precompute tournaments for n=5,6
    print("\n  Computing all tournaments for n=5...")
    tours5 = all_tournaments(5)
    print(f"  Done: {len(tours5)} tournaments")

    # ================================================================
    # IDEA 1: p-ADIC STRUCTURE OF H
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 1: p-ADIC STRUCTURE OF H")
    print(f"{'='*70}")

    print(f"""
  H is always odd (Redei's theorem). So v_2(H) = 0 always.
  But what about v_3, v_5, v_7?

  The 7-adic valuation is especially interesting:
  H=7 is FORBIDDEN, meaning v_7(H) >= 1 is impossible at n<=6.
  But at n=7, H can be divisible by 7 (e.g., H=35=5*7, H=91=7*13).

  QUESTION: Is there a pattern in the p-adic valuations of the H spectrum?
""")

    H_vals_5 = [H for _, H, _ in tours5]
    H_spectrum_5 = sorted(set(H_vals_5))

    print(f"  n=5 H-spectrum: {H_spectrum_5}")
    print(f"\n  p-adic valuations of achievable H values at n=5:")
    for H in H_spectrum_5:
        v3 = vp(H, 3)
        v5 = vp(H, 5)
        v7 = vp(H, 7)
        print(f"    H={H:3d}: v_3={v3}, v_5={v5}, v_7={v7}, "
              f"H mod 8 = {H%8}, H mod 6 = {H%6}")

    # 2-adic: H is always odd, but what about (H-1)/2?
    print(f"\n  2-adic valuation of (H-1)/2:")
    for H in H_spectrum_5:
        hh = (H-1)//2
        print(f"    H={H:3d}: (H-1)/2 = {hh}, v_2((H-1)/2) = {v2(hh)}")

    # ================================================================
    # IDEA 2: TOURNAMENT AS QUANTUM STATE
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 2: TOURNAMENT AS QUANTUM STATE")
    print(f"{'='*70}")

    print(f"""
  Consider the 2^m dimensional Hilbert space where each basis state
  corresponds to a tournament on n vertices (via arc orientations).

  Define |psi_H> = (1/sqrt(N_H)) * sum_{{T: H(T)=h}} |T>

  This is the uniform superposition over all tournaments with H = h.

  The DENSITY MATRIX rho_H = |psi_H><psi_H| encodes the tournament
  structure at fixed H value.

  Question: Are the H-eigenspaces ORTHOGONAL in some natural sense?
  The H operator: H_op |T> = H(T) |T> (diagonal in tournament basis).
  Yes, H-eigenspaces are automatically orthogonal since H_op is Hermitian.

  But more interesting: define the ARC FLIP operator F_e for each arc e.
  F_e |T> = |T xor e> (flip a single arc).

  Then H_op and F_e don't commute! The commutator [H_op, F_e] measures
  how much flipping arc e changes H.
""")

    n = 4
    m = n*(n-1)//2
    tours4 = all_tournaments(n)
    num_tours = len(tours4)

    # Build the H operator and arc flip operators for n=4
    H_vec = np.array([H for _, H, _ in tours4])

    # Compute <H^2> and <H>^2
    mean_H = np.mean(H_vec)
    var_H = np.var(H_vec)
    print(f"\n  n={n}: <H> = {mean_H:.4f}, Var(H) = {var_H:.4f}")

    # Transition matrix: M[i,j] = 1 if tournaments i and j differ by one arc
    # (This is the adjacency matrix of the tournament hypercube)
    print(f"\n  Building hypercube adjacency matrix for n={n}...")
    M = np.zeros((num_tours, num_tours))
    for i in range(num_tours):
        bits_i = tours4[i][0]
        for arc in range(m):
            bits_j = bits_i ^ (1 << arc)
            M[i][bits_j] = 1

    # The Laplacian L = m*I - M
    L = m * np.eye(num_tours) - M

    # H as diagonal operator
    H_diag = np.diag(H_vec.astype(float))

    # Commutator [L, H_diag]
    comm = L @ H_diag - H_diag @ L
    comm_norm = np.linalg.norm(comm, 'fro')
    print(f"  ||[L, H]||_F = {comm_norm:.4f}")
    print(f"  ||L||_F = {np.linalg.norm(L, 'fro'):.4f}")
    print(f"  ||H||_F = {np.linalg.norm(H_diag, 'fro'):.4f}")
    print(f"  Relative commutator = {comm_norm / (np.linalg.norm(L,'fro') * np.linalg.norm(H_diag,'fro')):.6f}")

    # Eigenvalues of L (hypercube graph Laplacian)
    eigs_L = np.linalg.eigvalsh(L)
    print(f"\n  Laplacian eigenvalues: {sorted(set(np.round(eigs_L, 4)))[:10]}")

    # The "quantum walk" propagator exp(-iLt) applied to an H-eigenstate
    # How fast does it spread to other H values?
    # Answer: the off-diagonal elements of exp(-iLt) H exp(iLt) give the
    # time evolution of H in the Heisenberg picture.

    print(f"\n  H-value entanglement under hypercube evolution:")
    # For each H value, compute the overlap with each arc flip
    H_groups = defaultdict(list)
    for i, (bits, H, adj) in enumerate(tours4):
        H_groups[H].append(i)

    for H_val in sorted(H_groups.keys()):
        indices = H_groups[H_val]
        # For each tournament in this group, where do arc flips take us?
        flip_H_dist = Counter()
        for i in indices:
            bits_i = tours4[i][0]
            for arc in range(m):
                bits_j = bits_i ^ (1 << arc)
                H_j = tours4[bits_j][1]
                flip_H_dist[H_j] += 1
        total_flips = len(indices) * m
        print(f"    H={H_val}: {len(indices)} tournaments, flip destinations:", end="")
        for H_target in sorted(flip_H_dist.keys()):
            pct = flip_H_dist[H_target] / total_flips * 100
            print(f" H={H_target}({pct:.0f}%)", end="")
        print()

    # ================================================================
    # IDEA 3: TROPICAL TOURNAMENT THEORY
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 3: TROPICAL TOURNAMENT THEORY")
    print(f"{'='*70}")

    print(f"""
  In TROPICAL (max-plus) algebra: addition = max, multiplication = +.

  The tropical permanent of the adjacency matrix A:
  tperm(A) = max_sigma sum_i A[i][sigma(i)]

  For a tournament: A[i][j] = 0 if i->j, -inf if j->i.
  tperm(A) = max_sigma |{{i: A[i][sigma(i)] = 0}}|
           = maximum number of "correctly oriented" arcs in a permutation.

  But we want Hamiltonian PATHS, not cycles.
  Tropical HP count: max_sigma sum_{{i=1}}^{{n-1}} A[sigma(i)][sigma(i+1)]
  = max "forward arc" count across all orderings.

  The TROPICAL H = maximum number of forward arcs in any single path.
  For a transitive tournament: tropical H = n-1 (all arcs forward).
  For any tournament: tropical H = n-1 (there's always a HP).

  So tropical H is trivial! Everyone gets n-1.

  MORE INTERESTING: The tropical SIGNED count.
  Define w(P) = (number of forward arcs) - (number of backward arcs).
  Then tropical w = max_P w(P) = ??? This is the MAX WEIGHT HP.
""")

    n = 5
    print(f"  n={n}: tropical max-weight HP for each tournament class:")
    H_to_tropical = defaultdict(list)
    for bits, H, adj in tours5:
        max_w = -n
        for perm in permutations(range(n)):
            w = 0
            valid = True
            for i in range(n-1):
                if adj[perm[i]][perm[i+1]]:
                    w += 1
                else:
                    w -= 1
                    # Still a path since tournament has arc either way
            # Actually every permutation gives an HP in a tournament
            # since there's always an arc between any two vertices
            fwd = sum(1 for i in range(n-1) if adj[perm[i]][perm[i+1]])
            bwd = (n-1) - fwd
            weight = fwd - bwd
            max_w = max(max_w, weight)
        H_to_tropical[H].append(max_w)

    for H_val in sorted(H_to_tropical.keys()):
        vals = H_to_tropical[H_val]
        distinct = sorted(set(vals))
        print(f"    H={H_val:3d}: tropical max weight in {distinct}")

    # ================================================================
    # IDEA 4: H AS A MARKOV CHAIN STATIONARY DISTRIBUTION
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 4: H AS MARKOV CHAIN STATIONARY DISTRIBUTION")
    print(f"{'='*70}")

    print(f"""
  Define a random walk on tournaments:
  At each step, choose a random arc and flip it.
  The stationary distribution is UNIFORM (symmetric random walk on hypercube).

  BUT: what if we weight the transition probabilities by H?
  P(T -> T') = H(T') / sum_{{T'' neighbors}} H(T'')

  This creates a METROPOLIS-like chain that favors high-H tournaments.
  The stationary distribution would be proportional to H(T).

  Alternative: The "H-weighted" random walk on the hypercube graph.
  Transition rates proportional to H at the destination.
  The spectral gap of this walk tells us about the "mixing time"
  of tournament optimization.
""")

    # For n=4, compute the H-weighted transition matrix
    n = 4
    m = n*(n-1)//2
    num_T = 1 << m  # 2^6 = 64

    # H values for all n=4 tournaments
    H_4 = np.array([tours4[i][1] for i in range(num_T)], dtype=float)

    # Build H-weighted transition matrix
    P = np.zeros((num_T, num_T))
    for i in range(num_T):
        bits_i = tours4[i][0]
        neighbors = []
        for arc in range(m):
            j = bits_i ^ (1 << arc)
            neighbors.append(j)
        total_H = sum(H_4[j] for j in neighbors)
        for j in neighbors:
            P[i][j] = H_4[j] / total_H

    # Check that P is stochastic
    row_sums = P.sum(axis=1)
    print(f"  n={n}: P is stochastic? max|row_sum - 1| = {max(abs(row_sums - 1)):.2e}")

    # Find stationary distribution
    eigs, vecs = np.linalg.eig(P.T)
    # Find eigenvector with eigenvalue 1
    idx = np.argmin(np.abs(eigs - 1.0))
    pi_stat = np.real(vecs[:, idx])
    pi_stat = pi_stat / pi_stat.sum()

    # Is pi proportional to H?
    pi_vs_H = pi_stat / (H_4 / H_4.sum())
    print(f"  Stationary distribution proportional to H? "
          f"max/min ratio = {max(pi_vs_H)/min(pi_vs_H):.4f}")

    # Actually the stationary distribution of a Metropolis chain is NOT proportional to H
    # for this construction. Let me check what it IS proportional to.
    # The detailed balance condition: pi(i) P(i,j) = pi(j) P(j,i)
    # P(i,j) = H(j) / sum_k H(k) where k are neighbors of i
    # pi(i) * H(j) / Z_i = pi(j) * H(i) / Z_j
    # where Z_i = sum of H over neighbors of i
    # So: pi(i)/pi(j) = H(i)*Z_j / (H(j)*Z_i)
    # This is NOT just H(i)/H(j) unless Z_i = Z_j for all i.

    # Check: is Z_i constant across tournaments?
    Z = np.array([sum(H_4[tours4[i][0] ^ (1 << arc)] for arc in range(m)) for i in range(num_T)])
    print(f"  Z_i (sum of neighbor H values): min={Z.min():.0f}, max={Z.max():.0f}, "
          f"std={Z.std():.2f}")
    print(f"  Z is NOT constant -> stationary dist is NOT proportional to H")

    # What IS Z proportional to? Check correlation with H
    corr_ZH = np.corrcoef(Z, H_4)[0,1]
    print(f"  Corr(Z, H) = {corr_ZH:.4f}")
    print(f"  Interpretation: neighbors of high-H tournaments also tend to have high H")

    # ================================================================
    # IDEA 5: THE JONES POLYNOMIAL ANALOG
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 5: TOURNAMENT BRACKET POLYNOMIAL")
    print(f"{'='*70}")

    print(f"""
  In knot theory, the BRACKET POLYNOMIAL uses the Kauffman state model:
  Each crossing has two smoothings (A-type and B-type).
  <K> = sum over states: A^a * B^b * (-A^2 - A^{{-2}})^{{|s|-1}}

  TOURNAMENT ANALOG:
  Each arc (i,j) has two orientations: i->j (A) and j->i (B).
  A "state" is a tournament.
  Define: <T> = sum over arcs: product of A^{{arc_forward}} * B^{{arc_backward}}
  where forward/backward is relative to some reference ordering.

  With A = q, B = q^{{-1}}:
  <T> = q^{{fwd - bwd}} where fwd = forward arcs, bwd = backward arcs.

  The SCORE SEQUENCE determines fwd and bwd:
  fwd = sum s_i, bwd = C(n,2) - sum s_i.
  So <T> = q^{{2*sum(s_i) - C(n,2)}}.

  This is TOO SIMPLE. The bracket polynomial gets interesting through
  the LOOP counting, which corresponds to CYCLE structure.

  Better analog: Weight each HP P of T by q^{{asc(P)}} where
  asc(P) is the number of ascents in P (viewed as a permutation).
  Then F(T, q) is the F-POLYNOMIAL, which we already know!

  F(T, q) = sum_P q^{{asc(P)}} where sum is over Hamiltonian paths.

  F(T, 1) = H(T), F(T, 0) = F_0(T) in {{0,1}}.
  F(T, -1) = "signed HP count" (alternating by ascent parity).

  NEW IDEA: Define J(T, q) = F(T, q) / F(T_trans, q)
  where T_trans is the transitive tournament.
  F(T_trans, q) = q^0 + q^1 + ... + q^{{n-1}}... no, T_trans has exactly 1 HP.
  So F(T_trans, q) = q^{{C(n,2)}} (all arcs ascending).
  Wait: F_0(T_trans) = 1, all others = 0. So F(T_trans, q) = 1.
  And F(T_trans^op, q) = q^{{n-1}} (one HP, all descents).
  Hmm, that's not right either.

  For the transitive tournament T_trans on [n] with i->j iff i<j:
  The unique HP is 1->2->...->n, which has n-1 ascents.
  So F(T_trans, q) = q^{{n-1}}.

  J(T, q) = F(T, q) / q^{{n-1}} = sum_k F_k q^{{k-(n-1)}}
           = sum_k F_k q^{{k-n+1}}
""")

    # Compute F(-1) for all n=5 tournaments
    n = 5
    F_neg1_vals = []
    for bits, H, adj in tours5:
        F = [0]*n
        for perm in permutations(range(n)):
            valid = True
            for i in range(n-1):
                if not adj[perm[i]][perm[i+1]]:
                    valid = False
                    break
            if valid:
                asc = sum(1 for i in range(n-1) if perm[i] < perm[i+1])
                F[asc] += 1
        F_neg1 = sum(F[k] * (-1)**k for k in range(n))
        F_neg1_vals.append((H, F_neg1, tuple(F)))

    print(f"\n  n=5: F(T, -1) statistics:")
    F_neg1_by_H = defaultdict(list)
    for H, f_neg1, F in F_neg1_vals:
        F_neg1_by_H[H].append(f_neg1)

    for H_val in sorted(F_neg1_by_H.keys()):
        vals = F_neg1_by_H[H_val]
        distinct = sorted(set(vals))
        print(f"    H={H_val:3d}: F(-1) in {distinct}")

    # Is F(-1) always +-1?
    all_F_neg1 = [f for _, f, _ in F_neg1_vals]
    print(f"\n  F(-1) values across all tournaments: {sorted(set(all_F_neg1))}")
    print(f"  Is |F(-1)| always 1? {all(abs(f) == 1 for f in all_F_neg1)}")

    # ================================================================
    # IDEA 6: ENTROPY OF THE H DISTRIBUTION
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 6: ENTROPY AND INFORMATION THEORY")
    print(f"{'='*70}")

    for n in [3, 4, 5]:
        tours = all_tournaments(n) if n < 5 else tours5
        H_counts = Counter(H for _, H, _ in tours)
        total = len(tours)
        probs = [count/total for count in H_counts.values()]
        entropy = -sum(p * log2(p) for p in probs if p > 0)
        m = n*(n-1)//2
        max_entropy = log2(len(H_counts))  # log of number of H values

        print(f"\n  n={n} (m={m}):")
        print(f"    |H-spectrum| = {len(H_counts)}")
        print(f"    Shannon entropy H(H) = {entropy:.4f} bits")
        print(f"    Max possible = log2({len(H_counts)}) = {max_entropy:.4f} bits")
        print(f"    Efficiency = {entropy/max_entropy:.4f}")
        print(f"    Entropy per arc = {entropy/m:.4f} bits/arc")
        print(f"    Fraction of hypercube entropy = {entropy/m:.4f}")

    # ================================================================
    # IDEA 7: DETERMINANT / PERMANENT DUALITY
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 7: DET vs PERM — THE SIGN STRUCTURE")
    print(f"{'='*70}")

    print(f"""
  H(T) = permanent of a certain matrix (the HP adjacency).
  But what about the DETERMINANT?

  For tournament T with adjacency matrix A (A[i][j] = 1 if i->j):
  perm(A) = sum_sigma prod A[sigma(i)][i] (counts directed cycles covering all vertices)
  det(A) = sum_sigma (-1)^sigma prod A[sigma(i)][i] (signed version)

  For the SKEW matrix B = A - A^T (B[i][j] = +1 if i->j, -1 if j->i):
  pf(B) = Pfaffian (exists when n is even)
  det(B) = pf(B)^2

  The Pfaffian counts PERFECT MATCHINGS of the complete graph
  with orientation signs. For tournaments:
  pf(B) = sum over perfect matchings M: prod_{{(i,j) in M}} B[i][j]
""")

    # Compute det(A), det(B), perm(A) for n=4 tournaments
    n = 4
    print(f"\n  n={n}: det(A), det(B=A-A^T), perm(A) for each H class:")
    det_by_H = defaultdict(list)
    for bits, H, adj in tours4:
        A = np.array(adj, dtype=float)
        B = A - A.T  # Skew-symmetric
        det_A = int(round(np.linalg.det(A)))
        det_B = int(round(np.linalg.det(B)))
        # Permanent (brute force for n=4)
        perm_A = 0
        for perm in permutations(range(n)):
            perm_A += int(np.prod([A[i][perm[i]] for i in range(n)]))
        det_by_H[H].append((det_A, det_B, perm_A))

    for H_val in sorted(det_by_H.keys()):
        vals = det_by_H[H_val]
        det_As = sorted(set(v[0] for v in vals))
        det_Bs = sorted(set(v[1] for v in vals))
        perm_As = sorted(set(v[2] for v in vals))
        print(f"    H={H_val}: det(A) in {det_As}, det(B) in {det_Bs}, perm(A) in {perm_As}")

    # n=5
    n = 5
    print(f"\n  n={n}: det(B=A-A^T) and Pfaffian structure")
    det_B_by_H = defaultdict(list)
    for bits, H, adj in tours5:
        A = np.array(adj, dtype=float)
        B = A - A.T
        det_B = round(np.linalg.det(B))
        det_B_by_H[H].append(det_B)

    for H_val in sorted(det_B_by_H.keys()):
        vals = det_B_by_H[H_val]
        distinct = sorted(set(int(v) for v in vals))
        print(f"    H={H_val:3d}: det(B) in {distinct}")

    # ================================================================
    # IDEA 8: TOURNAMENT PERSISTENCE DIAGRAM
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 8: TOPOLOGICAL DATA ANALYSIS — H FILTRATION")
    print(f"{'='*70}")

    print(f"""
  Define a filtration on the tournament hypercube by H value:
  X_h = {{T : H(T) <= h}} for h = 1, 3, 5, ...

  As h increases, more tournaments are "born."
  The persistence diagram tracks when connected components
  appear and merge.

  For n=4: H-spectrum = {{1, 3, 5}}
  X_1: 24 tournaments (transitive + isomorphic)
  X_3: 24 + 16 = 40 tournaments
  X_5: 40 + 24 = 64 = all tournaments

  Question: How many connected components does X_h have?
  (Connected via single arc flips)
""")

    # BFS to count connected components at each H threshold
    n = 4
    m = n*(n-1)//2

    for h_thresh in [1, 3, 5]:
        # Vertices in X_h
        vertices = set()
        for i, (bits, H, adj) in enumerate(tours4):
            if H <= h_thresh:
                vertices.add(bits)

        # BFS to find components
        visited = set()
        num_components = 0
        for v in vertices:
            if v not in visited:
                num_components += 1
                queue = [v]
                visited.add(v)
                while queue:
                    curr = queue.pop(0)
                    for arc in range(m):
                        nbr = curr ^ (1 << arc)
                        if nbr in vertices and nbr not in visited:
                            visited.add(nbr)
                            queue.append(nbr)

        print(f"    X_{h_thresh}: {len(vertices)} tournaments, "
              f"{num_components} connected components")

    # n=5
    n = 5
    m = n*(n-1)//2
    H_vals_sorted = sorted(set(H for _, H, _ in tours5))

    print(f"\n  n=5 filtration:")
    prev_components = 0
    for h_thresh in H_vals_sorted:
        vertices = set()
        for bits, H, adj in tours5:
            if H <= h_thresh:
                vertices.add(bits)

        visited = set()
        num_components = 0
        for v in vertices:
            if v not in visited:
                num_components += 1
                queue = [v]
                visited.add(v)
                while queue:
                    curr = queue.pop(0)
                    for arc in range(m):
                        nbr = curr ^ (1 << arc)
                        if nbr in vertices and nbr not in visited:
                            visited.add(nbr)
                            queue.append(nbr)

        delta = num_components - prev_components
        print(f"    X_{h_thresh:3d}: {len(vertices):5d} tournaments, "
              f"{num_components:3d} components (delta={delta:+d})")
        prev_components = num_components

    # ================================================================
    # IDEA 9: GOLDEN RATIO IN TOURNAMENT COUNTS
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 9: GOLDEN/TRIBONACCI RATIOS IN H-MULTIPLICITIES")
    print(f"{'='*70}")

    # Check if consecutive H-value multiplicities have golden-ratio-like ratios
    n = 5
    H_counts = Counter(H for _, H, _ in tours5)
    H_sorted = sorted(H_counts.keys())

    phi = (1 + sqrt(5)) / 2
    tau = 1.8392867552141612

    print(f"\n  n=5 H-value multiplicities and consecutive ratios:")
    for i, H in enumerate(H_sorted):
        count = H_counts[H]
        if i > 0:
            prev_count = H_counts[H_sorted[i-1]]
            ratio = count / prev_count if prev_count > 0 else 0
            print(f"    H={H:3d}: {count:5d} tournaments, "
                  f"ratio to prev = {ratio:.4f} "
                  f"(phi={phi:.4f}, tau={tau:.4f})")
        else:
            print(f"    H={H:3d}: {count:5d} tournaments")

    # Check: count / 8 (since GCD of n=5 counts is 8)
    gcd_all = H_counts[H_sorted[0]]
    for H in H_sorted[1:]:
        gcd_all = gcd(gcd_all, H_counts[H])
    print(f"\n  GCD of all multiplicities = {gcd_all}")
    print(f"  Normalized counts (divided by {gcd_all}):")
    normalized = []
    for H in H_sorted:
        norm = H_counts[H] // gcd_all
        normalized.append(norm)
        print(f"    H={H:3d}: {norm}")
    print(f"  Normalized sequence: {normalized}")
    print(f"  Sum = {sum(normalized)}")

    # Check if this sequence appears in OEIS or has a pattern
    # [15, 15, 30, 30, 15, 15, 8]
    # Palindromic? 15,15,30,30,15,15,8 — NO (8 breaks it)
    # But 15+8 = 23, 15+15 = 30...
    # Actually these should be symmetric under T -> T^op
    # H(T) = H(T^op), so the distribution should be symmetric
    # But the H VALUES are already symmetric: {1,3,5,9,11,13,15}
    # 1 <-> 15, 3 <-> 13, 5 <-> 11, 9 <-> 9
    # Counts: 120,120,240,240,120,120,64
    # So 120=120, 120=120, 240=240, but 64 is the self-complementary class.

    # ================================================================
    # IDEA 10: MODULAR ARITHMETIC PATTERNS
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 10: H mod p FOR VARIOUS PRIMES")
    print(f"{'='*70}")

    n = 5
    for p in [3, 5, 7, 8, 12]:
        residue_counts = Counter()
        for bits, H, adj in tours5:
            residue_counts[H % p] += 1
        total = len(tours5)
        print(f"\n  n=5, H mod {p}:")
        for r in sorted(residue_counts.keys()):
            count = residue_counts[r]
            pct = count / total * 100
            print(f"    H equiv {r:2d} mod {p}: {count:5d} ({pct:.1f}%)")

    # ================================================================
    # IDEA 11: THE IHARA ZETA OF Omega(T)
    # ================================================================
    print(f"\n{'='*70}")
    print("IDEA 11: COMPLEXITY / SPANNING TREE COUNT OF Omega(T)")
    print(f"{'='*70}")

    print(f"""
  The conflict graph Omega(T) is an undirected graph.
  Its COMPLEXITY (number of spanning trees) is given by
  the Matrix Tree Theorem: kappa = det(L_0) / n where
  L is the Laplacian and L_0 is any cofactor.

  For small Omega graphs:
  - Empty graph (alpha_1 = 0): kappa = 0
  - Single vertex (alpha_1 = 1): kappa = 1
  - Two disjoint vertices: kappa = 0
  - K_2 (edge): kappa = 1
  - K_3 (triangle): kappa = 3
  - Path P_3: kappa = 1

  QUESTION: Does the complexity of Omega(T) correlate with H(T)?
  H = I(Omega, 2) counts independent sets.
  kappa counts spanning trees.
  Are these related?
""")

    # For n=5, compute number of 3-cycles (= alpha_1 = vertices of Omega)
    # and check connectivity of Omega
    n = 5
    alpha_H = defaultdict(list)
    for bits, H, adj in tours5:
        t3 = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if adj[i][j] and adj[j][k] and adj[k][i]:
                        t3 += 1
                    elif adj[i][k] and adj[k][j] and adj[j][i]:
                        t3 += 1
        alpha_H[H].append(t3)

    print(f"\n  n=5: 3-cycle count (alpha_1 proxy) by H value:")
    for H_val in sorted(alpha_H.keys()):
        vals = alpha_H[H_val]
        distinct = sorted(set(vals))
        print(f"    H={H_val:3d}: t3 in {distinct}")

    print(f"\n{'='*70}")
    print("DONE — WILD IDEAS SESSION")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
