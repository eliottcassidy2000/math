"""
knot_tournament_bridge_v3.py -- kind-pasteur-2026-03-14-S69
FIXED contraction convention (THM-082) + deeper analysis of key findings.

v1/v2 findings:
- det(M) = 0 iff H is "too small" (threshold around 2n-1)
- |writhe| negatively correlates with H (regularity = high H)
- Alexander poly det(t*M - M^T) has palindromic structure when nontrivial
- F(-1) = signed ham path count is a finer tournament invariant

v3 focus:
1. CORRECT DC / skein relation verification
2. Deep writhe-H correlation and "knot crossing" analogy
3. det(t*A - A^T) for the ADJACENCY matrix (not transfer matrix)
4. "Reidemeister invariance" under local moves
5. Kauffman bracket state sum properties
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def tournament_from_bits(n, bits):
    T = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[j][i] = 1
            else:
                T[i][j] = 1
            idx += 1
    return T

def all_tournaments(n):
    m = n * (n - 1) // 2
    for bits in range(2**m):
        yield bits, tournament_from_bits(n, bits)

def ham_paths_digraph(D):
    """Count Hamiltonian paths in a general digraph."""
    n = len(D)
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if D[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def H_count(T):
    return ham_paths_digraph(T)

def score_seq(T):
    return tuple(sorted([sum(T[i]) for i in range(len(T))]))

def fwd_polynomial(T):
    n = len(T)
    poly = Counter()
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if T[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            fwd = sum(1 for i in range(n-1) if perm[i] < perm[i+1])
            poly[fwd] += 1
    return dict(sorted(poly.items()))

# ============================================================
# CORRECTED contraction (THM-082)
# ============================================================
def contract_arc(T, u, v):
    """
    Contract arc e = u->v in digraph T.
    THM-082 convention:
      x -> w iff x -> u  (IN from TAIL u)
      w -> x iff v -> x  (OUT from HEAD v)
    """
    n = len(T)
    remaining = [x for x in range(n) if x != u and x != v]
    new_n = n - 1
    D = np.zeros((new_n, new_n), dtype=int)

    # w = index 0
    idx = {remaining[i]: i+1 for i in range(len(remaining))}

    # Among remaining vertices: unchanged
    for x in remaining:
        for y in remaining:
            if x != y:
                D[idx[x]][idx[y]] = T[x][y]

    # w -> x iff v -> x (OUT from head)
    for x in remaining:
        D[0][idx[x]] = T[v][x]

    # x -> w iff x -> u (IN from tail)
    for x in remaining:
        D[idx[x]][0] = T[x][u]

    return D

def delete_arc(T, u, v):
    D = T.copy()
    D[u][v] = 0
    return D

# ============================================================
# PART 1: DC Verification (corrected)
# ============================================================
def part1():
    print("=" * 70)
    print("PART 1: DELETION-CONTRACTION VERIFICATION (CORRECTED)")
    print("  THM-082: H(D) = H(D\\e) + H(D/e)")
    print("  Convention: x->w iff x->u; w->x iff v->x")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        total = 0
        failures = 0
        for bits, T in all_tournaments(n):
            for u in range(n):
                for v in range(n):
                    if u != v and T[u][v] == 1:
                        H_T = H_count(T)
                        H_del = ham_paths_digraph(delete_arc(T, u, v))
                        H_con = ham_paths_digraph(contract_arc(T, u, v))
                        total += 1
                        if H_T != H_del + H_con:
                            failures += 1
                            if failures <= 3:
                                print(f"  FAIL: e=({u},{v}), H={H_T}, del={H_del}, con={H_con}")
        print(f"  {total} checks, {failures} failures {'<-- ALL PASS' if failures == 0 else ''}")

# ============================================================
# PART 2: The SKEIN RELATION for arc flips
# ============================================================
def part2():
    print("\n" + "=" * 70)
    print("PART 2: ARC-FLIP SKEIN RELATION (corrected DC)")
    print("  H(T) - H(T') = H(T/e) - H(T'/e')")
    print("  where T' flips arc e")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        total = 0
        failures = 0
        delta_distribution = Counter()

        for bits, T in all_tournaments(n):
            for u in range(n):
                for v in range(u+1, n):
                    total += 1
                    # T has either u->v or v->u
                    if T[u][v] == 1:
                        tail, head = u, v
                    else:
                        tail, head = v, u

                    # T_flip: reverse the arc
                    T_flip = T.copy()
                    T_flip[tail][head] = 0
                    T_flip[head][tail] = 1

                    H_T = H_count(T)
                    H_flip = H_count(T_flip)

                    # DC: H(T) = H(T\e) + H(T/e)
                    # Since T\e = T_flip\e' (same digraph without the u-v arc):
                    # H(T) - H(T') = H(T/e) - H(T'/e')
                    H_con = ham_paths_digraph(contract_arc(T, tail, head))
                    H_con_flip = ham_paths_digraph(contract_arc(T_flip, head, tail))

                    delta_H = H_T - H_flip
                    delta_con = H_con - H_con_flip

                    if delta_H != delta_con:
                        failures += 1
                        if failures <= 3:
                            print(f"  FAIL: arc ({tail},{head}), dH={delta_H}, dCon={delta_con}")

                    delta_distribution[delta_H] += 1

        print(f"  {total} flips, {failures} failures {'<-- ALL PASS' if failures == 0 else ''}")
        print(f"  Delta_H distribution: {dict(sorted(delta_distribution.items()))}")

# ============================================================
# PART 3: Writhe-H deep correlation
# ============================================================
def part3():
    print("\n" + "=" * 70)
    print("PART 3: WRITHE vs H -- DEEP CORRELATION ANALYSIS")
    print("  Writhe = #{ascending arcs} - #{descending arcs}")
    print("  Analogy: knot writhe measures 'twistedness'")
    print("=" * 70)

    for n in [4, 5, 6]:
        print(f"\n--- n = {n} ---")
        m = n * (n - 1) // 2
        data = []
        count = 0

        for bits, T in all_tournaments(n):
            count += 1
            if n >= 6 and count > 10000:
                break
            H = H_count(T)
            asc = sum(T[i][j] for i in range(n) for j in range(i+1, n))
            writhe = 2 * asc - m
            data.append((H, writhe))

        H_vals = np.array([d[0] for d in data], dtype=float)
        W_vals = np.array([d[1] for d in data], dtype=float)
        absW = np.abs(W_vals)
        W2 = W_vals**2

        print(f"  Correlation(H, writhe) = {np.corrcoef(H_vals, W_vals)[0,1]:.6f}")
        print(f"  Correlation(H, |writhe|) = {np.corrcoef(H_vals, absW)[0,1]:.6f}")
        print(f"  Correlation(H, writhe^2) = {np.corrcoef(H_vals, W2)[0,1]:.6f}")

        # Linear regression: H ~ a + b * writhe^2
        A_mat = np.column_stack([np.ones_like(W2), W2])
        result = np.linalg.lstsq(A_mat, H_vals, rcond=None)
        a, b = result[0]
        print(f"  Linear fit H ~ {a:.4f} + {b:.6f} * writhe^2")
        residuals = H_vals - (a + b * W2)
        R2 = 1 - np.var(residuals) / np.var(H_vals)
        print(f"  R^2 = {R2:.6f}")

        # KEY INSIGHT: H(T) + H(T^op) = const for same score? No...
        # But: writhe(T) + writhe(T^op) = 0 (reversing all arcs negates writhe)
        # And H(T^op) = H(T) (reversing arcs reverses all Ham paths, same count)
        # So the writhe-H relationship must be symmetric in writhe
        # => H depends on writhe^2, not writhe (confirmed by correlation)

        # At writhe = 0: H takes what values?
        zero_writhe = [d[0] for d in data if d[1] == 0]
        max_writhe = [d[0] for d in data if abs(d[1]) == m]
        if zero_writhe:
            print(f"  H at writhe=0: mean={np.mean(zero_writhe):.2f}, range=[{min(zero_writhe)},{max(zero_writhe)}]")
        if max_writhe:
            print(f"  H at |writhe|={m}: always H={set(max_writhe)}")

# ============================================================
# PART 4: det(t*A - A^T) for ADJACENCY matrix
# ============================================================
def part4():
    print("\n" + "=" * 70)
    print("PART 4: 'ALEXANDER POLYNOMIAL' OF TOURNAMENT ADJACENCY")
    print("  P_T(t) = det(t*A - A^T) where A is adjacency matrix")
    print("  For skew-symmetric S = A - A^T, det(S) = Pfaff(S)^2 (even n)")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        print(f"\n--- n = {n} ---")
        poly_by_H = defaultdict(list)
        count = 0

        for bits, T in all_tournaments(n):
            count += 1
            if n == 6 and count > 3000:
                break

            H = H_count(T)
            A = T.astype(float)

            # Evaluate det(t*A - A^T) at several points
            vals = {}
            for t in [-1, 0, 1, 2, 3]:
                mat = t * A - A.T
                vals[t] = round(np.linalg.det(mat))

            poly_by_H[H].append(vals)

        print(f"  H -> det(t*A - A^T) spectra:")
        for H in sorted(poly_by_H.keys()):
            items = poly_by_H[H]
            spectra = set(tuple(sorted(s.items())) for s in items)
            print(f"  H={H:3d}: {len(spectra)} distinct spectra (of {len(items)} tournaments)")

            # Show at most 3
            for spec in sorted(list(spectra))[:3]:
                d = dict(spec)
                print(f"    t -> det: {d}")

        # KEY: A - A^T has entries +1, -1, 0 (diagonal)
        # It's a skew-symmetric {-1,0,1} matrix
        # For tournaments: A[i][j] + A[j][i] = 1, so (A-A^T)[i][j] = 2*A[i][j] - 1 for i!=j
        # = a skew-symmetric sign matrix!
        # det(A-A^T) = Pfaffian^2 for even n; = 0 for odd n

        print(f"\n  NOTE: det(A-A^T) for n={n}:")
        all_det1 = set()
        count2 = 0
        for bits, T in all_tournaments(n):
            count2 += 1
            if count2 > 3000:
                break
            S = T.astype(float) - T.T.astype(float)
            d = round(np.linalg.det(S))
            all_det1.add(d)
        print(f"    Distinct values: {sorted(all_det1)}")
        if n % 2 == 1:
            print(f"    (n odd => always 0)")
        else:
            print(f"    (n even => Pfaffian^2, always >=0)")

# ============================================================
# PART 5: "Reidemeister moves" on tournaments
# ============================================================
def part5():
    print("\n" + "=" * 70)
    print("PART 5: 'REIDEMEISTER MOVES' ON TOURNAMENTS")
    print("  Knot: R1=loop, R2=bigon, R3=triangle")
    print("  Tournament: source/sink removal, 3-cycle reversal, etc.")
    print("=" * 70)

    # "R3 move" for tournaments: given a 3-cycle (a,b,c) with a->b->c->a,
    # reverse ALL three arcs to get a->c->b->a (still a 3-cycle, opposite orientation)
    # How does H change?

    print("\n  3-CYCLE REVERSAL: reverse all 3 arcs of a directed 3-cycle")
    for n in [4, 5, 6]:
        print(f"\n  --- n = {n} ---")
        delta_H_dist = Counter()
        total = 0
        count = 0

        for bits, T in all_tournaments(n):
            count += 1
            if n == 6 and count > 3000:
                break

            H_T = H_count(T)

            # Find all 3-cycles
            for a, b, c in combinations(range(n), 3):
                for perm in permutations([a, b, c]):
                    x, y, z = perm
                    if T[x][y] == 1 and T[y][z] == 1 and T[z][x] == 1:
                        # Found 3-cycle x->y->z->x
                        # Reverse: x->z->y->x
                        T_rev = T.copy()
                        T_rev[x][y] = 0; T_rev[y][x] = 1
                        T_rev[y][z] = 0; T_rev[z][y] = 1
                        T_rev[z][x] = 0; T_rev[x][z] = 1

                        H_rev = H_count(T_rev)
                        delta = H_T - H_rev
                        delta_H_dist[delta] += 1
                        total += 1
                        break  # one cycle per triple is enough to avoid overcounting orientations

        print(f"    {total} 3-cycle reversals")
        print(f"    Delta_H distribution: {dict(sorted(delta_H_dist.items()))}")

        # KEY: Is delta_H always 0 under 3-cycle reversal?
        # If yes, 3-cycle reversal is a "Reidemeister move" for H!
        if all(d == 0 for d in delta_H_dist.keys()):
            print(f"    *** H IS INVARIANT UNDER 3-CYCLE REVERSAL! ***")
        else:
            print(f"    H is NOT invariant under 3-cycle reversal")
            print(f"    Delta always even: {all(d % 2 == 0 for d in delta_H_dist.keys())}")

# ============================================================
# PART 6: F(T, -1) signed count — "signature" invariant
# ============================================================
def part6():
    print("\n" + "=" * 70)
    print("PART 6: SIGNED HAMILTONIAN PATH COUNT F(T, -1)")
    print("  F(T, -1) = sum_{paths P} (-1)^{fwd(P)}")
    print("  Analogous to: signature of knot, or Jones at t=-1")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        print(f"\n--- n = {n} ---")
        sig_dist = Counter()
        sig_by_H = defaultdict(list)
        count = 0

        for bits, T in all_tournaments(n):
            count += 1
            if n == 6 and count > 5000:
                break

            H = H_count(T)
            fp = fwd_polynomial(T)
            sig = sum((-1)**k * c for k, c in fp.items())

            sig_dist[sig] += 1
            sig_by_H[H].append(sig)

        print(f"  F(-1) distribution:")
        for val in sorted(sig_dist.keys()):
            print(f"    F(-1)={val:4d}: {sig_dist[val]:5d} tournaments")

        # Mod structure
        print(f"\n  F(-1) mod 2: always odd? {all(v % 2 != 0 for v in sig_dist.keys())}")
        print(f"  F(-1) mod 4 distribution: {Counter(v % 4 for v in sig_dist.keys() for _ in range(sig_dist[v]))}")

        # Does F(-1) determine H?
        print(f"\n  F(-1) vs H:")
        for H in sorted(sig_by_H.keys()):
            vals = sorted(set(sig_by_H[H]))
            print(f"    H={H:3d}: F(-1) in {vals}")

# ============================================================
# PART 7: Kauffman state sum analysis
# ============================================================
def part7():
    print("\n" + "=" * 70)
    print("PART 7: KAUFFMAN STATE SUM K_n(A) = sum_T A^{H(T)}")
    print("  Generating function over ALL 2^C(n,2) tournaments")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        H_dist = Counter()
        for bits, T in all_tournaments(n):
            H = H_count(T)
            H_dist[H] += 1

        # Print as polynomial
        terms = []
        for H in sorted(H_dist.keys()):
            terms.append(f"{H_dist[H]}*A^{H}")
        print(f"  K_{n}(A) = " + " + ".join(terms))

        # Key evaluations
        total = sum(H_dist.values())
        mean_H = sum(H * c for H, c in H_dist.items()) / total
        var_H = sum(H**2 * c for H, c in H_dist.items()) / total - mean_H**2

        print(f"  Mean H = {mean_H:.4f}")
        print(f"  n!/2^(n-1) = {math.factorial(n) / 2**(n-1):.4f}")
        print(f"  Var H = {var_H:.4f}")

        # Check: is mean_H = n!/2^(n-1)?
        expected_mean = math.factorial(n) / 2**(n-1)
        print(f"  Mean H = n!/2^(n-1)? {abs(mean_H - expected_mean) < 0.001}")

        # Logarithmic derivative structure
        print(f"  K'_n(1)/K_n(1) = mean log_H = {sum(math.log(H) * c for H, c in H_dist.items()) / total:.4f}")

        # Ratio of consecutive H counts
        H_list = sorted(H_dist.keys())
        print(f"  Count ratios:")
        for i in range(len(H_list) - 1):
            H1, H2 = H_list[i], H_list[i+1]
            ratio = H_dist[H2] / H_dist[H1]
            print(f"    K[H={H2}]/K[H={H1}] = {ratio:.4f}")

# ============================================================
# PART 8: Transfer matrix as "Seifert matrix" — linking numbers
# ============================================================
def part8():
    print("\n" + "=" * 70)
    print("PART 8: TRANSFER MATRIX AS 'SEIFERT MATRIX'")
    print("  In knot theory: det(V - V^T) = determinant of knot")
    print("  For tournaments: det(M - M^T) and det(M + M^T)")
    print("  where M[a,b] = #{Ham paths from a to b}")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        results = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = H_count(T)
            M = np.zeros((n, n), dtype=int)
            for perm in permutations(range(n)):
                valid = True
                for i in range(n-1):
                    if T[perm[i]][perm[i+1]] != 1:
                        valid = False
                        break
                if valid:
                    M[perm[0]][perm[-1]] += 1

            # Symmetric and antisymmetric parts
            S = M + M.T  # "symmetrized Seifert"
            A = M - M.T  # antisymmetric part

            det_M = round(np.linalg.det(M.astype(float)))
            det_S = round(np.linalg.det(S.astype(float)))
            det_A = round(np.linalg.det(A.astype(float)))
            tr_M = int(np.trace(M))

            # Signature (number of positive - negative eigenvalues of S)
            eigs_S = np.linalg.eigvals(S.astype(float))
            signature = sum(1 for e in eigs_S if e.real > 0.001) - sum(1 for e in eigs_S if e.real < -0.001)

            results[H].append({
                'det_M': det_M, 'det_S': det_S, 'det_A': det_A,
                'tr': tr_M, 'sig': signature,
                'M': M.tolist(),
            })

        for H in sorted(results.keys()):
            items = results[H]
            dets = sorted(set(d['det_M'] for d in items))
            det_S_vals = sorted(set(d['det_S'] for d in items))
            sigs = sorted(set(d['sig'] for d in items))
            print(f"  H={H:3d}: det(M) in {dets}, det(M+M^T) in {det_S_vals}, signature in {sigs}")

        # Show explicit M for interesting cases
        if n == 5:
            for H in sorted(results.keys()):
                if H in [9, 15]:  # Show regular and maximizer
                    ex = results[H][0]
                    print(f"\n  Example M for H={H}:")
                    for row in ex['M']:
                        print(f"    {row}")

# ============================================================
# PART 9: The DEEP ANALOGY — summary
# ============================================================
def part9():
    print("\n" + "=" * 70)
    print("PART 9: KNOT-TOURNAMENT DICTIONARY SUMMARY")
    print("=" * 70)

    dictionary = [
        ("Knot / Link", "Tournament T"),
        ("Crossings", "Arcs (directed edges)"),
        ("Writhe (sum of signs)", "|writhe| = |#asc - #desc arcs|"),
        ("Reidemeister moves", "Arc flips (single), 3-cycle reversals?"),
        ("Skein relation V(K+)-t*V(K-)=(..)*V(K0)", "H(T)-H(T')=H(T/e)-H(T'/e') [THM-082]"),
        ("Jones polynomial V(K,t)", "Forward polynomial F(T,x)"),
        ("Jones at t=-1 = det(knot)", "F(T,-1) = signed Ham path count"),
        ("Seifert matrix V", "Transfer matrix M[a,b]"),
        ("Alexander poly det(t*V-V^T)", "det(t*M-M^T)"),
        ("det(V-V^T) = knot determinant", "det(M-M^T)"),
        ("Kauffman bracket <K>", "State sum K_n(A) = sum_T A^{H(T)}"),
        ("Khovanov homology", "GLMY path homology"),
        ("Euler char of Khovanov = Jones", "Euler char of GLMY = ??"),
        ("Unknot (trivial knot)", "Transitive tournament (H=1)"),
        ("Alternating knots", "Regular tournaments?"),
    ]

    for knot_side, tournament_side in dictionary:
        print(f"  {knot_side:45s} <-> {tournament_side}")

# ============================================================
# PART 10: What does path homology Euler characteristic equal?
# ============================================================
def compute_boundary_matrix(T, p):
    """Compute GLMY boundary map d_p for tournament T."""
    n = len(T)

    # p-paths = allowed sequences (v0,...,vp) where T[vi][vi+1]=1 for all i
    # and all vi distinct (elementary paths)
    def get_paths(k):
        """Get all elementary k-paths (sequences of k+1 distinct vertices with arcs)."""
        if k == 0:
            return [tuple([v]) for v in range(n)]
        paths = []
        for perm in permutations(range(n), k+1):
            valid = True
            for i in range(k):
                if T[perm[i]][perm[i+1]] != 1:
                    valid = False
                    break
            if valid:
                paths.append(perm)
        return paths

    source_paths = get_paths(p)
    target_paths = get_paths(p-1)

    if not source_paths or not target_paths:
        return np.zeros((0, 0), dtype=int), source_paths, target_paths

    # Index maps
    target_idx = {path: i for i, path in enumerate(target_paths)}

    # Boundary: d_p(v0,...,vp) = sum_{i=0}^{p} (-1)^i (v0,...,v_{i-1},v_{i+1},...,vp)
    # Only include terms that are themselves valid elementary paths
    matrix = np.zeros((len(target_paths), len(source_paths)), dtype=int)

    for j, path in enumerate(source_paths):
        for i in range(p + 1):
            face = path[:i] + path[i+1:]
            if face in target_idx:
                sign = (-1) ** i
                matrix[target_idx[face]][j] += sign

    return matrix, source_paths, target_paths

def betti_numbers(T):
    """Compute GLMY path homology Betti numbers."""
    n = len(T)
    bettis = []

    for p in range(n):
        if p == 0:
            # H_0: kernel of d_0 (but d_0 doesn't exist; H_0 = connected components)
            # For tournaments, always connected, so beta_0 = 1
            bettis.append(1)
            continue

        # Get d_p and d_{p+1}
        d_p, _, _ = compute_boundary_matrix(T, p)
        d_p1, _, _ = compute_boundary_matrix(T, p + 1)

        if d_p.size == 0:
            # No p-paths
            bettis.append(0)
            continue

        # kernel of d_p
        if d_p.shape[1] == 0:
            ker_dim = 0
        else:
            rank_dp = np.linalg.matrix_rank(d_p.astype(float))
            ker_dim = d_p.shape[1] - rank_dp

        # image of d_{p+1}
        if d_p1.size == 0 or d_p1.shape[1] == 0:
            im_dim = 0
        else:
            im_dim = np.linalg.matrix_rank(d_p1.astype(float))

        bettis.append(ker_dim - im_dim)

    return bettis

def part10():
    print("\n" + "=" * 70)
    print("PART 10: PATH HOMOLOGY EULER CHARACTERISTIC")
    print("  chi = sum (-1)^k * beta_k")
    print("  Khovanov analogy: chi(Kh) = Jones polynomial")
    print("  What does chi equal for tournament path homology?")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        print(f"\n--- n = {n} ---")
        chi_by_H = defaultdict(list)
        betti_by_H = defaultdict(list)
        count = 0

        for bits, T in all_tournaments(n):
            count += 1
            if n == 6 and count > 1000:
                break

            H = H_count(T)
            bettis = betti_numbers(T)
            chi = sum((-1)**k * b for k, b in enumerate(bettis))

            chi_by_H[H].append(chi)
            betti_by_H[H].append(tuple(bettis))

        print(f"  H -> chi, Betti numbers:")
        for H in sorted(chi_by_H.keys()):
            chis = sorted(set(chi_by_H[H]))
            bettis_set = set(betti_by_H[H])
            print(f"    H={H:3d}: chi in {chis}, Betti patterns: {len(bettis_set)}")
            for b in sorted(bettis_set)[:3]:
                chi_b = sum((-1)**k * v for k, v in enumerate(b))
                print(f"      beta = {list(b)}, chi = {chi_b}")

        # KEY: Is chi constant? Does it equal anything nice?
        all_chis = []
        for chis in chi_by_H.values():
            all_chis.extend(chis)
        print(f"\n  All chi values: {sorted(set(all_chis))}")
        print(f"  Chi constant: {len(set(all_chis)) == 1}")
        if len(set(all_chis)) == 1:
            print(f"  *** chi = {all_chis[0]} for ALL tournaments at n={n} ***")

# ============================================================
def main():
    print("=" * 70)
    print("KNOT-TOURNAMENT BRIDGE v3 -- kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    part7()
    part8()
    part9()
    part10()

    print("\n" + "=" * 70)
    print("EXPLORATION COMPLETE")
    print("=" * 70)

if __name__ == '__main__':
    main()
