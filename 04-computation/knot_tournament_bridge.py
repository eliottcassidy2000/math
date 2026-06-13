"""
knot_tournament_bridge.py — kind-pasteur-2026-03-14-S69
Exploring deep connections between knot theory and tournament theory.

Key analogies:
1. THM-083 (F_T = F_{T\e} + (x-1)*F(T/e)) IS a skein relation
2. Transfer matrix M ↔ Seifert matrix V (Alexander polynomial)
3. Tournament → chord diagram → Vassiliev invariants
4. Writhe ↔ score sequence
5. Yang-Baxter equation ↔ transfer matrix symmetry
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
from functools import lru_cache
import sys

# ============================================================
# Part 1: Transfer Matrix and "Tournament Alexander Polynomial"
# ============================================================
# The Seifert matrix V of a knot gives Alexander poly: det(t*V - V^T)
# Our transfer matrix M[a,b] = #{Ham paths from a to b}
# Define: Alex_T(t) = det(t*M - M^T) / appropriate normalization
# For symmetric M (VT tournaments): det((t-1)*M) = (t-1)^n * det(M)
# For non-symmetric M: genuinely interesting!

def tournament_from_bits(n, bits):
    """Create tournament adjacency matrix from bit encoding."""
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
    """Generate all tournaments on n vertices."""
    m = n * (n - 1) // 2
    for bits in range(2**m):
        yield bits, tournament_from_bits(n, bits)

def ham_paths(T):
    """Find all Hamiltonian paths, return as list of permutations."""
    n = len(T)
    paths = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if T[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            paths.append(perm)
    return paths

def transfer_matrix(T):
    """Compute transfer matrix M[a,b] = #{Ham paths from a to b}."""
    n = len(T)
    paths = ham_paths(T)
    M = np.zeros((n, n), dtype=int)
    for p in paths:
        M[p[0]][p[-1]] += 1
    return M

def tournament_alexander(T):
    """
    Compute "Tournament Alexander Polynomial" coefficients.
    Alex_T(t) = det(t*M - M^T)
    Return as polynomial in t (list of coefficients).
    """
    n = len(T)
    M = transfer_matrix(T)

    # We compute det(t*M - M^T) symbolically
    # det(t*M - M^T) = sum over k: c_k * t^k
    # where c_k comes from expanding the determinant

    # For small n, use numpy with specific t values and interpolate
    from numpy.polynomial import polynomial as P

    # Evaluate at n+1 points to determine degree-n polynomial
    points = list(range(n + 1))
    values = []
    for t_val in points:
        mat = t_val * M - M.T
        det_val = round(np.linalg.det(mat.astype(float)))
        values.append(det_val)

    # Lagrange interpolation to get polynomial coefficients
    coeffs = np.zeros(n + 1)
    for i, (xi, yi) in enumerate(zip(points, values)):
        # Lagrange basis polynomial
        basis = np.array([1.0])
        for j, xj in enumerate(points):
            if j != i:
                basis = np.convolve(basis, [1.0, -xj]) / (xi - xj)
        # Pad to length n+1
        padded = np.zeros(n + 1)
        padded[:len(basis)] = basis
        coeffs += yi * padded

    return [round(c) for c in coeffs]

def antisymmetric_part(M):
    """Extract antisymmetric part: A = (M - M^T)/2"""
    return (M - M.T) / 2.0

# ============================================================
# Part 2: Tournament Skein Polynomial
# ============================================================
# Define S(T, q) recursively via deletion-contraction:
# S(T, q) = S(T\e, q) + q * S(T/e, q)
# where e is a specific arc, T\e deletes arc, T/e contracts
# Check: is this independent of arc choice?

def contract_arc(T, u, v):
    """
    Contract arc u→v: merge u,v into single vertex w.
    w inherits IN-arcs from v (head) and OUT-arcs from u (tail).
    For conflicts (x→u and x→v, or u→x and v→x), use tail (u) for OUT, head (v) for IN.
    """
    n = len(T)
    # Vertex w replaces both u and v
    # Remaining vertices: all except u and v, plus w
    remaining = [x for x in range(n) if x != u and x != v]
    new_n = n - 1
    T_new = np.zeros((new_n, new_n), dtype=int)

    # Map: w = 0, remaining vertices = 1, 2, ...
    idx_map = {remaining[i]: i+1 for i in range(len(remaining))}

    # Arcs among remaining vertices (unchanged)
    for x in remaining:
        for y in remaining:
            if x != y:
                T_new[idx_map[x]][idx_map[y]] = T[x][y]

    # Arcs involving w (= merged u,v)
    # Convention: w inherits IN from head (v), OUT from tail (u)
    for x in remaining:
        # w → x: u → x (tail's outgoing)
        T_new[0][idx_map[x]] = T[u][x]
        # x → w: x → v (head's incoming)
        T_new[idx_map[x]][0] = T[x][v]

    return T_new

def delete_arc(T, u, v):
    """Delete arc u→v from tournament T. Result is a digraph, not a tournament."""
    T_new = T.copy()
    T_new[u][v] = 0
    return T_new

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

def skein_check(T):
    """
    Check: F_T(x) = F_{T\e}(x) + (x-1)*F(T/e, x) at x=2 (i.e., H values)
    For each arc e = (u,v) of T.
    """
    n = len(T)
    H_T = len(ham_paths(T))

    results = []
    for u in range(n):
        for v in range(n):
            if u != v and T[u][v] == 1:
                D_del = delete_arc(T, u, v)
                T_con = contract_arc(T, u, v)

                H_del = ham_paths_digraph(D_del)
                H_con = len(ham_paths(T_con)) if T_con.shape[0] > 1 else 1

                # Check: H_T = H_del + H_con (x=2 means coefficient is 1)
                lhs = H_T
                rhs = H_del + H_con
                results.append({
                    'arc': (u, v),
                    'H_T': H_T,
                    'H_del': H_del,
                    'H_con': H_con,
                    'match': lhs == rhs
                })
    return results

# ============================================================
# Part 3: Chord Diagram / Vassiliev Invariant Analogy
# ============================================================
# A tournament = complete signed chord diagram on circle
# Place vertices on circle, each pair connected by signed chord
# Vassiliev invariants = sums over sub-diagrams with specific structure

def arrow_diagram_invariant(T, k):
    """
    Compute k-th order "tournament Vassiliev invariant":
    v_k(T) = sum over k-subsets of arcs: product of local weights
    where the weight of an arc (i→j) is +1 if i<j (ascending), -1 if i>j (descending)
    """
    n = len(T)
    arcs = [(i, j) for i in range(n) for j in range(n) if i != j and T[i][j] == 1]

    total = 0
    for combo in combinations(arcs, k):
        weight = 1
        for (i, j) in combo:
            weight *= (1 if i < j else -1)
        total += weight
    return total

def chord_intersection_number(T):
    """
    For tournament T, place vertices 1..n on circle.
    Two chords (i,j) and (k,l) "cross" if exactly one of k,l is between i,j on circle.
    Count the total crossing number of the tournament chord diagram.
    This is analogous to the writhe of a knot diagram.
    """
    n = len(T)
    crossings = 0
    arcs = [(i, j) for i in range(n) for j in range(n) if i != j and T[i][j] == 1]

    for a1_idx in range(len(arcs)):
        for a2_idx in range(a1_idx + 1, len(arcs)):
            i, j = arcs[a1_idx]
            k, l = arcs[a2_idx]
            if len({i, j, k, l}) < 4:
                continue  # share endpoint, skip
            # Check if chords (i,j) and (k,l) cross on circle
            # Arrange on circle as 0..n-1
            # Two chords cross iff one separates the other's endpoints
            def between(a, b, c, n):
                """Is c between a and b going clockwise on circle of size n?"""
                if a < b:
                    return a < c < b
                else:
                    return c > a or c < b

            cross = (between(i, j, k, n) != between(i, j, l, n))
            if cross:
                # Sign = product of arc signs
                sign1 = 1 if T[i][j] else -1
                sign2 = 1 if T[k][l] else -1
                crossings += sign1 * sign2

    return crossings

# ============================================================
# Part 4: Yang-Baxter Check
# ============================================================
# The Yang-Baxter equation: R12 R13 R23 = R23 R13 R12
# For tournaments, define R(i,j) based on the arc between i and j
# Check if the transfer matrix satisfies YB-like relations

def yang_baxter_check(T):
    """
    For a tournament T on n>=3 vertices, check if the 'R-matrix'
    R_{ij} = T[i][j] - T[j][i] (= ±1 for tournaments)
    satisfies any Yang-Baxter-like identity for triples (i,j,k).
    """
    n = len(T)
    results = []
    for i, j, k in combinations(range(n), 3):
        r_ij = T[i][j] - T[j][i]
        r_ik = T[i][k] - T[k][i]
        r_jk = T[j][k] - T[k][j]

        # Classic YB would be: R12*R13*R23 = R23*R13*R12
        # For scalars this is trivially multiplication commutativity
        # Instead check: r_ij + r_jk + r_ki = ±1 (triangle orientation)
        triangle_sum = r_ij + r_jk - r_ik  # = +3 or -1 (3-cycle) or +1 or -3 (transitive)

        # Actually: r_ij*r_jk*r_ki = -1 (3-cycle) or +1 (transitive triple)
        triple_product = r_ij * r_jk * (-r_ik)  # note sign for cyclic ordering

        results.append({
            'triple': (i, j, k),
            'arcs': (r_ij, r_ik, r_jk),
            'sum': triangle_sum,
            'product': triple_product,
            'is_3cycle': abs(triangle_sum) == 3
        })
    return results

# ============================================================
# Part 5: "Writhe" and Score — Tournament Linking Number
# ============================================================
def tournament_writhe(T):
    """
    Writhe = sum of arc signs, where sign(i→j) = +1 if i<j, -1 if i>j.
    This equals #{ascending arcs} - #{descending arcs}.
    For tournaments: total arcs = C(n,2), ascending = sum T[i][j] for i<j.
    """
    n = len(T)
    ascending = sum(T[i][j] for i in range(n) for j in range(i+1, n))
    descending = n*(n-1)//2 - ascending
    return ascending - descending

def fwd_polynomial(T):
    """
    Compute F(T, x) = sum over Ham paths P: x^{fwd(P)}
    where fwd(P) = #{i : P[i] < P[i+1]} (forward/ascending steps).
    Return as dict: degree -> coefficient.
    """
    n = len(T)
    paths = ham_paths(T)
    poly = Counter()
    for p in paths:
        fwd = sum(1 for i in range(n-1) if p[i] < p[i+1])
        poly[fwd] += 1
    return dict(sorted(poly.items()))

# ============================================================
# Part 6: Tournament → Braid → Knot attempt
# ============================================================
def tournament_to_braid_word(T):
    """
    Map tournament T to a braid word.
    Strategy: read arcs in a specific order and map to braid generators.

    For vertices {0,1,...,n-1}, consider pairs in lexicographic order.
    Arc i→j (i<j) maps to sigma_j (positive crossing).
    Arc j→i (i<j) maps to sigma_j^{-1} (negative crossing).

    The braid word encodes all C(n,2) comparisons.
    """
    n = len(T)
    word = []
    for i in range(n):
        for j in range(i+1, n):
            if T[i][j] == 1:
                word.append(j)      # sigma_j
            else:
                word.append(-j)     # sigma_j^{-1}
    return word

def braid_word_writhe(word):
    """Sum of signs of crossings in braid word."""
    return sum(1 if x > 0 else -1 for x in word)

# ============================================================
# Part 7: The key new invariant: det(t*M - M^T)
# ============================================================
def compute_tournament_det_spectrum(T):
    """
    Compute det(t*M - M^T) for t = -2, -1, 0, 1, 2, 3 and also
    det(M), tr(M), eigenvalues of M, eigenvalues of M - M^T.
    """
    M = transfer_matrix(T)
    n = len(T)
    H = np.trace(M) if n == 1 else sum(M[i][j] for i in range(n) for j in range(n))
    # Actually H = sum of all entries... no, H = total ham paths
    H_total = int(np.sum(M))  # = sum of all M[a,b] = H(T)

    results = {
        'M': M.tolist(),
        'H': H_total,
        'det_M': round(np.linalg.det(M.astype(float))),
        'tr_M': int(np.trace(M)),
        'symmetric': np.allclose(M, M.T),
    }

    # det(t*M - M^T) at various t
    for t in [-2, -1, 0, 1, 2, 3]:
        mat = t * M.astype(float) - M.T.astype(float)
        results[f'det_tM-Mt_t={t}'] = round(np.linalg.det(mat))

    # Eigenvalues of antisymmetric part
    A = (M - M.T).astype(float)
    if not np.allclose(A, 0):
        eigs = np.linalg.eigvals(A)
        results['antisym_eigenvalues'] = sorted([round(e.real, 4) + round(e.imag, 4)*1j for e in eigs], key=lambda x: abs(x))
    else:
        results['antisym_eigenvalues'] = 'all zero (M symmetric)'

    # Eigenvalues of M itself
    eigs_M = np.linalg.eigvals(M.astype(float))
    results['M_eigenvalues'] = sorted([round(e.real, 4) + round(e.imag, 4)*1j for e in eigs_M], key=lambda x: -abs(x))

    return results

# ============================================================
# Main exploration
# ============================================================
def main():
    print("=" * 70)
    print("KNOT-TOURNAMENT BRIDGE EXPLORATION")
    print("kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    # -----------------------------------
    # Part A: Transfer matrix and "Alexander polynomial" for n=3,4,5
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART A: TOURNAMENT ALEXANDER POLYNOMIAL det(t*M - M^T)")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        alex_counter = Counter()
        det_M_counter = Counter()
        all_det_spectra = []

        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            M = transfer_matrix(T)
            det_M = round(np.linalg.det(M.astype(float)))
            is_sym = np.allclose(M, M.T)

            # Alexander polynomial coefficients
            alex = tuple(tournament_alexander(T))
            alex_counter[alex] += 1
            det_M_counter[(H, det_M, is_sym)] += 1

            # Collect det spectrum for each H class
            spec = {}
            for t in [-2, -1, 0, 1, 2, 3]:
                mat = t * M.astype(float) - M.T.astype(float)
                spec[t] = round(np.linalg.det(mat))
            all_det_spectra.append((H, det_M, spec, is_sym))

        print(f"Total tournaments: {2**(n*(n-1)//2)}")
        print(f"Distinct Alexander polynomials: {len(alex_counter)}")
        print(f"\nAlexander polynomial distribution:")
        for alex, count in sorted(alex_counter.items(), key=lambda x: -x[1])[:15]:
            print(f"  {list(alex)}: {count} tournaments")

        print(f"\n(H, det(M), symmetric?) distribution:")
        for (H, det, sym), count in sorted(det_M_counter.items()):
            print(f"  H={H}, det(M)={det}, sym={sym}: {count}")

    # -----------------------------------
    # Part B: Skein relation verification
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART B: SKEIN RELATION F_T = F_{T\\e} + (x-1)*F(T/e) at x=2")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        total_checks = 0
        failures = 0

        count = 0
        for bits, T in all_tournaments(n):
            count += 1
            if count > 20:  # sample for speed
                break
            results = skein_check(T)
            for r in results:
                total_checks += 1
                if not r['match']:
                    failures += 1
                    print(f"  FAILURE: arc {r['arc']}, H={r['H_T']}, H_del={r['H_del']}, H_con={r['H_con']}")

        print(f"  Checked {total_checks} arcs, failures: {failures}")

    # -----------------------------------
    # Part C: Chord diagram crossing number
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART C: CHORD DIAGRAM CROSSING NUMBER vs H(T)")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n--- n = {n} ---")
        crossing_H = defaultdict(list)
        writhe_H = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            cn = chord_intersection_number(T)
            w = tournament_writhe(T)
            crossing_H[H].append(cn)
            writhe_H[H].append(w)

        print(f"  H -> crossing numbers (min, max, mean):")
        for H in sorted(crossing_H.keys()):
            vals = crossing_H[H]
            print(f"    H={H:3d}: cross in [{min(vals):4d}, {max(vals):4d}], mean={np.mean(vals):7.2f}, count={len(vals)}")

        print(f"\n  H -> writhe (min, max, mean):")
        for H in sorted(writhe_H.keys()):
            vals = writhe_H[H]
            print(f"    H={H:3d}: writhe in [{min(vals):3d}, {max(vals):3d}], mean={np.mean(vals):6.2f}")

    # -----------------------------------
    # Part D: Braid word analysis
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART D: TOURNAMENT → BRAID WORD")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        writhe_H_pairs = []

        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            braid = tournament_to_braid_word(T)
            w = braid_word_writhe(braid)
            writhe_H_pairs.append((H, w, len(braid)))

        # Analyze
        H_by_writhe = defaultdict(list)
        for H, w, _ in writhe_H_pairs:
            H_by_writhe[w].append(H)

        print(f"  Braid writhe -> H distribution:")
        for w in sorted(H_by_writhe.keys()):
            vals = H_by_writhe[w]
            print(f"    writhe={w:3d}: H in {sorted(set(vals))}, mean={np.mean(vals):.2f}, count={len(vals)}")

    # -----------------------------------
    # Part E: Forward polynomial F(T, x) — the "Jones polynomial of tournaments"
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART E: FORWARD POLYNOMIAL F(T,x) — 'JONES POLYNOMIAL OF TOURNAMENTS'")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        poly_counter = Counter()
        special_values = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            fpoly = fwd_polynomial(T)

            # Convert to tuple for counting
            max_deg = n - 1
            coeffs = tuple(fpoly.get(k, 0) for k in range(max_deg + 1))
            poly_counter[coeffs] += 1

            # Evaluate at special points
            # F(T, -1) = "signed" count
            F_neg1 = sum((-1)**k * c for k, c in fpoly.items())
            # F(T, 0) = constant term = #paths starting with descent
            F_0 = fpoly.get(0, 0)
            # F(T, 1) = sum of all coefficients = n!
            F_1 = sum(fpoly.values())
            # F(T, 2) = H(T)
            F_2 = sum(2**k * c for k, c in fpoly.items())
            # F(T, i) where i = sqrt(-1) — complex evaluation
            F_i_real = sum((1j**k).real * c for k, c in fpoly.items())
            F_i_imag = sum((1j**k).imag * c for k, c in fpoly.items())

            special_values[H].append({
                'F(-1)': F_neg1,
                'F(0)': F_0,
                'F(1)': F_1,
                'F(2)': F_2,
                'F(i)': complex(F_i_real, F_i_imag),
                'coeffs': coeffs
            })

        print(f"  Distinct F-polynomials: {len(poly_counter)}")
        print(f"\n  Special values by H:")
        for H in sorted(special_values.keys()):
            vals = special_values[H]
            F_neg1_set = set(v['F(-1)'] for v in vals)
            F_0_set = set(v['F(0)'] for v in vals)
            print(f"    H={H:3d}: F(-1) in {sorted(F_neg1_set)}, F(0) in {sorted(F_0_set)}, count={len(vals)}")

        # Check palindrome: F(T, x) = x^{n-1} * F(T, 1/x)?
        # i.e., coeffs are palindromic?
        palindromes = 0
        for coeffs, count in poly_counter.items():
            if coeffs == coeffs[::-1]:
                palindromes += count
        print(f"\n  Palindromic F-polynomials: {palindromes}/{2**(n*(n-1)//2)}")

    # -----------------------------------
    # Part F: Vassiliev-type invariants
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART F: VASSILIEV-TYPE INVARIANTS v_k(T)")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        vass_by_H = defaultdict(list)

        max_order = min(3, n*(n-1)//2)
        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            v_vals = []
            for k in range(1, max_order + 1):
                v_vals.append(arrow_diagram_invariant(T, k))
            vass_by_H[H].append(tuple(v_vals))

        print(f"  H -> Vassiliev invariants (v_1, v_2, v_3):")
        for H in sorted(vass_by_H.keys()):
            vals = vass_by_H[H]
            distinct = len(set(vals))
            # Show a few examples
            examples = list(set(vals))[:5]
            print(f"    H={H:3d}: {distinct} distinct tuples, examples: {examples[:3]}")

    # -----------------------------------
    # Part G: Deep dive — det(M) and H relationship
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART G: DEEP DIVE — det(M), eigenvalues, and H")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        det_H = defaultdict(list)
        eig_patterns = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            M = transfer_matrix(T)
            det_val = round(np.linalg.det(M.astype(float)))

            eigs = sorted(np.linalg.eigvals(M.astype(float)).real, reverse=True)
            eig_rounded = tuple(round(e, 2) for e in eigs)

            det_H[H].append(det_val)
            eig_patterns[H].append(eig_rounded)

        for H in sorted(det_H.keys()):
            dets = sorted(set(det_H[H]))
            eigs_set = set(eig_patterns[H])
            print(f"  H={H:3d}: det(M) in {dets}, #{len(eigs_set)} eig patterns")
            if len(eigs_set) <= 3:
                for e in sorted(eigs_set):
                    print(f"         eigenvalues: {e}")

    # -----------------------------------
    # Part H: THE KEY TEST — Does det(t*M - M^T) factor like Jones?
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART H: FACTORIZATION OF det(t*M - M^T)")
    print("=" * 70)

    for n in [3, 4]:
        print(f"\n--- n = {n} ---")
        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            M = transfer_matrix(T)

            # Compute det(t*M - M^T) as polynomial
            alex = tournament_alexander(T)

            # Also compute det(M + M^T) (symmetrized = Seifert-like)
            S = M + M.T  # "Seifert matrix" analog
            det_S = round(np.linalg.det(S.astype(float)))

            # Anti-symmetrized
            A = M - M.T
            det_A = round(np.linalg.det(A.astype(float)))

            if not np.allclose(M, M.T):  # only interesting for non-symmetric
                print(f"  bits={bits}, H={H}")
                print(f"    Alex poly coeffs: {alex}")
                print(f"    det(M+M^T) = {det_S}")
                print(f"    det(M-M^T) = {det_A}")
                # Check: does Alex poly have nice roots?
                if any(a != 0 for a in alex):
                    roots = np.roots(list(reversed(alex)))
                    print(f"    Alex roots: {[round(r.real, 4) + round(r.imag, 4)*1j for r in roots]}")

    # -----------------------------------
    # Part I: F(T, omega) and Jones at roots of unity
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART I: F(T, omega) — EVALUATION AT ROOTS OF UNITY")
    print("=" * 70)

    omega = np.exp(2j * np.pi / 3)  # cube root of unity
    omega4 = np.exp(2j * np.pi / 4)  # 4th root
    omega5 = np.exp(2j * np.pi / 5)  # 5th root

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        root_vals = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            fpoly = fwd_polynomial(T)

            # Evaluate at roots of unity
            F_om3 = sum(omega**k * c for k, c in fpoly.items())
            F_om4 = sum(omega4**k * c for k, c in fpoly.items())
            F_om5 = sum(omega5**k * c for k, c in fpoly.items())

            root_vals[H].append({
                'F(omega3)': F_om3,
                'F(omega4)': F_om4,
                'F(omega5)': F_om5,
                '|F(omega3)|': abs(F_om3),
                '|F(omega4)|': abs(F_om4),
            })

        print(f"  H -> |F(omega_3)|, |F(omega_4)| (modulus):")
        for H in sorted(root_vals.keys()):
            vals = root_vals[H]
            mod3 = [v['|F(omega3)|'] for v in vals]
            mod4 = [v['|F(omega4)|'] for v in vals]
            print(f"    H={H:3d}: |F(w3)| in [{min(mod3):.2f}, {max(mod3):.2f}], |F(w4)| in [{min(mod4):.2f}, {max(mod4):.2f}]")

    # -----------------------------------
    # Part J: NEW — "Kauffman state sum" for tournaments
    # -----------------------------------
    print("\n" + "=" * 70)
    print("PART J: KAUFFMAN STATE SUM — H-GENERATING FUNCTION OVER ALL TOURNAMENTS")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        H_distribution = Counter()
        for bits, T in all_tournaments(n):
            H = len(ham_paths(T))
            H_distribution[H] += 1

        print(f"  K_n(A) = sum_T A^{{H(T)}} = sum over all 2^C(n,2) tournaments")
        print(f"  = ", end="")
        terms = []
        for H in sorted(H_distribution.keys()):
            count = H_distribution[H]
            terms.append(f"{count}*A^{H}")
        print(" + ".join(terms))

        # Factor out common factors
        total = sum(H_distribution.values())
        print(f"  Total terms: {total}")

        # Check: K_n(1) = 2^C(n,2) (trivial)
        K_1 = sum(H_distribution.values())
        print(f"  K_n(1) = {K_1} = 2^C({n},2) = {2**(n*(n-1)//2)} ✓" if K_1 == 2**(n*(n-1)//2) else f"  K_n(1) = {K_1} ERROR")

        # K_n(2) would need 2^H which is huge, skip

        # Check if K_n has nice structure
        H_vals = sorted(H_distribution.keys())
        print(f"  H values: {H_vals}")
        print(f"  Gaps: {[H_vals[i+1] - H_vals[i] for i in range(len(H_vals)-1)]}")

        # All H are odd (Redei), check
        all_odd = all(H % 2 == 1 for H in H_vals)
        print(f"  All H odd: {all_odd}")

        # Substitute A = -1: alternating sum
        K_neg1 = sum((-1)**H * count for H, count in H_distribution.items())
        print(f"  K_n(-1) = sum (-1)^H * count = {K_neg1}")
        # Since all H odd, this is -sum(count) = -2^C(n,2)
        print(f"    = -{total} (since all H odd)")

    print("\n" + "=" * 70)
    print("EXPLORATION COMPLETE")
    print("=" * 70)

if __name__ == '__main__':
    main()
