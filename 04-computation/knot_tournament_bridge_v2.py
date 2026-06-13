"""
knot_tournament_bridge_v2.py -- kind-pasteur-2026-03-14-S69
Fixed version: correct DC convention, UTF-8 output, deeper analysis.

Key findings from v1:
- Tournament Alexander poly det(t*M - M^T) has palindromic structure
- Chord crossing number is CONSTANT (depends only on n, not tournament)
- H-maximizers have lower |writhe| range (more balanced)
- DC convention needed fixing
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
import sys

# Force UTF-8 output
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

def ham_paths(T):
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

def H_count(T):
    return len(ham_paths(T))

def transfer_matrix(T):
    n = len(T)
    paths = ham_paths(T)
    M = np.zeros((n, n), dtype=int)
    for p in paths:
        M[p[0]][p[-1]] += 1
    return M

def fwd_polynomial(T):
    n = len(T)
    paths = ham_paths(T)
    poly = Counter()
    for p in paths:
        fwd = sum(1 for i in range(n-1) if p[i] < p[i+1])
        poly[fwd] += 1
    return dict(sorted(poly.items()))

def score_sequence(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)]))

# ============================================================
# CORRECTED Deletion-Contraction (THM-082 convention)
# ============================================================
def ham_paths_digraph(D):
    """Count Hamiltonian paths in a general digraph (may have missing arcs)."""
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

def delete_arc_digraph(T, u, v):
    """Delete arc u->v. Returns digraph (NOT a tournament)."""
    D = T.copy()
    D[u][v] = 0
    return D

def contract_arc_correct(T, u, v):
    """
    Contract arc u->v in tournament T.
    THM-082 convention: merged vertex w inherits IN from v (head), OUT from u (tail).
    For x != u,v:
      w->x iff u->x (tail's outgoing)
      x->w iff x->v (head's incoming)
    """
    n = len(T)
    remaining = [x for x in range(n) if x != u and x != v]
    new_n = n - 1
    T_new = np.zeros((new_n, new_n), dtype=int)

    # w is index 0, remaining are 1..new_n-1
    idx_map = {}
    idx_map['w'] = 0
    for i, x in enumerate(remaining):
        idx_map[x] = i + 1

    # Arcs among remaining vertices
    for x in remaining:
        for y in remaining:
            if x != y:
                T_new[idx_map[x]][idx_map[y]] = T[x][y]

    # w -> x: iff u -> x
    for x in remaining:
        T_new[0][idx_map[x]] = T[u][x]

    # x -> w: iff x -> v
    for x in remaining:
        T_new[idx_map[x]][0] = T[x][v]

    return T_new

def verify_dc_formula(T, u, v):
    """
    THM-082: H(T) = H(T\\e) + H(T/e) where e = u->v
    T\\e = delete arc u->v from T (digraph with n vertices)
    T/e = contract arc u->v (tournament with n-1 vertices)
    """
    H_T = H_count(T)
    D_del = delete_arc_digraph(T, u, v)
    H_del = ham_paths_digraph(D_del)
    T_con = contract_arc_correct(T, u, v)
    H_con = ham_paths_digraph(T_con)  # Use digraph counter (contraction may not be tournament)
    return H_T, H_del, H_con, H_T == H_del + H_con

# ============================================================
# PART 1: Correct DC/Skein verification
# ============================================================
def part1_dc_skein():
    print("=" * 70)
    print("PART 1: DELETION-CONTRACTION (SKEIN RELATION) VERIFICATION")
    print("  THM-082: H(T) = H(T\\e) + H(T/e)")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        total = 0
        failures = 0
        for bits, T in all_tournaments(n):
            for u in range(n):
                for v in range(n):
                    if u != v and T[u][v] == 1:
                        H_T, H_del, H_con, ok = verify_dc_formula(T, u, v)
                        total += 1
                        if not ok:
                            failures += 1
                            if failures <= 5:
                                print(f"  FAIL: arc ({u},{v}), H={H_T}, H_del={H_del}, H_con={H_con}")
        print(f"  Total checks: {total}, failures: {failures}")

# ============================================================
# PART 2: Arc-flip skein relation
# ============================================================
def part2_arcflip_skein():
    print("\n" + "=" * 70)
    print("PART 2: ARC-FLIP SKEIN RELATION")
    print("  F(T+,x) - F(T-,x) = (x-1) * D(x)")
    print("  where D(x) = F(T+/e, x) - F(T-/e', x)")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        anti_palindrome_count = 0
        total_flips = 0

        for bits, T in all_tournaments(n):
            if bits > 50:  # sample
                break
            for u in range(n):
                for v in range(u+1, n):
                    total_flips += 1

                    # T+ has u->v, T- has v->u
                    T_plus = T.copy()
                    T_minus = T.copy()
                    if T[u][v] == 1:
                        # T is already T+
                        T_minus[u][v] = 0
                        T_minus[v][u] = 1
                    else:
                        # T is T-
                        T_plus[v][u] = 0
                        T_plus[u][v] = 1
                        T_plus, T_minus = T_minus, T_plus  # swap so T_plus has u->v

                    # Forward polynomials
                    fp_plus = fwd_polynomial(T_plus)
                    fp_minus = fwd_polynomial(T_minus)

                    # Difference F(T+) - F(T-)
                    max_deg = n - 1
                    diff = {}
                    for k in range(max_deg + 1):
                        d = fp_plus.get(k, 0) - fp_minus.get(k, 0)
                        if d != 0:
                            diff[k] = d

                    # Contraction polynomials
                    T_plus_con = contract_arc_correct(T_plus, u, v)
                    T_minus_con = contract_arc_correct(T_minus, v, u)

                    fp_con_plus = fwd_polynomial(T_plus_con)
                    fp_con_minus = fwd_polynomial(T_minus_con)

                    # D = F(T+/e) - F(T-/e')
                    D = {}
                    max_deg_con = n - 2
                    for k in range(max_deg_con + 1):
                        d = fp_con_plus.get(k, 0) - fp_con_minus.get(k, 0)
                        if d != 0:
                            D[k] = d

                    # Check: diff = (x-1) * D
                    # (x-1)*D coefficients: D[k-1] - D[k]
                    expected = {}
                    for k in range(max_deg + 1):
                        val = D.get(k-1, 0) - D.get(k, 0)
                        if val != 0:
                            expected[k] = val

                    match = (diff == expected)

                    # Check anti-palindrome of D
                    D_coeffs = [D.get(k, 0) for k in range(max_deg_con + 1)]
                    D_rev = list(reversed(D_coeffs))
                    D_neg_rev = [-c for c in D_rev]
                    is_anti_palindrome = (D_coeffs == D_neg_rev)
                    if is_anti_palindrome:
                        anti_palindrome_count += 1

                    if total_flips <= 5:
                        print(f"  arc ({u},{v}): diff={diff}, D={D}")
                        print(f"    skein match: {match}, D anti-palindrome: {is_anti_palindrome}")

        print(f"  Total flips: {total_flips}")
        print(f"  D anti-palindromic: {anti_palindrome_count}/{total_flips}")

# ============================================================
# PART 3: Transfer matrix eigenvalue structure
# ============================================================
def part3_eigenvalues():
    print("\n" + "=" * 70)
    print("PART 3: TRANSFER MATRIX EIGENVALUE STRUCTURE")
    print("  M = transfer matrix, eigenvalues, det, trace")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        eig_data = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = H_count(T)
            M = transfer_matrix(T)
            eigs = sorted(np.linalg.eigvals(M.astype(float)), key=lambda x: -abs(x))
            eigs_rounded = tuple(round(e.real, 3) + round(e.imag, 3)*1j for e in eigs)

            det_val = round(np.linalg.det(M.astype(float)))
            tr_val = int(np.trace(M))
            is_sym = np.allclose(M, M.T)

            eig_data[H].append({
                'det': det_val,
                'tr': tr_val,
                'sym': is_sym,
                'eigs': eigs_rounded,
                'H/n': H / n if n > 0 else 0,
                'top_eig': eigs_rounded[0].real if eigs_rounded else 0,
            })

        for H in sorted(eig_data.keys()):
            items = eig_data[H]
            dets = set(d['det'] for d in items)
            trs = set(d['tr'] for d in items)
            syms = set(d['sym'] for d in items)
            top_eigs = sorted(set(d['top_eig'] for d in items))

            print(f"  H={H:3d}: det(M) in {sorted(dets)}, tr(M) in {sorted(trs)}, "
                  f"sym={syms}, top_eig in {top_eigs[:3]}, count={len(items)}")

# ============================================================
# PART 4: "Tournament Alexander Polynomial" deep analysis
# ============================================================
def tournament_alexander_poly(T):
    """det(t*M - M^T) as polynomial coefficients [c_0, c_1, ..., c_n]."""
    n = len(T)
    M = transfer_matrix(T)

    # Use n+1 evaluation points
    points = list(range(n + 2))
    values = []
    for t_val in points:
        mat = t_val * M.astype(float) - M.T.astype(float)
        det_val = np.linalg.det(mat)
        values.append(det_val)

    # Vandermonde interpolation
    V = np.vander(points, increasing=True)
    coeffs = np.linalg.solve(V, values)
    return [round(c) for c in coeffs]

def part4_alexander():
    print("\n" + "=" * 70)
    print("PART 4: TOURNAMENT ALEXANDER POLYNOMIAL det(t*M - M^T)")
    print("  Deeper analysis: palindrome, roots, relation to H")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        alex_by_H = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = H_count(T)
            alex = tournament_alexander_poly(T)
            alex_by_H[H].append(tuple(alex))

        for H in sorted(alex_by_H.keys()):
            distinct = set(alex_by_H[H])
            print(f"  H={H:3d}: {len(distinct)} distinct Alex polys")
            for a in sorted(distinct):
                # Check palindrome / anti-palindrome
                a_list = list(a)
                a_rev = list(reversed(a_list))
                is_palin = a_list == a_rev
                is_anti = a_list == [-c for c in a_rev]

                # Check: a(1) and a(-1)
                a_at_1 = sum(a_list)
                a_at_neg1 = sum((-1)**k * c for k, c in enumerate(a_list))

                symm = "PALINDROME" if is_palin else ("ANTI-PALINDROME" if is_anti else "neither")
                print(f"      {a_list} [{symm}] a(1)={a_at_1}, a(-1)={a_at_neg1}")

# ============================================================
# PART 5: The Key New Idea — "Tournament Jones Polynomial"
# ============================================================
def part5_jones():
    print("\n" + "=" * 70)
    print("PART 5: TOURNAMENT JONES POLYNOMIAL")
    print("  J_T(q) = det(q*M + M^T) / det(M+M^T)")
    print("  Analogy: Jones poly = det(t*V - V^T) / det(V - V^T)")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        jones_by_H = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = H_count(T)
            M = transfer_matrix(T)

            # Evaluate det(q*M + M^T) at q = -2, -1, 0, 1, 2, 3
            spec = {}
            for q in [-2, -1, 0, 1, 2, 3]:
                mat = q * M.astype(float) + M.T.astype(float)
                spec[q] = round(np.linalg.det(mat))

            jones_by_H[H].append(spec)

        for H in sorted(jones_by_H.keys()):
            items = jones_by_H[H]
            # Show all distinct spectra
            distinct = set(tuple(sorted(s.items())) for s in items)
            print(f"  H={H:3d}: {len(distinct)} distinct spectra")
            for d in sorted(list(distinct))[:3]:
                d_dict = dict(d)
                print(f"      q->det: {d_dict}")

# ============================================================
# PART 6: Braid word and "knot type" exploration
# ============================================================
def part6_braid():
    print("\n" + "=" * 70)
    print("PART 6: TOURNAMENT -> BRAID WORD ANALYSIS")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        # Compute writhe (= sum of signs) for each tournament
        writhe_by_H = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = H_count(T)
            # Writhe = #ascending arcs - #descending arcs
            asc = sum(T[i][j] for i in range(n) for j in range(i+1, n))
            desc = n*(n-1)//2 - asc
            writhe = asc - desc
            writhe_by_H[H].append(writhe)

        print(f"  H -> writhe distribution:")
        for H in sorted(writhe_by_H.keys()):
            vals = writhe_by_H[H]
            dist = Counter(vals)
            print(f"    H={H:3d}: writhe values = {dict(sorted(dist.items()))}")

        # KEY: Does |writhe| predict H?
        print(f"\n  |writhe| -> H statistics:")
        abs_writhe_H = defaultdict(list)
        for bits, T in all_tournaments(n):
            H = H_count(T)
            asc = sum(T[i][j] for i in range(n) for j in range(i+1, n))
            desc = n*(n-1)//2 - asc
            abs_writhe_H[abs(asc - desc)].append(H)

        for w in sorted(abs_writhe_H.keys()):
            vals = abs_writhe_H[w]
            print(f"    |w|={w:2d}: H_mean={np.mean(vals):6.2f}, H_range=[{min(vals)},{max(vals)}], count={len(vals)}")

# ============================================================
# PART 7: F(T, roots of unity) -- the Jones evaluation analogy
# ============================================================
def part7_roots_of_unity():
    print("\n" + "=" * 70)
    print("PART 7: F(T, roots of unity)")
    print("  Jones poly at specific roots gives linking/coloring info")
    print("  What does F(T, omega_k) tell us about tournaments?")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n--- n = {n} ---")

        omega3 = np.exp(2j * np.pi / 3)
        omega4 = np.exp(2j * np.pi / 4)
        omega5 = np.exp(2j * np.pi / 5)
        omega6 = np.exp(2j * np.pi / 6)

        root_by_H = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = H_count(T)
            fp = fwd_polynomial(T)

            vals = {}
            for name, omega in [('w3', omega3), ('w4', omega4), ('w5', omega5), ('w6', omega6)]:
                F_omega = sum(omega**k * c for k, c in fp.items())
                vals[f'|F({name})|'] = abs(F_omega)
                vals[f'F({name})'] = F_omega

            # F(-1) = signed count
            F_neg1 = sum((-1)**k * c for k, c in fp.items())
            vals['F(-1)'] = F_neg1

            root_by_H[H].append(vals)

        print(f"  H -> |F(w_k)| and F(-1):")
        for H in sorted(root_by_H.keys()):
            items = root_by_H[H]
            F_neg1_vals = sorted(set(int(v['F(-1)']) for v in items))
            mod3_vals = sorted(set(round(v['|F(w3)|'], 2) for v in items))
            mod4_vals = sorted(set(round(v['|F(w4)|'], 2) for v in items))
            print(f"    H={H:3d}: F(-1) in {F_neg1_vals}, |F(w3)| in {mod3_vals[:5]}, |F(w4)| in {mod4_vals[:5]}")

        # KEY CHECK: Is F(-1) always = F(1) = n! mod something? Or is it a finer invariant?
        print(f"\n  F(-1) distribution:")
        F_neg1_counter = Counter()
        for bits, T in all_tournaments(n):
            fp = fwd_polynomial(T)
            F_neg1 = sum((-1)**k * c for k, c in fp.items())
            F_neg1_counter[int(F_neg1)] += 1
        for val in sorted(F_neg1_counter.keys()):
            print(f"    F(-1)={val:4d}: {F_neg1_counter[val]} tournaments")

# ============================================================
# PART 8: The Kauffman State Sum = H-generating function
# ============================================================
def part8_kauffman():
    print("\n" + "=" * 70)
    print("PART 8: KAUFFMAN STATE SUM (H-GENERATING FUNCTION)")
    print("  K_n(A) = sum over all 2^C(n,2) tournaments: A^H(T)")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        print(f"\n--- n = {n} ---")
        H_dist = Counter()
        total = 2**(n*(n-1)//2)

        if n <= 5:
            for bits, T in all_tournaments(n):
                H = H_count(T)
                H_dist[H] += 1
        else:
            # Sample for n=6
            import random
            for _ in range(10000):
                bits = random.randint(0, total - 1)
                T = tournament_from_bits(n, bits)
                H = H_count(T)
                H_dist[H] += 1

        print(f"  H distribution (n={n}){'  [SAMPLED 10k]' if n > 5 else ''}:")
        for H in sorted(H_dist.keys()):
            count = H_dist[H]
            pct = 100 * count / sum(H_dist.values())
            print(f"    H={H:3d}: {count:6d} ({pct:5.2f}%)")

        # Key evaluations of K_n
        # K_n(1) = total tournaments = 2^C(n,2)
        K_1 = sum(H_dist.values())
        # K_n(-1) = sum (-1)^H * count = -K_1 (since all H odd)
        K_neg1 = sum((-1)**H * count for H, count in H_dist.items())

        print(f"\n  K_n(1) = {K_1}")
        print(f"  K_n(-1) = {K_neg1}")
        print(f"  All H odd: {all(H % 2 == 1 for H in H_dist.keys())}")

        # Mean H
        if n <= 5:
            total_H = sum(H * count for H, count in H_dist.items())
            print(f"  Mean H = {total_H / K_1:.4f}")
            print(f"  n! / 2^(n-1) = {np.math.factorial(n) / 2**(n-1):.4f}")

# ============================================================
# PART 9: NEW IDEA — Determinant of adjacency vs Jones determinant
# ============================================================
def part9_adj_det():
    print("\n" + "=" * 70)
    print("PART 9: ADJACENCY MATRIX DETERMINANT AND H")
    print("  det(A) where A is tournament adjacency matrix")
    print("  Jones: det(V-V^T) = det of knot")
    print("  Our analogy: det(A-A^T) = det(skew-symmetric part)")
    print("=" * 70)

    for n in [3, 4, 5, 6, 7]:
        print(f"\n--- n = {n} ---")
        det_by_H = defaultdict(list)
        pfaff_by_H = defaultdict(list)

        count = 0
        for bits, T in all_tournaments(n):
            count += 1
            if n >= 6 and count > 5000:
                break

            H = H_count(T)

            # Tournament adjacency: A[i][j] = 1 if i->j, 0 otherwise
            A = T.astype(float)

            # det(A)
            det_A = round(np.linalg.det(A))
            det_by_H[H].append(det_A)

            # Skew-symmetric part: S = A - A^T (has entries +1, -1, 0)
            S = A - A.T
            det_S = round(np.linalg.det(S))
            pfaff_by_H[H].append(det_S)

        print(f"  H -> det(A), det(A-A^T) (Pfaffian^2 for even n):")
        for H in sorted(det_by_H.keys()):
            det_vals = sorted(set(det_by_H[H]))
            pfaff_vals = sorted(set(pfaff_by_H[H]))
            # Truncate display
            det_str = str(det_vals[:5])
            pfaff_str = str(pfaff_vals[:5])
            if len(det_vals) > 5:
                det_str += f"... ({len(det_vals)} distinct)"
            if len(pfaff_vals) > 5:
                pfaff_str += f"... ({len(pfaff_vals)} distinct)"
            print(f"    H={H:3d}: det(A) in {det_str}, det(A-A^T) in {pfaff_str}")

        # KEY: Is det(A-A^T) constant for tournaments?
        # A-A^T is the signed adjacency: entries +-1
        # For odd n, det of skew-symmetric matrix = 0 (Pfaffian doesn't exist)
        # For even n, det = Pfaffian^2
        all_pfaffs = []
        for vals in pfaff_by_H.values():
            all_pfaffs.extend(vals)
        distinct_pfaffs = set(all_pfaffs)
        print(f"\n  ALL det(A-A^T) values: {sorted(distinct_pfaffs)[:20]}")
        if n % 2 == 1:
            print(f"  (n={n} is odd, so det of skew-symmetric = 0 always)")

# ============================================================
# PART 10: DEEP — Can we define a "Jones polynomial" J_T(t)?
# ============================================================
def part10_jones_attempt():
    print("\n" + "=" * 70)
    print("PART 10: DEFINING A 'JONES POLYNOMIAL' FOR TOURNAMENTS")
    print("  Attempt: J_T(t) = det(t*A - A^T) / (t-1)^floor(n/2)")
    print("  where A is adjacency matrix, A^T is transpose")
    print("  Note: A-A^T is skew, t*A - A^T at t=1 is A-A^T")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        jones_by_H = defaultdict(list)

        for bits, T in all_tournaments(n):
            H = H_count(T)
            A = T.astype(float)

            # Evaluate det(t*A - A^T) at t = -3, -2, -1, 0, 1, 2, 3, 4, 5
            vals = {}
            for t in [-3, -2, -1, 0, 1, 2, 3, 4, 5]:
                mat = t * A - A.T
                vals[t] = round(np.linalg.det(mat))

            jones_by_H[H].append(vals)

        for H in sorted(jones_by_H.keys()):
            items = jones_by_H[H]
            distinct = set(tuple(sorted(s.items())) for s in items)
            print(f"  H={H:3d}: {len(distinct)} distinct det spectra")
            for d in sorted(list(distinct))[:3]:
                d_dict = dict(d)
                print(f"    t -> det(t*A - A^T): {d_dict}")

    # Now check: does det(t*A - A^T) determine H?
    print("\n  QUESTION: Does det(t*A - A^T) determine H(T)?")
    for n in [4, 5]:
        print(f"\n  n = {n}:")
        det_to_H = defaultdict(set)
        for bits, T in all_tournaments(n):
            H = H_count(T)
            A = T.astype(float)

            # Use t=2 and t=3 as discriminants
            spec = []
            for t in [0, 2, 3]:
                mat = t * A - A.T
                spec.append(round(np.linalg.det(mat)))
            det_to_H[tuple(spec)].add(H)

        ambiguous = sum(1 for s in det_to_H.values() if len(s) > 1)
        print(f"    Distinct det spectra: {len(det_to_H)}")
        print(f"    Ambiguous (multiple H): {ambiguous}")
        if ambiguous > 0:
            for spec, H_set in sorted(det_to_H.items()):
                if len(H_set) > 1:
                    print(f"      spec={spec} -> H in {sorted(H_set)}")

# ============================================================
# PART 11: WRITHE vs H — Deep correlation analysis
# ============================================================
def part11_writhe_deep():
    print("\n" + "=" * 70)
    print("PART 11: WRITHE vs H — DEEP ANALYSIS")
    print("  Writhe_T = sum_{i<j} (2*T[i][j] - 1)")
    print("  = #{ascending arcs} - #{descending arcs}")
    print("  Analogy: writhe of knot diagram")
    print("=" * 70)

    for n in [5, 6]:
        print(f"\n--- n = {n} ---")
        data = []
        count = 0
        for bits, T in all_tournaments(n):
            count += 1
            if n == 6 and count > 10000:
                break
            H = H_count(T)
            asc = sum(T[i][j] for i in range(n) for j in range(i+1, n))
            writhe = 2 * asc - n*(n-1)//2
            scores = score_sequence(T)
            data.append((H, writhe, scores))

        # Correlation
        H_vals = [d[0] for d in data]
        W_vals = [d[1] for d in data]
        corr = np.corrcoef(H_vals, W_vals)[0, 1] if len(data) > 1 else 0

        print(f"  Pearson correlation(H, writhe) = {corr:.6f}")
        print(f"  Pearson correlation(H, |writhe|) = {np.corrcoef(H_vals, [abs(w) for w in W_vals])[0,1]:.6f}")
        print(f"  Pearson correlation(H, writhe^2) = {np.corrcoef(H_vals, [w**2 for w in W_vals])[0,1]:.6f}")

        # H vs score variance (score = out-degree sequence)
        var_vals = []
        for H, w, scores in data:
            var_vals.append(np.var(scores))
        print(f"  Pearson correlation(H, score_var) = {np.corrcoef(H_vals, var_vals)[0,1]:.6f}")

        # Key: H-maximizers have what writhe?
        max_H = max(H_vals)
        max_writhes = [d[1] for d in data if d[0] == max_H]
        max_scores = [d[2] for d in data if d[0] == max_H]
        print(f"\n  H-maximizers (H={max_H}):")
        print(f"    writhes: {sorted(set(max_writhes))}")
        print(f"    scores: {sorted(set(max_scores))[:5]}")

# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("KNOT-TOURNAMENT BRIDGE v2 — kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    part1_dc_skein()
    part2_arcflip_skein()
    part3_eigenvalues()
    part4_alexander()
    part5_jones()
    part6_braid()
    part7_roots_of_unity()
    part8_kauffman()
    part9_adj_det()
    part10_jones_attempt()
    part11_writhe_deep()

    print("\n" + "=" * 70)
    print("EXPLORATION COMPLETE")
    print("=" * 70)

if __name__ == '__main__':
    main()
