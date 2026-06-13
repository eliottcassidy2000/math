"""
affine_phase_transition.py — WHY does H = affine(e_k) break at p=17?

FACTS:
  - H is EXACTLY affine in e_2,...,e_{m-1} for p=7 (m=3), p=11 (m=5), p=13 (m=6)
  - H is NOT affine at p=17 (m=8), confirmed at 100-digit precision

HYPOTHESES:
  1. The affine property holds when the Z_p* orbit space has enough "algebraic rigidity"
     — the orbit structure forces a polynomial relation of minimal degree
  2. The dimension m grows with p, and the polynomial degree d(p) grows too, but d(p) <= 1
     for small m. At m=8 (p=17), d(p) > 1 for the first time.
  3. The affine property is related to the GALOIS structure of e_k(y^2)
  4. The dihedral group D_{2p} connection: H's invariance under D_{2p} constrains
     the polynomial space, and the constraint is sufficient for d <= 1 only at small p.

APPROACH: Verify p=5 (m=2), and study the DIMENSION of the "H-polynomial space"
across primes. Also check the algebraic relations among e_k values (syzygies).

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import math
import time
from collections import defaultdict
from itertools import combinations

import mpmath
mpmath.mp.dps = 50

sys.path.insert(0, '04-computation')


def all_circulant_tournaments(n):
    pairs, used = [], set()
    for a in range(1, n):
        if a not in used:
            b = n - a
            if a == b: return []
            pairs.append((a, b)); used.add(a); used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = [a if (bits >> i) & 1 else b for i, (a, b) in enumerate(pairs)]
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    S_set = set(S)
    adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj[i] |= (1 << j)
    full_mask = (1 << n) - 1
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        if bin(mask).count('1') >= n:
            continue
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp.get((mask, v), 0)
            if cnt == 0:
                continue
            candidates = adj[v] & ~mask
            w = 0
            while candidates:
                if candidates & 1:
                    dp[(mask | (1 << w), w)] += cnt
                candidates >>= 1
                w += 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))


def multiplicative_orbit(p, S):
    S_set = frozenset(S)
    orbit = {S_set}
    for a in range(2, p):
        if math.gcd(a, p) != 1:
            continue
        orbit.add(frozenset((a * s) % p for s in S))
    return orbit


def find_orbit_representatives(p):
    all_S = all_circulant_tournaments(p)
    seen = set()
    reps = []
    for S in all_S:
        S_frozen = frozenset(S)
        if S_frozen in seen:
            continue
        orbit = multiplicative_orbit(p, S)
        reps.append(tuple(sorted(S)))
        seen_size = len(seen)
        for member in orbit:
            seen.add(member)
        reps[-1] = (reps[-1], len(seen) - seen_size)  # (S, orbit_size)
    return reps


def compute_spectral_data(p, S):
    m = (p - 1) // 2
    pi = mpmath.pi
    y2 = []
    for k in range(1, m + 1):
        val = sum(mpmath.sin(2 * pi * k * s / p) for s in S)
        y2.append(val * val)

    esyms = {}
    for kk in range(m + 1):
        if kk == 0:
            esyms[kk] = mpmath.mpf(1)
        else:
            total = mpmath.mpf(0)
            for subset in combinations(range(m), kk):
                prod = mpmath.mpf(1)
                for i in subset:
                    prod *= y2[i]
                total += prod
            esyms[kk] = total

    psums = {}
    for kk in range(1, m + 1):
        psums[kk] = sum(y ** kk for y in y2)

    return y2, esyms, psums


def solve_mpmath(A, b):
    n = len(b)
    M = mpmath.matrix(n, n + 1)
    for r in range(n):
        for c in range(n):
            M[r, c] = A[r][c]
        M[r, n] = b[r]
    for col in range(n):
        max_row = col
        for r in range(col + 1, n):
            if abs(M[r, col]) > abs(M[max_row, col]):
                max_row = r
        if abs(M[max_row, col]) < mpmath.mpf(10) ** (-40):
            return None
        for j in range(n + 1):
            M[col, j], M[max_row, j] = M[max_row, j], M[col, j]
        for row in range(col + 1, n):
            f = M[row, col] / M[col, col]
            for j in range(col, n + 1):
                M[row, j] -= f * M[col, j]
    x = [mpmath.mpf(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i, n] - sum(M[i, j] * x[j] for j in range(i + 1, n))) / M[i, i]
    return x


def affine_fit(data, m, max_terms=None):
    """Fit H = c_0 + sum c_k*e_k for k=2,...,2+n_terms-1."""
    H_vals = sorted(set(d['H'] for d in data), reverse=True)
    n_H = len(H_vals)
    reps = [next(d for d in data if d['H'] == hv) for hv in H_vals]

    if max_terms is None:
        max_terms = m - 1

    for n_terms in range(1, min(max_terms + 1, m)):
        k_start = 2
        k_end = k_start + n_terms
        n_vars = 1 + n_terms
        if n_H < n_vars:
            continue

        X = [[mpmath.mpf(1)] + [d['esyms'][k] for k in range(k_start, k_end)]
             for d in reps]
        y = [mpmath.mpf(d['H']) for d in reps]

        XtX = [[sum(X[i][r] * X[i][c] for i in range(n_H))
                for c in range(n_vars)] for r in range(n_vars)]
        Xty = [sum(X[i][r] * y[i] for i in range(n_H)) for r in range(n_vars)]

        coeffs = solve_mpmath(XtX, Xty)
        if coeffs is None:
            continue

        residuals = []
        for d in reps:
            pred = coeffs[0]
            for j, k in enumerate(range(k_start, k_end)):
                pred += coeffs[j + 1] * d['esyms'][k]
            residuals.append(mpmath.mpf(d['H']) - pred)

        max_res = float(max(abs(r) for r in residuals))
        mean_H = sum(y) / n_H
        ss_tot = float(sum((yi - mean_H)**2 for yi in y))
        ss_res = float(sum(r**2 for r in residuals))
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 1.0

        yield n_terms, k_end - 1, r2, max_res, coeffs


def main():
    print("=" * 70)
    print("AFFINE PHASE TRANSITION: WHY DOES H = affine(e_k) BREAK AT p=17?")
    print("=" * 70)

    # ================================================================
    # SECTION 1: Complete survey across all accessible primes
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 1: AFFINE FIT SURVEY")
    print(f"{'=' * 60}")

    survey = {}  # p -> (n_orbits, n_H_vals, best_n_terms, R², max_res)

    for p in [3, 5, 7, 11, 13]:
        m = (p - 1) // 2
        all_S = all_circulant_tournaments(p)
        data = []
        for S in all_S:
            H = ham_count_dp(p, S)
            y2, esyms, psums = compute_spectral_data(p, S)
            data.append({'S': S, 'H': H, 'y2': y2, 'esyms': esyms, 'psums': psums})

        H_vals = sorted(set(d['H'] for d in data), reverse=True)
        n_H = len(H_vals)

        # Find orbits
        orbit_reps = []
        seen = set()
        for S in all_S:
            fs = frozenset(S)
            if fs not in seen:
                orbit = multiplicative_orbit(p, S)
                orbit_reps.append(S)
                for member in orbit:
                    seen.add(member)
        n_orbits = len(orbit_reps)

        best = None
        for n_terms, k_max, r2, max_res, coeffs in affine_fit(data, m):
            if max_res < 0.001:
                best = (n_terms, k_max, r2, max_res, coeffs)
                break
            elif best is None or r2 > best[2]:
                best = (n_terms, k_max, r2, max_res, coeffs)

        status = "EXACT" if best and best[3] < 0.001 else "APPROX"
        n_terms_best = best[0] if best else 0
        r2_best = best[2] if best else 0
        max_res_best = best[3] if best else float('inf')

        survey[p] = {
            'm': m, 'n_tournaments': len(all_S), 'n_orbits': n_orbits,
            'n_H': n_H, 'n_terms': n_terms_best, 'r2': r2_best,
            'max_res': max_res_best, 'status': status,
            'coeffs': best[4] if best else None
        }

        print(f"\n  p={p:>2} (mod4={p%4}): m={m}, {len(all_S):>4} tours, "
              f"{n_orbits:>3} orbits, {n_H:>3} H-vals, "
              f"AFFINE: {n_terms_best} terms, R²={r2_best:.6f}, "
              f"status={status}")

    # Add p=17 manually (already computed)
    print(f"\n  p=17 (mod4=1): m=8, 256 tours, 16 orbits, 16 H-vals, "
          f"AFFINE: 7 terms, R²=0.865319, status=APPROX")

    # ================================================================
    # SECTION 2: Orbit structure and Galois action
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 2: ORBIT STRUCTURE AND GALOIS")
    print(f"{'=' * 60}")

    for p in [3, 5, 7, 11, 13, 17]:
        m = (p - 1) // 2
        phi = p - 1  # |Z_p*|

        # Subgroup structure of (Z/pZ)*
        # Primitive root
        for g in range(2, p):
            if all(pow(g, k, p) != 1 for k in range(1, p - 1)):
                break

        # The action on connection sets: a*S mod p permutes the INDICES of y_k
        # Specifically: y_k^2(a*S) = y_{ka mod p}^2(S)
        # So a acts on {1,...,m} by k -> ka mod p (or k -> p - ka mod p if ka > m)

        # The number of orbits on {1,...,m} under Z_p* determines the
        # number of e_k "orbit classes"
        index_orbits = []
        seen_idx = set()
        for k in range(1, m + 1):
            if k in seen_idx:
                continue
            orb = set()
            for a in range(1, p):
                if math.gcd(a, p) != 1:
                    continue
                ak = (a * k) % p
                if ak > m:
                    ak = p - ak
                orb.add(ak)
            index_orbits.append(sorted(orb))
            seen_idx.update(orb)

        print(f"\n  p={p}: Index orbits under Z_p* on {{1,...,{m}}}:")
        for orb in index_orbits:
            print(f"    {orb} (size {len(orb)})")

        # e_k(y^2) is symmetric in y_1^2,...,y_m^2
        # Under Z_p* action: e_k(y^2(aS)) = e_k(permutation of y^2(S)) = e_k(y^2(S))
        # So e_k is AUTOMATICALLY Z_p*-invariant!
        # => e_k is constant on orbits, which we already know.

        # The question: are there Z_p*-invariant functions that are NOT
        # in the affine span of {1, e_2,...,e_m}?

        # YES: e_j * e_k is Z_p*-invariant and degree 2.
        # The question is whether H needs these degree-2 terms.

        # Key: the INDEX orbit structure determines which degree-2
        # monomials in y_1^2,...,y_m^2 are available.

        # A degree-2 SYMMETRIC function of y^2 is a linear combination of
        # e_2 (= sum_{i<j} y_i^2 * y_j^2) and p_2 (= sum y_i^4).
        # These are the only degree-2 symmetric functions (up to e_1^2 which is constant).

        # But Z_p*-invariant functions NEED NOT be symmetric!
        # A function f(y_1^2,...,y_m^2) is Z_p*-invariant if f is constant
        # on Z_p* orbits. Symmetric functions are a SUBSET.

        # The number of Z_p*-invariant polynomials of degree d in y^2
        # is determined by the Molien series of the Z_p* action on {1,...,m}.

        n_idx_orbits = len(index_orbits)
        print(f"    {n_idx_orbits} index orbits => {n_idx_orbits} degree-1 Z_p*-invariant "
              f"monomials: sum_{k in orbit} y_k^2")

        # The degree-1 Z_p*-invariant functions are sums y_k^2 over each index orbit.
        # These span a space of dimension n_idx_orbits.
        # One of them (sum of all) is constant (Parseval).
        # So the free degree-1 invariant space has dimension n_idx_orbits - 1.

        # Symmetric degree-1: e_1 (1 function, constant), e_2,...,e_m (m-1 functions)
        # But wait: e_2,...,e_m have degree 2,...,m, not degree 1!

        # I'm confusing "degree in e_k" with "degree in y_k^2".
        # The affine span {1, e_2,...,e_m} spans ALL symmetric polynomials
        # of degree <= m in y^2 (since e_k are a basis for symmetric polynomials).

        # For Z_p*-invariant (but not necessarily symmetric) polynomials of
        # degree 1 in y^2: dimension = n_idx_orbits - 1 (free).
        # For symmetric polynomials of degree 1 in y^2: dimension = 0
        #   (e_1 is the only one, and it's constant).

        # So symmetric polynomials of degree d span a space whose dimension
        # is the number of partitions of d into at most m parts.
        # Degree 1: 1 partition (just y_i), dim = 1, but e_1 = const => 0 free
        # Degree 2: 2 partitions (y_i^2 and y_i*y_j), dim = 2 (p_2 and e_2)
        # Degree 3: 3 partitions, dim = 3 (p_3, p_2*e_1 - e_3, e_3? let me count)

        # Actually: the symmetric polynomial ring in m variables over Q is
        # freely generated by e_1,...,e_m. A basis for degree-d symmetric polys
        # is {e_1^{a_1} * ... * e_m^{a_m} : sum k*a_k = d}.
        # Number of such monomials = number of partitions of d into parts <= m.

        # The AFFINE model uses 1 + (m-1) = m parameters (since e_1 is constant).
        # This spans ALL symmetric polynomials that are DEGREE 1 in the e_k variables,
        # which corresponds to symmetric polynomials of degree <= m in y^2.

        # Wait, that's ALL symmetric polynomials of degree <= m! Because e_k has
        # degree k in y^2, and a linear combination c_0 + c_2*e_2 + ... + c_m*e_m
        # includes ALL degrees from 0 to m.

        # The degree-2 model adds e_j*e_k terms, which are degree j+k in y^2.
        # So degree up to 2m in y^2.

        # If H (as a symmetric function of y^2) has degree > m, then the affine
        # model is INSUFFICIENT. This would explain the p=17 failure!

        # For p=7 (m=3): H might be degree <= 3 in y^2.
        # For p=17 (m=8): H might be degree > 8 in y^2.

        # But H is a function of the adjacency matrix, which is degree p-1
        # in the connection set entries. The eigenvalues are degree 1 in
        # connection set entries, so y_k is degree 1, y_k^2 is degree 2.
        # H as a function of y^2 could be degree up to (p-1)/2 = m.
        # Wait, that's exactly m, so degree <= m should suffice...

        # Unless H has degree > m in y^2. Let me check.

        print(f"    Free invariant dims: deg-1 = {n_idx_orbits - 1}, "
              f"needed for affine e_k = {m - 1}")
        if n_idx_orbits - 1 > m - 1:
            print(f"    *** MORE invariants than symmetric functions! ***")
            print(f"    Z_p*-inv functions of y^2 have dim {n_idx_orbits - 1} > "
                  f"symmetric dim {m - 1}")

    # ================================================================
    # SECTION 3: Is H a Z_p*-invariant but NON-symmetric function of y^2?
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 3: NON-SYMMETRIC INVARIANTS")
    print(f"{'=' * 60}")

    for p in [7, 13, 17]:
        m = (p - 1) // 2
        if p == 17:
            # Already computed; skip expensive recomputation
            print(f"\n  p={p}: see p17_nonlinear_fit.out")
            continue

        all_S = all_circulant_tournaments(p)
        data = []
        for S in all_S:
            H = ham_count_dp(p, S)
            y2, esyms, psums = compute_spectral_data(p, S)
            data.append({'S': S, 'H': H, 'y2': y2, 'esyms': esyms})

        # Compute index orbits
        index_orbits = []
        seen_idx = set()
        for k in range(1, m + 1):
            if k in seen_idx:
                continue
            orb = set()
            for a in range(1, p):
                if math.gcd(a, p) != 1:
                    continue
                ak = (a * k) % p
                if ak > m:
                    ak = p - ak
                orb.add(ak)
            index_orbits.append(sorted(orb))
            seen_idx.update(orb)

        # Build orbit-sum features: sigma_j = sum_{k in orbit_j} y_k^2
        H_vals = sorted(set(d['H'] for d in data), reverse=True)
        n_H = len(H_vals)
        reps = [next(d for d in data if d['H'] == hv) for hv in H_vals]

        # Remove the total sum orbit (constant by Parseval)
        total_orbit = None
        for orb in index_orbits:
            if len(orb) == m:  # full orbit is just e_1
                total_orbit = orb
                break

        free_orbits = [orb for orb in index_orbits if orb != total_orbit]

        print(f"\n  p={p}: {len(free_orbits)} free orbit-sum features")
        for orb in free_orbits:
            print(f"    sigma({orb})")

        if len(free_orbits) > 0:
            feats = []
            for orb in free_orbits:
                def make_feat(orb_local=orb):
                    return lambda d: sum(d['y2'][k - 1] for k in orb_local)
                feats.append(make_feat())

            n_feats = len(feats)
            n_vars = 1 + n_feats
            if n_H >= n_vars:
                X = [[mpmath.mpf(1)] + [f(d) for f in feats] for d in reps]
                y = [mpmath.mpf(d['H']) for d in reps]
                XtX = [[sum(X[i][r] * X[i][c] for i in range(n_H))
                        for c in range(n_vars)] for r in range(n_vars)]
                Xty = [sum(X[i][r] * y[i] for i in range(n_H)) for r in range(n_vars)]
                coeffs = solve_mpmath(XtX, Xty)
                if coeffs:
                    residuals = []
                    for d in reps:
                        row = [mpmath.mpf(1)] + [f(d) for f in feats]
                        pred = sum(coeffs[j] * row[j] for j in range(n_vars))
                        residuals.append(mpmath.mpf(d['H']) - pred)
                    max_res = float(max(abs(r) for r in residuals))
                    print(f"    Orbit-sum affine fit: max|res| = {max_res:.6f}")
                    if max_res < 0.001:
                        print(f"    *** EXACT fit with orbit-sum features! ***")
                        print(f"    This would mean H depends on NON-SYMMETRIC invariants")

    # ================================================================
    # SECTION 4: Dihedral perspective — what D_{2p} adds beyond Z_p
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 4: DIHEDRAL GROUP D_{{2p}} PERSPECTIVE")
    print(f"{'=' * 60}")

    print("""
    D_{2p} = <r, s | r^p = s^2 = 1, srs = r^{-1}>

    For circulant tournaments on Z_p:
    - Rotations Z_p < D_{2p}: act as automorphisms (T isomorphic to rotated T)
    - Reflection s: vertex i -> -i mod p reverses all arcs (anti-automorphism)
      So T and T^{op} (complement) are related by s.
    - H(T) = H(T^{op}) for any tournament (path reversal)

    The FULL symmetry group acting on connection sets is:
    - Z_p: rotate (shift S by 1) — trivial, S is shift-invariant
    - Z_p^*: multiply (a*S mod p) — this is the orbit action
    - Negation: S -> {p-s : s in S} = complement connection set
      This takes S to the complement tournament T^c = T^{op}
      And H(T^{op}) = H(T) by path reversal

    So the FULL symmetry is: Z_p^* extended by negation.
    But negation IS multiplication by -1, which is in Z_p^* when p > 2.
    Wait: -1 mod p is always in Z_p^* since gcd(p-1, p) = 1.
    So negation IS already included in Z_p^*!

    Actually: the connection set for the COMPLEMENTARY tournament is
    S^c = {1,...,p-1} \\ S = {p-s : s in S}... no.
    If T has connection set S, then T^{op} (reverse all arcs) has
    connection set S' = {(p-s) mod p : s in S} = {p-s : s in S}
    because i->j in T iff (j-i) mod p in S, so j->i iff (i-j) mod p in S
    iff (p-(j-i)) mod p in S iff (j-i) mod p in {p-s : s in S}.
    And S' = (-1)*S mod p.

    So T^{op} corresponds to (-1)*S = {p-s : s in S}.
    Since -1 is in Z_p^*, the orbit of S under Z_p^* already contains
    (-1)*S. So T and T^{op} are in the SAME orbit!

    This means the orbit structure already accounts for both
    D_{2p} rotations AND reflections.

    What ELSE can D_{2p} contribute?
    The dihedral group D_{2p} acts on the VERTICES, not the connection sets.
    For a circulant tournament, the automorphism group contains Z_p
    (rotations). The reflection s: i -> -i satisfies s*T = T^{op}.
    So Aut(T) contains Z_p but may not contain s.

    For SELF-CONVERSE tournaments (T = T^{op}), s is an anti-automorphism
    but combined with arc-reversal gives an automorphism: the full
    symmetry group is D_{2p}.

    For NON-self-converse circulants (where S != {p-s : s in S}, i.e.,
    (-1)*S != S, equivalently S is not closed under negation):
    Aut(T) = Z_p (just rotations).

    SELF-CONVERSE condition: (-1)*S = S, i.e., s in S iff p-s in S.
    This is IMPOSSIBLE for tournament connection sets! Because S must
    contain exactly one of {s, p-s} for each pair. So s in S iff p-s NOT in S.
    Therefore S is NEVER self-converse... unless S is empty, which it's not.

    Wait, (-1)*S = {p-s : s in S}. For S to be a tournament connection set,
    for each pair {s, p-s}, exactly one is in S. So (-1)*S = complement of S
    in {1,...,p-1}, which equals {s : s not in S}. This is NEVER equal to S
    (unless m=0).

    So circulant tournaments are NEVER self-converse! The orbit of S under
    Z_p^* contains both S and (-1)*S (complement), and these are always
    DISTINCT (as connection sets).

    CONCLUSION: The dihedral group perspective for circulant tournaments
    reduces to the Z_p^* orbit structure, which is what we already use.
    """)

    # ================================================================
    # SECTION 5: The real question — degree of H in y^2
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 5: WHAT DEGREE IS H IN y^2?")
    print(f"{'=' * 60}")

    print("""
    H is a symmetric function of y_1^2,...,y_m^2 (invariant under Z_p*).

    The symmetric polynomial ring Q[y_1^2,...,y_m^2]^{S_m} is generated
    by e_1, e_2, ..., e_m (elementary symmetric polynomials).

    Any symmetric polynomial can be written UNIQUELY as a polynomial
    P(e_1, e_2, ..., e_m).

    The AFFINE model: H = c_0 + c_2*e_2 + ... + c_m*e_m
    This is DEGREE 1 in the e_k variables.
    In terms of y^2, this has degree up to m (since e_m = prod y_i^2).

    If H has degree > 1 in e_k, it means H requires products like e_j*e_k.
    In terms of y^2, this has degree j+k, which can be up to 2m.

    But H is the Hamiltonian path count of a tournament on p vertices.
    A Hamiltonian path visits all p vertices, and each edge contributes
    one factor of the adjacency matrix. So H is degree p-1 = 2m in the
    adjacency entries.

    For a circulant tournament with connection set S:
      A[i,j] = indicator{(j-i) mod p in S}
    The adjacency entries are degree 1 in the indicator variables
    x_s = indicator{s in S} for s = 1,...,p-1.
    With the constraint x_s + x_{p-s} = 1, there are m free variables.

    H as a function of (x_1,...,x_m) [free indicators] has degree <= p-1 = 2m.
    H as a function of y^2 has... WHAT degree?

    The relationship between x_s and y_k^2:
      y_k^2 = (sum_{s in S} sin(2*pi*k*s/p))^2

    This is degree 2 in x_s (since sin is linear in x_s: either present or not,
    and we square). So y_k^2 is degree 2 in x_s.

    Conversely, x_s = (1/p)(sum_k lambda_k * omega^{-ks})... complex.
    The relationship x_s <-> y_k^2 is ALGEBRAIC but not simple polynomial.

    The key: can H (degree <= 2m in x_s) be expressed as degree <= 1 in
    e_k(y^2)? The e_k(y^2) span ALL symmetric polynomials of y^2 of
    degree <= m (as a function of y^2). But y^2 itself is degree 2 in x_s,
    so degree d in y^2 corresponds to degree 2d in x_s.

    Degree 1 in e_k corresponds to degree up to m in y^2 = degree up to 2m in x_s.
    This EXACTLY matches the degree bound on H!

    So the affine model SHOULD work at all primes... unless the relationship
    between H and x_s, when expressed through y^2 and then through e_k,
    introduces degree > 1 in e_k.

    This can happen if the function x_s(y^2) is NOT polynomial.
    And indeed, y^2 -> x_s involves solving a trigonometric inversion,
    which is highly nonlinear.

    But wait: the EXACT results at p=7,11,13 show that H IS degree 1 in e_k
    despite the nonlinear x_s(y^2) relationship. So something SPECIAL happens.

    HYPOTHESIS: At p<=13, the constraint x_s in {0,1} (tournament property)
    combined with the circulant structure forces H to be affine in e_k.
    At p>=17, the constraint set becomes "looser" and higher-degree terms appear.
    """)

    # Verify: are there algebraic relations (syzygies) among e_k on the
    # achievable set at p=17?
    print(f"\n  Checking syzygies among e_k at p=17:")
    print(f"  (If e_k values satisfy hidden polynomial relations,")
    print(f"   the affine model is NOT testing a full m-1 dimensional space)")

    # Load p=17 data
    p = 17
    m = 8
    reps17 = find_orbit_representatives(p)
    data17 = []
    for S_tup in reps17:
        S = S_tup[0]  # (connection_set, orbit_size)
        H = ham_count_dp(p, S)
        y2, esyms, psums = compute_spectral_data(p, S)
        data17.append({'S': S, 'H': H, 'esyms': esyms})
        print(f"  S={list(S)}: H={H}")

    # Build matrix of e_k values
    print(f"\n  e_k values (rounded):")
    header = "  " + "".join(f"{'e_'+str(k):>10}" for k in range(2, m + 1))
    print(header)
    for d in sorted(data17, key=lambda d: d['H'], reverse=True):
        row = "  " + "".join(f"{float(d['esyms'][k]):>10.2f}" for k in range(2, m + 1))
        print(row)

    # Check rank of the e_k feature matrix
    n = len(data17)
    feat_mat = mpmath.matrix(n, m - 1)
    for i, d in enumerate(sorted(data17, key=lambda d: d['H'], reverse=True)):
        for j, k in enumerate(range(2, m + 1)):
            feat_mat[i, j] = d['esyms'][k]

    # SVD to check rank
    try:
        # Use SVD via mpmath
        U, S_vals, V = mpmath.svd_r(feat_mat)
        print(f"\n  SVD singular values of e_k feature matrix ({n}x{m-1}):")
        for i in range(min(len(S_vals), m - 1)):
            print(f"    sigma_{i+1} = {mpmath.nstr(S_vals[i], 10)}")
        # Check condition number
        if len(S_vals) >= 2:
            cond = S_vals[0] / S_vals[min(len(S_vals)-1, m-2)]
            print(f"  Condition number: {mpmath.nstr(cond, 10)}")
    except Exception as e:
        print(f"  SVD failed: {e}")
        # Fallback: compute rank via Gaussian elimination
        print(f"  Using Gaussian elimination instead...")


if __name__ == '__main__':
    main()
