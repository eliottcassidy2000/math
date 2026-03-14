"""
degree_drop_n6.py -- kind-pasteur-2026-03-14-S68
Compute H(T) as a SYMBOLIC multilinear polynomial in C(6,2)=15 arc variables
for n=6 tournaments.

Key question: is the maximum degree n-2=4 (degree drop!) or n-1=5 (no drop)?

At n=4: degree = 2 = n-2 (degree drop by 1)
At n=5: degree = 4 = n-1 (no drop)
What happens at n=6?

Approach: expand all 720 permutations, each contributing a product of 5 arc
indicators (each x_{ij} or (1-x_{ij})), and collect coefficients by monomial.

Also analyze:
- Which monomials survive at each degree
- Triangle triples at degree 3 (always vanish) vs non-triangle
- If degree drops, how many quintic monomials would exist before cancellation
"""

import sys
from itertools import permutations, combinations
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("DEGREE DROP ANALYSIS: H(T) AS SYMBOLIC POLYNOMIAL AT n=6")
    print("kind-pasteur-2026-03-14-S68")
    print("=" * 70)

    n = 6
    # C(6,2) = 15 arc variables: x_{ij} for i < j
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    arc_idx = {a: i for i, a in enumerate(arcs)}
    num_arcs = len(arcs)

    print(f"\n  n = {n}")
    print(f"  Number of arc variables: {num_arcs}")
    print(f"  Arcs: {arcs}")
    print(f"  n! = {720} permutations to expand")
    print(f"  Each path uses {n-1} = 5 arcs")
    print(f"  Each product expands to 2^5 = 32 terms")
    print(f"  Total terms before collection: 720 * 32 = 23040")

    # For each permutation, build the list of factors
    # Each factor is (arc_index, positive?)
    # positive=True means the factor is x_{arc}, False means (1 - x_{arc})

    print(f"\n  Expanding H(T) as symbolic multilinear polynomial...")

    # Polynomial: maps frozenset of arc indices -> integer coefficient
    poly = Counter()
    num_perms = 0

    for perm in permutations(range(n)):
        num_perms += 1
        factors = []
        for step in range(n - 1):
            i, j = perm[step], perm[step + 1]
            if i < j:
                factors.append((arc_idx[(i, j)], True))   # x_{ij}
            else:
                factors.append((arc_idx[(j, i)], False))   # (1 - x_{ji})

        # Expand product of 5 terms, each either x_k or (1 - x_k)
        # Iterate over 2^5 = 32 subsets
        num_factors = len(factors)
        for mask in range(1 << num_factors):
            coeff = 1
            arc_set = set()
            valid = True
            for bit_pos in range(num_factors):
                arc_id, positive = factors[bit_pos]
                if mask & (1 << bit_pos):
                    # Choose the x_{arc_id} part
                    if positive:
                        # factor is x, choosing x: coeff *= 1
                        pass
                    else:
                        # factor is (1-x), choosing -x: coeff *= -1
                        coeff *= -1
                    if arc_id in arc_set:
                        # Multilinear: x_i^2 = x_i, but in our setting
                        # two steps of the same path can't use the same arc
                        # Actually in a Hamiltonian path, can two consecutive
                        # edges use the same unordered pair? No! Because
                        # perm[step] != perm[step+2] (all distinct vertices).
                        # But (i,j) and (j,i) map to the same arc variable!
                        # E.g., path 0->1->0 impossible (vertices distinct in perm).
                        # So this shouldn't happen in a permutation path.
                        # But let's handle it correctly anyway.
                        pass  # arc already in set, multilinear collapse
                    arc_set.add(arc_id)
                else:
                    # Choose the constant 1 part
                    if positive:
                        # factor is x, choosing constant 1: but x has no constant
                        # x = 0 + 1*x, so choosing "constant" = 0
                        coeff = 0
                        valid = False
                        break
                    else:
                        # factor is (1-x), choosing 1: coeff *= 1
                        pass

            if valid and coeff != 0:
                poly[frozenset(arc_set)] += coeff

    print(f"  Expanded {num_perms} permutations.")

    # Remove zero coefficients
    poly = {S: c for S, c in poly.items() if c != 0}

    # ===== DEGREE ANALYSIS =====
    print(f"\n{'=' * 70}")
    print("DEGREE ANALYSIS")
    print(f"{'=' * 70}")

    degree_counts = Counter()
    degree_coeffs = defaultdict(list)
    for S, c in poly.items():
        d = len(S)
        degree_counts[d] += 1
        degree_coeffs[d].append((S, c))

    max_deg = max(len(S) for S in poly.keys()) if poly else 0
    min_deg = min(len(S) for S in poly.keys()) if poly else 0

    print(f"\n  Maximum degree: {max_deg}")
    print(f"  Minimum degree: {min_deg}")
    print(f"\n  Nonzero coefficient counts by degree:")
    for d in range(max_deg + 1):
        possible = len(list(combinations(range(num_arcs), d)))
        actual = degree_counts.get(d, 0)
        print(f"    Degree {d}: {actual:5d} nonzero out of C({num_arcs},{d}) = {possible} possible")

    # KEY RESULT
    print(f"\n  {'*' * 50}")
    if max_deg == n - 2:
        print(f"  *** DEGREE DROP CONFIRMED! Max degree = {max_deg} = n-2 ***")
        print(f"  *** All degree-{n-1} monomials cancel! ***")
    elif max_deg == n - 1:
        print(f"  *** NO DEGREE DROP: Max degree = {max_deg} = n-1 ***")
    else:
        print(f"  *** UNEXPECTED: Max degree = {max_deg} ***")
    print(f"  {'*' * 50}")

    # ===== PRINT MONOMIALS BY DEGREE =====
    print(f"\n{'=' * 70}")
    print("MONOMIAL DETAILS BY DEGREE")
    print(f"{'=' * 70}")

    for d in range(max_deg + 1):
        terms = degree_coeffs.get(d, [])
        if not terms:
            print(f"\n  Degree {d}: NO nonzero monomials")
            continue

        print(f"\n  Degree {d}: {len(terms)} nonzero monomials")
        coeff_vals = [c for _, c in terms]
        print(f"    Coefficient range: [{min(coeff_vals)}, {max(coeff_vals)}]")
        print(f"    Sum of coefficients: {sum(coeff_vals)}")

        if len(terms) <= 30:
            for S, c in sorted(terms, key=lambda x: sorted(x[0])):
                arc_names = [f"x{arcs[i][0]}{arcs[i][1]}" for i in sorted(S)]
                label = "*".join(arc_names) if arc_names else "1"
                print(f"      {c:+d} * {label}")

    # ===== DEGREE-5 CANCELLATION ANALYSIS =====
    print(f"\n{'=' * 70}")
    print("QUINTIC (DEGREE-5) CANCELLATION ANALYSIS")
    print(f"{'=' * 70}")

    # Count how many degree-5 monomials appear BEFORE cancellation
    # i.e., how many distinct 5-arc subsets appear in at least one path
    pre_cancel_quintic = Counter()
    for perm in permutations(range(n)):
        factors = []
        for step in range(n - 1):
            i, j = perm[step], perm[step + 1]
            if i < j:
                factors.append((arc_idx[(i, j)], True))
            else:
                factors.append((arc_idx[(j, i)], False))

        # The degree-5 term (all x chosen) has coefficient (-1)^(number of False factors)
        arc_set = set()
        sign = 1
        for arc_id, positive in factors:
            arc_set.add(arc_id)
            if not positive:
                sign *= -1

        # Only contributes to degree-5 if all 5 arcs are distinct
        if len(arc_set) == 5:
            pre_cancel_quintic[frozenset(arc_set)] += sign

    total_quintic_pre = len(pre_cancel_quintic)
    surviving_quintic = sum(1 for v in pre_cancel_quintic.values() if v != 0)
    cancelled_quintic = total_quintic_pre - surviving_quintic

    print(f"\n  Quintic monomials appearing (before cancellation): {total_quintic_pre}")
    print(f"  Quintic monomials surviving (nonzero after sum): {surviving_quintic}")
    print(f"  Quintic monomials that cancelled to zero: {cancelled_quintic}")
    print(f"  Total possible quintic monomials C(15,5) = {len(list(combinations(range(15), 5)))}")

    if surviving_quintic == 0:
        print(f"\n  *** ALL quintic monomials cancel! ***")
        # Show some examples of cancellation
        print(f"\n  Examples of quintic cancellation (first 10):")
        count = 0
        for S in sorted(pre_cancel_quintic.keys(), key=lambda s: sorted(s)):
            if pre_cancel_quintic[S] == 0 and count < 10:
                arc_names = [f"x{arcs[i][0]}{arcs[i][1]}" for i in sorted(S)]
                # Count positive and negative contributions
                pos = 0
                neg = 0
                for perm in permutations(range(n)):
                    factors = []
                    for step in range(n - 1):
                        vi, vj = perm[step], perm[step + 1]
                        if vi < vj:
                            factors.append((arc_idx[(vi, vj)], True))
                        else:
                            factors.append((arc_idx[(vj, vi)], False))
                    aset = set()
                    sgn = 1
                    for aid, p in factors:
                        aset.add(aid)
                        if not p:
                            sgn *= -1
                    if len(aset) == 5 and frozenset(aset) == S:
                        if sgn > 0:
                            pos += sgn
                        else:
                            neg += sgn
                print(f"      {'*'.join(arc_names)}: +{pos} and {neg} => net {pos+neg}")
                count += 1
    else:
        print(f"\n  Surviving quintic monomials:")
        for S in sorted([s for s in pre_cancel_quintic if pre_cancel_quintic[s] != 0],
                        key=lambda s: sorted(s)):
            c = pre_cancel_quintic[S]
            arc_names = [f"x{arcs[i][0]}{arcs[i][1]}" for i in sorted(S)]
            print(f"      {c:+d} * {'*'.join(arc_names)}")

    # ===== DEGREE-3 ANALYSIS: TRIANGLE VS NON-TRIANGLE =====
    print(f"\n{'=' * 70}")
    print("DEGREE-3 ANALYSIS: TRIANGLE VS NON-TRIANGLE TRIPLES")
    print(f"{'=' * 70}")

    # A triangle triple is {(a,b), (b,c), (a,c)} for some {a,b,c}
    triangle_triples = set()
    for triple in combinations(range(n), 3):
        a, b, c = triple
        arc_set = frozenset([arc_idx[(a,b)], arc_idx[(a,c)], arc_idx[(b,c)]])
        triangle_triples.add(arc_set)

    print(f"\n  Total vertex triples: C({n},3) = {len(list(combinations(range(n), 3)))}")
    print(f"  Triangle arc-triples: {len(triangle_triples)}")
    print(f"  Total arc-triples C(15,3): {len(list(combinations(range(15), 3)))}")
    print(f"  Non-triangle arc-triples: {len(list(combinations(range(15), 3))) - len(triangle_triples)}")

    cubic_terms = degree_coeffs.get(3, [])
    cubic_dict = {S: c for S, c in cubic_terms}

    triangle_nonzero = 0
    triangle_zero = 0
    nontriangle_nonzero = 0
    nontriangle_zero = 0

    for triple_arcs in combinations(range(num_arcs), 3):
        S = frozenset(triple_arcs)
        c = cubic_dict.get(S, 0)
        is_tri = S in triangle_triples
        if is_tri:
            if c != 0:
                triangle_nonzero += 1
            else:
                triangle_zero += 1
        else:
            if c != 0:
                nontriangle_nonzero += 1
            else:
                nontriangle_zero += 1

    print(f"\n  Triangle triples with nonzero cubic coeff: {triangle_nonzero}")
    print(f"  Triangle triples with zero cubic coeff: {triangle_zero}")
    print(f"  Non-triangle triples with nonzero cubic coeff: {nontriangle_nonzero}")
    print(f"  Non-triangle triples with zero cubic coeff: {nontriangle_zero}")

    if triangle_nonzero == 0:
        print(f"\n  *** ALL triangle triples have zero cubic coefficient (as expected) ***")

    # Show the non-triangle cubic monomials
    if nontriangle_nonzero > 0 and nontriangle_nonzero <= 50:
        print(f"\n  Non-triangle cubic monomials with nonzero coefficient:")
        for triple_arcs in combinations(range(num_arcs), 3):
            S = frozenset(triple_arcs)
            if S not in triangle_triples:
                c = cubic_dict.get(S, 0)
                if c != 0:
                    arc_names = [f"x{arcs[i][0]}{arcs[i][1]}" for i in sorted(S)]
                    # Describe the arc triple
                    verts = set()
                    for i in sorted(S):
                        verts.add(arcs[i][0])
                        verts.add(arcs[i][1])
                    print(f"      {c:+d} * {'*'.join(arc_names)}  (vertices: {sorted(verts)})")

    # ===== DEGREE-4 ANALYSIS =====
    print(f"\n{'=' * 70}")
    print("DEGREE-4 ANALYSIS")
    print(f"{'=' * 70}")

    quartic_terms = degree_coeffs.get(4, [])
    print(f"\n  Number of nonzero quartic monomials: {len(quartic_terms)}")
    print(f"  Total possible quartic monomials C(15,4) = {len(list(combinations(range(15), 4)))}")

    if quartic_terms:
        coeff_vals = [c for _, c in quartic_terms]
        print(f"  Coefficient range: [{min(coeff_vals)}, {max(coeff_vals)}]")
        print(f"  Sum of all quartic coefficients: {sum(coeff_vals)}")

        # Classify by number of vertices spanned
        vertex_span_counts = Counter()
        for S, c in quartic_terms:
            verts = set()
            for i in S:
                verts.add(arcs[i][0])
                verts.add(arcs[i][1])
            vertex_span_counts[len(verts)] += 1

        print(f"\n  Quartic monomials by number of vertices spanned:")
        for v_count in sorted(vertex_span_counts.keys()):
            print(f"    {v_count} vertices: {vertex_span_counts[v_count]} monomials")

        if len(quartic_terms) <= 40:
            print(f"\n  All quartic monomials:")
            for S, c in sorted(quartic_terms, key=lambda x: sorted(x[0])):
                arc_names = [f"x{arcs[i][0]}{arcs[i][1]}" for i in sorted(S)]
                verts = set()
                for i in S:
                    verts.add(arcs[i][0])
                    verts.add(arcs[i][1])
                print(f"      {c:+d} * {'*'.join(arc_names)}  (vertices: {sorted(verts)})")

    # ===== VERIFICATION: evaluate polynomial at specific tournaments =====
    print(f"\n{'=' * 70}")
    print("VERIFICATION: CHECK POLYNOMIAL AGAINST BRUTE-FORCE H(T)")
    print(f"{'=' * 70}")

    import random
    random.seed(42)

    # Generate a few random tournaments and verify
    for trial in range(5):
        # Random tournament
        x = {}
        for a in arcs:
            x[a] = random.randint(0, 1)

        # Evaluate polynomial
        h_poly = 0
        for S, c in poly.items():
            monomial_val = 1
            for i in S:
                monomial_val *= x[arcs[i]]
            h_poly += c * monomial_val

        # Brute force: count Hamiltonian paths
        h_brute = 0
        for perm in permutations(range(n)):
            valid = True
            for step in range(n - 1):
                i, j = perm[step], perm[step + 1]
                if i < j:
                    if x[(i, j)] != 1:
                        valid = False
                        break
                else:
                    if x[(j, i)] != 0:
                        valid = False
                        break
            if valid:
                h_brute += 1

        match = "OK" if h_poly == h_brute else "MISMATCH!"
        print(f"  Trial {trial}: H(poly)={h_poly}, H(brute)={h_brute}  [{match}]")

    # ===== COMPARISON TABLE =====
    print(f"\n{'=' * 70}")
    print("COMPARISON TABLE: DEGREE DROP ACROSS n")
    print(f"{'=' * 70}")

    print(f"\n  {'n':>3}  {'n-1':>4}  {'max_deg':>8}  {'drop':>5}  {'status':>12}")
    print(f"  {'---':>3}  {'----':>4}  {'--------':>8}  {'-----':>5}  {'------------':>12}")
    print(f"  {'4':>3}  {'3':>4}  {'2':>8}  {'1':>5}  {'CONFIRMED':>12}")
    print(f"  {'5':>3}  {'4':>4}  {'4':>8}  {'0':>5}  {'NO DROP':>12}")
    print(f"  {'6':>3}  {'5':>4}  {max_deg:>8}  {n-1-max_deg:>5}  {'CONFIRMED' if max_deg < n-1 else 'NO DROP':>12}")

    # ===== CONSTANT TERM CHECK =====
    print(f"\n{'=' * 70}")
    print("CONSTANT TERM AND LINEAR ANALYSIS")
    print(f"{'=' * 70}")

    const = poly.get(frozenset(), 0)
    print(f"\n  Constant term (H when all arcs point upward): {const}")

    linear_terms = degree_coeffs.get(1, [])
    print(f"  Number of nonzero linear monomials: {len(linear_terms)}")
    if linear_terms:
        coeff_vals_lin = [c for _, c in linear_terms]
        print(f"  Linear coefficient range: [{min(coeff_vals_lin)}, {max(coeff_vals_lin)}]")
        print(f"  Sum of linear coefficients: {sum(coeff_vals_lin)}")
        if len(linear_terms) <= 20:
            print(f"  Linear monomials:")
            for S, c in sorted(linear_terms, key=lambda x: sorted(x[0])):
                arc_names = [f"x{arcs[i][0]}{arcs[i][1]}" for i in sorted(S)]
                print(f"    {c:+d} * {'*'.join(arc_names)}")

    # ===== QUADRATIC ANALYSIS =====
    print(f"\n{'=' * 70}")
    print("QUADRATIC ANALYSIS")
    print(f"{'=' * 70}")

    quad_terms = degree_coeffs.get(2, [])
    print(f"\n  Number of nonzero quadratic monomials: {len(quad_terms)}")
    print(f"  Total possible C(15,2) = {len(list(combinations(range(15), 2)))}")
    if quad_terms:
        coeff_vals_q = [c for _, c in quad_terms]
        print(f"  Coefficient range: [{min(coeff_vals_q)}, {max(coeff_vals_q)}]")
        print(f"  Sum of quadratic coefficients: {sum(coeff_vals_q)}")

    # ===== FINAL SUMMARY =====
    print(f"\n{'=' * 70}")
    print("FINAL SUMMARY")
    print(f"{'=' * 70}")

    print(f"\n  H(T) at n=6 as multilinear polynomial in 15 arc variables:")
    print(f"  Maximum degree: {max_deg}")
    if max_deg == n - 2:
        print(f"  DEGREE DROP: from n-1={n-1} to n-2={n-2}")
        print(f"  All degree-{n-1} terms cancel despite 720 permutations contributing them.")
        print(f"  Pattern: n=4 drops, n=5 does NOT, n=6 drops?")
        print(f"  Conjecture: EVEN n has degree drop, ODD n does not?")
    elif max_deg == n - 1:
        print(f"  NO DEGREE DROP: max degree = n-1 = {n-1}")
        print(f"  Degree drop at n=4 was a SMALL CASE PHENOMENON")
    else:
        print(f"  UNEXPECTED degree = {max_deg}")

    total_monomials = sum(degree_counts.values())
    print(f"\n  Total nonzero monomials: {total_monomials}")
    for d in range(max_deg + 1):
        print(f"    Degree {d}: {degree_counts.get(d, 0)}")

    print(f"\n{'=' * 70}")
    print("DONE")
    print(f"{'=' * 70}")

if __name__ == '__main__':
    main()
