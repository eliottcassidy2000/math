#!/usr/bin/env python3
"""
PFAFFIAN APPROACH to M[a,b] symmetry.

Key idea: For skew-symmetric matrices, Pf(A)^2 = det(A), and Pfaffians
are sums over perfect matchings. Can M[a,b] be expressed via a Pfaffian
or ratio of Pfaffians?

Motivation:
- M[a,b] has "skew-symmetric" flavor (s_ij = -s_ji)
- M[a,b] = M[b,a] is a SYMMETRY property
- Pfaffians of skew-symmetric matrices naturally produce symmetric expressions

At n=4: M[0,1] is degree 2 in arc weights, same as Pf of a 4x4 matrix.
At n=5: M[0,1] is degree 3, but Pf only defined for even-dim matrices.
  However, we could take Pf of a 6x6 matrix (bordering with a,b).

ALSO: Check if M[a,b] relates to the ADJUGATE of the skew-symmetric part S.
For skew-symmetric S: adj(S) = Pf(S) * S^{-1} * Pf(S)... no, that's not right.
For 2m x 2m skew-symmetric S: adj(S)_{ij} = Pf(S_{ij}) where S_{ij} is S
with rows/cols i,j removed, times a sign.

Actually: (adj S)_{ij} = (-1)^{i+j} Pf(S with rows i,j and cols i,j removed) * ...
No, for skew-symmetric: if we remove row i and col j (i != j), the result
is NOT skew-symmetric. But removing BOTH row/col i AND row/col j gives
a skew-symmetric (n-2) x (n-2) matrix, whose Pfaffian is well-defined.

KEY FORMULA: For n x n skew-symmetric A (n even):
  Pf(A) * (A^{-1})_{ij} = (-1)^{i+j+1+[i<j]} * Pf(A_{ij})
where A_{ij} is A with rows/cols i,j removed.

So Pf(A_{ab}) could relate to M[a,b]!

kind-pasteur-2026-03-06-S23b
"""
from itertools import permutations, combinations
from sympy import symbols, expand, Poly, Matrix, sqrt, Rational
from collections import defaultdict

def setup(n):
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')
    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]
    def t(i, j):
        if i == j: return 0
        return r + s(i, j)
    return r, sv, s, t

def pfaffian(M):
    """Compute Pfaffian of a skew-symmetric matrix."""
    n = M.shape[0]
    if n % 2 == 1:
        return 0
    if n == 0:
        return 1
    if n == 2:
        return M[0, 1]
    # Expansion along first row
    result = 0
    for j in range(1, n):
        # Remove rows/cols 0 and j
        rows_cols = [i for i in range(n) if i != 0 and i != j]
        submatrix = M[rows_cols, :][:, rows_cols]
        sign = (-1)**(j - 1)
        result += sign * M[0, j] * pfaffian(submatrix)
    return expand(result)

def transfer_M(t_fn, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)
        S_set = set(S) | {a}
        R_set = set(R) | {b}

        ea = 0
        for p in permutations(sorted(S_set)):
            if p[-1] != a: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t_fn(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1

        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t_fn(p[i], p[i+1])
            bb += prod
        if len(R_set) == 1: bb = 1

        result += sign * ea * bb
    return expand(result)

print("=" * 70)
print("PFAFFIAN APPROACH TO M[a,b]")
print("=" * 70)

# ============================================================
# Part 1: Build the skew-symmetric matrix S and compute Pfaffians
# ============================================================
for n in [4, 5, 6]:
    r, sv, s, t = setup(n)
    a, b = 0, 1

    # Build the n x n skew-symmetric matrix S
    S_mat = Matrix(n, n, lambda i, j: s(i, j))

    # Build the full tournament matrix T = rJ' + S
    T_mat = Matrix(n, n, lambda i, j: t(i, j) if i != j else 0)

    M_ab = transfer_M(t, n, a, b)

    print(f"\nn={n}:")
    print(f"  M[0,1] = {M_ab}")

    # For even n: compute Pf(S) and Pf(S with rows/cols a,b removed)
    if n % 2 == 0:
        pf_S = pfaffian(S_mat)
        print(f"  Pf(S) = {pf_S}")

        # S with rows/cols a,b removed
        remaining = [i for i in range(n) if i != a and i != b]
        S_ab = S_mat[remaining, :][:, remaining]
        pf_S_ab = pfaffian(S_ab)
        print(f"  Pf(S_{{0,1}}) = {pf_S_ab}")

        # Check: is M related to Pf(S_ab)?
        diff = expand(M_ab - pf_S_ab)
        print(f"  M - Pf(S_ab) = {diff}")

        # Also try Pf(T) for even n where T is skew-symmetric at r=0
        # T is NOT skew-symmetric (it has r on off-diagonal), but T - rJ' = S is.
        # Actually T_ij + T_ji = 2r for i!=j, so T is NOT skew-symmetric.
        # But T - rJ' = S IS skew-symmetric.

        # Try: Pf(S) with modified entries involving r
        # E.g., Pf of S + r*something?

        # Build S + r*(something skew-symmetric)
        # The simplest skew-symmetric matrix with constant entries is... there's none
        # that's "natural" for tournaments.
        # Actually: J' - J'^T = 0 (J' is symmetric for complete graph).
        # So adding r to both s_ij and s_ji changes nothing (they cancel in skew part).

        # KEY INSIGHT: The tournament matrix T is NOT skew-symmetric.
        # T_ij = r + s_ij, T_ji = r - s_ij.
        # T_ij + T_ji = 2r != 0.
        # So Pfaffian of T doesn't apply directly.

        # But we can decompose: T = r*J' + S where J'_ij = 1 for i!=j.
        # And J' = J - I where J is the all-ones matrix.

    # For odd n: Pf is 0 for the full matrix, but we can border it
    if n % 2 == 1:
        # Build (n+1) x (n+1) bordered matrix
        # Border with row/col for a "phantom" vertex
        # Pf of bordered matrix relates to Hamiltonian paths?
        pass

    # ============================================================
    # NEW IDEA: M[a,b] at r=0 (pure s-part)
    # ============================================================
    M_at_0 = expand(M_ab.subs(r, 0))
    print(f"  M[0,1](r=0) = {M_at_0}")

    if n % 2 == 0:
        # For even n: is M(0) related to Pf(S_ab)?
        diff_at_0 = expand(M_at_0 - pf_S_ab)
        print(f"  M(0) - Pf(S_ab) = {diff_at_0}")

    # ============================================================
    # Check: does M[a,b] at r=0 equal some determinant?
    # For n=4: M(0) is degree 2 in s. Could be det of 2x2.
    # For n=5: M(0) is degree 3. Could be Pf of 6x6... no.
    # For n=6: M(0) is degree 4. Could be det of 4x4 or Pf of 4x4.
    # ============================================================

    if n == 4:
        # M(0) = s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23
        # Try: det of [[A, B], [C, D]] = AD - BC
        # Need degree 1 entries: linear in s_ij
        #
        # Actually, Pf of 4x4 skew-symmetric with s-entries would have degree 2.
        # Pf([[0,a,b,c],[-a,0,d,e],[-b,-d,0,f],[-c,-e,-f,0]])
        # = a*f - b*e + c*d
        # This has 3 terms. M(0) has 6 terms. So Pf of 4x4 skew won't work.
        #
        # What about Pf of a 4x4 matrix built from s and extra parameters?
        # Or: determinant of a 2x2 matrix?
        # det [[A,B],[C,D]] = AD - BC has 2 terms (each of degree 1).
        # But M(0) has 6 terms of degree 2. Can we find 2x2 with degree-2 entries?
        # det [[Q1, Q2],[Q3, Q4]] where Qi are degree 1 in s
        # = Q1*Q4 - Q2*Q3, which has at most 2*k^2 terms if Qi has k terms.
        # With k=2: at most 8 terms. Possible!
        #
        # Let me try: Q1 = s02 + s12, Q4 = s23 + s13, Q2 = s03, Q3 = -s23
        # Q1*Q4 = (s02+s12)(s23+s13) = s02*s23 + s02*s13 + s12*s23 + s12*s13
        # Q2*Q3 = s03*(-s23) = -s03*s23
        # det = s02*s23 + s02*s13 + s12*s23 + s12*s13 + s03*s23
        # M(0) = s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23
        # Not matching.
        #
        # Let me try systematic search for Q1,Q2,Q3,Q4 linear in s
        print("\n  --- Searching for 2x2 determinant representation at n=4 ---")
        s_vars_list = list(sv.values())
        # Q_i = sum c_ij * s_j, with c_ij in {-1, 0, 1}
        # Too many combos (7^6 ~ 100k for each Q). Let me parameterize.
        # Actually: use symbolic coefficients and solve.
        from sympy import symbols as sym
        cs = {}
        for k in range(4):
            for var in s_vars_list:
                cs[(k, var)] = sym(f'c{k}_{var.name}')

        Q = [sum(cs[(k, v)] * v for v in s_vars_list) for k in range(4)]
        det_Q = expand(Q[0]*Q[3] - Q[1]*Q[2])

        # Match coefficients with M(0)
        # M(0) = s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23
        from sympy import Poly as SPoly, solve
        all_s = s_vars_list
        target = M_at_0

        # This is a quadratic system in the c's — hard to solve in general.
        # Let me just check a few natural candidates.

        # Candidate: inspired by the structure a=0, b=1, U={2,3}
        # Try: [[t(2,0), t(3,0)], [t(2,1), t(3,1)]] at r=0
        # = [[s(2,0), s(3,0)], [s(2,1), s(3,1)]]
        # = [[-s02, -s03], [s12, s13]]  (using s(i,j) = -s(j,i) for i>j)
        # Wait: s(2,0) = -s02 (since 2>0), s(3,0) = -s03
        # s(2,1) = s12... no, s(2,1) = -s12 (since 2>1)
        # s(3,1) = -s13

        cand_mat = Matrix([
            [s(2,0), s(3,0)],
            [s(2,1), s(3,1)]
        ])
        det_cand = expand(cand_mat.det())
        print(f"    det [[s(2,0),s(3,0)],[s(2,1),s(3,1)]] = {det_cand}")
        print(f"    = {expand(det_cand)}")
        print(f"    M(0) = {M_at_0}")
        print(f"    Match: {expand(M_at_0 - det_cand) == 0}")

        # Try with t instead of s
        cand_mat_t = Matrix([
            [t(2,0), t(3,0)],
            [t(2,1), t(3,1)]
        ])
        det_cand_t = expand(cand_mat_t.det())
        print(f"    det [[t(2,0),t(3,0)],[t(2,1),t(3,1)]] = {det_cand_t}")
        print(f"    Match with M: {expand(M_ab - det_cand_t) == 0}")

        # Another candidate: sum over u in U of t(u,a)*t(b,...) type things
        # Or: the "transfer determinant"
        # For n=4, U={2,3}. Define L_{ij} for i,j in U:
        # L_{22} = sum of paths a -> 2 -> b = t(a,2)*t(2,b) = t(0,2)*t(2,1)
        # L_{23} = path a -> 2 -> 3 -> ... wait, this overcounts.

        # Actually: L_{ij} = t(i, j) for i,j in U. Then det(L_U) involves only
        # internal edges. But M involves edges to a and b too.

        # Try: L_{ij} = t(a,i)*t(j,b) for i,j in U (rank 1, det=0). No good.

        # Try: L_{ij} = t(a,i)*delta_{ij} + t(j,b) ... getting speculative.

        # Let me try ALL 2x2 matrices with entries from {t(i,j): 0<=i,j<=3, i!=j}
        # and check which determinant matches M_ab.
        print("\n  --- Exhaustive 2x2 det search ---")
        all_t_entries = [(i,j) for i in range(4) for j in range(4) if i != j]
        found = False
        for (a1,b1) in all_t_entries:
            for (a2,b2) in all_t_entries:
                for (a3,b3) in all_t_entries:
                    for (a4,b4) in all_t_entries:
                        d = expand(t(a1,b1)*t(a4,b4) - t(a2,b2)*t(a3,b3))
                        if expand(d - M_ab) == 0:
                            print(f"    FOUND: det [[t({a1},{b1}), t({a2},{b2})], "
                                  f"[t({a3},{b3}), t({a4},{b4})]] = M[0,1]")
                            found = True
                            break
                    if found: break
                if found: break
            if found: break
        if not found:
            print("    No single 2x2 det of t-entries matches M[0,1]")

        # Try: det of sum/difference combinations
        # det [[t(i,j)+t(k,l), ...], [..., ...]]
        # This is 12^4 * many combos... too many.

        # Instead: try SUMS of determinants
        # M = sum of det(...) type expressions
        # At n=4, M = sum over S of (-1)^|S| * E_a * B_b
        # Each E_a * B_b is a product of at most 2 t-values
        # So M = sum of products of t-values, not naturally a single determinant.

print("\n" + "=" * 70)
print("Part 2: EVEN POWERS via skew-symmetry structure")
print("=" * 70)

# The key observation: t(i,j) = r + s(i,j) and t(j,i) = r - s(i,j).
# So t(i,j) * t(j,i) = r^2 - s(i,j)^2 (EVEN in r!)
#
# If every arc in M[a,b] appeared paired with its reverse, we'd be done.
# But in a Hamiltonian PATH, arcs are NOT paired.
#
# However, the INCLUSION-EXCLUSION might create the pairing!
# When we sum over all S with signs (-1)^|S|, the arc (i,j) in one
# cover might pair with (j,i) in another cover (with opposite sign).
#
# Let's check: for n=4, list all arcs used and see if they pair up.

n = 4
r, sv, s, t = setup(n)
a, b = 0, 1
U = [v for v in range(n) if v != a and v != b]

print(f"\nn=4: Arc analysis")
arc_contributions = defaultdict(list)

for mask in range(1 << len(U)):
    S = tuple(sorted([U[i] for i in range(len(U)) if mask & (1 << i)]))
    R = tuple(sorted([U[i] for i in range(len(U)) if not (mask & (1 << i))]))
    sign = (-1)**len(S)

    S_set = set(S) | {a}
    R_set = set(R) | {b}

    for p1 in permutations(sorted(S_set)):
        if p1[-1] != a: continue
        for p2 in permutations(sorted(R_set)):
            if p2[0] != b: continue
            arcs = []
            for i in range(len(p1)-1):
                arcs.append((p1[i], p1[i+1]))
            for i in range(len(p2)-1):
                arcs.append((p2[i], p2[i+1]))

            arc_tuple = tuple(sorted(arcs))
            arc_contributions[arc_tuple].append((sign, S, p1, p2))

for arc_set, covers in sorted(arc_contributions.items()):
    total_sign = sum(sign for sign, _, _, _ in covers)
    reversed_arcs = tuple(sorted([(j,i) for (i,j) in arc_set]))
    has_reverse = reversed_arcs in arc_contributions
    rev_sign = sum(sign for sign, _, _, _ in arc_contributions.get(reversed_arcs, []))
    print(f"  Arcs {arc_set}: total_sign={total_sign}, "
          f"reverse={reversed_arcs}: reverse_sign={rev_sign}")

# Check: for each arc set, does its reverse arc set appear with the same total sign?
print("\n  Arc-reversal pairing check:")
checked = set()
for arc_set in arc_contributions:
    if arc_set in checked: continue
    reversed_arcs = tuple(sorted([(j,i) for (i,j) in arc_set]))
    sign_fwd = sum(sign for sign, _, _, _ in arc_contributions[arc_set])
    sign_rev = sum(sign for sign, _, _, _ in arc_contributions.get(reversed_arcs, []))
    checked.add(arc_set)
    checked.add(reversed_arcs)
    if arc_set == reversed_arcs:
        print(f"  SELF-REVERSE: {arc_set}, sign={sign_fwd}")
    else:
        print(f"  PAIR: {arc_set} (sign {sign_fwd}) <-> {reversed_arcs} (sign {sign_rev})")
        # If they pair, the product t(i,j)*t(j,i) = r^2 - s^2 appears
        # Actually no — the arcs in a cover form a SEQUENCE, not a single product.
        # The contribution is sign * prod_{(i,j) in arcs} t(i,j).
        # If arcs_fwd and arcs_rev have the same product structure, great.
        # But reversing all arcs changes each t(i,j) -> t(j,i) = r - s(i,j).
        # prod t(j,i) = prod (r - s(i,j)) which is M evaluated at -s (NOT -r).
        # Hmm, that's different.

print("\n" + "=" * 70)
print("Part 3: Substitution r -> -r at the arc level")
print("=" * 70)

# Under r -> -r: t(i,j) = r + s(i,j) -> -r + s(i,j) = -(r - s(i,j)) = -t(j,i)
# So t(i,j) at -r = -t(j,i) at r.
# Product of m arcs at -r = (-1)^m * product of REVERSED arcs at r.
# m = n-2.

# For a cover (P1: ...->a, P2: b->...):
# P1 reversed becomes a->... (a path STARTING at a through same vertices)
# P2 reversed becomes ...->b (a path ENDING at b through same vertices)
# So reversed cover has: path a->... through S+a starting at a,
#                         path ...->b through R+b ending at b.
# This is B_a(S+a) * E_b(R+b) — the SWAPPED version!

# So: contribution of (S, P1, P2) at -r
# = (-1)^m * sign * B_a(S+a at r) * E_b(R+b at r)
# And summing over all S:
# M[a,b](-r) = (-1)^m * sum_S (-1)^|S| B_a(S+a) E_b(R+b)

# Now relabel S <-> R (and |S| -> |R| = m-|S|, sign flips by (-1)^m):
# = (-1)^m * (-1)^m * sum_R (-1)^|R| E_b(R+b) B_a(S+a)
# Wait, let me be more careful.

# sum_S (-1)^|S| B_a(S+a) E_b(R+b)   where R = U\S
# Let S' = R = U\S. Then |S'| = m - |S|, (-1)^|S| = (-1)^m * (-1)^|S'|.
# sum_{S'} (-1)^m * (-1)^|S'| B_a((U\S')+a) E_b(S'+b)
# = (-1)^m * sum_{S'} (-1)^|S'| E_b(S'+b) B_a((U\S')+a)
# = (-1)^m * M[b,a](r)

# So M[a,b](-r) = (-1)^m * (-1)^m * M[b,a](r) = (-1)^{2m} M[b,a](r) = M[b,a](r).
# This confirms M[a,b](-r) = M[b,a](r). Already known.

# The question is: WHY does M[a,b](r) = M[b,a](r)?

# APPROACH: Direct comparison.
# M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
# M[b,a] = sum_S (-1)^|S| E_b(S+b) B_a(R+a)
# These are DIFFERENT sums. For them to be equal, we need a bijection
# between the terms.

# At n=4: let's compare term by term.
print(f"\nn=4: Term-by-term comparison M[0,1] vs M[1,0]")
for M_type, aa, bb in [("M[0,1]", 0, 1), ("M[1,0]", 1, 0)]:
    print(f"\n  {M_type}:")
    for mask in range(1 << len(U)):
        S = tuple(sorted([U[i] for i in range(len(U)) if mask & (1 << i)]))
        R = tuple(sorted([U[i] for i in range(len(U)) if not (mask & (1 << i))]))
        sign = (-1)**len(S)
        S_set = set(S) | {aa}
        R_set = set(R) | {bb}

        ea = 0
        for p in permutations(sorted(S_set)):
            if p[-1] != aa: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1

        bbb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != bb: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t(p[i], p[i+1])
            bbb += prod
        if len(R_set) == 1: bbb = 1

        contrib = expand(sign * ea * bbb)
        print(f"    S={set(S)}: ({sign:+d}) * E_{aa}({sorted(S_set)}) * B_{bb}({sorted(R_set)}) = {contrib}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
