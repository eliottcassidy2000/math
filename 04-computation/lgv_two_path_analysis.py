#!/usr/bin/env python3
"""
LGV / TWO-PATH-COVER ANALYSIS for M[a,b]

M[a,b] = sum_{S subset U} (-1)^|S| E_a(S+{a}) B_b(R+{b})

This is an INCLUSION-EXCLUSION formula. Each term is a product of two
Hamiltonian path sums through complementary vertex sets.

KEY QUESTION: Can we interpret this as a DETERMINANT via Lindstrom-Gessel-Viennot?

LGV says: for a DAG, det of path-weight matrix = signed sum of non-crossing
path systems. But our graph is complete (not a DAG), and we have
crossing paths by construction.

ALTERNATIVE: The Inclusion-Exclusion structure itself might be a determinant.
Let's test: define matrix L where L_{u,v} depends on the internal vertices u,v in U.
Does M[a,b] = det(L) for some natural L?

Also: analyze the PAIRING structure of the r^1 coefficient at the
individual 2-path-cover level.

kind-pasteur-2026-03-06-S23b
"""
from itertools import permutations
from sympy import symbols, expand, Poly, Matrix, det, zeros
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
print("LGV / TWO-PATH-COVER ANALYSIS")
print("=" * 70)

# ============================================================
# Part 1: Test if M[a,b] = det(I + something)
# ============================================================
# At n=4, U={2,3}, M is degree 2.
# Try: M[0,1] = det([[f(2,2), f(2,3)], [f(3,2), f(3,3)]]) for some f?

for n in [4, 5]:
    r, sv, s, t = setup(n)
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m = len(U)

    M_ab = transfer_M(t, n, a, b)
    print(f"\nn={n}: M[0,1] = {M_ab}")

    # For n=4 (m=2): try det of 2x2 matrix with entries f(u,v)
    # for u,v in U = {2,3}
    if n == 4:
        # What 2x2 matrix has determinant = M[0,1]?
        # M = ad - bc for entries a,b,c,d
        # M = r^2 + r^0 terms (since r^1 = 0)
        # At r=0: M(0) = s(0,2)*s(3,1) + s(0,3)*s(2,1) + ...
        # Let's compute M(0) explicitly
        M_at_0 = M_ab.subs(r, 0)
        print(f"  M(0) = {expand(M_at_0)}")

        # Try: L[u,v] = t(u,a)*t(b,v) - t(u,b)*t(a,v) ... some LGV-like thing
        # Or simpler: L[u,v] = t(a,u) * t(v,b)... (path a->u->...->v->b weights)
        # Hmm, let me think about what matrix L gives det(L) = M

        # Brute force: parameterize L = [[p, q], [s_var, t_var]] with linear entries
        # in arc weights. Too many unknowns.

        # Better: L_{u,v} = "contribution of vertex u being in position ... and v in position ..."
        # For n=4: a 2-path-cover is: path ?->a through {?,a} and path b->? through {b,?}
        # where ? ranges over U = {2,3}.
        # If u in S: path u->a, and 3-u goes to b-path
        # Actually: S={}: both 2,3 go with b. S={2}: 2 with a, 3 with b. etc.
        # E_a({a}) = 1, B_b({b,2,3}) = t(b,2)*t(2,3) + t(b,3)*t(3,2)
        # E_a({2,a}) = t(2,a) = t(2,0), B_b({b,3}) = t(b,3) = t(1,3)
        # E_a({3,a}) = t(3,a) = t(3,0), B_b({b,2}) = t(b,2) = t(1,2)
        # E_a({2,3,a}) = t(2,3)*t(3,a) + t(3,2)*t(2,a) = t(2,3)*t(3,0)+t(3,2)*t(2,0)
        # B_b({b}) = 1

        # M = 1*B_b({1,2,3}) - t(2,0)*t(1,3) - t(3,0)*t(1,2) + (t(2,3)*t(3,0)+t(3,2)*t(2,0))*1
        # = t(1,2)*t(2,3) + t(1,3)*t(3,2) - t(2,0)*t(1,3) - t(3,0)*t(1,2) + t(2,3)*t(3,0) + t(3,2)*t(2,0)

        # Let me check: this is
        # t(1,2)[t(2,3) - t(3,0)] + t(1,3)[t(3,2) - t(2,0)] + t(2,3)*t(3,0) + t(3,2)*t(2,0)
        # Hmm, let me just try det of various natural matrices

        # Try: L = [[t(a,u) + t(v,b), ...]] — nah
        # Let me try: det of the matrix with L_{ij} = delta_{ij} - t(U[i], U[j])
        L1 = Matrix(m, m, lambda i, j: (1 if i==j else 0) - t(U[i], U[j]) if i != j else 1)
        print(f"  det(I - T_UU) = {expand(det(L1))}")

        # Try: det with L_{ij} = t(a, U[i]) * t(U[j], b) for i != j, something on diag
        # This would be rank 1, so det = 0. No good.

        # Try the "path matrix" approach:
        # P[u] = t(a,u) for paths a->u
        # Q[v] = t(v,b) for paths v->b
        # M might relate to det of a matrix involving P and Q
        # Actually at n=4 with just 2 internal vertices:
        # det [[t(0,2)*t(2,1) + something, ...], [...]] — too speculative

        # Let me instead check: is M[a,b] the PERMANENT of something minus
        # the DETERMINANT?
        # per(A) - det(A) = 2*a12*a21 for 2x2
        # det(A) = a11*a22 - a12*a21
        # per(A) = a11*a22 + a12*a21

        pass

    # For any n: check the r^k coefficients
    p_r = Poly(M_ab, r)
    print(f"  r-degree: {p_r.degree()}")
    for k in range(p_r.degree() + 1):
        coeff = expand(p_r.nth(k))
        if coeff == 0:
            print(f"  r^{k}: 0")
        else:
            print(f"  r^{k}: {coeff}")

# ============================================================
# Part 2: The "complementary path" involution
# ============================================================
# For each 2-path-cover (P1: ...->a through S+a, P2: b->... through R+b),
# consider the "complementary" cover obtained by REVERSING the role:
# P1': ...->b through S+b (same vertices, different endpoint)
# P2': a->... through R+a (same vertices, different start)
# This is M[b,a], not M[a,b]. So M[a,b] = M[b,a] iff the signed sums match.

# But we proved M[a,b](-r) = M[b,a](r). So M[a,b](r) = M[b,a](r) iff
# M[a,b](r) = M[a,b](-r), i.e., only even powers of r.

# Let's look at the r^1 coefficient more carefully.
# At r^1: each 2-path-cover contributes e_{m-1}(s_1,...,s_m)
# = sum_i prod_{j!=i} s_j
# = (prod all s_j) * sum_i (1/s_i)   [when no s_i = 0]

# But with the inclusion-exclusion sign (-1)^|S|, the sum vanishes.

# KEY IDEA: Can we pair covers (S, P1, P2) with (S', P1', P2') such that
# their r^1 contributions cancel?
# The pairing should flip the sign (-1)^|S| -> -(-1)^|S|, i.e., |S'| = |S| +/- 1.

print("\n" + "=" * 70)
print("Part 2: r^1 contributions by individual cover")
print("=" * 70)

n = 4
r, sv, s, t = setup(n)
a, b = 0, 1
U = [v for v in range(n) if v != a and v != b]

covers = []
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

            # Product of (r + s(arc))
            prod_expr = 1
            for (u, v) in arcs:
                prod_expr *= (r + s(u, v))
            prod_expr = expand(prod_expr)

            # Extract r^1 coefficient
            p_cover = Poly(prod_expr, r)
            r1 = expand(p_cover.nth(1)) if p_cover.degree() >= 1 else 0
            r0 = expand(p_cover.nth(0))

            covers.append({
                'S': S, 'sign': sign,
                'p1': p1, 'p2': p2,
                'arcs': arcs,
                'r1_unsigned': r1,
                'r1_signed': expand(sign * r1),
                'r0_signed': expand(sign * r0)
            })

print(f"\nn=4: {len(covers)} individual covers")
for c in covers:
    print(f"  S={set(c['S'])}, sign={c['sign']:+d}, "
          f"P1={'->'.join(map(str,c['p1']))}, P2={'->'.join(map(str,c['p2']))}")
    print(f"    arcs={c['arcs']}, r^1(unsigned)={c['r1_unsigned']}, "
          f"r^1(signed)={c['r1_signed']}")

total_r1 = expand(sum(c['r1_signed'] for c in covers))
print(f"\n  Total r^1 = {total_r1}")

# ============================================================
# Part 3: The MATRIX FORMULATION
# ============================================================
# Define the "augmented" transfer matrix T as the nxn matrix with T_{ij} = t(i,j).
# M[a,b] relates to paths in this matrix.
#
# IDEA: The Hamiltonian path from a through V ending at b is the (a,b) entry
# of the "Hamiltonian adjacency matrix" which we can write as:
# H(T)_{a,b} = sum over permutations sigma with sigma(a)=b and sigma is a
#              single cycle containing all vertices, of product of t(v, sigma(v)).
# But M[a,b] is NOT H(T)_{a,b}. M involves the INCLUSION-EXCLUSION of two-path-covers.
#
# Actually, M[a,b] counts something specific:
# It's the (a,b) cofactor of (I - T) or similar? Let's check.

print("\n" + "=" * 70)
print("Part 3: Matrix formulations")
print("=" * 70)

for n in [4, 5]:
    r, sv, s, t = setup(n)

    # Build T matrix (n x n)
    T = Matrix(n, n, lambda i, j: t(i, j) if i != j else 0)
    print(f"\nn={n}: T matrix built")

    # Check: is M[a,b] related to cofactors of (I - T)?
    I_minus_T = Matrix(n, n, lambda i, j: (1 if i==j else 0) - (t(i,j) if i!=j else 0))
    a_idx, b_idx = 0, 1

    # Cofactor C_{a,b} of (I-T)
    # = (-1)^{a+b} * det of (I-T) with row a, col b removed
    minor = I_minus_T.minor_submatrix(a_idx, b_idx)
    cofactor = expand((-1)**(a_idx + b_idx) * det(Matrix(minor)))
    M_ab = transfer_M(t, n, 0, 1)

    print(f"  M[0,1] = {expand(M_ab)}")
    print(f"  Cofactor(I-T, 0, 1) = {cofactor}")
    print(f"  Match: {expand(M_ab - cofactor) == 0}")

    # Try adj(I-T)
    adj_entry = expand((-1)**(a_idx + b_idx) * det(Matrix(I_minus_T.minor_submatrix(b_idx, a_idx))))
    print(f"  adj(I-T)[0,1] = {adj_entry}")
    print(f"  Match adj: {expand(M_ab - adj_entry) == 0}")

    # Try (I+T) instead
    I_plus_T = Matrix(n, n, lambda i, j: (1 if i==j else 0) + (t(i,j) if i!=j else 0))
    cof_plus = expand((-1)**(a_idx + b_idx) * det(Matrix(I_plus_T.minor_submatrix(a_idx, b_idx))))
    print(f"  Cofactor(I+T, 0, 1) = {cof_plus}")
    print(f"  Match I+T: {expand(M_ab - cof_plus) == 0}")

    # Try (-1)^{n-1} * cofactor of (T - I)
    T_minus_I = Matrix(n, n, lambda i, j: (t(i,j) if i!=j else 0) - (1 if i==j else 0))
    cof_TmI = expand((-1)**(a_idx + b_idx) * det(Matrix(T_minus_I.minor_submatrix(a_idx, b_idx))))
    print(f"  Cofactor(T-I, 0, 1) = {cof_TmI}")
    print(f"  Match T-I: {expand(M_ab - cof_TmI) == 0}")
    print(f"  Match -cof(T-I): {expand(M_ab + cof_TmI) == 0}")

    if n <= 5:
        # What IS det(I-T)?
        det_I_T = expand(det(I_minus_T))
        print(f"  det(I-T) = {det_I_T}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
