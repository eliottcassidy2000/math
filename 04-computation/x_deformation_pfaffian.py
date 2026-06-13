#!/usr/bin/env python3
"""
x-DEFORMATION OF THE PFAFFIAN IDENTITY
opus-2026-03-14-S68

THM-174 says: det(I + 2A) = Pf(A - A^T)²

QUESTION: What happens when we replace 2 with a general parameter x?

det(I + xA) = ???

The decomposition I + xA = J + (x-1)A + (A - A^T) doesn't obviously simplify.
But: I + xA = (1-x/2)I + (x/2)(I + 2A) = ... hmm, that's not clean.

Alternative: I + xA with S = A - A^T (skew), A = (J-I+S)/2.
So I + xA = I + x(J-I+S)/2 = (1-x/2)I + (x/2)J + (x/2)S
           = (1-x/2)I + (x/2)(J+S)

When x=2: (0)I + (1)(J+S) = J+S, and det(J+S) = Pf(S)² ✓.

For general x: det((1-x/2)I + (x/2)(J+S))
             = det(J+S) · det((1-x/2)(J+S)^{-1} + (x/2)I)  ... circular.

Let's just COMPUTE det(I+xA) for all tournaments and see what happens.
"""

import numpy as np
from itertools import combinations, permutations
from sympy import symbols, Matrix, det, factor, sqrt, expand, Rational
from sympy import Poly, Symbol

x = Symbol('x')

print("=" * 78)
print("  x-DEFORMATION OF THE PFAFFIAN IDENTITY — det(I + xA)")
print("=" * 78)

# ============================================================================
# PART 1: SYMBOLIC COMPUTATION FOR n=3
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 1: SYMBOLIC det(I + xA) FOR n=3                                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")

def adj_matrix_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

# n=3: 3 edges, 8 tournaments (4 iso classes: transitive + 2 rotations of 3-cycle)
print("n=3 tournaments — det(I + x·A) as polynomial in x:")
print()

n = 3
m = n*(n-1)//2
seen_polys = {}

for bits in range(1 << m):
    A_list = adj_matrix_bits(bits, n)
    A_sym = Matrix(A_list)
    M = Matrix.eye(n) + x * A_sym
    d = det(M)
    d_expanded = expand(d)
    d_factored = factor(d)

    # Also compute at x=2
    d_at_2 = d_expanded.subs(x, 2)

    key = str(d_factored)
    if key not in seen_polys:
        seen_polys[key] = 0
        print(f"  T(bits={bits:03b}): det(I+xA) = {d_factored}")
        print(f"    expanded: {d_expanded}")
        print(f"    at x=2: {d_at_2}")

        # Check if perfect square as polynomial
        p = Poly(d_expanded, x)
        print(f"    coefficients: {p.all_coeffs()}")
        print()
    seen_polys[key] += 1

print(f"  {len(seen_polys)} distinct polynomials from {1<<m} tournaments")

# ============================================================================
# PART 2: SYMBOLIC COMPUTATION FOR n=4
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 2: SYMBOLIC det(I + xA) FOR n=4                                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")

n = 4
m = n*(n-1)//2
seen_polys = {}

for bits in range(1 << m):
    A_list = adj_matrix_bits(bits, n)
    A_sym = Matrix(A_list)
    M = Matrix.eye(n) + x * A_sym
    d = det(M)
    d_expanded = expand(d)
    d_factored = factor(d)

    key = str(d_expanded)
    if key not in seen_polys:
        seen_polys[key] = (d_expanded, d_factored, bits)
    else:
        pass

print(f"n=4: {len(seen_polys)} distinct determinant polynomials from {1<<m} tournaments\n")

for key, (dexp, dfact, bits) in sorted(seen_polys.items(), key=lambda kv: kv[1][0].subs(x, 2)):
    d_at_2 = dexp.subs(x, 2)
    p = Poly(dexp, x)

    # Check: is it a perfect square?
    # Try sqrt
    try:
        sq = sqrt(dexp)
        sq_expanded = expand(sq**2)
        is_square = (expand(sq_expanded - dexp) == 0)
    except:
        is_square = False

    print(f"  bits={bits:06b}: det = {dfact}")
    print(f"    expanded: {dexp}")
    print(f"    at x=2: {d_at_2} = √{d_at_2}² ({'✓ square' if int(d_at_2)**0.5 == int(int(d_at_2)**0.5) else '?'})")
    print()

# ============================================================================
# PART 3: IS det(I+xA) ALWAYS A PERFECT SQUARE IN x?
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 3: IS det(I+xA) A PERFECT SQUARE POLYNOMIAL IN x?                    ║
╚══════════════════════════════════════════════════════════════════════════════╝

At x=2, det(I+2A) = Pf(S)² is always a perfect square (THM-174).
QUESTION: Is det(I+xA) a perfect square as a polynomial in x?
""")

# For n=3: check each polynomial
n = 3
m = n*(n-1)//2
print(f"n={n}:")
seen = set()
for bits in range(1 << m):
    A_list = adj_matrix_bits(bits, n)
    A_sym = Matrix(A_list)
    M = Matrix.eye(n) + x * A_sym
    d = expand(det(M))
    key = str(d)
    if key in seen: continue
    seen.add(key)

    # Try to find sqrt as polynomial
    p = Poly(d, x)
    coeffs = p.all_coeffs()
    deg = p.degree()

    # For a perfect square, degree must be even
    if deg % 2 != 0:
        print(f"  det = {d}: degree {deg} (ODD) → NOT a perfect square polynomial")
    else:
        # Try to compute square root
        # For degree 2: ax²+bx+c = (αx+β)² requires b²=4ac
        if deg == 2:
            a, b, c = coeffs
            disc = b**2 - 4*a*c
            print(f"  det = {d}: degree 2, discriminant = {disc} → {'SQUARE' if disc == 0 else 'NOT square'}")
        elif deg == 0:
            print(f"  det = {d}: degree 0 → {'SQUARE' if d >= 0 else 'NOT'}")

# For n=4
n = 4
m = n*(n-1)//2
print(f"\nn={n}:")
seen = set()
for bits in range(1 << m):
    A_list = adj_matrix_bits(bits, n)
    A_sym = Matrix(A_list)
    M = Matrix.eye(n) + x * A_sym
    d = expand(det(M))
    key = str(d)
    if key in seen: continue
    seen.add(key)

    p = Poly(d, x)
    coeffs = p.all_coeffs()
    deg = p.degree()

    if deg % 2 != 0:
        print(f"  det = {factor(d)}: degree {deg} → NOT a perfect square polynomial")
    elif deg == 4:
        # Check if a⁴x⁴+...+e = (αx²+βx+γ)²
        a4, a3, a2, a1, a0 = coeffs
        # (αx²+βx+γ)² = α²x⁴ + 2αβx³ + (β²+2αγ)x² + 2βγx + γ²
        # Need: a4 = α², a0 = γ² → α = √a4, γ = √a0
        # Then: a3 = 2αβ → β = a3/(2α)
        # Check: a2 = β²+2αγ, a1 = 2βγ
        from sympy import sqrt as ssqrt, Rational as R
        if a4 > 0 and a0 > 0:
            try:
                alpha = ssqrt(a4)
                gamma = ssqrt(a0)
                beta = R(a3, 2*alpha)
                check_a2 = beta**2 + 2*alpha*gamma
                check_a1 = 2*beta*gamma
                is_sq = (check_a2 == a2 and check_a1 == a1)
                if is_sq:
                    root_poly = alpha*x**2 + beta*x + gamma
                    print(f"  det = {factor(d)}: PERFECT SQUARE = ({root_poly})²")
                else:
                    print(f"  det = {factor(d)}: degree 4, NOT a perfect square (a2 mismatch: {check_a2} vs {a2})")
            except:
                print(f"  det = {factor(d)}: degree 4, could not extract square root")
        else:
            print(f"  det = {factor(d)}: degree 4, leading/trailing coeff not positive square")
    elif deg == 2:
        a, b, c = coeffs
        disc = b**2 - 4*a*c
        print(f"  det = {factor(d)}: degree 2, disc = {disc} → {'SQUARE' if disc == 0 else 'NOT square'}")
    else:
        print(f"  det = {factor(d)}: degree {deg}")

# ============================================================================
# PART 4: THE x-DEFORMED H AND PFAFFIAN — NUMERICAL
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 4: det(I+xA) vs H(x) = I(CG, x) — NUMERICAL COMPARISON              ║
╚══════════════════════════════════════════════════════════════════════════════╝

We know H(2) = I(CG, 2) and det(I+2A) = Pf(S)² ≤ H(2)².
What about det(I+xA) vs I(CG, x)² for other x?
""")

def count_hp_numeric(A_np, n):
    """Count Hamiltonian paths (= H(T))."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A_np[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

# For n=4, compare at various x values
n = 4
m = n*(n-1)//2

print(f"n={n}: det(I+xA) vs I(CG,x)² at various x values")
print(f"(showing one representative per H-value)\n")

# Group by H value
h_groups = {}
for bits in range(1 << m):
    b = [(bits >> i) & 1 for i in range(m)]
    A_np = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if b[idx]: A_np[i][j] = 1
            else: A_np[j][i] = 1
            idx += 1
    H = count_hp_numeric(A_np, n)
    if H not in h_groups:
        h_groups[H] = (bits, A_np)

for H_val in sorted(h_groups.keys()):
    bits, A_np = h_groups[H_val]
    A_list = adj_matrix_bits(bits, n)
    A_sym = Matrix(A_list)
    det_poly = expand(det(Matrix.eye(n) + x * A_sym))

    print(f"  H={H_val} (bits={bits:06b}):")
    print(f"    det(I+xA) = {factor(det_poly)}")

    for x_val in [0, 1, 2, 3, 4, -1]:
        d = int(det_poly.subs(x, x_val))
        h_sq = H_val**2 if x_val == 2 else "?"
        print(f"      x={x_val:2d}: det = {d:5d}", end="")
        if x_val == 2:
            print(f"  H² = {H_val**2}, gap = {H_val**2 - d}")
        else:
            print()
    print()

# ============================================================================
# PART 5: THE I(CG, x) POLYNOMIAL — WHAT IT LOOKS LIKE
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 5: INDEPENDENCE POLYNOMIAL I(CG, x) — THE H(x) GENERALIZATION        ║
╚══════════════════════════════════════════════════════════════════════════════╝

H(x) = I(CG(T), x) = Σ_k α_k(T) · x^k

At x=2 this gives H(T).
What is the polynomial itself?
""")

def find_odd_cycles(A, n):
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    mi = perm.index(min(perm))
                    canon = perm[mi:] + perm[:mi]
                    cycles.append(canon)
    return list(set(cycles))

def build_cg(cycles):
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = True
    return adj

def indep_poly(adj, nc):
    """Return alpha_k for k=0,...,nc."""
    alpha = [0]*(nc+1)
    for mask in range(1 << nc):
        verts = [i for i in range(nc) if mask & (1 << i)]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False; break
            if not ok: break
        if ok:
            alpha[len(verts)] += 1
    return alpha

for n in [3, 4, 5]:
    m = n*(n-1)//2
    print(f"n={n}:")
    seen = {}

    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A_np = np.zeros((n,n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if b[idx]: A_np[i][j] = 1
                else: A_np[j][i] = 1
                idx += 1

        H = count_hp_numeric(A_np, n)

        cycles = find_odd_cycles(A_np, n)
        if len(cycles) > 12: continue

        cg = build_cg(cycles)
        alpha = indep_poly(cg, len(cycles))

        poly_str = " + ".join(f"{alpha[k]}x^{k}" for k in range(len(alpha)) if alpha[k])
        key = tuple(alpha[:max((k for k,a in enumerate(alpha) if a), default=0)+1])

        if key not in seen:
            seen[key] = True

            # Build I(CG, x) and det(I+xA) symbolically
            Hx = sum(alpha[k] * x**k for k in range(len(alpha)))
            A_list = adj_matrix_bits(bits, n)
            A_sym = Matrix(A_list)
            det_x = expand(det(Matrix.eye(n) + x * A_sym))

            # Compare: Hx² vs det_x at several points
            diff = expand(Hx**2 - det_x)

            print(f"  H={H}, α={list(key)}")
            print(f"    I(CG,x) = {Hx}")
            print(f"    det(I+xA) = {factor(det_x)}")
            print(f"    I(CG,x)² - det(I+xA) = {factor(diff)}")

            # Check: is the difference divisible by x²(x+2)²/4 or similar?
            if diff != 0:
                # Evaluate at a few points
                for xv in [0, 1, 2, 3, -1, -2]:
                    val = diff.subs(x, xv)
                    print(f"      diff at x={xv}: {val}")
            print()
    print()

# ============================================================================
# PART 6: THE KEY FORMULA — I(CG,x)² - det(I+xA)
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 6: THE GAP FORMULA I(CG,x)² - det(I+xA)                              ║
╚══════════════════════════════════════════════════════════════════════════════╝

For each tournament T, define:
  Q(x) = I(CG(T), x)² - det(I + xA)

At x=2: Q(2) = H² - Pf² = 8·q (always divisible by 8).
At x=0: Q(0) = 1² - 1 = 0 (trivially).
At x=-2: Q(-2) = I(CG,-2)² - det(I-2A).

QUESTION: Does Q(x) have a nice factorization? Is it always ≥ 0 for x ≥ 0?
""")

# Detailed analysis for n=5
n = 5
m = n*(n-1)//2
print(f"n={n}: Analyzing Q(x) = I(CG,x)² - det(I+xA)")
print()

seen_Q = {}
for bits in range(1 << m):
    b = [(bits >> i) & 1 for i in range(m)]
    A_np = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if b[idx]: A_np[i][j] = 1
            else: A_np[j][i] = 1
            idx += 1

    H = count_hp_numeric(A_np, n)

    cycles = find_odd_cycles(A_np, n)
    if len(cycles) > 12: continue

    cg = build_cg(cycles)
    alpha = indep_poly(cg, len(cycles))

    Hx = sum(alpha[k] * x**k for k in range(len(alpha)))
    A_list = adj_matrix_bits(bits, n)
    A_sym = Matrix(A_list)
    det_x = expand(det(Matrix.eye(n) + x * A_sym))
    Q = expand(Hx**2 - det_x)

    Q_key = str(factor(Q))
    if Q_key not in seen_Q:
        seen_Q[Q_key] = (H, Q, factor(Q), Hx, det_x)

print(f"Found {len(seen_Q)} distinct Q(x) polynomials:\n")
for Q_key, (H, Q_poly, Q_fact, Hx, det_x) in sorted(seen_Q.items(), key=lambda kv: kv[1][0]):
    print(f"  H={H}: Q(x) = {Q_fact}")
    # Check divisibility by x²
    Q_poly_expanded = expand(Q_poly)
    if Q_poly_expanded != 0:
        p = Poly(Q_poly_expanded, x)
        coeffs = p.all_coeffs()
        min_degree = p.degree() - len(coeffs) + 1
        # Evaluate
        print(f"    Q(0)={Q_poly_expanded.subs(x,0)}, Q(1)={Q_poly_expanded.subs(x,1)}, Q(2)={Q_poly_expanded.subs(x,2)}, Q(-1)={Q_poly_expanded.subs(x,-1)}")
    else:
        print(f"    Q(x) = 0 identically!")
    print()

print("""
════════════════════════════════════════════════════════════════════════════════
KEY FINDING SUMMARY
════════════════════════════════════════════════════════════════════════════════

If Q(x) = 0 identically for some tournaments, those are the cases where
I(CG,x)² = det(I+xA) as POLYNOMIALS, not just at x=2.

If Q(x) factors as x² · P(x) with P(x) ≥ 0 for x ≥ 0, this would explain
why Q(2) = 8q with q ≥ 0.

The x-deformation reveals the algebraic structure behind the H ≥ |Pf| inequality.
""")
