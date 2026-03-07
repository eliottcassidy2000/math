#!/usr/bin/env python3
"""
Reduced polynomial P(u, x) analysis at n=5.

For n=5 (odd), the generating function G_T(t, x) is palindromic of degree 4 in t.
Writing u = t + 1/t, we get G_T(t, x) = t^2 * P(u, x) where P is quadratic in u:

  P(u, x) = p_0(x) + p_1(x)*u + p_2(x)*u^2

The invariants at n=5 are:
  - t3: number of directed 3-cycles (f=2, parts=1)
  - t5: number of directed 5-cycles (f=0, parts=1)
  - No bc (no room for disjoint 3-cycle pairs at n=5)

We compute P(u, x) for ALL tournaments on 5 vertices and express p_0, p_1, p_2
as functions of (t3, t5), comparing with the n=7 structure.

At n=7, the results were:
  p_3(x) = I(Omega, x) = independence polynomial
  p_2(x) = A(7,1) + (24*t3 - 6*t7)*x
  p_1(x) = (A(7,2)-3) + 12*(t3-t5+t7)*x - 12*bc*x^2
  p_0(x) = (A(7,3)-2*A(7,1)) + (-128*t3+16*t5-8*t7)*x + 16*bc*x^2
"""

from itertools import permutations, combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction


# ====================================================================
# Core tournament utilities
# ====================================================================

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def count_t3(A, n):
    """Count directed 3-cycles."""
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3


def count_directed_cycles(A, n, cl):
    """Count directed cycles of length cl."""
    if n < cl: return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cl):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total


def forward_edge_dist_dp(A, n):
    """Compute the distribution of Hamiltonian paths by number of forward edges.
    Returns dict: k -> number of H-paths with exactly k forward edges."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    new_fwd = fwd + A[v][u]
                    key = (mask | (1 << u), u, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)


def eulerian_number(n, k):
    """Eulerian number A(n,k): number of permutations of [n] with exactly k ascents."""
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))


def find_odd_cycles_vertex_sets(A, n, v):
    """Find all directed odd cycles through vertex v. Returns list of vertex sets (as frozensets)."""
    cycles = []
    other = [u for u in range(n) if u != v]
    for L in range(3, n+1, 2):
        for subset in combinations(other, L-1):
            for perm in permutations(subset):
                if A[v][perm[0]] != 1: continue
                valid = True
                for i in range(len(perm)-1):
                    if A[perm[i]][perm[i+1]] != 1:
                        valid = False; break
                if valid and A[perm[-1]][v] == 1:
                    cycles.append(frozenset([v] + list(perm)))
    return cycles


def independence_poly(A, n, v):
    """Compute independence polynomial I(Omega(T-v), x) of the conflict graph.
    Omega(T-v) has odd cycles of T-v as vertices, edges when they share a vertex."""
    other = [u for u in range(n) if u != v]
    sub = [[A[other[i]][other[j]] for j in range(len(other))] for i in range(len(other))]

    # Find odd cycles of T-v
    m = len(other)
    odd_cycles = []
    for L in range(3, m+1, 2):
        for verts in combinations(range(m), L):
            ssub = [[sub[verts[i]][verts[j]] for j in range(len(verts))] for i in range(len(verts))]
            dp = [[0]*len(verts) for _ in range(1 << len(verts))]
            dp[1][0] = 1
            for mask in range(1, 1 << len(verts)):
                for vi in range(len(verts)):
                    if not (mask & (1 << vi)) or dp[mask][vi] == 0: continue
                    for ui in range(len(verts)):
                        if mask & (1 << ui): continue
                        if ssub[vi][ui]: dp[mask | (1 << ui)][ui] += dp[mask][vi]
            full = (1 << len(verts)) - 1
            cyc_count = sum(dp[full][vi] for vi in range(1, len(verts)) if ssub[vi][0])
            if cyc_count > 0:
                odd_cycles.append(frozenset(verts))

    # Build conflict graph
    nc = len(odd_cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if odd_cycles[i] & odd_cycles[j]:
                adj[i][j] = adj[j][i] = True

    # Independence polynomial by inclusion-exclusion over independent sets
    # I(G, x) = sum_{k} alpha_k * x^k where alpha_k = # independent sets of size k
    coeffs = [0] * (nc + 1)
    for mask in range(1 << nc):
        verts_in = [i for i in range(nc) if mask & (1 << i)]
        # Check independence
        indep = True
        for i in range(len(verts_in)):
            for j in range(i+1, len(verts_in)):
                if adj[verts_in[i]][verts_in[j]]:
                    indep = False; break
            if not indep: break
        if indep:
            coeffs[len(verts_in)] += 1
    return coeffs


# ====================================================================
# G_T(t, x) as a polynomial in t, with coefficients depending on x
# ====================================================================

def inflated_eulerian(f, d, k):
    """Coefficient of t^k in A_{f+1}(t) * (t-1)^{d-f}, which is the contribution
    to [t^k] of a cycle family with f forward edges in a tournament of degree d=n-1."""
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total


def G_T_coefficients(A, n):
    """Compute G_T(t, x) as array of (n-1+1) lists, where entry k is [x^0, x^1, ...] coeff of t^k.
    Uses the forward-edge distribution for x=0 (Eulerian) and cycle counts for higher x."""
    d = n - 1

    # Forward edge distribution gives the t-distribution at x=2 (full generating function)
    dist = forward_edge_dist_dp(A, n)

    # For the polynomial approach: G_T(t, x) = A_n(t) + sum over cycle families I of
    # x^{parts(I)} * count(I) * 2^{parts} * A_{f+1}(t) * (t-1)^{d-f}

    # Invariants at n=5: t3 (f=2, parts=1), t5 (f=0, parts=1)
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)

    # Coefficients of t^k for each x-power
    # x^0: Eulerian numbers A(n, k)
    # x^1: 2 * [t3 * inflated_eulerian(2, d, k) + t5 * inflated_eulerian(0, d, k)]
    # x^2: nothing (no bc at n=5)

    result = []
    for k in range(n):
        x0 = eulerian_number(n, k)
        x1 = 2 * (t3 * inflated_eulerian(2, d, k) + t5 * inflated_eulerian(0, d, k))
        result.append([x0, x1])

    return result, t3, t5


# ====================================================================
# Main analysis
# ====================================================================

print("=" * 70)
print("REDUCED POLYNOMIAL P(u, x) AT n=5")
print("=" * 70)

n = 5
d = n - 1  # = 4
m = d // 2  # = 2, so P is quadratic in u

print(f"""
n = {n}, degree in t = {d} (palindromic), degree in u = {m}
G_T(t, x) = t^{m} * P(u, x)  where  u = t + 1/t
P(u, x) = p_0(x) + p_1(x)*u + p_2(x)*u^2

Invariants: t3 (directed 3-cycles), t5 (directed 5-cycles)
No disjoint cycle pairs (bc) at n=5.
""")

# For the u-substitution, we need to convert from t-coefficients to u-coefficients.
# t^2 * P(u) = t^2 * (p0 + p1*(t+1/t) + p2*(t+1/t)^2)
#            = t^2 * (p0 + p1*t + p1/t + p2*t^2 + 2*p2 + p2/t^2)
#            = p2*t^4 + p1*t^3 + (p0 + 2*p2)*t^2 + p1*t + p2
#
# So: a_4 = p2, a_3 = p1, a_2 = p0 + 2*p2, a_1 = p1, a_0 = p2
# Palindromic check: a_0 = a_4 = p2, a_1 = a_3 = p1. OK!
#
# Inversion: p2 = a_4, p1 = a_3, p0 = a_2 - 2*a_4

print("Converting t-coefficients to u-coefficients:")
print("  a_4 = p_2")
print("  a_3 = p_1")
print("  a_2 = p_0 + 2*p_2")
print("  a_1 = p_1  (palindrome)")
print("  a_0 = p_2  (palindrome)")
print()
print("  => p_2 = a_4,  p_1 = a_3,  p_0 = a_2 - 2*a_4")
print()

# Collect data for all tournaments
print("-" * 70)
print("Scanning all tournaments on 5 vertices...")

# Group by isomorphism class (score sequence)
data_by_inv = defaultdict(list)

count = 0
for A in all_tournaments(n):
    coeffs, t3, t5 = G_T_coefficients(A, n)
    count += 1

    # coeffs[k] = [x^0 coeff, x^1 coeff] of t^k
    # Extract a_k for each x-power
    a_x0 = [coeffs[k][0] for k in range(n)]  # Eulerian distribution
    a_x1 = [coeffs[k][1] for k in range(n)]  # x^1 coefficient of t^k

    # Verify palindromic symmetry
    for k in range(n):
        assert a_x0[k] == a_x0[d-k], f"Eulerian not palindromic at k={k}"
        assert a_x1[k] == a_x1[d-k], f"x^1 not palindromic at k={k}"

    # Convert to u-coefficients
    # p_2 = a_4, p_1 = a_3, p_0 = a_2 - 2*a_4
    p2_x0 = a_x0[4]
    p1_x0 = a_x0[3]
    p0_x0 = a_x0[2] - 2 * a_x0[4]

    p2_x1 = a_x1[4]
    p1_x1 = a_x1[3]
    p0_x1 = a_x1[2] - 2 * a_x1[4]

    data_by_inv[(t3, t5)].append({
        'a_x0': a_x0, 'a_x1': a_x1,
        'p0': (p0_x0, p0_x1), 'p1': (p1_x0, p1_x1), 'p2': (p2_x0, p2_x1),
    })

print(f"Total tournaments: {count}")
print(f"Distinct (t3, t5) classes: {len(data_by_inv)}")
print()

# ====================================================================
# Display results by invariant class
# ====================================================================
print("=" * 70)
print("P-COEFFICIENTS BY INVARIANT CLASS (t3, t5)")
print("=" * 70)
print()
print(f"{'t3':>4s} {'t5':>4s} | {'p2(x)':>24s} | {'p1(x)':>24s} | {'p0(x)':>24s} | {'#T':>4s}")
print("-" * 90)

all_results = []
for (t3, t5) in sorted(data_by_inv.keys()):
    entries = data_by_inv[(t3, t5)]
    # Check all tournaments with same invariants give same P
    p0s = set(e['p0'] for e in entries)
    p1s = set(e['p1'] for e in entries)
    p2s = set(e['p2'] for e in entries)
    assert len(p0s) == 1, f"p0 not unique for (t3,t5)=({t3},{t5}): {p0s}"
    assert len(p1s) == 1, f"p1 not unique for (t3,t5)=({t3},{t5}): {p1s}"
    assert len(p2s) == 1, f"p2 not unique for (t3,t5)=({t3},{t5}): {p2s}"

    p0 = entries[0]['p0']
    p1 = entries[0]['p1']
    p2 = entries[0]['p2']

    def fmt_poly(c0, c1):
        if c1 == 0:
            return f"{c0}"
        elif c1 > 0:
            return f"{c0} + {c1}x"
        else:
            return f"{c0} - {-c1}x"

    p2_str = fmt_poly(p2[0], p2[1])
    p1_str = fmt_poly(p1[0], p1[1])
    p0_str = fmt_poly(p0[0], p0[1])

    print(f"{t3:4d} {t5:4d} | {p2_str:>24s} | {p1_str:>24s} | {p0_str:>24s} | {len(entries):4d}")

    all_results.append((t3, t5, p0, p1, p2, len(entries)))

# ====================================================================
# Analyze the pattern
# ====================================================================
print()
print("=" * 70)
print("PATTERN ANALYSIS")
print("=" * 70)

# Extract the constant terms (x^0)
print("\n--- Constant terms p_j(0) ---")
p2_consts = set(r[4][0] for r in all_results)
p1_consts = set(r[3][0] for r in all_results)
p0_consts = set(r[2][0] for r in all_results)

print(f"  p_2(0) values: {sorted(p2_consts)}")
print(f"  p_1(0) values: {sorted(p1_consts)}")
print(f"  p_0(0) values: {sorted(p0_consts)}")

print(f"\n  Eulerian numbers A(5,k): ", end="")
for k in range(n):
    print(f"A(5,{k})={eulerian_number(5,k)} ", end="")
print()

# Check: are p_j(0) Eulerian-related?
A5 = [eulerian_number(5, k) for k in range(n)]
print(f"\n  A(5,0)={A5[0]}, A(5,1)={A5[1]}, A(5,2)={A5[2]}, A(5,3)={A5[3]}, A(5,4)={A5[4]}")

# From the t-expansion: a_k(x=0) = A(5,k)
# p_2(0) = a_4(0) = A(5,4) = 1
# p_1(0) = a_3(0) = A(5,3) = 26
# p_0(0) = a_2(0) - 2*a_4(0) = A(5,2) - 2*A(5,4) = 66 - 2 = 64
print(f"\n  Expected from Eulerian:")
print(f"    p_2(0) = A(5,4) = {A5[4]}")
print(f"    p_1(0) = A(5,3) = {A5[3]}")
print(f"    p_0(0) = A(5,2) - 2*A(5,4) = {A5[2]} - 2*{A5[4]} = {A5[2] - 2*A5[4]}")
print(f"    (for n=7: p_0(0) = A(7,3) - 2*A(7,1) = 2416 - 240 = 2176)")

# Linear coefficients (x^1)
print("\n--- Linear coefficients (coefficient of x in p_j) ---")
print("  Fitting: p_j(x) = p_j(0) + (a*t3 + b*t5)*x")
print()

for j, label in [(4, 'p_2'), (3, 'p_1'), (2, 'p_0')]:
    # Collect (t3, t5, x1_coeff)
    if j == 4:
        data_pts = [(r[0], r[1], r[4][1]) for r in all_results]
    elif j == 3:
        data_pts = [(r[0], r[1], r[3][1]) for r in all_results]
    else:
        data_pts = [(r[0], r[1], r[2][1]) for r in all_results]

    # Solve a*t3 + b*t5 = x1_coeff for all data points
    # Use first two with different (t3, t5) to solve, then verify
    if len(data_pts) >= 2:
        # Try to fit: coeff_x1 = a * t3 + b * t5
        # Use least squares via two data points
        from fractions import Fraction
        solved = False
        for i in range(len(data_pts)):
            for jj in range(i+1, len(data_pts)):
                t3_i, t5_i, c_i = data_pts[i]
                t3_j, t5_j, c_j = data_pts[jj]
                det = t3_i * t5_j - t3_j * t5_i
                if det != 0:
                    a = Fraction(c_i * t5_j - c_j * t5_i, det)
                    b = Fraction(t3_i * c_j - t3_j * c_i, det)
                    # Verify against all
                    ok = True
                    for t3_k, t5_k, c_k in data_pts:
                        if a * t3_k + b * t5_k != c_k:
                            ok = False; break
                    if ok:
                        print(f"  {label}_1 = {a}*t3 + {b}*t5")
                        solved = True
                        break
            if solved: break
        if not solved:
            # Maybe there's a constant term too
            # coeff_x1 = c + a*t3 + b*t5
            for i in range(len(data_pts)):
                for jj in range(i+1, len(data_pts)):
                    for kk in range(jj+1, len(data_pts)):
                        t3_i, t5_i, c_i = data_pts[i]
                        t3_j, t5_j, c_j = data_pts[jj]
                        t3_k, t5_k, c_k = data_pts[kk]
                        # c + a*t3 + b*t5 = val
                        # [1 t3_i t5_i] [c]   [c_i]
                        # [1 t3_j t5_j] [a] = [c_j]
                        # [1 t3_k t5_k] [b]   [c_k]
                        det = (1*(t3_j*t5_k - t3_k*t5_j) -
                               t3_i*(1*t5_k - 1*t5_j) +
                               t5_i*(1*t3_k - 1*t3_j))
                        # Simpler: just use matrix
                        det = (t3_j*t5_k - t3_k*t5_j - t3_i*t5_k + t3_i*t5_j +
                               t5_i*t3_k - t5_i*t3_j)
                        if det != 0:
                            # Cramer's rule
                            det_a = (c_i*(t3_j*t5_k - t3_k*t5_j) -
                                     t3_i*(c_j*t5_k - c_k*t5_j) +
                                     t5_i*(c_j*t3_k - c_k*t3_j))
                            # Too messy, use numpy
                            import numpy as np
                            M = np.array([[1, t3_i, t5_i],
                                          [1, t3_j, t5_j],
                                          [1, t3_k, t5_k]], dtype=float)
                            rhs = np.array([c_i, c_j, c_k], dtype=float)
                            try:
                                sol = np.linalg.solve(M, rhs)
                                c_val, a_val, b_val = sol
                                # Verify
                                ok = True
                                for t3_m, t5_m, c_m in data_pts:
                                    pred = c_val + a_val * t3_m + b_val * t5_m
                                    if abs(pred - c_m) > 0.01:
                                        ok = False; break
                                if ok:
                                    # Convert to fractions
                                    c_fr = Fraction(round(c_val)).limit_denominator(1000)
                                    a_fr = Fraction(round(a_val)).limit_denominator(1000)
                                    b_fr = Fraction(round(b_val)).limit_denominator(1000)
                                    print(f"  {label}_1 = {c_fr} + {a_fr}*t3 + {b_fr}*t5")
                                    solved = True
                                    break
                            except: pass
                    if solved: break
                if solved: break
            if not solved:
                print(f"  {label}_1: could not find clean linear relation")
                print(f"    data: {data_pts}")

# ====================================================================
# Verification: evaluate P(u, x) at specific u values
# ====================================================================
print()
print("=" * 70)
print("VERIFICATION: P(u, x) AT SPECIAL VALUES")
print("=" * 70)

for idx, (t3, t5, p0, p1, p2, cnt) in enumerate(all_results):
    if idx > 3: break
    p0_x0, p0_x1 = p0
    p1_x0, p1_x1 = p1
    p2_x0, p2_x1 = p2

    # At x=0: P(u, 0) = p0_x0 + p1_x0*u + p2_x0*u^2
    # This should be related to Eulerian
    print(f"\n  (t3={t3}, t5={t5}):")
    print(f"    P(u, 0) = {p0_x0} + {p1_x0}*u + {p2_x0}*u^2")

    # Check: P(2, 0) should give A_5(t=1) = 5! = 120
    P_at_u2 = p0_x0 + p1_x0 * 2 + p2_x0 * 4
    print(f"    P(2, 0) = {P_at_u2}  (should be {factorial(n)} = n!)")

    # At x=2: P(u, 2)
    P_at_u2_x2_const = p0_x0 + 2*p0_x1
    P_at_u2_x2_lin = p1_x0 + 2*p1_x1
    P_at_u2_x2_quad = p2_x0 + 2*p2_x1
    H = P_at_u2_x2_const + P_at_u2_x2_lin * 2 + P_at_u2_x2_quad * 4
    print(f"    P(2, 2) = {H}  (= H(T) = {factorial(n)} + 2*2*({t3}+{t5}) = {factorial(n) + 4*(t3+t5)})")

    # At u = -2 (t = -1): gives the alternating sum
    P_at_um2 = p0_x0 - p1_x0 * 2 + p2_x0 * 4
    print(f"    P(-2, 0) = {P_at_um2}  (alternating Eulerian)")

    # At u = 0 (t = i): gives P(0, 0)
    print(f"    P(0, 0) = {p0_x0}")


# ====================================================================
# Compare with n=7 structure
# ====================================================================
print()
print("=" * 70)
print("COMPARISON WITH n=7 STRUCTURE")
print("=" * 70)

# n=7 results (from the problem statement):
# p_3(x) = I(Omega, x)                  [p_3(0) = 1 = A(7,6)]
# p_2(x) = 120 + (24*t3 - 6*t7)*x       [p_2(0) = 120 = A(7,1) = A(7,5)]
# p_1(x) = 1188 + 12*(t3-t5+t7)*x - 12*bc*x^2  [p_1(0) = 1188 = A(7,2)-3]
# p_0(x) = 2176 + (-128*t3+16*t5-8*t7)*x + 16*bc*x^2  [p_0(0) = 2176 = A(7,3)-2*A(7,1)]

print(f"""
n=5 structure (this computation):
  P(u, x) = p_0(x) + p_1(x)*u + p_2(x)*u^2
  degree in u = 2 = (n-1)/2

n=7 structure (given):
  P(u, x) = p_0(x) + p_1(x)*u + p_2(x)*u^2 + p_3(x)*u^3
  degree in u = 3 = (n-1)/2

Constant terms (x=0):
  n=5: p_2(0) = {A5[4]} = A(5,4),  p_1(0) = {A5[3]} = A(5,3),  p_0(0) = {A5[2] - 2*A5[4]} = A(5,2) - 2*A(5,4)
  n=7: p_3(0) = 1 = A(7,6),  p_2(0) = 120 = A(7,5),  p_1(0) = 1188,  p_0(0) = 2176

  For general n, the pattern for top two:
    p_m(0) = A(n, n-1) = 1
    p_{{m-1}}(0) = A(n, n-2) = 2^n - n - 1  [A(5,3)={A5[3]}, A(7,5)={eulerian_number(7,5)}]
""")

# Check if A(5,3) = 2^5 - 5 - 1 = 26
print(f"  A(5,3) = {A5[3]},  2^5 - 6 = {2**5 - 6}")
print(f"  A(7,5) = {eulerian_number(7,5)},  2^7 - 8 = {2**7 - 8}")

# Check the Chebyshev-like formula for p_0(0)
# n=5: p_0(0) = A(5,2) - 2*A(5,4) = 66 - 2 = 64
# n=7: p_0(0) = A(7,3) - 2*A(7,1) = 2416 - 240 = 2176
# Pattern? For degree m, p_0(0) = A(n, m) - 2*A(n, m+2) + ...?
# Actually: this is from the Chebyshev relation.
# t^m * u^j has t^k coefficient = C(j, (k-m+j)/2) if k-m+j is even.
# For u^0: only t^m. For u^2: t^{m+2} + 2*t^m + t^{m-2}. Etc.
# So a_m = p_0 + 2*p_2 + ... and p_0 = a_m - 2*p_2 + ... (Mobius inversion).

print(f"\n  Chebyshev decomposition of A_n(t):")
print(f"  n=5: A_5(t) = {A5[0]} + {A5[1]}t + {A5[2]}t^2 + {A5[3]}t^3 + {A5[4]}t^4")
print(f"  t^2 * (p_0 + p_1*u + p_2*u^2):")
print(f"    = p_2*t^4 + p_1*t^3 + (p_0+2p_2)*t^2 + p_1*t + p_2")
print(f"  So: A(5,4) = p_2 = {A5[4]}")
print(f"      A(5,3) = p_1 = {A5[3]}")
print(f"      A(5,2) = p_0 + 2*p_2 => p_0 = A(5,2) - 2 = {A5[2] - 2}")


# ====================================================================
# Summary table
# ====================================================================
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()

for t3, t5, p0, p1, p2, cnt in all_results:
    print(f"  (t3={t3:2d}, t5={t5:3d}): p_2 = {p2[0]:4d} + {p2[1]:5d}x,  "
          f"p_1 = {p1[0]:4d} + {p1[1]:5d}x,  "
          f"p_0 = {p0[0]:4d} + {p0[1]:5d}x  [{cnt} tournaments]")

print()
print("DONE")
