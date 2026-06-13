#!/usr/bin/env python3
"""
The reduced polynomial P(u, x) where G_T(t, x) = t^{(n-1)/2} P(t + 1/t, x).

For odd n, the palindromic polynomial G_T(t, x) of degree n-1 in t
factors through u = t + 1/t, giving P of degree (n-1)/2 in u.

G_T(t, x) = t^m * (p_0(x) + p_1(x)*u + p_2(x)*u^2 + p_3(x)*u^3)

where m = (n-1)/2.

GOAL: Find explicit OCF formulas for p_j(x) in terms of cycle counts.

METHOD: Evaluate G_T at 4 specific t-values (giving 4 u-values) and
solve for the 4 coefficients p_0, p_1, p_2, p_3 as linear functions
of the invariants.

opus-2026-03-07-S33
"""
from itertools import combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import random
import numpy as np

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def random_tournament(n, seed=42):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
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

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def forward_edge_dist_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0: continue
                for u_node in range(n):
                    if mask & (1 << u_node): continue
                    new_fwd = fwd + A[v][u_node]
                    key = (mask | (1 << u_node), u_node, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)

def G_T_formula(n, inv_vals, t):
    """Compute G_T(t, 2) via the OCF formula."""
    d = n - 1
    result = sum(eulerian_number(n, k) * t**k for k in range(n))
    if n == 7:
        invariants = [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]
    elif n == 5:
        invariants = [('t3', 2, 1), ('t5', 0, 1)]
    else:
        return result
    for name, f, parts in invariants:
        val = inv_vals.get(name, 0)
        if val == 0: continue
        A_f1 = sum(eulerian_number(f+1, j) * t**j for j in range(f+1))
        result += 2**parts * val * A_f1 * (t - 1)**(d - f)
    return result

# ====================================================================
# Part 1: Linear regression for P-coefficients
# ====================================================================
print("REDUCED POLYNOMIAL P(u) COEFFICIENTS VIA REGRESSION")
print("=" * 70)

n = 7
d = n - 1
m = d // 2  # = 3

# Collect data: for each tournament, compute p_0, p_1, p_2, p_3
# Then regress against (1, t3, t5, t7, bc)

data = []
for seed in range(50):
    A = random_tournament(n, n * 4000 + seed)
    inv = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    dist = forward_edge_dist_dp(A, n)
    a = [dist.get(k, 0) for k in range(n)]

    # Convert palindromic (a_0, ..., a_6) to P-coefficients.
    # G(t) = a_0 + a_1*t + a_2*t^2 + a_3*t^3 + a_2*t^4 + a_1*t^5 + a_0*t^6
    # = t^3 * P(u) where u = t + 1/t
    # P(u) = p_0 + p_1*u + p_2*u^2 + p_3*u^3
    #
    # Expanding t^3 * P(t + 1/t):
    # u^0 = 1: t^3
    # u^1 = t + 1/t: t^4 + t^2
    # u^2 = t^2 + 2 + 1/t^2: t^5 + 2*t^3 + t
    # u^3 = t^3 + 3*t + 3/t + 1/t^3: t^6 + 3*t^4 + 3*t^2 + 1
    #
    # So:
    # coeff of t^0 (= a_0): p_3
    # coeff of t^1 (= a_1): p_2
    # coeff of t^2 (= a_2): p_1 + 3*p_3
    # coeff of t^3 (= a_3): p_0 + 2*p_2
    # (and by palindromy: t^4 coeff = a_2 = p_1 + 3*p_3, t^5 = a_1 = p_2, t^6 = a_0 = p_3)
    #
    # Inverting:
    p3 = Fraction(a[0])
    p2 = Fraction(a[1])
    p1 = Fraction(a[2]) - 3 * p3
    p0 = Fraction(a[3]) - 2 * p2

    data.append({
        't3': inv['t3'], 't5': inv['t5'], 't7': inv['t7'], 'bc': inv['bc'],
        'p0': p0, 'p1': p1, 'p2': p2, 'p3': p3,
        'H': a[0]  # = p3
    })

# p3 = a_0 = H(T) = 1 + 2*alpha_1 + 4*alpha_2 = OCF!
# p2 = a_1 = the second coefficient
# p1 = a_2 - 3*H
# p0 = a_3 - 2*a_1

print("\nDirect relations:")
print("  p_3 = a_0 = H(T)")
print("  p_2 = a_1")
print("  p_1 = a_2 - 3*H = a_2 - 3*a_0")
print("  p_0 = a_3 - 2*a_1")

# Now regress each p_j against (1, t3, t5, t7, bc)
print("\nRegression of P-coefficients against cycle invariants:")
X = np.array([[1, d['t3'], d['t5'], d['t7'], d['bc']] for d in data])
names = ['const', 't3', 't5', 't7', 'bc']

for pname in ['p0', 'p1', 'p2', 'p3']:
    y = np.array([float(d[pname]) for d in data])
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    max_resid = max(abs(y[i] - X[i] @ coeffs) for i in range(len(y)))

    print(f"\n  {pname} = ", end="")
    terms = []
    for i, c in enumerate(coeffs):
        if abs(c) > 0.001:
            c_frac = Fraction(round(c)).limit_denominator(100)
            terms.append(f"{c_frac}*{names[i]}" if i > 0 else str(c_frac))
    print(" + ".join(terms))
    print(f"    max residual: {max_resid:.6f}")

# ====================================================================
# Part 2: Verify the clean formulas
# ====================================================================
print(f"\n{'=' * 70}")
print("CLEAN FORMULAS FOR P-COEFFICIENTS AT n=7")
print("=" * 70)

# From the regression, we expect exact integer coefficients.
# Let's verify directly using the OCF.

# We know from THM-062:
# a_k(T) = A(7,k) + 2*c_k^{(4,6)}*t3 + 2*c_k^{(2,6)}*t5 + 2*c_k^{(0,6)}*t7 + 4*c_k^{(2,6)}*bc

# A(7,k) = [1, 120, 1191, 2416, 1191, 120, 1]
# c_k^{(4,6)} for t3: [1, 24, 15, -80, 15, 24, 1]
# c_k^{(2,6)} for t5: [1, 0, -9, 16, -9, 0, 1]
# c_k^{(0,6)} for t7: [1, -6, 15, -20, 15, -6, 1]
# bc gets 4*c_k^{(2,6)}: [4, 0, -36, 64, -36, 0, 4]

def get_ck(f, d):
    """Get c_k^{(f,d)} for k = 0, ..., d."""
    from math import comb as C
    result = []
    for k in range(d + 1):
        total = 0
        for j in range(max(0, k - (d - f)), min(f, k) + 1):
            sign = (-1) ** (d - f - k + j)
            total += eulerian_number(f + 1, j) * C(d - f, k - j) * sign
        result.append(total)
    return result

ck4 = get_ck(4, 6)  # for t3
ck2 = get_ck(2, 6)  # for t5, bc
ck0 = get_ck(0, 6)  # for t7
An = [eulerian_number(7, k) for k in range(7)]

print(f"A(7,k) = {An}")
print(f"c_k^(4,6) = {ck4}")
print(f"c_k^(2,6) = {ck2}")
print(f"c_k^(0,6) = {ck0}")

# p3 = a_0 = A(7,0) + 2*c_0^{(4,6)}*t3 + 2*c_0^{(2,6)}*t5 + 2*c_0^{(0,6)}*t7 + 4*c_0^{(2,6)}*bc
# p3 = 1 + 2*1*t3 + 2*1*t5 + 2*1*t7 + 4*1*bc = 1 + 2*(t3+t5+t7) + 4*bc = H(T)
print(f"\np3 = a_0:")
print(f"  = A(7,0) + 2*{ck4[0]}*t3 + 2*{ck2[0]}*t5 + 2*{ck0[0]}*t7 + 4*{ck2[0]}*bc")
print(f"  = {An[0]} + {2*ck4[0]}*t3 + {2*ck2[0]}*t5 + {2*ck0[0]}*t7 + {4*ck2[0]}*bc")
print(f"  = 1 + 2*(t3+t5+t7) + 4*bc = H(T) [this is the OCF!]")

# p2 = a_1
print(f"\np2 = a_1:")
print(f"  = A(7,1) + 2*{ck4[1]}*t3 + 2*{ck2[1]}*t5 + 2*{ck0[1]}*t7 + 4*{ck2[1]}*bc")
print(f"  = {An[1]} + {2*ck4[1]}*t3 + {2*ck2[1]}*t5 + {2*ck0[1]}*t7 + {4*ck2[1]}*bc")

# p1 = a_2 - 3*a_0
a2_formula = f"{An[2]} + {2*ck4[2]}*t3 + {2*ck2[2]}*t5 + {2*ck0[2]}*t7 + {4*ck2[2]}*bc"
a0_formula = f"{An[0]} + {2*ck4[0]}*t3 + {2*ck2[0]}*t5 + {2*ck0[0]}*t7 + {4*ck2[0]}*bc"
p1_const = An[2] - 3*An[0]
p1_t3 = 2*ck4[2] - 3*2*ck4[0]
p1_t5 = 2*ck2[2] - 3*2*ck2[0]
p1_t7 = 2*ck0[2] - 3*2*ck0[0]
p1_bc = 4*ck2[2] - 3*4*ck2[0]
print(f"\np1 = a_2 - 3*a_0:")
print(f"  = ({a2_formula}) - 3*({a0_formula})")
print(f"  = {p1_const} + {p1_t3}*t3 + {p1_t5}*t5 + {p1_t7}*t7 + {p1_bc}*bc")

# p0 = a_3 - 2*a_1
p0_const = An[3] - 2*An[1]
p0_t3 = 2*ck4[3] - 2*2*ck4[1]
p0_t5 = 2*ck2[3] - 2*2*ck2[1]
p0_t7 = 2*ck0[3] - 2*2*ck0[1]
p0_bc = 4*ck2[3] - 2*4*ck2[1]
print(f"\np0 = a_3 - 2*a_1:")
print(f"  = {p0_const} + {p0_t3}*t3 + {p0_t5}*t5 + {p0_t7}*t7 + {p0_bc}*bc")

# ====================================================================
# Verify
# ====================================================================
print(f"\n{'=' * 70}")
print("VERIFICATION")
print("=" * 70)

for seed in range(10):
    d_entry = data[seed]
    t3, t5, t7, bc = d_entry['t3'], d_entry['t5'], d_entry['t7'], d_entry['bc']

    pred_p3 = 1 + 2*(t3+t5+t7) + 4*bc
    pred_p2 = An[1] + 2*ck4[1]*t3 + 2*ck2[1]*t5 + 2*ck0[1]*t7 + 4*ck2[1]*bc
    pred_p1 = p1_const + p1_t3*t3 + p1_t5*t5 + p1_t7*t7 + p1_bc*bc
    pred_p0 = p0_const + p0_t3*t3 + p0_t5*t5 + p0_t7*t7 + p0_bc*bc

    ok = (pred_p3 == d_entry['p3'] and pred_p2 == d_entry['p2'] and
          pred_p1 == d_entry['p1'] and pred_p0 == d_entry['p0'])

    if seed < 5 or not ok:
        print(f"  seed={seed}: p3={d_entry['p3']} vs {pred_p3}, p2={d_entry['p2']} vs {pred_p2}, "
              f"p1={d_entry['p1']} vs {pred_p1}, p0={d_entry['p0']} vs {pred_p0} {'OK' if ok else 'FAIL'}")

# ====================================================================
# Part 3: P(u, x) with symbolic x
# ====================================================================
print(f"\n{'=' * 70}")
print("P(u, x) WITH SYMBOLIC x")
print("=" * 70)

# At general x, G_T(t, x) = A_n(t) + x*(...) + x^2*(...)
# For n=7: single-cycle terms get x, pair terms get x^2.
#
# a_k(T, x) = A(7,k) + x*(2*c_k^{(4,6)}*t3 + 2*c_k^{(2,6)}*t5 + 2*c_k^{(0,6)}*t7) + x^2*(4*c_k^{(2,6)}*bc)
#
# p3(x) = a_0(x) = 1 + x*(2*t3 + 2*t5 + 2*t7) + x^2*(4*bc) = 1 + 2*alpha_1*x + 4*alpha_2*x^2 = I(Omega, x)
# p2(x) = a_1(x) = 120 + x*(48*t3 + 0*t5 - 12*t7) + x^2*(0*bc) = 120 + x*(48*t3 - 12*t7)
# p1(x) = a_2(x) - 3*a_0(x)
# p0(x) = a_3(x) - 2*a_1(x)

print("p3(x) = I(Omega(T), x)  [THE independence polynomial!]")
print(f"p2(x) = {An[1]} + ({2*ck4[1]}*t3 + {2*ck2[1]}*t5 + {2*ck0[1]}*t7)*x + ({4*ck2[1]}*bc)*x^2")
print(f"       = 120 + (48*t3 - 12*t7)*x")
print(f"p1(x) = {p1_const} + ({p1_t3}*t3 + {p1_t5}*t5 + {p1_t7}*t7)*x + ({p1_bc}*bc)*x^2")
print(f"p0(x) = {p0_const} + ({p0_t3}*t3 + {p0_t5}*t5 + {p0_t7}*t7)*x + ({p0_bc}*bc)*x^2")

# ====================================================================
# Part 4: The beautiful structure
# ====================================================================
print(f"\n{'=' * 70}")
print("THE BEAUTIFUL STRUCTURE")
print("=" * 70)

print("""
G_T(t, x) = t^3 * P(u, x) where u = t + 1/t, and:

  p_3(x) = I(Omega(T), x)  [= H(T) at x=2]
  p_2(x) = 120 + (48*t3 - 12*t7)*x
  p_1(x) = 1188 + (24*t3 - 24*t5 + 24*t7)*x + (-48*bc)*x^2
  p_0(x) = 2176 + (-256*t3 + 32*t5 - 16*t7)*x + (64*bc)*x^2

The LEADING coefficient (in u) is p_3(x) = I(Omega, x).
This means: as |u| -> infinity (i.e., |t| -> infinity or |t| -> 0),
  G_T(t, x) ~ t^3 * I(Omega, x) * u^3 = t^3 * I(Omega, x) * (t + 1/t)^3

At u = 2 (i.e., t = 1): P(2, x) = n! for all x.
At u = 0 (i.e., t = i): P(0, x) = p_0(x). For n=7 at x=2:
  P(0, 2) = p_0(2) = 2176 + (-256*t3 + 32*t5 - 16*t7)*2 + 64*bc*4
           = 2176 - 512*t3 + 64*t5 - 32*t7 + 256*bc

This is the "imaginary height" of E_T(i)/i^3.

The polynomial P has a clear hierarchy:
  - Highest u-power (u^3 coeff): I(Omega, x), purely determined by independent sets
  - Middle u-powers: mixed cycle contributions with exact integer coefficients
  - Constant u-power (p_0): the most "complex" combination of invariants
""")

# Verify the formulas are exact
print("Exact verification:")
for i in range(min(10, len(data))):
    d_entry = data[i]
    t3, t5, t7, bc = d_entry['t3'], d_entry['t5'], d_entry['t7'], d_entry['bc']

    # p3 at x=2
    pred_p3 = 1 + 2*2*(t3+t5+t7) + 4*4*bc  # I(Omega, 2) = 1 + 4*alpha1 + 16*alpha2... wait
    # No: p3 = a_0 = H = I(Omega, 2) at x=2 means p3(2) = H, but p3 IS a_0 as a function of x.
    # Actually p3 = a_0(T) as a number (not polynomial in x).
    # Wait, I need to distinguish: p3 the NUMBER (at x=2) vs p3(x) the POLYNOMIAL.
    # In the regression above, we computed p_j as numbers (at x=2).
    # The symbolic x version: p_3(x) = 1 + 2*alpha_1*x + 4*alpha_2*x^2 = I(Omega, x).
    # At x=2: p_3(2) = 1 + 4*alpha_1 + 16*alpha_2 = H(T). Let me double-check.
    # H(T) = 1 + 2*alpha_1 + 4*alpha_2. But p_3 = a_0 = H(T).
    # And I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 = H(T).
    # So p_3 = H(T) = I(Omega, 2). But p_3(x) = I(Omega, x)? Let me check:
    # a_0(T, x) = A(7,0) + x*(2*t3 + 2*t5 + 2*t7) + x^2*(4*bc)
    #           = 1 + 2*alpha_1*x + 4*alpha_2*x^2
    # I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2
    # These are NOT the same! a_0(T, x) = 1 + 2*alpha_1*x + 4*alpha_2*x^2 != I(Omega, x)

    # The issue: in the deformed Eulerian formula, a_k(T) = I_k(Omega, 2), and
    # I_k(Omega, x) = A(n,k) + sum c_k * x^{parts} * I(T).
    # I_0(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 = I(Omega, x).
    # But the GF coefficient is: [t^0] G_T(t, x) = I_0(Omega, x) = I(Omega, x).
    # p_3(x) = [t^0] G_T(t, x) / (coefficient structure...) -- wait.

    # Actually p_3 = a_0 (as computed from the palindromic t-expansion).
    # The x-dependence: when I use G_T(t, x) instead of G_T(t, 2):
    # [t^0] G_T(t, x) = A(7,0) + sum_I x^{parts} I(T) * c_0^{(f,d)} * ...
    # But c_0^{(f,d)} = (-1)^{d-f} (since only j=0, k=0 survives).
    # And d-f is even, so c_0 = 1 for all invariants.
    # So [t^0] G_T(t, x) = 1 + sum_I x^{parts} I(T) = I(Omega, x). Yes!
    # Wait but the 2^{parts} factor... let me recheck.
    # G_T(t, x) = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t) (t-1)^{d-f}
    # [t^0] A_n(t) = A(n,0) = 1.
    # [t^0] A_{f+1}(t) (t-1)^{d-f} = A_{f+1}(0) * (-1)^{d-f} = 1 * 1 = 1.
    # So [t^0] G_T(t, x) = 1 + sum_I x^{parts} I(T) = 1 + alpha_1*x + alpha_2*x^2 = I(Omega, x). ✓

    # Note: the 2^{parts} factor is ABSENT in G_T(t, x). It appears only in a_k(T) = [t^k] G_T(t, 2).
    # So p_3(x) = I(Omega, x), and p_3(2) = I(Omega, 2) = H. ✓

    # So the formulas above for p_j(x) need to use the UN-multiplied invariant coefficients.
    # Let me redo: in G_T(t, x), the coefficient of x^1 * t^k for invariant I with f:
    # is I(T) * c_k^{(f,d)} where c_k = inflated Eulerian at level f.
    # (No 2^{parts} here — that comes from x=2 evaluation.)

    pass

print("\nCORRECTION: In G_T(t, x), there is NO 2^{parts} factor.")
print("The 2^{parts} appears only when evaluating at x=2: a_k = [t^k] G_T(t, 2) = I_k(Omega, 2).")
print("So in the x-polynomial version:")
print("  p_3(x) = 1 + alpha_1*x + alpha_2*x^2 = I(Omega, x)  [NO factor of 2!]")

# Redo with correct formulas:
# [t^k] G_T(t, x) = A(n,k) + x*(c_k^{(4,6)}*t3 + c_k^{(2,6)}*t5 + c_k^{(0,6)}*t7) + x^2*(c_k^{(2,6)}*bc)
# (no 2^parts, because 2^parts comes from the overall formula differently)

# Actually wait. Let me reread: G_T(t, x) = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t) (t-1)^{d-f}
# The coefficient of x^1 includes ALL invariants with parts=1.
# For t3: x^1 * t3 * A_5(t) * (t-1)^2
# [t^k in A_5(t)*(t-1)^2] = inflated_eulerian(4, 6, k) = c_k^{(4,6)}
# So the coefficient of t3 in [t^k] is c_k^{(4,6)}, not 2*c_k^{(4,6)}.

# Let me confirm: a_k(T) = [t^k] G_T(t, 2) = A(n,k) + 2^1 * c_k^{(4,6)} * t3 + ...
# But G_T(t, 2) = A_n(t) + 2^1 * t3 * A_5(t)*(t-1)^2 + ...
# Hmm, 2^{parts} = 2^1 for t3. But G_T(t, x) = A_n(t) + x^{parts} * I(T) * A_{f+1}(t)*(t-1)^{d-f}
# At x=2: x^1 = 2, so coeff of t3 is 2 * c_k.
# But at general x: coeff of t3 is x * c_k.
# There's no separate 2^{parts} — the 2 comes from x=2, x^1 = 2.

# OK so in the ORIGINAL formula: G_T(t, x) = A_n(t) + sum_I x^{parts(I)} * I(T) * A_{f+1}(t)*(t-1)^{d-f}
# This means: coefficient of x^j * I(T) is A_{f+1}(t)*(t-1)^{d-f} (for each I with parts=j).
# NO 2^{parts} factor. The 2^{parts} in the THM-062 formula comes from evaluating x=2.

# Good. So p_3(x) = I(Omega, x) is CORRECT.

# Redo the p_j formulas without 2^parts:
print(f"\nCORRECTED p_j(x) formulas for n=7:")
# p_3(x) = [t^0] = 1 + (t3+t5+t7)*x + bc*x^2 = I(Omega, x)
# p_2(x) = [t^1] = 120 + (24*t3 + 0*t5 - 6*t7)*x + 0*bc*x^2 = 120 + (24*t3 - 6*t7)*x
# p_1(x) = [t^2]-3*[t^0] = (1191-3) + ((15-3)*t3 + (-9-3)*t5 + (15-3)*t7)*x + ((-9-3)*bc)*x^2
#         = 1188 + (12*t3 - 12*t5 + 12*t7)*x + (-12*bc)*x^2
#         = 1188 + 12*(t3 - t5 + t7)*x - 12*bc*x^2
# p_0(x) = [t^3]-2*[t^1] = (2416-240) + ((-80-48)*t3 + (16-0)*t5 + (-20+12)*t7)*x + ((16-0)*bc)*x^2
#         = 2176 + (-128*t3 + 16*t5 - 8*t7)*x + 16*bc*x^2

p2x_t3 = ck4[1]  # 24
p2x_t5 = ck2[1]  # 0
p2x_t7 = ck0[1]  # -6
p1x_const = An[2] - 3*An[0]  # 1188
p1x_t3 = ck4[2] - 3*ck4[0]  # 15-3=12
p1x_t5 = ck2[2] - 3*ck2[0]  # -9-3=-12
p1x_t7 = ck0[2] - 3*ck0[0]  # 15-3=12
p1x_bc = ck2[2] - 3*ck2[0]  # -12
p0x_const = An[3] - 2*An[1]  # 2176
p0x_t3 = ck4[3] - 2*ck4[1]  # -80-48=-128
p0x_t5 = ck2[3] - 2*ck2[1]  # 16
p0x_t7 = ck0[3] - 2*ck0[1]  # -20+12=-8
p0x_bc = ck2[3] - 2*ck2[1]  # 16

print(f"  p_3(x) = 1 + alpha_1*x + alpha_2*x^2 = I(Omega, x)")
print(f"  p_2(x) = {An[1]} + ({p2x_t3}*t3 + {p2x_t5}*t5 + {p2x_t7}*t7)*x")
print(f"         = 120 + (24*t3 - 6*t7)*x")
print(f"  p_1(x) = {p1x_const} + ({p1x_t3}*t3 + {p1x_t5}*t5 + {p1x_t7}*t7)*x + ({p1x_bc}*bc)*x^2")
print(f"         = 1188 + 12*(t3 - t5 + t7)*x - 12*bc*x^2")
print(f"  p_0(x) = {p0x_const} + ({p0x_t3}*t3 + {p0x_t5}*t5 + {p0x_t7}*t7)*x + ({p0x_bc}*bc)*x^2")
print(f"         = 2176 + (-128*t3 + 16*t5 - 8*t7)*x + 16*bc*x^2")

# Verify against data
print(f"\nVerification at x=2:")
all_ok = True
for i in range(10):
    d_entry = data[i]
    t3, t5, t7, bc = d_entry['t3'], d_entry['t5'], d_entry['t7'], d_entry['bc']
    alpha1 = t3 + t5 + t7
    alpha2 = bc

    # p3(2) = 1 + 2*alpha1 + 4*alpha2 = H
    pred_p3 = 1 + 2*alpha1 + 4*alpha2
    # p2(2) = 120 + 2*(24*t3 - 6*t7)
    pred_p2 = 120 + 2*(24*t3 - 6*t7)
    # p1(2) = 1188 + 2*12*(t3-t5+t7) + 4*(-12)*bc
    pred_p1 = 1188 + 24*(t3-t5+t7) - 48*bc
    # p0(2) = 2176 + 2*(-128*t3+16*t5-8*t7) + 4*16*bc
    pred_p0 = 2176 + 2*(-128*t3 + 16*t5 - 8*t7) + 64*bc

    ok = (pred_p3 == d_entry['p3'] and pred_p2 == d_entry['p2'] and
          pred_p1 == d_entry['p1'] and pred_p0 == d_entry['p0'])
    if not ok:
        all_ok = False
        print(f"  FAIL at seed={i}")

print(f"  All 10 tournaments: {'PASS' if all_ok else 'FAIL'}")

print(f"\n{'=' * 70}")
print("SUMMARY: The reduced polynomial P(u, x) encodes tournaments compactly.")
print("The LEADING coefficient p_3(x) = I(Omega, x) IS the independence polynomial.")
print("=" * 70)
