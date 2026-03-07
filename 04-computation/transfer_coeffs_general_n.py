#!/usr/bin/env python3
"""
Derive transfer matrix coefficients at general n from moment formulas.

Key identity (THM-055): tr(c_{n-1-2k}) = sum_P e_{2k}(s_P)
where s_i = A[p_i, p_{i+1}] - 1/2, and e_{2k} are elementary symmetric polynomials.

Since s_i in {-1/2, +1/2}, all power sums reduce to functions of g = f - (n-1)/2:
  p_{2l} = (n-1)/4^l
  p_{2l+1} = g/4^l

Newton's identities give e_k as polynomial in g (degree k), hence tr(c_{n-1-2k})
is a linear combination of m_0, ..., m_{2k}.

opus-2026-03-06-S28
"""
from math import factorial, comb
from fractions import Fraction
import random
from itertools import permutations, combinations

# =====================================================================
# Symbolic computation of e_{2k} in terms of g = f - (n-1)/2
# =====================================================================

def compute_e2k_coefficients(n, max_k):
    """
    Compute e_{2k} as polynomial in g = f - (n-1)/2.
    Returns list of coefficients: e_{2k}(g) = sum_j c_j * g^(2j) for j=0,...,k.
    (Only even powers of g appear in e_{2k}.)

    Uses Newton's identities: k*e_k = sum_{i=1}^{k} (-1)^{i+1} * e_{k-i} * p_i
    """
    # Power sums as functions of g:
    # p_{2l} = (n-1)/4^l
    # p_{2l+1} = g/4^l

    # Represent e_k as a polynomial in g with Fraction coefficients.
    # e_k[j] = coefficient of g^j.
    N = n - 1  # number of edges = n-1
    half = Fraction(1, 2)

    def p(k):
        """Power sum p_k as polynomial in g."""
        if k % 2 == 0:
            # Constant: (n-1)/4^(k/2)
            poly = {0: Fraction(N, 4**(k//2))}
        else:
            # g/4^((k-1)/2)
            poly = {1: Fraction(1, 4**((k-1)//2))}
        return poly

    def poly_mult(a, b):
        """Multiply two polynomials (dict: degree -> coeff)."""
        result = {}
        for da, ca in a.items():
            for db, cb in b.items():
                d = da + db
                result[d] = result.get(d, Fraction(0)) + ca * cb
        return {d: c for d, c in result.items() if c != 0}

    def poly_add(a, b, sign=1):
        result = dict(a)
        for d, c in b.items():
            result[d] = result.get(d, Fraction(0)) + sign * c
        return {d: c for d, c in result.items() if c != 0}

    def poly_scale(a, s):
        return {d: c * s for d, c in a.items() if c * s != 0}

    # e_0 = 1
    e = [{0: Fraction(1)}]  # e[0] = 1

    for k_val in range(1, 2*max_k + 1):
        # k*e_k = sum_{i=1}^{k} (-1)^{i+1} * e_{k-i} * p_i
        total = {}
        for i in range(1, k_val + 1):
            term = poly_mult(e[k_val - i], p(i))
            sign = (-1)**(i + 1)
            total = poly_add(total, term, sign)
        e_k = poly_scale(total, Fraction(1, k_val))
        e.append(e_k)

    # Return e_{2k} for k=0,...,max_k
    return [e[2*k] for k in range(max_k + 1)]

# =====================================================================
# Convert e_{2k}(g) to sum_P e_{2k} using moment formulas
# =====================================================================

def moment_formula(j, n):
    """
    Return m_j = sum_P f^j as (const, t3_coeff, t5_coeff, bc_coeff).
    Uses proved closed forms for j=0,...,5.
    """
    if j == 0:
        return (factorial(n), 0, 0, 0)
    if j == 1:
        return (factorial(n) * (n-1) // 2, 0, 0, 0)
    if j == 2:
        return (factorial(n) * (3*n*n - 5*n + 4) // 12,
                4 * factorial(n-2), 0, 0)
    if j == 3:
        return (factorial(n) * (n-1) * (n*n - n + 2) // 8,
                6 * factorial(n-1), 0, 0)
    if j == 4:
        return (factorial(n) * (15*n**4 - 30*n**3 + 65*n**2 - 82*n + 48) // 240,
                2 * factorial(n-2) * (3*n*n - 5*n + 4),
                48 * factorial(n-4),
                96 * factorial(n-4))
    if j == 5:
        return (factorial(n) * (n-1) * (3*n**4 - 2*n**3 + 13*n**2 - 14*n + 16) // 96,
                5 * (n-1) * (n*n - n + 2) * factorial(n-2),
                120 * (n-1) * factorial(n-4),
                240 * (n-1) * factorial(n-4))
    return None

def compute_trace_coeff(k, n):
    """
    Compute tr(c_{n-1-2k}) in terms of (const, t3, t5, bc) at general n.

    e_{2k}(g) is a polynomial in g of degree 2k.
    sum_P e_{2k}(g_P) = sum_P e_{2k}(f_P - (n-1)/2)
    """
    e_coeffs = compute_e2k_coefficients(n, k)
    e2k = e_coeffs[k]  # polynomial in g

    # sum_P g^j = sum_P (f - (n-1)/2)^j = sum_{l=0}^{j} C(j,l) * m_l * (-(n-1)/2)^{j-l}
    # where m_l = sum_P f^l

    # For each power g^j in e_{2k}, expand and sum
    half_n1 = Fraction(n-1, 2)

    result = [Fraction(0)] * 4  # const, t3, t5, bc

    for deg, coeff in e2k.items():
        # sum_P g^deg = sum_P (f - half_n1)^deg
        # = sum_{l=0}^{deg} C(deg, l) * (-half_n1)^{deg-l} * m_l
        for l in range(deg + 1):
            m = moment_formula(l, n)
            if m is None:
                print(f"  WARNING: m_{l} not available at n={n}")
                return None
            binom = comb(deg, l)
            factor = coeff * binom * (-half_n1)**(deg - l)
            for i in range(4):
                result[i] += factor * m[i]

    return result

# =====================================================================
# Compute and verify
# =====================================================================
print("=" * 70)
print("Transfer matrix coefficients at general n")
print("=" * 70)

# k=0: tr(c_{n-1})
print("\nk=0: tr(c_{n-1})")
for n in [5, 7, 9, 11, 13]:
    r = compute_trace_coeff(0, n)
    print(f"  n={n}: {[float(x) for x in r]}")

# k=1: tr(c_{n-3})
print("\nk=1: tr(c_{n-3})")
for n in [5, 7, 9, 11, 13]:
    r = compute_trace_coeff(1, n)
    if r is None: continue
    print(f"  n={n}: const={r[0]}, t3={r[1]}")
    # Check against known: 2*(n-2)!*(t3 - C(n,3)/4)
    pred_t3 = 2 * factorial(n-2)
    pred_const = -2 * factorial(n-2) * Fraction(comb(n,3), 4)
    ok = (r[1] == pred_t3 and r[0] == pred_const)
    print(f"         expected: {float(pred_const)} + {float(pred_t3)}*t3  {'OK' if ok else 'FAIL'}")

# k=2: tr(c_{n-5})
print("\nk=2: tr(c_{n-5})")
for n in [7, 9, 11, 13]:
    r = compute_trace_coeff(2, n)
    if r is None: continue
    print(f"  n={n}: const={float(r[0]):.1f}, t3={float(r[1]):.1f}, "
          f"t5={float(r[2]):.1f}, bc={float(r[3]):.1f}")

# Check n=7 against known: 24*bc - 60*t3 + 12*t5 + 231
print("\n  n=7 check: expected 231 - 60*t3 + 12*t5 + 24*bc")
r7 = compute_trace_coeff(2, 7)
print(f"  Got: {float(r7[0]):.1f} + {float(r7[1]):.1f}*t3 + {float(r7[2]):.1f}*t5 + {float(r7[3]):.1f}*bc")

# Check n=9 against known: -4200*t3 + 240*t5 + 480*bc + 40320
print("\n  n=9 check: expected 40320 - 4200*t3 + 240*t5 + 480*bc")
r9 = compute_trace_coeff(2, 9)
print(f"  Got: {float(r9[0]):.1f} + {float(r9[1]):.1f}*t3 + {float(r9[2]):.1f}*t5 + {float(r9[3]):.1f}*bc")

# =====================================================================
# Extract closed forms
# =====================================================================
print(f"\n{'='*70}")
print("Closed forms for tr(c_{n-5})")
print(f"{'='*70}")

for n in [7, 9, 11, 13, 15, 17]:
    r = compute_trace_coeff(2, n)
    if r is None: continue

    # Normalize
    fac_n4 = factorial(n-4)
    print(f"\n  n={n}:")
    for i, name in enumerate(["const", "t3", "t5", "bc"]):
        val = r[i]
        ratio = val / fac_n4 if fac_n4 > 0 else 0
        print(f"    {name}: {float(val):.1f}, /(n-4)! = {float(ratio):.6f}")

# =====================================================================
# Brute-force verify at n=7
# =====================================================================
print(f"\n{'='*70}")
print("Brute-force verification of tr(c_{n-5}) at n=7")
print(f"{'='*70}")

def compute_transfer_trace_bf(A, n, r_val):
    """Compute tr(M(r)) using DP."""
    dp_fwd = {}
    for v in range(n):
        dp_fwd[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp_fwd.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[v][u] - 0.5)
                key = (mask | (1 << u), u)
                dp_fwd[key] = dp_fwd.get(key, 0) + val * wt
    dp_bwd = {}
    for v in range(n):
        dp_bwd[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp_bwd.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[u][v] - 0.5)
                key = (mask | (1 << u), u)
                dp_bwd[key] = dp_bwd.get(key, 0) + val * wt
    full = (1 << n) - 1
    tr_val = 0.0
    for v in range(n):
        for mask_before in range(1 << n):
            if mask_before & (1 << v): continue
            mask_with_v = mask_before | (1 << v)
            fwd = dp_fwd.get((mask_with_v, v), 0)
            if fwd == 0: continue
            mask_after = full ^ mask_before
            if not (mask_after & (1 << v)): continue
            bwd = dp_bwd.get((mask_after, v), 0)
            if bwd == 0: continue
            k = bin(mask_before).count('1')
            tr_val += ((-1)**k) * fwd * bwd
    return tr_val

import numpy as np

n = 7
r_formula = compute_trace_coeff(2, n)

for trial in range(5):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    # Compute t3, t5, bc
    t3 = sum(1 for a,b,c in combinations(range(n),3)
             if (A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]))
    t5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%5]] for i in range(5)):
                t5 += 1
    t5 //= 5
    cyc_triples = [set(t) for t in combinations(range(n), 3)
                   if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
                      A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    bc = sum(1 for i in range(len(cyc_triples)) for j in range(i+1, len(cyc_triples))
             if cyc_triples[i].isdisjoint(cyc_triples[j]))

    # Compute actual c_{n-5} = c_2 from transfer matrix
    r_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
    traces = [compute_transfer_trace_bf(A, n, r) for r in r_vals]
    # tr(M(r)) = c_0 + c_2*r^2 + c_4*r^4 + c_6*r^6
    V = np.array([[r**(2*k) for k in range(4)] for r in r_vals[:4]])
    c_actual = np.linalg.solve(V, traces[:4])
    c2_actual = c_actual[1]

    c2_pred = float(r_formula[0] + r_formula[1]*t3 + r_formula[2]*t5 + r_formula[3]*bc)

    print(f"  T{trial}: t3={t3}, t5={t5}, bc={bc}, c2_actual={c2_actual:.1f}, pred={c2_pred:.1f}, "
          f"{'OK' if abs(c2_actual-c2_pred)<0.1 else 'FAIL'}")

# =====================================================================
# General n closed form for tr(c_{n-5})
# =====================================================================
print(f"\n{'='*70}")
print("Closed form extraction for tr(c_{n-5})")
print(f"{'='*70}")

print("\n  tr(c_{n-5}) coefficients:")
for n in range(7, 22, 2):
    r = compute_trace_coeff(2, n)
    if r is None: continue
    fac = factorial(n-4)
    ratios = [r[i] / fac for i in range(4)]
    print(f"  n={n:2d}: const/(n-4)! = {float(ratios[0]):12.4f}, "
          f"t3/(n-4)! = {float(ratios[1]):10.4f}, "
          f"t5/(n-4)! = {float(ratios[2]):8.4f}, "
          f"bc/(n-4)! = {float(ratios[3]):8.4f}")

# Check bc/t5 ratio
print("\n  bc/t5 ratio for tr(c_{n-5}):")
for n in range(7, 22, 2):
    r = compute_trace_coeff(2, n)
    if r is None or r[2] == 0: continue
    ratio = r[3] / r[2]
    print(f"  n={n}: {float(ratio):.4f}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
