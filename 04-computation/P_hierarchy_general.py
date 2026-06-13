#!/usr/bin/env python3
"""
P(u,x) coefficient hierarchy at general n.

G_T(t,x) = t^m * P(u,x), u = t+1/t, m = (n-1)/2.
P(u,x) = p_0(x) + p_1(x)*u + ... + p_m(x)*u^m

Key facts (verified n=7):
  p_m(x) = I(Omega(T), x)   [independence polynomial]
  P(2,x) = n!                [t=1 evaluation]
  p_0(2) = (-2)^m * W(i/2)  [imaginary evaluation]

Questions for general n:
  1. What is p_{m-1}(x)? Can it be expressed in terms of cycle invariants?
  2. Is there a "depth" hierarchy: p_j involves invariants up to "depth" m-j?
  3. Does p_{m-1}(x) always lack the highest-order alpha term?

kind-pasteur-2026-03-07-S28
"""
from itertools import combinations, permutations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import numpy as np

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def tournament_from_bits(n, bits_int):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits_int >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

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
    return [dist.get(k, 0) for k in range(n)]

def count_directed_cycles(A, n, cl):
    if n < cl: return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp_c = [[0]*cl for _ in range(1 << cl)]
        dp_c[1][0] = 1
        for mask in range(1, 1 << cl):
            for v in range(cl):
                if not (mask & (1 << v)) or dp_c[mask][v] == 0: continue
                for u in range(cl):
                    if mask & (1 << u): continue
                    if sub[v][u]: dp_c[mask | (1 << u)][u] += dp_c[mask][v]
        full = (1 << cl) - 1
        total += sum(dp_c[full][v] for v in range(1, cl) if sub[v][0])
    return total

def get_ck(f, d):
    """Get c_k^{(f,d)} coefficients for k = 0, ..., d."""
    result = []
    for k in range(d + 1):
        total = 0
        for j in range(max(0, k - (d - f)), min(f, k) + 1):
            sign = (-1) ** (d - f - k + j)
            total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
        result.append(total)
    return result


# =============================================
# n = 5 (m = 2)
# =============================================
print("=" * 70)
print("n = 5: P(u,x) hierarchy")
print("=" * 70)

n = 5
m_half = (n-1)//2  # = 2
m_edges = n*(n-1)//2  # = 10
d = n - 1  # = 4

An5 = [eulerian_number(5, k) for k in range(5)]
print(f"A(5,k) = {An5}")

# n=5 invariants: t3 (f=2, parts=1), t5 (f=0, parts=1)
ck_t3 = get_ck(2, 4)  # f=2 for 3-cycles at n=5
ck_t5 = get_ck(0, 4)  # f=0 for 5-cycles at n=5

print(f"c_k^(2,4) for t3 = {ck_t3}")
print(f"c_k^(0,4) for t5 = {ck_t5}")

# a_k = A(5,k) + 2*ck_t3[k]*t3 + 2*ck_t5[k]*t5
# p2 = a_0, p1 = a_1, p0 = a_2 - 2*a_0

# p2(x) = A(5,0) + (2*ck_t3[0]*t3 + 2*ck_t5[0]*t5)*x = 1 + 2*(t3+t5)*x = I(Omega,x)
print(f"\np2(x) = 1 + 2*(t3+t5)*x = I(Omega(T), x)")

# p1(x) = A(5,1) + (2*ck_t3[1]*t3 + 2*ck_t5[1]*t5)*x
p1_c = An5[1]
p1_t3 = 2*ck_t3[1]
p1_t5 = 2*ck_t5[1]
print(f"p1(x) = {p1_c} + ({p1_t3}*t3 + {p1_t5}*t5)*x")

# p0(x) = a_2(x) - 2*a_0(x) = [A(5,2) + (2*ck_t3[2]*t3 + 2*ck_t5[2]*t5)*x] - 2*[1 + 2*(t3+t5)*x]
p0_c = An5[2] - 2*An5[0]
p0_t3 = 2*ck_t3[2] - 2*2*ck_t3[0]
p0_t5 = 2*ck_t5[2] - 2*2*ck_t5[0]
print(f"p0(x) = {p0_c} + ({p0_t3}*t3 + {p0_t5}*t5)*x")

# Verify exhaustively
print(f"\nExhaustive verification ({2**m_edges} tournaments):")
all_ok = True
for bits in range(2**m_edges):
    A = tournament_from_bits(n, bits)
    a = forward_edge_dist_dp(A, n)
    t3 = count_directed_cycles(A, n, 3)
    t5 = count_directed_cycles(A, n, 5)

    # p2 = a[0], p1 = a[1], p0 = a[2] - 2*a[0]
    actual_p2 = a[0]
    actual_p1 = a[1]
    actual_p0 = a[2] - 2*a[0]

    pred_p2 = 1 + 2*t3 + 2*t5  # a_0 at x=2 = H(T)
    pred_p1 = p1_c + p1_t3*t3 + p1_t5*t5  # a_1 at x=2 (c_k already has factor 2)
    pred_p0 = p0_c + p0_t3*t3 + p0_t5*t5  # a_2 - 2*a_0 at x=2

    if actual_p2 != pred_p2 or actual_p1 != pred_p1 or actual_p0 != pred_p0:
        all_ok = False
        print(f"  FAIL at bits={bits}: t3={t3}, t5={t5}")
        print(f"    a={a}")
        print(f"    p2: {actual_p2} vs {pred_p2}")
        print(f"    p1: {actual_p1} vs {pred_p1}")
        print(f"    p0: {actual_p0} vs {pred_p0}")
        break

if all_ok:
    print(f"  ALL {2**m_edges} tournaments MATCH!")

# =============================================
# n = 7 (m = 3) — from opus's results
# =============================================
print(f"\n{'=' * 70}")
print("n = 7: P(u,x) hierarchy")
print("=" * 70)

n = 7
d = n - 1  # = 6
An7 = [eulerian_number(7, k) for k in range(7)]
ck4 = get_ck(4, 6)  # t3: f=4
ck2 = get_ck(2, 6)  # t5: f=2, bc: f=2
ck0 = get_ck(0, 6)  # t7: f=0

print(f"A(7,k) = {An7}")
print(f"c_k^(4,6) for t3 = {ck4}")
print(f"c_k^(2,6) for t5,bc = {ck2}")
print(f"c_k^(0,6) for t7 = {ck0}")

# p3(x) = I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2
# where alpha_1 = t3+t5+t7, alpha_2 = bc
print(f"\np3(x) = 1 + (t3+t5+t7)*x + bc*x^2 = I(Omega(T), x)")

# p2(x) = a_1(x) = A(7,1) + (2*ck4[1]*t3 + 2*ck2[1]*t5 + 2*ck0[1]*t7)*x + (4*ck2[1]*bc)*x^2
p2_c = An7[1]
p2_t3 = 2*ck4[1]
p2_t5 = 2*ck2[1]
p2_t7 = 2*ck0[1]
p2_bc = 4*ck2[1]
print(f"p2(x) = {p2_c} + ({p2_t3}*t3 + {p2_t5}*t5 + {p2_t7}*t7)*x + ({p2_bc}*bc)*x^2")
print(f"       = 120 + (48*t3 - 12*t7)*x + 0*bc*x^2")
print(f"       => p2(x) is LINEAR in x! (bc coefficient = 0)")

# p1(x)
p1_c = An7[2] - 3*An7[0]
p1_t3 = 2*ck4[2] - 3*2*ck4[0]
p1_t5 = 2*ck2[2] - 3*2*ck2[0]
p1_t7 = 2*ck0[2] - 3*2*ck0[0]
p1_bc = 4*ck2[2] - 3*4*ck2[0]
print(f"\np1(x) = {p1_c} + ({p1_t3}*t3 + {p1_t5}*t5 + {p1_t7}*t7)*x + ({p1_bc}*bc)*x^2")

# p0(x)
p0_c = An7[3] - 2*An7[1]
p0_t3 = 2*ck4[3] - 2*2*ck4[1]
p0_t5 = 2*ck2[3] - 2*2*ck2[1]
p0_t7 = 2*ck0[3] - 2*2*ck0[1]
p0_bc = 4*ck2[3] - 2*4*ck2[1]
print(f"p0(x) = {p0_c} + ({p0_t3}*t3 + {p0_t5}*t5 + {p0_t7}*t7)*x + ({p0_bc}*bc)*x^2")

# =============================================
# KEY OBSERVATION: p_{m-1}(x) linearity
# =============================================
print(f"\n{'=' * 70}")
print("KEY OBSERVATION: p_{m-1}(x) linearity")
print("=" * 70)

print(f"""
n=5 (m=2): p_{m_half-1}(x) = p_1(x) = {An5[1]} + ({p1_t3}*t3 + {p1_t5}*t5)*x
  -> LINEAR in x (no alpha_2 term, because alpha_2=0 at n=5 trivially)

n=7 (m=3): p_{{m-1}}(x) = p_2(x) = 120 + (48*t3 - 12*t7)*x
  -> LINEAR in x! The bc (alpha_2) coefficient is EXACTLY ZERO.

WHY? Because c_1^(2,6) = 0 for the f=2 invariants (t5 and bc).
The c_k^(f,d) coefficients at k=1 (which determine p_{{m-1}}) involve:
  c_1^(f,d) = sum_j (-1)^{{d-f-1+j}} * A(f+1,j) * C(d-f, 1-j)
For 1-j >= 0 we need j=0 or j=1.
  j=0: (-1)^{{d-f-1}} * 1 * C(d-f, 1) = (-1)^{{d-f-1}} * (d-f)
  j=1: (-1)^{{d-f}} * A(f+1, 1) * 1 = (-1)^{{d-f}} * (2^{{f+1}} - f - 2)

For n=7, d=6:
  t3 (f=4): c_1 = (-1)^1 * 2 + (-1)^2 * (2^5-6) = -2 + 26 = 24
  t5 (f=2): c_1 = (-1)^3 * 4 + (-1)^4 * (2^3-4) = -4 + 4 = 0  [!]
  t7 (f=0): c_1 = (-1)^5 * 6 + (-1)^6 * (2^1-2) = -6 + 0 = -6

So c_1^(2,6) = 0 for f=2. Since bc also has f=2, its p_2 coefficient is 4*0 = 0.
This is NOT a coincidence — it happens because A(3,1) = 2^3-4 = 4,
and C(4,1) = 4, so the two terms cancel exactly!

GENERAL PATTERN: c_1^(f,d) = 0 when d-f is even and 2^{{f+1}}-f-2 = d-f.
At n=7: d-f=4, need 2^3-4 = 4 = d-f. YES! So f=2, d=6 gives c_1=0.

Does this generalize? At n=9 (d=8, m=4):
""")

# Check at n=9
d9 = 8
for f in range(0, 9, 2):  # even f values
    c1 = 0
    for j in range(2):
        sign = (-1) ** (d9 - f - 1 + j)
        a_val = eulerian_number(f+1, j) if j <= f else 0
        binom_val = comb(d9-f, 1-j)
        c1 += sign * a_val * binom_val
    print(f"  n=9: c_1^({f},{d9}) = {c1}")

print(f"\nAt n=9, the f-values for invariants are:")
print(f"  t3: f=6, t5: f=4, t7: f=2, t9: f=0")
print(f"  bc33: f=4, bc35: f=2, a3: f=2")
print(f"  So p_3(x) (= p_{{m-1}} at n=9) has:")
print(f"    t3 coeff: 2*c_1^(6,8)")
print(f"    t5 coeff: 2*c_1^(4,8)")
print(f"    t7 coeff: 2*c_1^(2,8)")
print(f"    t9 coeff: 2*c_1^(0,8)")
print(f"    bc33 coeff: 4*c_1^(4,8)")
print(f"    bc35 coeff: 4*c_1^(2,8)")
print(f"    a3 coeff: 8*c_1^(2,8)")

print(f"\n  If c_1^(2,8) = 0, then t7, bc35, a3 all DROP OUT of p_{{m-1}}(x)!")
print(f"  The f=2 group is invisible to p_{{m-1}} just like it was at n=7.")
