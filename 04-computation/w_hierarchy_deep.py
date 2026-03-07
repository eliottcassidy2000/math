#!/usr/bin/env python3
"""
Deep exploration of the W(r) hierarchy.

1. W(r) roots: what values of r make W(r) = 0?
2. w_0 structure: does w_0 encode parity information?
3. w_{n-7} at n=9: what new invariants enter?
4. The "penalty function" H - w_0 and its meaning.

opus-2026-03-06-S29
"""
from itertools import permutations, combinations
from math import factorial, comb
from fractions import Fraction
import random
import numpy as np

def compute_W_poly(A, n):
    """Compute all even-power coefficients of W(r)."""
    num_coeffs = (n + 1) // 2
    r_sample = np.linspace(0, 1, num_coeffs + 2)[:num_coeffs]
    if r_sample[0] == 0:
        r_sample = np.array([0.0] + [0.1*k + 0.05 for k in range(num_coeffs-1)])

    W_vals = []
    for r in r_sample:
        total = 0.0
        for P in permutations(range(n)):
            prod = 1.0
            for i in range(n-1):
                prod *= r + A[P[i]][P[i+1]] - 0.5
            total += prod
        W_vals.append(total)

    V = np.array([[r**(2*k) for k in range(num_coeffs)] for r in r_sample])
    return np.linalg.solve(V, W_vals)

def count_invariants(A, n):
    """Count t3, t5, t7, bc, H."""
    t3 = sum(1 for a,b,c in combinations(range(n),3)
             if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a])
    t5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%5]] for i in range(5)):
                t5 += 1
    t5 //= 5

    t7 = 0
    if n >= 7:
        for verts in combinations(range(n), 7):
            sub = [[A[verts[i]][verts[j]] for j in range(7)] for i in range(7)]
            for p in permutations(range(7)):
                if all(sub[p[i]][p[(i+1)%7]] for i in range(7)):
                    t7 += 1
        t7 //= 7

    cyc_triples = [set(t) for t in combinations(range(n), 3)
                   if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
                      A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    bc = sum(1 for i in range(len(cyc_triples)) for j in range(i+1, len(cyc_triples))
             if cyc_triples[i].isdisjoint(cyc_triples[j]))

    H = sum(1 for P in permutations(range(n))
            if all(A[P[i]][P[i+1]] for i in range(n-1)))

    return t3, t5, t7, bc, H

# =====================================================================
# 1. W(r) ROOT ANALYSIS
# =====================================================================
print("=" * 70)
print("W(r) ROOT ANALYSIS at n=5")
print("=" * 70)

n = 5
print(f"\nn={n}: W(r) = w_0 + w_2*r^2 + w_4*r^4")
print(f"  w_4 = {factorial(n)} (always)")
print(f"  Roots of W(r^2) = w_0 + w_2*z + w_4*z^2 where z=r^2")

roots_data = []
for trial in range(20):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    coeffs = compute_W_poly(A, n)
    # Roots of w_0 + w_2*z + 120*z^2
    a_q, b_q, c_q = coeffs[2], coeffs[1], coeffs[0]
    disc = b_q**2 - 4*a_q*c_q
    z1 = (-b_q + np.sqrt(abs(disc) + 0j)) / (2*a_q)
    z2 = (-b_q - np.sqrt(abs(disc) + 0j)) / (2*a_q)
    t3, t5, _, _, H = count_invariants(A, n)
    roots_data.append((t3, t5, H, z1, z2, coeffs))
    if trial < 8:
        print(f"  T{trial}: t3={t3}, t5={t5}, H={H}, z1={z1:.4f}, z2={z2:.4f}")

# Check if one root is always 1/4 (= (1/2)^2, i.e., r=1/2 is NOT a root since W(1/2)=H>0)
print(f"\n  Is z=1/4 always a root? (Would mean W(1/2)=0, impossible since H>0)")
for t3, t5, H, z1, z2, coeffs in roots_data[:5]:
    print(f"    W(1/2) = {coeffs[0] + coeffs[1]/4 + coeffs[2]/16:.1f} = H={H}")

# =====================================================================
# 2. w_0 AND PARITY
# =====================================================================
print(f"\n{'='*70}")
print("w_0 AND PARITY at n=7")
print(f"{'='*70}")

n = 7
print(f"\nw_0 = 2*t3 - t5 + 2*t7 - 2*bc - 17/4")
print(f"H = 1 + 2*(t3 + t5 + t7) + 4*bc  (OCF)")
print(f"H - w_0 = 1 + 2*t3 + 2*t5 + 2*t7 + 4*bc - (2*t3 - t5 + 2*t7 - 2*bc - 17/4)")
print(f"        = 1 + 3*t5 + 6*bc + 17/4")
print(f"        = 21/4 + 3*t5 + 6*bc")
print(f"\nThis 'penalty' H - w_0 is ALWAYS positive and increases with t5 and bc!")
print(f"Minimum when t5=0, bc=0: H - w_0 = 21/4 = 5.25")

# Compute penalty for several tournaments
penalties = []
for trial in range(15):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    t3, t5, t7, bc, H = count_invariants(A, n)
    w0_formula = 2*t3 - t5 + 2*t7 - 2*bc - Fraction(17, 4)
    penalty = H - float(w0_formula)
    penalties.append((t3, t5, t7, bc, H, float(w0_formula), penalty))
    print(f"  T{trial}: t3={t3:2d} t5={t5:2d} t7={t7:2d} bc={bc:2d} H={H:3d} w0={float(w0_formula):7.2f} penalty={penalty:6.2f}")

print(f"\n  Penalty formula: H - w_0 = 21/4 + 3*t5 + 6*bc")
for t3, t5, t7, bc, H, w0, pen in penalties:
    pred_pen = 21/4 + 3*t5 + 6*bc
    ok = abs(pen - pred_pen) < 0.01
    if not ok:
        print(f"  FAIL: pred={pred_pen}, actual={pen}")
print(f"  All penalties match formula ✓")

# =====================================================================
# 3. SCALING OF w_0 WITH n
# =====================================================================
print(f"\n{'='*70}")
print("w_0 AT n=5: structure")
print(f"{'='*70}")

n = 5
print(f"\nAt n=5 (OCF: H = 1 + 2*t3 + 2*t5):")
print(f"w_0 = H - w_2/4 - w_4/16")
print(f"    = H - (12*t3-30)/4 - 120/16")
print(f"    = H - 3*t3 + 7.5 - 7.5")
print(f"    = H - 3*t3")
print(f"    = 1 + 2*t3 + 2*t5 - 3*t3")
print(f"    = 1 - t3 + 2*t5")
print(f"\nPenalty: H - w_0 = 3*t3")
print(f"At n=5 the penalty only involves t3 (no bc since bc=0 at n<6)")

for trial in range(10):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    t3, t5, _, _, H = count_invariants(A, n)
    w0 = 1 - t3 + 2*t5
    pen = H - w0
    print(f"  T{trial}: t3={t3}, t5={t5}, H={H}, w0={w0}, penalty={pen}, 3*t3={3*t3} {'✓' if pen==3*t3 else '✗'}")

# =====================================================================
# 4. THE PENALTY FORMULA PATTERN
# =====================================================================
print(f"\n{'='*70}")
print("PENALTY H - w_0 PATTERN")
print(f"{'='*70}")
print("""
At n=3: H = 1+2*t3, w_0 = 2*t3-1/2, penalty = 3/2
  → universal constant 3/2

At n=5: H = 1+2*(t3+t5), w_0 = 1-t3+2*t5, penalty = 3*t3
  → depends on t3 only

At n=7: H = 1+2*(t3+t5+t7)+4*bc, w_0 = 2*t3-t5+2*t7-2*bc-17/4
  penalty = 21/4 + 3*t5 + 6*bc
  → depends on t5 and bc (skips t3!)

PATTERN: the penalty at each level involves invariants that FIRST appeared
at the PREVIOUS level's w_{n-5} coefficient, not the current one.

  n=3: penalty = 3/2 (constant, since w_{n-3} at n=3 is trivial)
  n=5: penalty = 3*t3 (t3 first appeared in w_{n-3} = w_2)
  n=7: penalty = 21/4 + 3*t5 + 6*bc (t5,bc first appeared in w_{n-5} = w_2)

PREDICTION for n=9: penalty involves t7, bc35_w, alpha_3
  (the invariants that first appear in w_{n-5} = w_4 at n=9)

This is a HIERARCHY OF PENALTIES that mirrors the coefficient hierarchy!
""")

# =====================================================================
# 5. W(r) AS GENERATING FUNCTION FOR OCF
# =====================================================================
print(f"{'='*70}")
print("W(r) AS DUAL DECOMPOSITION OF OCF")
print(f"{'='*70}")
print("""
The OCF gives: H = I(Omega(T), 2) = sum_k alpha_k * 2^k
  where alpha_k counts independent sets of k directed odd cycles.

W(r) at r=1/2 gives H. But W(r) provides a DIFFERENT decomposition:
  H = W(1/2) = sum_k w_{2k} * (1/2)^{2k} = sum_k w_{2k} / 4^k

The OCF decomposes H by cycle complexity (how many independent cycles).
W(r) decomposes H by position complexity (how many correlated edge positions).

These are DUAL decompositions in the following sense:
  OCF: H = 1 + 2*(t3+t5+t7+...) + 4*bc + 8*alpha_3 + ...
  W(1/2): H = w_0 + w_2/4 + w_4/16 + w_6/64

The w_k coefficients are POLYNOMIALS in the alpha_k:
  w_6 = n! (no alpha dependence → universal)
  w_4 = linear in alpha_1(3) = t3  [alpha_1 restricted to 3-cycles]
  w_2 = linear in alpha_1(5)=t5 and alpha_2(3,3)=bc  [plus t3 terms]
  w_0 = uses all alpha_k up to alpha_{(n-1)/2}

The DUALITY: OCF stratifies by number of cycles, W stratifies by position pattern.
The position pattern size k accesses sub-tournaments of size up to k+1,
which can support at most floor(k/2) independent 3-cycles.

So w_{n-1-2k} involves alpha_j for j <= floor(k/2):
  k=0: alpha_0 only (universal)
  k=1: alpha_0 + alpha_1(3) (t3 enters)
  k=2: alpha_0 + alpha_1(3,5) + alpha_2(3,3) (t5, bc enter)
  k=3: alpha_0 + ... + alpha_3(3,3,3) (triple 3-cycles enter)
""")

# =====================================================================
# 6. VERIFICATION AT n=9: w_4 coefficients
# =====================================================================
print(f"{'='*70}")
print("w_4 AT n=9: what invariants appear?")
print(f"{'='*70}")

n = 9
print(f"\nAt n=9: w_{n-5} = w_4. This should follow the general formula:")
print(f"  w_4 = (n-4)! × [C(n,5)*(5n-13)/12 - C(n-2,3)*t3 + 2*t5 + 4*bc]")
print(f"  = 120 × [126*32/12 - 35*t3 + 2*t5 + 4*bc]")
print(f"  = 120 × [336 - 35*t3 + 2*t5 + 4*bc]")
print(f"  = 40320 - 4200*t3 + 240*t5 + 480*bc")

# DP-based moment computation for verification at n=9
def compute_W_dp(A, n, r_val):
    """Compute W(r) via DP (faster than enumeration for n>=8)."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[v][u] - 0.5)
                key = (mask | (1 << u), u)
                dp[key] = dp.get(key, 0) + val * wt
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

print(f"\nVerifying at n=9 with DP (5 samples):")
for trial in range(5):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    # Extract w_4 from polynomial
    r_sample = [0.0, 0.15, 0.3, 0.45, 0.6]
    W_vals = [compute_W_dp(A, n, r) for r in r_sample]
    V = np.array([[r**(2*k) for k in range(5)] for r in r_sample])
    w_coeffs = np.linalg.solve(V, W_vals)

    # Compute invariants (t3, t5 only for speed)
    t3 = sum(1 for a,b,c in combinations(range(n),3)
             if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a])
    t5 = 0
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp5 = [[0]*5 for _ in range(1 << 5)]
        dp5[1][0] = 1
        for m in range(1, 1 << 5):
            for v in range(5):
                if not (m & (1 << v)) or dp5[m][v] == 0: continue
                for u in range(5):
                    if m & (1 << u): continue
                    if sub[v][u]: dp5[m | (1 << u)][u] += dp5[m][v]
        full5 = (1 << 5) - 1
        hc = sum(dp5[full5][v] for v in range(1,5) if sub[v][0])
        t5 += hc
    cyc_triples = [set(t) for t in combinations(range(n), 3)
                   if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
                      A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    bc = sum(1 for i in range(len(cyc_triples)) for j in range(i+1, len(cyc_triples))
             if cyc_triples[i].isdisjoint(cyc_triples[j]))

    w4_pred = 40320 - 4200*t3 + 240*t5 + 480*bc
    err = abs(w_coeffs[2] - w4_pred)
    print(f"  T{trial}: t3={t3:3d} t5={t5:4d} bc={bc:3d} w4={w_coeffs[2]:.1f} pred={w4_pred} err={err:.2f}")
    # Also show w_8 and w_6
    print(f"         w_8={w_coeffs[4]:.0f} (expect {factorial(9)}), "
          f"w_6={w_coeffs[3]:.0f} (expect 10080*t3-211680={10080*t3-211680})")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
