#!/usr/bin/env python3
"""
gk_extend_s112.py — Extend g_k analysis: recurrences, GF, OEIS connections
kind-pasteur-2026-03-15-S112

Key structure discovered:
- g_k(m) has parity (-1)^k in m
- Even k: g_k(m) = P_{k/2}(m^2), P(0)=0
- Odd k: g_k(m) = m * Q_{(k-1)/2}(m^2), Q_k(0) = 1/(2k+1)

Goals:
1. Find recurrence in k for fixed m
2. Identify coefficient sequences in OEIS
3. Connect to known polynomial families
4. Find generating function
"""

from fractions import Fraction
from math import factorial, gcd

def transfer_gk_values(k, m_max):
    results = {}
    for m in range(0, m_max + 1):
        n = m + 2*k
        num_edges = n - 2
        k_max = (n - 1) // 2
        if k > k_max:
            continue
        state = [[Fraction(1)] + [Fraction(0)] * k_max,
                 [Fraction(0)] * (k_max + 1),
                 [Fraction(0)] * (k_max + 1)]
        for step in range(num_edges):
            A, B, C = state
            nA = [A[i] + C[i] for i in range(k_max + 1)]
            nB = [Fraction(0)] + [2*A[i] + C[i] for i in range(k_max)]
            nC = list(B)
            state = [nA, nB, nC]
        total = [state[0][i] + state[1][i] + state[2][i] for i in range(k_max + 1)]
        if k < len(total) and total[k] != 0:
            results[m] = total[k] / 2
        else:
            results[m] = Fraction(0)
    return results

# ============================================================
# PART 1: Recurrence in k for fixed m
# ============================================================
print("="*60)
print("PART 1: g_k(m) for fixed m, varying k")
print("="*60)

for m in range(1, 8):
    print(f"\nm={m}: g_k({m}) for k=1..{min(m+5, 12)}")
    seq = []
    for k in range(1, min(m + 8, 15)):
        tm = transfer_gk_values(k, m)
        if m in tm:
            seq.append(tm[m])
            print(f"  k={k}: g_{k}({m}) = {tm[m]}")
        else:
            break

    # Check for recurrence: a*g_{k+2} + b*g_{k+1} + c*g_k = 0?
    if len(seq) >= 5:
        # Try: g_{k+2}(m) = A*g_{k+1}(m) + B*g_k(m) + C
        # From 3 consecutive: solve for A,B,C
        for start in range(len(seq) - 4):
            a, b, c, d, e = seq[start:start+5]
            # c = A*b + B*a + C
            # d = A*c + B*b + C
            # e = A*d + B*c + C
            # (d-c) = A*(c-b) + B*(b-a)
            # (e-d) = A*(d-c) + B*(c-b)
            num1 = d - c
            rhs1_a = c - b
            rhs1_b = b - a
            num2 = e - d
            rhs2_a = d - c
            rhs2_b = c - b

            det = rhs1_a * rhs2_b - rhs1_b * rhs2_a
            if det != 0:
                A = (num1 * rhs2_b - num2 * rhs1_b) / det
                B = (rhs1_a * num2 - rhs2_a * num1) / det
                C = c - A*b - B*a
                # Verify
                ok = True
                for i in range(len(seq) - 2):
                    pred = A * seq[i+1] + B * seq[i] + C
                    if pred != seq[i+2]:
                        ok = False
                        break
                if ok and start == 0:
                    print(f"  RECURRENCE: g_{{k+2}}({m}) = {A}*g_{{k+1}}({m}) + {B}*g_k({m}) + {C}")
                    break

# ============================================================
# PART 2: Q_k(0) = 1/(2k+1) pattern
# ============================================================
print("\n" + "="*60)
print("PART 2: Q_k(0) = g_{2k+1}(m)/m as m->0")
print("="*60)

for k_half in range(0, 8):
    k = 2*k_half + 1
    tm = transfer_gk_values(k, 5)
    # g_k(m)/m at m=1 gives Q_k_half(1)
    # Q_k_half(0) = lim m->0 g_k(m)/m
    # Since g_k(m) = c_1*m + c_3*m^3 + ..., Q(0) = c_1
    # c_1 = (coefficient of m in D_k*g_k) / D_k

    # From the polynomial, the m coefficient is known
    # g_k(1) = 1, so Q(1) = 1
    # g_k(2) = 2k, so Q(4) = 2k/2 = k
    # Let me compute Q(0) from the polynomial

    # Use L'Hopital: g_k(m)/m -> g_k'(0).
    # For the polynomial D*g_k = a_1*m + a_3*m^3 + ..., g_k(m)/m = (a_1 + a_3*m^2 + ...)/D
    # Q(0) = a_1/D

    # We already know the constants from the previous script
    print(f"  k={k}: Q_{k_half}(0) = g_{k}(m)/m at m->0 = 1/{2*k_half+1} (conjectured)")

# ============================================================
# PART 3: Eigenvalue connection
# ============================================================
print("\n" + "="*60)
print("PART 3: F(N,x) = h(N,x) generating function")
print("="*60)

# Compute F(N,x) for several N
# F(N,x) = [1,0,0] * M(x)^N * [1,1,1]^T

def compute_F(N, max_k):
    """Compute F(N,x) = sum_{k>=0} f_k * x^k as a list of coefficients."""
    state = [[Fraction(1)] + [Fraction(0)] * max_k,
             [Fraction(0)] * (max_k + 1),
             [Fraction(0)] * (max_k + 1)]

    for step in range(N):
        A, B, C = state
        nA = [A[i] + C[i] for i in range(max_k + 1)]
        nB = [Fraction(0)] + [2*A[i] + C[i] for i in range(max_k)]
        nC = list(B)
        state = [nA, nB, nC]

    total = [state[0][i] + state[1][i] + state[2][i] for i in range(max_k + 1)]
    return total

print("F(N,x) coefficients [x^k] for N=0..15:")
for N in range(0, 16):
    F = compute_F(N, min(N//2+1, 8))
    coeffs = [int(f) for f in F if f != 0]
    # g_k(m) = F[k]/2 where m = N-2k+2
    print(f"  N={N:2d}: {coeffs}")

# ============================================================
# PART 4: The g_k table as a 2D array
# ============================================================
print("\n" + "="*60)
print("PART 4: g_k(m) table (k=1..8, m=1..8)")
print("="*60)

print(f"{'':>5}", end="")
for m in range(1, 9):
    print(f"{'m='+str(m):>10}", end="")
print()

for k in range(1, 9):
    tm = transfer_gk_values(k, 8)
    print(f"k={k:2d}:", end="")
    for m in range(1, 9):
        v = tm.get(m, Fraction(0))
        print(f"{str(v):>10}", end="")
    print()

# ============================================================
# PART 5: Check column recurrence (fixed m)
# ============================================================
print("\n" + "="*60)
print("PART 5: Column recurrence g_{k+1}(m) in terms of g_k(m), g_{k-1}(m)")
print("="*60)

# For fixed m, check if g_{k+1} = (am + b) * g_k + (cm^2 + d) * g_{k-1}
for m in [1, 2, 3, 4, 5]:
    g = {}
    for k in range(1, 12):
        tm = transfer_gk_values(k, m)
        if m in tm:
            g[k] = tm[m]

    print(f"\nm={m}: ratios g_{'{k+1}'}/g_k:")
    for k in range(1, min(9, len(g))):
        if k in g and k+1 in g and g[k] != 0:
            ratio = g[k+1] / g[k]
            print(f"  g_{k+1}/g_{k} = {g[k+1]}/{g[k]} = {float(ratio):.4f}")

# ============================================================
# PART 6: OEIS-ready sequences
# ============================================================
print("\n" + "="*60)
print("PART 6: OEIS-ready sequences")
print("="*60)

print("\ng_k(1) = 1 for all k (trivial)")

print("\ng_k(2) sequence (should be 2k):")
for k in range(1, 12):
    tm = transfer_gk_values(k, 2)
    if 2 in tm:
        print(f"  g_{k}(2) = {tm[2]}")

print("\ng_k(3) sequence:")
seq3 = []
for k in range(1, 12):
    tm = transfer_gk_values(k, 3)
    if 3 in tm:
        seq3.append(int(tm[3]))
print(f"  {seq3}")

print("\ng_k(4) sequence:")
seq4 = []
for k in range(1, 12):
    tm = transfer_gk_values(k, 4)
    if 4 in tm:
        seq4.append(int(tm[4]))
print(f"  {seq4}")

print("\ng_k(5) sequence:")
seq5 = []
for k in range(1, 10):
    tm = transfer_gk_values(k, 5)
    if 5 in tm:
        seq5.append(int(tm[5]))
print(f"  {seq5}")

# Row sums: sum_k g_k(m) for fixed m
print("\nRow sums S(m) = sum_{k>=1} g_k(m):")
for m in range(1, 8):
    total = Fraction(0)
    for k in range(1, m + 5):
        tm = transfer_gk_values(k, m)
        if m in tm:
            total += tm[m]
    print(f"  S({m}) = {total}")

# Diagonal: g_k(k) sequence
print("\nDiagonal g_k(k):")
diag = []
for k in range(1, 10):
    tm = transfer_gk_values(k, k)
    if k in tm:
        diag.append(int(tm[k]))
print(f"  {diag}")

print("\nDone!")
