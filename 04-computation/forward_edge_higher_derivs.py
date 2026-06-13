#!/usr/bin/env python3
"""
Higher derivatives of W(r) at r=1/2 — extracting a_{n-k} formulas.

KEY RESULTS:
  W'(1/2) = (n-1)*H + a_{n-2}     [1 backward edge allowed]
  W''(1/2) = 2*[C(n-1,2)*H + (n-2)*a_{n-2} + a_{n-3}]   [factor of 2!]

  F'_f(1/2) = 2^{f+1} - 2
  F''_f(1/2) = f(f-1) + 2(f-1)(2^{f+1}-f-2) + 2*A(f+1,2)

At n=7:
  a_{n-2} = 120 + 48*t3 - 12*t7
  a_{n-3} = 1191 + 30*t3 - 18*t5 + 30*t7 - 36*bc  [PREDICTED]

opus-2026-03-07-S32
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb, factorial
import random

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

def count_H(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            c = dp.get((mask, v), 0)
            if c == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_forward_edge_dist(A, n):
    """DP: count permutations by number of forward edges."""
    # dp[(mask, v, k)] = # partial perms ending at v with k forward edges
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for k in range(n):
                c = dp.get((mask, v, k), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    new_k = k + A[v][u]
                    key = (mask | (1 << u), u, new_k)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for k in range(n):
            dist[k] += dp.get((full, v, k), 0)
    return dict(dist)

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cycle_len):
    """Count directed cycles of given length."""
    if n < cycle_len: return 0
    total = 0
    for verts in combinations(range(n), cycle_len):
        sub = [[A[verts[i]][verts[j]] for j in range(cycle_len)] for i in range(cycle_len)]
        dp = [[0]*cycle_len for _ in range(1 << cycle_len)]
        dp[1][0] = 1
        for m in range(1, 1 << cycle_len):
            for v in range(cycle_len):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cycle_len):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cycle_len) - 1
        total += sum(dp[full][v] for v in range(1, cycle_len) if sub[v][0])
    return total

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def eulerian_number(n, k):
    """A(n,k) = number of permutations of [n] with k descents."""
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

# ==================================================================
# Part 1: Verify F''_f(1/2) formula
# ==================================================================
print("=" * 70)
print("F''_f(1/2) = f(f-1) + 2(f-1)(2^{f+1}-f-2) + 2*A(f+1,2)")
print("=" * 70)

for f in range(8):
    # Theoretical
    if f >= 2:
        theory = f*(f-1) + 2*(f-1)*(2**(f+1)-f-2) + 2*eulerian_number(f+1, 2)
    elif f == 1:
        theory = 0  # f(f-1)=0, A(2,2)=0
    else:
        theory = 0

    # Numerical
    def eval_Ff(r, f=f):
        total = 0
        for sigma in permutations(range(f+1)):
            prod = 1.0
            for i in range(f):
                eps = 1 if sigma[i+1] > sigma[i] else -1
                prod *= r + eps*0.5
            total += prod
        return total

    eps_h = 1e-5
    d2_num = (eval_Ff(0.5+eps_h) - 2*eval_Ff(0.5) + eval_Ff(0.5-eps_h)) / eps_h**2

    print(f"  f={f}: theory={theory}, numerical={d2_num:.1f}, A(f+1,2)={eulerian_number(f+1,2)}")

# ==================================================================
# Part 2: Theoretical a_{n-3} formula at n=7
# ==================================================================
print(f"\n{'=' * 70}")
print("THEORETICAL a_{n-3} at n=7")
print("=" * 70)

n = 7
# W''(1/2) = F''_6(1/2) + 2*F''_4(1/2)*t3 + 2*F''_2(1/2)*t5 + 2*F''_0(1/2)*t7 + 4*F''_2(1/2)*bc
# = 3612 + 600*t3 + 24*t5 + 0*t7 + 48*bc

F6pp = 6*5 + 2*5*(2**7-8) + 2*eulerian_number(7, 2)
F4pp = 4*3 + 2*3*(2**5-6) + 2*eulerian_number(5, 2)
F2pp = 2*1 + 2*1*(2**3-4) + 2*eulerian_number(3, 2)
F0pp = 0

print(f"  F''_6(1/2) = {F6pp}")
print(f"  F''_4(1/2) = {F4pp}")
print(f"  F''_2(1/2) = {F2pp}")
print(f"  F''_0(1/2) = {F0pp}")

print(f"\n  W''(1/2) = {F6pp} + {2*F4pp}*t3 + {2*F2pp}*t5 + {2*F0pp}*t7 + {4*F2pp}*bc")
print(f"  W''(1/2)/2 = {F6pp//2} + {F4pp}*t3 + {F2pp}*t5 + {F0pp}*t7 + {2*F2pp}*bc")

# a_{n-3} = W''(1/2)/2 - C(n-1,2)*H - (n-2)*a_{n-2}
# C(6,2) = 15, n-2 = 5
# H = 1 + 2t3 + 2t5 + 2t7 + 4bc
# a_{n-2} = 120 + 48t3 - 12t7

print(f"\n  15*H = 15 + 30*t3 + 30*t5 + 30*t7 + 60*bc")
print(f"  5*a_{{n-2}} = 600 + 240*t3 - 60*t7")
print(f"\n  a_{{n-3}} = W''(1/2)/2 - 15*H - 5*a_{{n-2}}")

# Constant: 1806 - 15 - 600 = 1191 = A(7,2)!
const = F6pp//2 - comb(n-1, 2) - (n-2)*120
t3_coeff = F4pp - 30 - (n-2)*48
t5_coeff = F2pp - 30
t7_coeff = F0pp - 30 + (n-2)*12
bc_coeff = 2*F2pp - 60

print(f"  = {const} + {t3_coeff}*t3 + {t5_coeff}*t5 + {t7_coeff}*t7 + {bc_coeff}*bc")
print(f"\n  Note: constant = {const} = A(7,2) = {eulerian_number(7,2)}  (Eulerian number!)")

# ==================================================================
# Part 3: Numerical verification at n=7
# ==================================================================
print(f"\n{'=' * 70}")
print("NUMERICAL VERIFICATION at n=7")
print("=" * 70)

num_samples = 30
all_match = True
for trial in range(num_samples):
    A = random_tournament(n, n*300 + trial)
    dist = count_forward_edge_dist(A, n)

    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)

    a_n2_actual = dist.get(n-2, 0)
    a_n3_actual = dist.get(n-3, 0)

    a_n2_pred = 120 + 48*t3 - 12*t7
    a_n3_pred = const + t3_coeff*t3 + t5_coeff*t5 + t7_coeff*t7 + bc_coeff*bc

    match2 = (a_n2_actual == a_n2_pred)
    match3 = (a_n3_actual == a_n3_pred)

    if trial < 10 or not match2 or not match3:
        print(f"  T{trial:2d}: t3={t3:3d} t5={t5:4d} t7={t7:5d} bc={bc:3d}  "
              f"a_{{n-2}}={a_n2_actual:5d}({'OK' if match2 else 'FAIL'})  "
              f"a_{{n-3}}={a_n3_actual:5d}({'OK' if match3 else 'FAIL':>4s}) pred={a_n3_pred}")

    if not match2 or not match3:
        all_match = False

print(f"\n  All {num_samples} samples: {'PASS' if all_match else 'FAIL'}")

# ==================================================================
# Part 4: General formula for a_{n-2} and a_{n-3}
# ==================================================================
print(f"\n{'=' * 70}")
print("GENERAL FORMULA STRUCTURE")
print("=" * 70)

print("""
a_{n-1}(T) = H(T)  (Hamiltonian paths)

a_{n-2}(T) = F'_{n-1}(1/2) + sum_I 2^parts * F'_{f_I}(1/2) * I(T) - (n-1)*H(T)

  where F'_f(1/2) = 2^{f+1} - 2.

a_{n-3}(T) = W''(1/2)/2 - C(n-1,2)*H(T) - (n-2)*a_{n-2}(T)

  where W''(1/2) = F''_{n-1}(1/2) + sum_I 2^parts * F''_{f_I}(1/2) * I(T)
  and   F''_f(1/2) = f(f-1) + 2(f-1)(2^{f+1}-f-2) + 2*A(f+1,2).

PATTERN: The k-th derivative of W at r=1/2 involves F^{(k)}_f(1/2),
which depends on Eulerian numbers A(f+1, j) for j <= k.

Each a_{n-1-k}(T) is expressible as an OCF polynomial in the invariants
{t3, t5, t7, ..., bc, bc35, ...}.
""")

# ==================================================================
# Part 5: a_{n-2} at n=5 and n=9 for comparison
# ==================================================================
print(f"{'=' * 70}")
print("a_{n-2} AT DIFFERENT n VALUES")
print("=" * 70)

for n in [5, 9]:
    print(f"\nn={n}:")
    # At n=5: invariants are just t3, t5
    # f_t3 = n-1-2 = n-3, f_t5 = n-1-4 = n-5
    # W'(1/2) = F'_{n-1}(1/2) + 2*F'_{n-3}(1/2)*t3 + [if n>=6: 2*F'_{n-5}(1/2)*t5 + ...]
    # a_{n-2} = W'(1/2) - (n-1)*H

    Fp_n1 = 2**n - 2
    print(f"  F'_{n-1}(1/2) = 2^{n} - 2 = {Fp_n1}")

    if n == 5:
        # H = 1 + 2*t3 + 2*t5 at n=5
        # F'_2(1/2) = 6, F'_0(1/2) = 0
        # W'(1/2) = 30 + 12*t3 + 0*t5
        # a_{n-2} = 30 + 12*t3 - 4*(1 + 2*t3 + 2*t5) = 26 + 4*t3 - 8*t5
        Fp_2 = 2**3 - 2  # F'_2(1/2) = 6
        Fp_0 = 2**1 - 2  # F'_0(1/2) = 0
        print(f"  H = 1 + 2*t3 + 2*t5")
        print(f"  W'(1/2) = {Fp_n1} + {2*Fp_2}*t3 + {2*Fp_0}*t5")
        print(f"  a_{{n-2}} = {Fp_n1} + {2*Fp_2}*t3 + {2*Fp_0}*t5 - {n-1}*(1 + 2*t3 + 2*t5)")
        pred_const = Fp_n1 - (n-1)
        pred_t3 = 2*Fp_2 - 2*(n-1)
        pred_t5 = 2*Fp_0 - 2*(n-1)
        print(f"         = {pred_const} + {pred_t3}*t3 + {pred_t5}*t5")

    if n == 9:
        # Invariants: t3 (f=6), t5 (f=4), t7 (f=2), t9 (f=0), bc (f=4)
        # bc35 (f=2), bc37 (f=0), a3 (f=2)
        Fp = {6: 2**7-2, 4: 2**5-2, 2: 2**3-2, 0: 2**1-2}
        print(f"  F' values: {Fp}")
        print(f"  W'(1/2) = {Fp_n1} + 2*{Fp[6]}*t3 + 2*{Fp[4]}*t5 + 2*{Fp[2]}*t7 + 2*{Fp[0]}*t9 + 4*{Fp[4]}*bc + 4*{Fp[2]}*bc35 + 4*{Fp[0]}*bc37 + 4*{Fp[2]}*a3")
        print(f"  H = 1 + 2*t3 + 2*t5 + 2*t7 + 2*t9 + 4*bc + 4*bc35 + 4*bc37 + 8*a3")
        # a_{n-2} = W'(1/2) - 8*H
        # Each term: coeff_W' - 8*coeff_H
        c0 = Fp_n1 - 8*1
        ct3 = 2*Fp[6] - 8*2
        ct5 = 2*Fp[4] - 8*2
        ct7 = 2*Fp[2] - 8*2
        ct9 = 2*Fp[0] - 8*2
        cbc = 4*Fp[4] - 8*4
        cbc35 = 4*Fp[2] - 8*4
        cbc37 = 4*Fp[0] - 8*4
        ca3 = 4*Fp[2] - 8*8
        print(f"  a_{{n-2}} = {c0} + {ct3}*t3 + {ct5}*t5 + {ct7}*t7 + {ct9}*t9 + {cbc}*bc + {cbc35}*bc35 + {cbc37}*bc37 + {ca3}*a3")

    # Verify numerically
    num_check = 10
    all_ok = True
    for trial in range(num_check):
        A = random_tournament(n, n*500 + trial)
        dist = count_forward_edge_dist(A, n)
        a_n2 = dist.get(n-2, 0)
        t3 = count_t3(A, n)

        if n == 5:
            t5 = count_directed_cycles(A, n, 5)
            pred = pred_const + pred_t3 * t3 + pred_t5 * t5
            ok = (a_n2 == pred)
            if trial < 3 or not ok:
                print(f"    T{trial}: t3={t3}, t5={t5}, a_{{n-2}}={a_n2}, pred={pred}, {'OK' if ok else 'FAIL'}")
            if not ok: all_ok = False
        elif n == 9:
            t5 = count_directed_cycles(A, n, 5)
            t7 = count_directed_cycles(A, n, 7)
            # t9 and bc at n=9 would need more code; skip for now
            if trial < 3:
                print(f"    T{trial}: t3={t3}, t5={t5}, t7={t7}, a_{{n-2}}={a_n2}")

    if n == 5:
        print(f"  All {num_check}: {'PASS' if all_ok else 'FAIL'}")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
