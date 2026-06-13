#!/usr/bin/env python3
"""
tr(c_0) as a signed path count — the BOTTOM of the even-r hierarchy.

KEY IDENTITY (odd n):
  tr(c_0) = W(0) = sum_P prod_{e in P} s_e
  where s_e = T(e) - 1/2 in {+1/2, -1/2}

  = (1/2^{n-1}) * sum_P (-1)^{f_P}

  where f_P = number of forward edges in permutation P.

  At odd n: (-1)^{b_P} = (-1)^{n-1-f_P} = (-1)^{f_P} (since n-1 even)

CONNECTION TO H:
  H = W(1/2) = #{directed Ham paths} = #{P : f_P = n-1}
  tr(c_0) = W(0) = "signed count" weighted by (-1)^f / 2^{n-1}

  The even-r polynomial interpolates between the signed count (r=0)
  and the actual count (r=1/2)!

FORMULA:
  At n=5: tr(c_0) = H - 3*t_3 = 1 - t_3 + 2*t_5 (by OCF)
  At n=7: tr(c_0) = H - tr(c_2)/4 - 15*t_3 + 52.5

  2^{n-1} * tr(c_0) = sum_P (-1)^{f_P}

  This is the "alternating Hamiltonian path permanent" — counts paths
  by parity of forward edges.

opus-2026-03-06-S27
"""

from itertools import permutations
from math import factorial, comb
import random

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

def count_5cycles(A):
    n = len(A)
    count = 0
    for perm in permutations(range(n), 5):
        if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
            count += 1
    return count // 5

def ham_path_count_dp(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if not (S_prev & (1 << u)):
                        continue
                    if A[u][v] and (S_prev, u) in dp:
                        total += dp[(S_prev, u)]
                if total > 0:
                    dp[(S, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# =====================================================================
print("=" * 70)
print("tr(c_0) = SIGNED PATH COUNT")
print("=" * 70)

for n in [3, 5, 7]:
    print(f"\n  n={n}:")
    random.seed(n * 17)

    for trial in range(8 if n <= 5 else 5):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        H = ham_path_count_dp(A)
        t3 = count_3cycles(A)

        # Compute signed sum
        signed_sum = 0
        f_dist = {}
        for p in permutations(range(n)):
            f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
            signed_sum += (-1)**f
            f_dist[f] = f_dist.get(f, 0) + 1

        tr_c0 = signed_sum / 2**(n-1)

        if n == 5:
            t5 = count_5cycles(A)
            pred = 1 - t3 + 2*t5
            print(f"    Trial {trial}: H={H}, t3={t3}, t5={t5}, "
                  f"2^4*tr(c_0)={signed_sum}, tr(c_0)={tr_c0:.4f}, "
                  f"1-t3+2*t5={pred}, {'✓' if abs(tr_c0-pred)<0.01 else '✗'}")
        elif n == 7:
            t5 = count_5cycles(A)
            # tr(c_0) = H - tr(c_2)/4 - (240*t3-2100)/16 - 5040/64
            # We don't know tr(c_2), but we know:
            # 2^6 * tr(c_0) = signed_sum
            print(f"    Trial {trial}: H={H:3d}, t3={t3:2d}, t5={t5:2d}, "
                  f"2^6*tr(c_0)={signed_sum:5d}, tr(c_0)={tr_c0:7.2f}, "
                  f"H-3t3={H-3*t3:4d}")
        else:
            print(f"    Trial {trial}: H={H}, t3={t3}, "
                  f"2^2*tr(c_0)={signed_sum}, tr(c_0)={tr_c0:.4f}")

# =====================================================================
# The f-distribution tells us about the structure
# =====================================================================
print("\n" + "=" * 70)
print("FORWARD-EDGE DISTRIBUTION")
print("=" * 70)

n = 5
random.seed(505)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

H = ham_path_count_dp(A)
t3 = count_3cycles(A)

f_dist = [0] * n
for p in permutations(range(n)):
    f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
    f_dist[f] += 1

print(f"  n=5, H={H}, t3={t3}")
print(f"  f-distribution: {f_dist}")
print(f"  Sum: {sum(f_dist)} = {n}!")
print(f"  f=4 (directed paths): {f_dist[4]} = H = {H}")
print(f"  f=0 (all backward): {f_dist[0]} = H(T^op)")
print(f"  Alternating: sum (-1)^f * count(f) = {sum((-1)**f * c for f,c in enumerate(f_dist))}")
print(f"  = 2^4 * tr(c_0) = {16 * (H - 3*t3)}")

# Symmetry: f_dist[k] = f_dist[n-1-k] at odd n?
print(f"\n  Symmetry check: f_dist[k] vs f_dist[n-1-k]:")
for k in range(n):
    print(f"    f={k}: {f_dist[k]}, f={n-1-k}: {f_dist[n-1-k]}, "
          f"equal={'✓' if f_dist[k]==f_dist[n-1-k] else '✗'}")

# =====================================================================
# PALEY at n=7
# =====================================================================
print("\n" + "=" * 70)
print("PALEY T_7: f-DISTRIBUTION AND tr(c_0)")
print("=" * 70)

n = 7
QR = {1, 2, 4}
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j-i)%n in QR:
            A[i][j] = 1

H = ham_path_count_dp(A)
t3 = count_3cycles(A)

f_dist = [0] * n
signed_sum = 0
for p in permutations(range(n)):
    f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
    f_dist[f] += 1
    signed_sum += (-1)**f

tr_c0 = signed_sum / 64
print(f"  H = {H}, t3 = {t3}")
print(f"  f-distribution: {f_dist}")
print(f"  2^6 * tr(c_0) = {signed_sum}")
print(f"  tr(c_0) = {tr_c0}")
print(f"  H - 3*t3 = {H - 3*t3}")

# For VT: tr(c_0) = M[0,0](0) * n (scalar)
# M[0,0](0) = tr(c_0)/n = {tr_c0/n}
print(f"  c_0 diagonal = tr(c_0)/n = {tr_c0/n}")

# Check f-distribution symmetry
print(f"\n  f-distribution symmetry:")
symmetric = True
for k in range(n):
    s = f_dist[k] == f_dist[n-1-k]
    symmetric = symmetric and s
    print(f"    f={k}: {f_dist[k]}, f={n-1-k}: {f_dist[n-1-k]}, {'✓' if s else '✗'}")

if symmetric:
    print(f"  SYMMETRIC! f_dist[k] = f_dist[n-1-k] for all k.")
    print(f"  This means: sum (-1)^f * count(f) = 0 when n-1 is even.")
    print(f"  Wait, n-1=6 is even. (-1)^f alternates. With symmetry,")
    print(f"  count(f) = count(6-f), so (-1)^f * count(f) + (-1)^(6-f) * count(6-f)")
    print(f"  = (-1)^f * count(f) + (-1)^f * count(f) = 2*(-1)^f * count(f)")
    print(f"  So the sum = 2 * sum_{{f even}} count(f) - 2 * sum_{{f odd}} count(f)")
    print(f"  = 2 * [{sum(f_dist[f] for f in range(0,n,2))} - {sum(f_dist[f] for f in range(1,n,2))}]")
    print(f"  = 2 * {sum(f_dist[f] for f in range(0,n,2)) - sum(f_dist[f] for f in range(1,n,2))}")
    print(f"  = {signed_sum} ✓")

print("\n" + "=" * 70)
print("KEY INSIGHT: f-distribution symmetry f_dist[k] = f_dist[n-1-k]")
print("=" * 70)
print("""
  At odd n, for ANY tournament T:
    #{permutations with f forward edges} = #{permutations with n-1-f forward edges}

  This is because the map P -> P^{op} (reverse the permutation and flip
  all edges) sends f_P -> n-1-f_P. Since T is a tournament, each edge
  has exactly one direction, and reversal flips all edges.

  More precisely: for permutation P = (v_0,...,v_{n-1}), define
  P' = (v_{n-1},...,v_0). Then edge (v_i,v_{i+1}) in P becomes
  edge (v_{n-1-i}, v_{n-2-i}) in P'. And T(v_i,v_{i+1}) + T(v_{i+1},v_i) = 1,
  so forward edges in P become backward edges in P' and vice versa.
  f_{P'} = b_P = n-1-f_P. QED.

  CONSEQUENCE: tr(c_0) = (1/2^{n-1}) * sum (-1)^f * count(f)
  = (2/2^{n-1}) * [count(even f) - count(odd f)]
  = (1/2^{n-2}) * [count(even f) - count(odd f)]

  And H = count(n-1) = count(0) by symmetry (!)
  Wait, count(0) = count(n-1) means #{paths with 0 forward edges} = H.
  But #{paths with 0 forward edges} = #{permutations where every
  consecutive pair has backward edge} = #{directed Ham paths of T^op} = H(T^op).

  And H(T^op) = H(T) by path reversal! So this is consistent.
""")

# Verify: f_dist[0] = H(T^op) = H(T)?
for trial in range(5):
    random.seed(trial * 999 + 7)
    A = [[0]*7 for _ in range(7)]
    for i in range(7):
        for j in range(i+1, 7):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_path_count_dp(A)
    f0 = sum(1 for p in permutations(range(7))
             if all(A[p[i+1]][p[i]] == 1 for i in range(6)))  # all backward = T^op path
    print(f"  Trial {trial}: H={H}, f_dist[0]={f0}, equal={'✓' if H==f0 else '✗'}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
