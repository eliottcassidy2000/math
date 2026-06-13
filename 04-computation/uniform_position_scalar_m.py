#!/usr/bin/env python3
"""
CONJECTURE: At odd n, uniform position matrix <==> M = (H/n)*I.

P[v,k] = # Ham paths where vertex v is at position k.
"Uniform" means P[v,k] = H/n for all v,k.

This is equivalent to: every vertex appears equally often at every position.

If true, this gives a clean characterization of when M is scalar:
  M = (H/n)*I iff the tournament has "position-uniform" Hamiltonian paths.

PROOF DIRECTION: If P[v,k] = H/n for all v,k, then:
  M[a,a] = sum_k (-1)^k P[a,k] = (H/n) sum_k (-1)^k = (H/n) * 1 (odd n) or 0 (even n)
  This matches tr(M) = H for odd n.

  For M[a,b] = 0 (a != b), we need: the PAIR position distribution is also uniform.
  Define P2[a,j,b,k] = # Ham paths where a is at pos j and b is at pos k.
  Then M[a,b] = ... (some signed combination of P2 values).

  If P2[a,j,b,k] = H/(n*(n-1)) for all j != k, then the signed sum might vanish.

Let's test this at n=5.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

def pair_position_matrix(T, n):
    """P2[a,j,b,k] = # Ham paths where a at pos j and b at pos k."""
    P2 = np.zeros((n, n, n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for i in range(n):
                for j in range(n):
                    if i != j:
                        pi = list(perm).index(i)
                        pj = list(perm).index(j)
                        P2[i, pi, j, pj] += 1
    return P2

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        from itertools import combinations as comb
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            ea = 0
            if len(S_set) == 1:
                ea = 1
            else:
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    ea += prod
            bb = 0
            if len(R_set) == 1:
                bb = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    bb += prod
            val += sign * ea * bb
        return val

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T


# ============================================================
# n=5: Exhaustive test of conjecture
# ============================================================
print("=" * 70)
print("n=5: Uniform position matrix <==> M scalar?")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

uniform_pos_scalar = 0
uniform_pos_not_scalar = 0
not_uniform_scalar = 0
not_uniform_not_scalar = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)

    is_uniform_pos = all(P[v,k] == H // n for v in range(n) for k in range(n)) if H % n == 0 else False

    # Check scalar M: just check a few off-diagonal entries
    is_scalar = True
    for a in range(n):
        for b in range(a+1, n):
            if compute_M_entry(T, n, a, b) != 0:
                is_scalar = False
                break
        if not is_scalar:
            break

    if is_uniform_pos and is_scalar:
        uniform_pos_scalar += 1
    elif is_uniform_pos and not is_scalar:
        uniform_pos_not_scalar += 1
    elif not is_uniform_pos and is_scalar:
        not_uniform_scalar += 1
    else:
        not_uniform_not_scalar += 1

print(f"\n  Results:")
print(f"    Uniform pos AND scalar M:     {uniform_pos_scalar}")
print(f"    Uniform pos but NOT scalar M: {uniform_pos_not_scalar}")
print(f"    NOT uniform pos but scalar M: {not_uniform_scalar}")
print(f"    Neither:                       {not_uniform_not_scalar}")

if uniform_pos_not_scalar == 0 and not_uniform_scalar == 0:
    print(f"\n  PERFECT EQUIVALENCE: Uniform position <==> M scalar at n=5!")
elif uniform_pos_not_scalar == 0:
    print(f"\n  Uniform position => M scalar (but not converse)")
elif not_uniform_scalar == 0:
    print(f"\n  M scalar => uniform position (but not converse)")


# ============================================================
# Pair position analysis for uniform-position tournaments
# ============================================================
print("\n" + "=" * 70)
print("Pair position analysis for uniform-pos tournaments")
print("=" * 70)

# Take a representative with uniform positions
for bits in [40, 76]:  # bits=40 is non-regular H=15, bits=76 is regular H=15
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))

    P2 = pair_position_matrix(T, n)

    print(f"\n  bits={bits}: H={H}, scores={scores}")

    # Check if P2 is "uniform": P2[a,j,b,k] = H/(n*(n-1)) for all a!=b, j!=k
    target = H / (n * (n-1))
    max_dev = 0
    for a in range(n):
        for b in range(n):
            if a == b: continue
            for j in range(n):
                for k in range(n):
                    if j == k: continue
                    dev = abs(P2[a,j,b,k] - target)
                    max_dev = max(max_dev, dev)

    print(f"    Pair target H/(n*(n-1)) = {target:.3f}, max deviation = {max_dev:.3f}")

    # Show P2 for specific pair
    a, b = 0, 1
    print(f"    P2[{a},j,{b},k]:")
    for j in range(n):
        row = [P2[a,j,b,k] for k in range(n)]
        print(f"      j={j}: {row}")

    # Also check: is sum_k (-1)^k P2[a,j,b,k] = something simple?
    print(f"    sum_k (-1)^k P2[{a},j,{b},k]:")
    for j in range(n):
        val = sum((-1)**k * P2[a,j,b,k] for k in range(n))
        print(f"      j={j}: {val}")


# ============================================================
# PROOF IDEA: M[a,b] expressed via P2
# ============================================================
print("\n" + "=" * 70)
print("M[a,b] as function of P2[a,j,b,k]")
print("=" * 70)

# Recall M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
# E_a(S+a) counts paths in T[S+a] ending at a.
# B_b(R+b) counts paths in T[R+b] starting at b.
# The pair (path in S+a, path in R+b) is NOT the same as a global path.
# But the FULL inclusion-exclusion sum gives M.

# However, let's try to express M[a,b] as sum_{j,k} c(j,k) P2[a,j,b,k]
# for some coefficients c(j,k).

# If M[a,b] = sum_P weight(P), then weight(P) depends on (pos(a,P), pos(b,P)).
# But M[a,b] is NOT a simple sum over Hamiltonian paths of T!

# For diagonal: M[a,a] = sum_P (-1)^{pos(a,P)} — a sum over paths with position-dependent sign.
# For off-diagonal: M[a,b] is defined by inclusion-exclusion, not by summing over paths.

# However, there might be a FORMULA expressing M[a,b] in terms of P2.
# Let me search for coefficients c(j,k) such that M[a,b] = sum c(j,k) P2[a,j,b,k].

for bits in [40, 76]:
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P2 = pair_position_matrix(T, n)

    print(f"\n  bits={bits}: H={H}")

    for a in range(n):
        for b in range(a+1, min(a+3, n)):
            M_val = compute_M_entry(T, n, a, b)

            # Build system: find c(j,k) for j,k in range(n), j!=k
            # such that M[a,b] = sum c(j,k) P2[a,j,b,k]
            # This is underdetermined — but check if (-1)^j works
            check_1 = sum((-1)**j * P2[a,j,b,k] for j in range(n) for k in range(n))
            check_2 = sum((-1)**(j+k) * P2[a,j,b,k] for j in range(n) for k in range(n))
            check_3 = sum((-1)**j * (1 if j < k else -1 if j > k else 0) * P2[a,j,b,k]
                         for j in range(n) for k in range(n))

            # Most promising: weighted by relative position
            check_4 = 0
            for j in range(n):
                for k in range(n):
                    if j == k: continue
                    w = (-1)**j * (1 if k > j else -1) / 2
                    check_4 += w * P2[a,j,b,k]

            print(f"    M[{a},{b}]={M_val}: (-1)^j={check_1}, (-1)^(j+k)={check_2}, "
                  f"(-1)^j*sgn(k-j)={check_3}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
