#!/usr/bin/env python3
"""
THM-073: Master Decomposition of P(u,x).

G_T(t,x) = A_n(t) + sum_I x^|I| * A_{f+1}(t) * (t-1)^{S_I}

where S_I = sum(l_i - 1) for cycles in independent set I, f = (n-1) - S_I.

Converting to u = t + 1/t and P = G/t^m:

P(u,x) = P_n(u,0) + sum_{S>0} C_S(x) * P_{n-S}(u,0) * (u-2)^{S/2}

where C_S(x) = sum_{I: S_I=S} x^{|I|} * c_I (weighted count of independent sets with given S).

PROVED algebraically. Verified:
- n=5: exhaustive (1024 tournaments x 3 t-values = 3072 checks, 0 failures)
- n=7: random sample (200 tournaments x 4 t-values)

kind-pasteur-2026-03-07-S29
"""
from itertools import combinations
from collections import defaultdict
from math import comb
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def P_base(n_val, u_val):
    """Compute P_n(u, 0) = A_n(t)/t^m in u = t + 1/t."""
    d = n_val - 1
    m = d // 2
    coeffs = [eulerian_number(n_val, k) for k in range(n_val)]
    U = [0] * (m + 1)
    U[0] = 2
    if m >= 1: U[1] = u_val
    for j in range(2, m+1): U[j] = u_val * U[j-1] - U[j-2]
    result = coeffs[m]
    for j in range(1, m+1): result += coeffs[m+j] * U[j]
    return result

def tournament_from_bits(n, bits_int):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits_int >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def find_directed_odd_cycles(A, n):
    all_cycles = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            sub = [[A[subset[i]][subset[j]] for j in range(length)] for i in range(length)]
            dp = [[0]*length for _ in range(1 << length)]
            dp[1][0] = 1
            for mask in range(1, 1 << length):
                for v in range(length):
                    if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                    for w in range(length):
                        if mask & (1 << w): continue
                        if sub[v][w]: dp[mask | (1 << w)][w] += dp[mask][v]
            full = (1 << length) - 1
            count = sum(dp[full][v] for v in range(1, length) if sub[v][0])
            for _ in range(count):
                all_cycles.append((set(subset), length))
    return all_cycles

def enumerate_indep_by_S(cycles, x_val):
    """Return dict S -> sum of x^|I| for independent sets with that S value."""
    nc = len(cycles)
    if nc == 0:
        return {0: 1}
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][0] & cycles[j][0]:
                adj[i][j] = adj[j][i] = True
    S_sums = defaultdict(float)
    S_sums[0] = 1
    for mask in range(1, 1 << nc):
        nodes = [i for i in range(nc) if (mask >> i) & 1]
        ok = True
        for a in range(len(nodes)):
            for b in range(a+1, len(nodes)):
                if adj[nodes[a]][nodes[b]]:
                    ok = False
                    break
            if not ok: break
        if ok:
            parts = len(nodes)
            S = sum(cycles[i][1] - 1 for i in nodes)
            S_sums[S] += x_val ** parts
    return dict(S_sums)

def verify_n5_exhaustive():
    n = 5
    m = 2
    m_edges = 10
    x_val = 2
    print(f"n={n}: Exhaustive verification ({2**m_edges} tournaments)")

    errors = 0
    for bits in range(2**m_edges):
        A = tournament_from_bits(n, bits)
        cycles = find_directed_odd_cycles(A, n)
        S_sums = enumerate_indep_by_S(cycles, x_val)

        for t_val in [2, -1, 3]:
            u_val = t_val + 1.0/t_val

            # Direct from definition
            G_direct = sum(eulerian_number(n, k) * t_val**k for k in range(n))
            for S, coeff in S_sums.items():
                if S == 0: continue
                f = (n-1) - S
                A_f1 = sum(eulerian_number(f+1, k) * t_val**k for k in range(f+1))
                G_direct += coeff * A_f1 * (t_val - 1)**S
            P_direct = G_direct / t_val**m

            # Master formula
            P_master = P_base(n, u_val)
            for S, coeff in S_sums.items():
                if S == 0: continue
                P_master += coeff * P_base(n - S, u_val) * (u_val - 2)**(S//2)

            if abs(P_direct - P_master) > 0.01:
                errors += 1

    if errors == 0:
        print(f"  ALL {2**m_edges * 3} checks PASS!")
    else:
        print(f"  {errors} FAILURES")

def verify_n7_random(num_samples=200):
    n = 7
    m = 3
    m_edges = 21
    x_val = 2
    print(f"\nn={n}: Random sample verification ({num_samples} tournaments)")

    rng = random.Random(42)
    errors = 0

    for trial in range(num_samples):
        bits = rng.getrandbits(m_edges)
        A = tournament_from_bits(n, bits)
        cycles = find_directed_odd_cycles(A, n)
        S_sums = enumerate_indep_by_S(cycles, x_val)

        for t_val in [2, -1, 3, 0.5]:
            u_val = t_val + 1.0/t_val

            G_direct = sum(eulerian_number(n, k) * t_val**k for k in range(n))
            for S, coeff in S_sums.items():
                if S == 0: continue
                f = (n-1) - S
                A_f1 = sum(eulerian_number(f+1, k) * t_val**k for k in range(f+1))
                G_direct += coeff * A_f1 * (t_val - 1)**S
            P_direct = G_direct / t_val**m

            P_master = P_base(n, u_val)
            for S, coeff in S_sums.items():
                if S == 0: continue
                P_master += coeff * P_base(n - S, u_val) * (u_val - 2)**(S//2)

            if abs(P_direct - P_master) > 0.01:
                errors += 1
                if errors <= 3:
                    print(f"  FAIL trial={trial}, t={t_val}: {P_direct:.4f} vs {P_master:.4f}")

        if trial < 3:
            H = int(round(sum(S_sums.values())))  # at x=2, this is H(T)
            print(f"  trial={trial}: S_groups={dict(S_sums)}, H={H}")

    if errors == 0:
        print(f"  ALL {num_samples * 4} checks PASS!")
    else:
        print(f"  {errors} FAILURES")

def main():
    print("=" * 70)
    print("THM-073: MASTER DECOMPOSITION OF P(u,x)")
    print("=" * 70)
    print()
    print("P(u,x) = P_n(u,0) + sum_{S>0} C_S(x) * P_{n-S}(u,0) * (u-2)^{S/2}")
    print()
    print("Base polynomials P_k(u,0):")
    for k in [1, 3, 5, 7, 9]:
        print(f"  P_{k}(2,0) = {P_base(k, 2)} = {k}!")
    print()

    verify_n5_exhaustive()
    verify_n7_random()

    print()
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print()
    print("G_T(t, x) is NOT the HP descent polynomial E_T(t).")
    print("It is the 'inflated independence polynomial':")
    print("  G_T(0, x) = I(Omega(T), x)  [independence polynomial]")
    print("  G_T(0, 2) = H(T)            [Hamiltonian path count]")
    print("  G_T(1, x) = n!              [universal, corrections vanish]")
    print()
    print("The master decomposition groups cycle types by S = sum(l_i-1):")
    print("  Different types with same S get same u-polynomial P_{n-S}(u,0)*(u-2)^{S/2}")
    print("  Only the x-weighting (number of parts) distinguishes them")

if __name__ == "__main__":
    main()
