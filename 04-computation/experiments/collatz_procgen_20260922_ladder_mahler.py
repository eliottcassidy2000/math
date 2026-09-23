#!/usr/bin/env python3
"""collatz_procgen_20260922_ladder_mahler.py -- the Mahler 3/2 rung of the sibling ladder (re-check only).

Re-verifies, by three independent counts, the safe-prefix language of THM-3848 section 6:
  P_m = binary words c_0..c_{m-1} with 2*sum_{j<m-i} c_{i+j} 2^j 3^{m-i-1-j} < 3^{m-i} for all i  (S2);
  (a) incremental exact DFS over (S2);
  (b) the renewal law a_m = 1 + sum_{k<=m} d_k a_{m-k}, d = greedy (3/2)-expansion of 1          (S6),(S8);
  (c) Parry's lexicographic condition: every suffix <=_lex the equal-length prefix of d
      (i.e. K is the closed beta-shift for beta = 3/2; standard identification, checked here finitely);
and the asymptotics a_m ~ C (3/2)^m, C = 9/(2 D'(2/3)) = 1.5510451884 (S10), so the prefix-count
dimension is log_2(3/2) in the binary ultrametric, transported isometrically to Z_2 by THM-2228's Phi
(checked: r_m is a bijection from length-m words onto Z/2^m for m <= 12).
"""
import math
from fractions import Fraction


def greedy_digits(n):
    z, d = Fraction(1), [0]
    for _ in range(n):
        v = Fraction(3, 2) * z
        dk = math.floor(v)
        d.append(dk)
        z = v - dk
        assert 0 < z < 1
    return d  # d[1..n]


def renewal(d, M):
    a = [1]
    for m in range(1, M + 1):
        a.append(1 + sum(d[k] * a[m - k] for k in range(1, m + 1)))
    return a


def dfs_S2(M):
    """count words of each length m<=M satisfying all strict suffix inequalities (S2), exactly."""
    counts = [0] * (M + 1)
    # state: tuple of V_i for i<m, where V_i = sum_j c_{i+j} 2^j 3^{m-i-1-j}; condition 2 V_i < 3^{m-i}
    stack = [(0, ())]
    while stack:
        m, V = stack.pop()
        counts[m] += 1
        if m == M:
            continue
        for c in (0, 1):
            # append c at position m: new V_i = 3 V_i + c 2^{m-i}; new V_m = c
            W = tuple(3 * v + c * (1 << (m - i)) for i, v in enumerate(V)) + (c,)
            if all(2 * W[i] < 3 ** (m + 1 - i) for i in range(m + 1)):
                stack.append((m + 1, W))
    return counts


def dfs_parry(d, M):
    counts = [0] * (M + 1)
    stack = [(0, "")]
    ds = "".join(str(x) for x in d[1:M + 2])
    while stack:
        m, w = stack.pop()
        counts[m] += 1
        if m == M:
            continue
        for c in "01":
            u = w + c
            if all(u[i:] <= ds[:len(u) - i] for i in range(len(u))):
                stack.append((m + 1, u))
    return counts


def phi_bijection(m):
    seen = set()
    for w in range(1 << m):
        C = sum(((w >> j) & 1) * (1 << j) * 3 ** (m - 1 - j) for j in range(m))
        r = (-C * pow(3, -m, 1 << m)) % (1 << m)
        seen.add(r)
    return len(seen) == (1 << m)


def main():
    print("=" * 78)
    print("LADDER-MAHLER. THM-3848 safe-prefix language re-counted three ways")
    print("=" * 78)
    M = 30
    d = greedy_digits(400)
    print("greedy word d (first 40): " + "".join(str(x) for x in d[1:41]))
    a = renewal(d, 400)
    s2 = dfs_S2(M)
    pa = dfs_parry(d, M)
    ok = a[:M + 1] == s2 == pa
    print(f"counts a_m, m=0..20: {a[:21]}")
    print(f"(a) DFS over (S2) == (b) renewal == (c) Parry-lexicographic for m<=30: {'YES' if ok else 'NO'}")
    Dp = sum(k * d[k] * Fraction(2, 3) ** (k - 1) for k in range(1, 401))
    C = 9 / (2 * float(Dp))
    print(f"C = 9/(2 D'(2/3)) = {C:.10f} (THM-3848 states 1.5510451884)")
    for m in (10, 20, 50, 100, 200, 400):
        print(f"  m={m:>3}: a_m/(3/2)^m = {a[m] / 1.5 ** m:.10f}   log2(a_m)/m = {math.log2(a[m]) / m:.6f}   "
              f"(log2 a_m - log2 C)/m = {(math.log2(a[m]) - math.log2(C)) / m:.8f}")
    print(f"log2(3/2) = {math.log2(1.5):.8f}  (entropy in bits = binary-ultrametric / 2-adic dimension)")
    print(f"Haar mass of the level-m safe cylinders a_m/2^m: m=20 {a[20]/2**20:.3e}, m=100 {a[100]/2**100:.3e}")
    ok2 = all(phi_bijection(m) for m in range(1, 13))
    print(f"THM-2228 Phi: length-m words -> Z/2^m bijective for m<=12 (isometry onto Z_2): {'YES' if ok2 else 'NO'}")
    print(f"LADDER-MAHLER TOTAL: {'ALL CHECKS PASS' if ok and ok2 and abs(C - 1.5510451884) < 1e-9 else 'CHECK FAILED'}")


if __name__ == "__main__":
    main()
