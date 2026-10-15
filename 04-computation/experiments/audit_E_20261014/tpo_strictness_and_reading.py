#!/usr/bin/env python3
"""Audit E (2026-10-14), THM-4655 statement 1, two supplements (exact, via the path <-> A bijection, which
tpo_dfs_bruteforce.c verified against DFS-enumerated Hamiltonian paths for every word with m <= 13):
(1) strictness: at the unique equality word w = 1 0 1^(m-3) 0 (m >= 5), P(u_3 before z) = 2^(m-4) / (1 + 2^(m-3)) lies strictly
    inside (1/3, 1/2); at m = 4 (w = 1010) it equals 1/3 and the best pair is reported; so delta(T) > 1/3 for EVERY strong
    transitive-plus-one tournament with m >= 3 (n >= 4);
    also P(u_2 before u_1) = 1/3 exactly there (u_2's two best partners u_1, z are both at balance exactly 1/3);
(2) the Reading's claim: for w = 1^a 0^b, most pairs i < j have P(u_i before u_j) close to 3/4."""
from fractions import Fraction as F

def tallies(w):
    m = len(w); n = m + 1; Z = m
    H = 0; bef = [[0] * n for _ in range(n)]
    for A in range(1 << m):
        B = ((1 << m) - 1) ^ A
        if A and w[A.bit_length() - 1] != 0: continue
        if B and w[(B & -B).bit_length() - 1] != 1: continue
        H += 1
        order = [i for i in range(m) if A >> i & 1] + [Z] + [i for i in range(m) if not A >> i & 1]
        for s in range(n):
            for t in range(s + 1, n):
                bef[order[s]][order[t]] += 1
    return H, bef

print("(1) equality word w = 1 0 1^(m-3) 0")
for m in range(4, 15):
    w = (1, 0) + (1,) * (m - 3) + (0,)
    H, bef = tallies(w)
    p3 = F(bef[2][m], H); pred = F(2 ** (m - 4), 1 + 2 ** (m - 3))
    p2z = F(bef[1][m], H); p21 = F(bef[1][0], H)
    print(f"  m={m}: H={H} (pred 3+3*2^(m-3) = {3 + 3 * 2 ** (m - 3)}), P(u2<z)={p2z}, P(u2<u1)={p21}, "
          f"P(u3<z)={p3} (formula 2^(m-4)/(1+2^(m-3)) = {pred}, {'ok' if p3 == pred else 'MISMATCH'}), strictly inside (1/3,2/3): {F(1,3) < p3 < F(2,3)}")
    if m == 4:
        nm = lambda v: 'z' if v == m else 'u%d' % (v + 1)
        best = max((min(F(bef[x][y], H), 1 - F(bef[x][y], H)), nm(x), nm(y), F(bef[x][y], H)) for x in range(m + 1) for y in range(m + 1) if x < y)
        allbest = [(nm(x), nm(y), str(F(bef[x][y], H))) for x in range(m + 1) for y in range(m + 1) if x < y and min(F(bef[x][y], H), 1 - F(bef[x][y], H)) == best[0]]
        print(f"     m=4: delta = {best[0]}, attained by pairs (x, y, P(x before y)): {allbest}")

print("(2) w = 1^a 0^b: distribution of P(u_i before u_j), i < j")
for a, b in [(6, 7), (7, 7)]:
    w = (1,) * a + (0,) * b
    H, bef = tallies(w)
    m = a + b
    ps = [F(bef[i][j], H) for i in range(m) for j in range(i + 1, m)]
    near = sum(1 for p in ps if abs(p - F(3, 4)) < F(1, 20))
    unbal = sum(1 for p in ps if p > F(2, 3) or p < F(1, 3))
    print(f"  a={a} b={b}: {len(ps)} pairs; |P - 3/4| < 0.05 for {near}; unbalanced (outside [1/3,2/3]) {unbal}; "
          f"min {float(min(ps)):.4f} max {float(max(ps)):.4f} mean {float(sum(ps) / len(ps)):.4f}")
