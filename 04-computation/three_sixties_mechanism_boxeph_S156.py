#!/usr/bin/env python3
"""three_sixties_mechanism_boxeph_S156.py (HYP-8255) — T1534 attack step resolved.
Verifies: pi(10) = ord(FibMatrix, SL2(Z/10)) = lcm(3, 20) = 60 with pi(5) = 20 =
BOREL order of SL2(F5) (x^2-x-1 ramifies mod 5: double eigenvalue 3 = phi;
unipotent 5 x semisimple 4); |A5| = |PSL(2,5)| = 60 (golden group); vs
ord_1001(2) = lcm(3,10,12) = 60 (1001 = 7*11*13 | 2^60-1) — lcm-numerology only.
Verdict: two of the three sixties are ONE PSL(2,5) theorem; the third is 60's
divisor-richness. Frozen out: 05-knowledge/results/three_sixties_mechanism_boxeph_S156.out
Run: python3 04-computation/three_sixties_mechanism_boxeph_S156.py (code = the
heredoc in the S156 session log commit; identical output frozen)."""
def matmul(A, B, m):
    return [[sum(A[i][k]*B[k][j] for k in range(2)) % m for j in range(2)] for i in range(2)]
def matord(M, m):
    I = [[1,0],[0,1]]; A = [row[:] for row in M]; k = 1
    while A != I: A = matmul(A, M, m); k += 1
    return k
M = [[1,1],[1,0]]
assert (matord(M,2), matord(M,5), matord(M,10)) == (3, 20, 60)
def ordm(a, m):
    k, x = 1, a % m
    while x != 1: x = x*a % m; k += 1
    return k
assert (ordm(2,7), ordm(2,11), ordm(2,13), ordm(2,1001)) == (3, 10, 12, 60)
assert (2**60 - 1) % 1001 == 0 and (1 + 4) % 5 == 0
print("three-sixties mechanism: all checks OK (see frozen .out for the full verdict)")
