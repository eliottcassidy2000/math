#!/usr/bin/env python3
"""Pisano periods pi(n), n <= NMAX: the ratio pi(n)/n, its maximum 6 exactly at n = 2*5^k; the owner's examples
pi(11)=10, pi(19)=18, pi(12)=24, pi(20)=60, pi(10)=60, pi(50)=300; split primes pi | p-1, inert pi | 2(p+1), pi(5)=20;
the Fibonacci matrix mod 5 is projectively of order 5 in PGL(2,5) and lies in PSL(2,5) = A5 (icosahedral rotations)."""
from sympy import primerange
NMAX = 20000
def pisano(n):
    a, b, k = 0, 1, 0
    while True:
        a, b = b, (a + b) % n; k += 1
        if a == 0 and b == 1: return k
P = {n: pisano(n) for n in range(2, NMAX + 1)}
mx = max(P[n] / n for n in P)
print("examples:", {n: P[n] for n in (10, 11, 12, 19, 20, 50, 5, 2, 3, 4, 25, 250)})
print(f"max pi(n)/n for n <= {NMAX}: {mx}; attained at", [n for n in P if P[n] == 6 * n])
bad = [p for p in primerange(7, NMAX) if (p % 5 in (1, 4) and (p - 1) % P[p]) or (p % 5 in (2, 3) and (2 * (p + 1)) % P[p])]
print("split pi | p-1 and inert pi | 2(p+1), primes 7..NMAX: exceptions", bad)
eq_split = [p for p in primerange(7, 400) if p % 5 in (1, 4) and P[p] == p - 1]
eq_inert = [p for p in primerange(7, 400) if p % 5 in (2, 3) and P[p] == 2 * (p + 1)]
print("primes < 400 with maximal period: split pi = p-1:", eq_split, "; inert pi = 2(p+1):", eq_inert)
# Fibonacci matrix mod 5
Q = ((1, 1), (1, 0))
def mm(A, B, m): return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2)) % m for j in range(2)) for i in range(2))
M = ((1, 0), (0, 1)); powers = []
for k in range(1, 21):
    M = mm(M, Q, 5); powers.append((k, M))
scal = [k for k, Mk in powers if Mk[0][1] == 0 and Mk[1][0] == 0 and Mk[0][0] == Mk[1][1]]
print("Q^k mod 5 scalar for k =", scal, "; Q^5 mod 5 =", powers[4][1], "; det Q = -1 = 4 = 2^2 is a square mod 5, so Q lies in PSL(2,5) = A5")
