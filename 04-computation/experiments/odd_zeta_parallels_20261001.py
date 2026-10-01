#!/usr/bin/env python3
"""Odd zeta values: the arithmetic gain in arXiv:2609.22316, Brown's dinner parties, and the Petersen graph
inside M_{0,5} (opus S15, 2026-10-01).

Checks:
  A. the preprint's exponent v_p (its eq. (15)) for h0 = 3n+2, h_1..h_11 = n+1: v_p = 3 if 3 (n mod p) >= 2p - 1, else 0
     (all primes sqrt(h0) < p <= n, many n; a carry-free split exists iff some residue r has 2m-p < r < p-m);
     hence phi(x) = 3 * 1_[2/3,1); the arithmetic gain
     lim (1/n) log Phi_n = 3 (psi(1) - psi(2/3)) - 3/2 = 0.72306, confirmed by direct prime sums; the preprint's
     C2' = 5.27694 equals 9 - 2.22306 - 1.5 (wrong sign); correct C2' = 8.27694 > C0' = 5.75349, and even the trivial
     bound log Phi_n <= 3 theta(n) gives C2' >= 6
  B. Brown's convergent configurations (dinner parties): Hamiltonian cycles of the complement of C_N
     (1, 3, 23, 177, 1553 for N = 5..9, OEIS A002493) and their classes under the dihedral symmetry (1, 1, 5, 19, 112)
  C. Petersen = Kneser KG(5,2) = intersection graph of the 10 boundary lines D_ij of M_{0,5}-bar; for the
     dihedral order delta = (12345) and its convergent partner delta' = (13524), the consecutive pairs of delta and
     of delta' are two disjoint 5-cycles of Petersen joined by a perfect matching (the standard drawing): Brown's
     zeta(2) configuration is the Petersen graph split into its two pentagons
  D. boundary divisors D_S of M_{0,8}-bar with |S| = 4 (each = M_{0,5}-bar x M_{0,5}-bar): 35 of them
Reproduce: python 04-computation/experiments/odd_zeta_parallels_20261001.py   (a few minutes)
"""
import math
from itertools import combinations

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def primes_upto(N):
    s = bytearray([1]) * (N + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]


print("=== A. the arithmetic gain of arXiv:2609.22316 ===")


def v_kp(n, k, p):
    h0, h = 3 * n + 2, n + 1
    f = lambda x: x // p
    A = f(k - 1) + f(h0 - k - 1) - f(k - h) - f(h0 - h - k) - 2 * f(h - 1)
    B = f(h0 - 2 * h) - f(k - h) - f(h0 - h - k)
    return 3 * A + 8 * B


ok = True
tested = 0
for n in list(range(10, 200)) + [997, 1500, 2023]:
    h0 = 3 * n + 2
    for p in primes_upto(n):
        if p * p <= h0:
            continue
        vp = min(v_kp(n, k, p) for k in range(n + 1, 2 * n + 2))
        m = n % p
        ok &= vp == (3 if 3 * m >= 2 * p - 1 else 0)
        tested += 1
check(ok, f"v_p = 3 [3 (n mod p) >= 2p - 1] (else 0) for all {tested} (n, p) pairs tested (sqrt(h0) < p <= n); "
      "asymptotically phi(x) = 3 * 1_[2/3,1)(frac x)")
gamma = 0.5772156649015329
psi1 = -gamma
psi23 = -gamma - 1.5 * math.log(3) + math.pi / (2 * math.sqrt(3))
gain = 3 * (psi1 - psi23) - 1.5
print(f"   int phi dpsi = {3 * (psi1 - psi23):.7f}, int phi dx/x^2 = 1.5, gain = {gain:.7f}")
N = 3 * 10 ** 6
P = primes_upto(N)
lg = sum(math.log(p) for p in P if p * p > 3 * N + 2 and (N % p) / p >= 2 / 3) * 3 / N
print(f"   direct: (1/n) log Phi_n at n = {N}: {lg:.5f}")
check(abs(lg - gain) < 0.01, "the prime sum confirms the gain 0.7231 (not 3.7231)")
paper_C2 = 9 - 3 * (psi1 - psi23) - 1.5
right_C2 = 9 - gain
C0 = 5.75349395
check(abs(paper_C2 - 5.27694374) < 1e-6, f"the preprint's C2' = 5.27694374 is 9 - 2.22306 - 1.5 (second integral subtracted)")
check(right_C2 > C0 and 9 - 3 > C0,
      f"correct C2' = {right_C2:.5f} > C0' = {C0}; even the trivial bound gives C2' >= 6 > C0': the criterion fails")

print("=== B. Brown's dinner parties ===")


def ham_cycles_complement(n):
    adj = lambda a, b: (a - b) % n in (1, n - 1)
    res = set()

    def rec(path, used):
        if len(path) == n:
            if not adj(path[-1], path[0]):
                c = tuple(path)
                res.add(min(c, (c[0],) + tuple(reversed(c[1:]))))
            return
        for v in range(n):
            if not used >> v & 1 and not adj(path[-1], v):
                rec(path + [v], used | 1 << v)
    rec([0], 1)
    return res


def canon_cycle(c):
    n = len(c)
    best = None
    for s in range(n):
        r = c[s:] + c[:s]
        for t in (r, (r[0],) + tuple(reversed(r[1:]))):
            if best is None or t < best:
                best = t
    return best


raw, orb = [], []
for n in range(5, 10):
    cyc = ham_cycles_complement(n)
    o = set()
    for c in cyc:
        o.add(min(canon_cycle(tuple(((-x if g >= n else x) + g % n) % n for x in c)) for g in range(2 * n)))
    raw.append(len(cyc))
    orb.append(len(o))
print("   raw:", raw, " up to the dihedral symmetry of delta:", orb)
check(raw == [1, 3, 23, 177, 1553] and orb == [1, 1, 5, 19, 112],
      "Hamiltonian cycles of complement(C_N): 1, 3, 23, 177, 1553 (A002493); dihedral classes 1, 1, 5, 19, 112")

print("=== C. the Petersen graph inside M_{0,5}-bar ===")
pairs = [frozenset(c) for c in combinations(range(1, 6), 2)]
pet = {a: {b for b in pairs if not a & b} for a in pairs}
check(len(pairs) == 10 and all(len(pet[a]) == 3 for a in pairs) and sum(len(v) for v in pet.values()) == 30,
      "Kneser KG(5,2): 10 vertices (the lines D_ij), 15 edges (D_ij meets D_kl iff {i,j}, {k,l} disjoint)")
delta = [1, 2, 3, 4, 5]
deltap = [1, 3, 5, 2, 4]
cons = lambda c: [frozenset((c[i], c[(i + 1) % 5])) for i in range(5)]
outer, inner = cons(delta), cons(deltap)
is_cycle = lambda V: all(sum(1 for b in V if b in pet[a]) == 2 for a in V)
matching = all(sum(1 for b in inner if b in pet[a]) == 1 for a in outer)
check(not set(outer) & set(inner) and is_cycle(outer) and is_cycle(inner) and matching,
      "consecutive pairs of (12345) and of its convergent partner (13524) are two disjoint 5-cycles of Petersen "
      "joined by a perfect matching: Brown's zeta(2) configuration = Petersen's pentagon + pentagram + spokes")
print("=== D. M_{0,8}-bar ===")
nS = sum(1 for S in combinations(range(8), 4)) // 2
check(nS == 35, "35 boundary divisors D_S of M_{0,8}-bar with |S| = 4, each isomorphic to M_{0,5}-bar x M_{0,5}-bar")

print()
print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
