#!/usr/bin/env python3
"""Hypothesis battery 2.
(a) River span law: a river (equal-time partner class on the Mersenne line) spans O(sqrt K) exponents, because
    K and K-D are partners only if the lag-D debt walk (start D) returns to 0 within the ~7.6K available Terras steps.
(b) Jump rate of the staircase sigma(M_K) by window; river count growth.
(c) Real shadow of the +1 barrier: the child y* = 2/3^D - 1 sheds one factor 3 per odd step and lands on an integer N_D
    after exactly D odd steps (landing time s_0, landing debt ceil(s_0/2)). Kesten-Goldie prediction: the real recursion
    v -> v/2, (3v+1)/2 with fair parities is a perpetuity with E[A^kappa] = 1 at kappa = 1, so P(N_D > n) ~ C/n."""
import math
from collections import defaultdict, Counter
data = {}
for line in open('../runcompress_20261007/mersenne_sigma_12800.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
Ks = sorted(data)
riv = defaultdict(list)
for K in Ks: riv[data[K][0]].append(K)          # rivers = level sets of the odd count (sigma_T - K is then forced)
# sanity: key (o, sigma_T - K) is a function of o
assert all(len({data[K][1] - K for K in v}) == 1 for o, v in riv.items() if min(v) >= 5)
print("(a) rivers = level sets of the odd count: confirmed (sigma_T - K is a function of o for K >= 5)")
rows = sorted((max(v), min(v), len(v)) for v in riv.values())
print("    river span (max K - min K) / sqrt(max K), by max-K window:")
for lo, hi in ((100, 400), (400, 1600), (1600, 6400), (6400, 12801)):
    sp = [(a - b) / math.sqrt(a) for a, b, n in rows if lo <= a < hi and n > 1]
    if sp: print(f"      maxK in [{lo},{hi}): {len(sp)} rivers, span/sqrt(K): mean {sum(sp)/len(sp):.2f}, max {max(sp):.2f}")
big = sorted(rows, key=lambda r: -(r[0] - r[1]))[:6]
print("    widest rivers (maxK, minK, size, span/sqrt(maxK)):", [(a, b, n, round((a - b)/math.sqrt(a), 2)) for a, b, n in big])
# (b) jumps by window
jumpK = [K for K in Ks if K - 1 in data and data[K][0] != data[K-1][0]]
print("(b) staircase jumps by window (rate, rate*sqrt(Kmid)):")
for lo, hi in ((100, 400), (400, 1600), (1600, 6400), (6400, 12801)):
    n = sum(1 for K in jumpK if lo <= K < hi); r = n / (hi - lo)
    print(f"      K in [{lo},{hi}): {n} jumps, rate {r:.4f}, rate*sqrt(Kmid) {r*math.sqrt((lo+hi)/2):.2f}")
first = {}
for K in Ks: first.setdefault(data[K][0], K)
starts = sorted(first.values())
for X in (400, 1600, 6400, 12800):
    R = sum(1 for s in starts if s <= X)
    print(f"      rivers started by K <= {X}: {R}   (R/X^0.41 = {R/X**0.41:.2f}, R/sqrt X = {R/math.sqrt(X):.2f})")
# (c) landing integers
def landing(D):
    a, m, s, o = 2 - 3**D, D, 0, 0
    while m > 0:
        if a & 1:
            a = (a + 3**(m-1)) // 2; m -= 1; o += 1
        else:
            a //= 2
        s += 1
    return a, s
Ns, S0 = [], []
DMAX = 4000
for D in range(3, DMAX + 1):
    N, s0 = landing(D); Ns.append(N); S0.append(s0 / D)
neg = sum(1 for N in Ns if N <= 0)
print(f"(c) landing integers N_D for 3 <= D <= {DMAX}: nonpositive {neg}; max {max(Ns)}; median {sorted(Ns)[len(Ns)//2]}; s_0/D in [{min(S0):.3f}, {max(S0):.3f}]")
for n in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024):
    p = sum(1 for N in Ns if N > n) / len(Ns)
    print(f"      P(N_D > {n:5d}) = {p:.4f}   n*P = {n*p:.3f}")
