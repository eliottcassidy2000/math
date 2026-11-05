#!/usr/bin/env python3
"""signed_clock_recipe_20261105.py -- opus-2026-11-05-S80

Exact constructive verification of signed-response clocks on the one-apex O/Z orbit (Session198's calculus, which is
the G(1,0) automaton in profile coordinates (D, E, U, d) = (H, 1 - last bit, 1 + binary value, 1)).

  O(D,E,U,d) = (D, 0, 2U, d),   Z(D,E,U,d) = (D - E + U, d, 2U - d, d).

For a slope a/b (coprime, 1/2 < a/b < 1) with b odd: B = 2^k = 1 (mod b), A = a(B-1)/b, W = the letters of the bits of A
(O for 0, Z for 1, least significant bit first).  Session198 (a even, b odd: A even, W starts with O): the response
aU - bD drops by K = aR - bC per period.  Extension (a, b both odd: A odd, W starts with Z): the drop is K_cyc = K_lin - b
once the input is strong (E = d = 1), where K_lin uses Session198's C and K_cyc adds the cyclic term eps_{k-1} eps_0 to C;
the first period after an O (E = 0) drops by K_lin.  Every strong orbit point is a one-apex tournament, so each
verified response is a physical signed count.
Checks: (1) the per-period drops on random strong inputs (exact); (2) the recipe seed -> O -> W^h hits random targets
exactly, with the final point strong (E = 1); (3) the mod-4 confinement for b even (class -(a+b) on the positive side
when a/b <= 3/4, class 2 - a or -a on the negative side when a/b > 3/4).
"""
import random, sys
from math import gcd

def g(n):
    if n == 1: return 1
    if n % 2 == 0: return g(n // 2)
    m = (n + 1) // 2
    return g(m) + m - (m % 2)

def O(p): D, E, U, d = p; return (D, 0, 2 * U, d)
def Z(p): D, E, U, d = p; return (D - E + U, d, 2 * U - d, d)
def apply_word(p, word):
    for ch in word: p = O(p) if ch == 'O' else Z(p)
    return p

def clock(a, b):
    k = 1
    while pow(2, k, b) != 1: k += 1
    B = 1 << k; A = a * (B - 1) // b; assert b * A == a * (B - 1)
    eps = [(A >> i) & 1 for i in range(k)]
    word = ''.join('Z' if e else 'O' for e in eps)
    R = sum(eps[i] << (k - 1 - i) for i in range(k))
    C_lin = sum(eps[j] * eps[i] << (i - j - 1) for i in range(k) for j in range(i)) + sum(eps[i - 1] * eps[i] for i in range(1, k))
    C_cyc = C_lin + eps[k - 1] * eps[0]
    return k, B, A, eps, word, R, a * R - b * C_lin, a * R - b * C_cyc

def resp(p, a, b): return a * p[2] - b * p[0]

def main():
    rng = random.Random(20261105)
    slopes = [(2, 3), (4, 5), (6, 7), (4, 7), (8, 11), (3, 5), (5, 7), (5, 9), (7, 9), (7, 11), (9, 13)]
    for (a, b) in slopes:
        k, B, A, eps, word, R, K_lin, K_cyc = clock(a, b)
        # (1) per-period drops
        ok = True
        for _ in range(200):
            n = rng.randrange(3, 1 << 20) | 1
            p = (g(n), 1, n, 1)                                   # strong one-apex source (E = 1 since n odd)
            q = apply_word(p, word)
            if resp(q, a, b) != resp(p, a, b) - K_cyc or q[1] != 1: ok = False
            p0 = O(p); q0 = apply_word(p0, word)                   # input with E = 0
            if resp(q0, a, b) != resp(p0, a, b) - K_lin: ok = False
        # (2) recipe for random targets of the right parity
        par = (a - b) % 2; hits = 0; tried = 0; maxh = 0
        for _ in range(40):
            r = rng.randrange(-10000, 10000)
            if r % 2 != par: r += 1
            tried += 1
            target_S = r + b * eps[0]            # need S = 2aU - bD = target_S (mod K_cyc), S - target_S >= K_cyc... S - b eps_0 - r >= K_cyc
            found = None
            for n in range(3, 1 << 22, 2):
                S = 2 * a * n - b * g(n)
                if (S - target_S) % K_cyc == 0 and S - target_S >= K_cyc:
                    found = n; break
            if found is None: continue
            n = found; S = 2 * a * n - b * g(n); h = (S - target_S) // K_cyc
            p = apply_word(O((g(n), 1, n, 1)), word * h)
            if resp(p, a, b) == r and p[1] == 1 and p[3] == 1: hits += 1; maxh = max(maxh, h)
        print(f"slope {a}/{b}: k={k} B={B} A={A} W={word} R={R} K_lin={K_lin} K_cyc={K_cyc} | drops verified: {ok} | recipe: {hits}/{tried} random targets of parity {par} in [-10000, 10000] hit exactly (max repetitions {maxh})")
    # (3) mod-4 confinement for b even
    print("b even: residues mod 4 of positive / negative values over odd n < 2^18")
    for (a, b) in [(3, 4), (5, 8), (7, 10), (7, 12), (9, 14), (11, 16), (9, 16), (7, 8), (5, 6), (9, 10), (11, 12), (13, 16), (11, 14), (13, 14), (15, 16)]:
        pos = set(); neg = set(); negsmall = set()
        for n in range(3, 1 << 18, 2):
            v = a * n - b * g(n)
            if v > 0: pos.add(v % 4)
            elif v < 0:
                neg.add(v % 4)
                if v < -200: negsmall.add(v % 4)
        pred = (-(a + b)) % 4 if 4 * a <= 3 * b else None
        predneg = ((2 - a) % 4 if b % 4 == 0 else (-a) % 4) if 4 * a > 3 * b else None
        print(f"  {a}/{b} (a/b {'<=' if 4*a <= 3*b else '>'} 3/4): positive classes {sorted(pos)} (predicted single class {pred}); negative classes {sorted(neg)}, below -200: {sorted(negsmall)} (predicted single class {predneg})")

main()
