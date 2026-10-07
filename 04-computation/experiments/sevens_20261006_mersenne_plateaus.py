#!/usr/bin/env python3
"""The Mersenne line under U(x) = oddpart(3x+1): ternary repunits, the reset pairing, odd shift distances,
2-adic periodicity of switches, and plateaus of the odd-step stopping time (mac-mini, 2026-10-06).

M_a = 2^a - 1;  R_a = (3^a - 1)/2 = 3 R_(a-1) + 1 (ternary repunit, the undivided 3x+1 chain from 0);
sigma(n) = number of odd steps (U-steps) from n to 1 (OEIS A390816 for n = M_a; orbits stopped at 1).

  1. U^a(M_a) = oddpart(3^a - 1) = oddpart(R_a)  (a >= 2): the binary repunit passes through the ternary repunit.
  2. Pair gluing: for even a = 2k >= 4, U(R_(a-1)) = oddpart(R_a), so sigma(M_(2k)) = sigma(M_(2k-1))
     (the reset switch of checked_switch_phase19 (3) on the Mersenne line; standard-step form: A193688,
     a(2n) = 1 + a(2n-1), Herrera 2023).  For odd a, R_a = 3 * 2^v * oddpart(R_(a-1)) + 1 with v = 1 + v_2(a-1):
     a debt state (e, h) = (1, v) of checked_switch (4).
  3. Odd shift distance: if M_a (a odd) merges at equal odd-step time with M_b, b < a odd, it also merges with
     M_(b+1); hence the least shift D = a - b is odd, and the partner M_b (b even) = 3 * (4^(b/2) - 1)/3 is a
     multiple of 3 (a leaf of the backward tree) whose cofactor (4^(b/2) - 1)/3 maps to 1 in one step.
  4. 2-adic periodicity: for a >= 2 and K >= 3 the reduced 2-adic word of M_a (after its run of a - 1 ones; words
     continued past 1 with U(1) = 1, exponent 2) up to total K depends only on a mod 2^(K-2) (it is the word of
     2*3^(a-1) - 1 mod 2^(K+1)).  Hence a collision switch
     (THM-4555) M_a0 -> M_(a0-D) with colliding prefixes of total K propagates to every a = a0 (mod 2^(K-2)),
     a > D, with merge step (a - 1) + |u|.
  5. Certified density: the classes mod 2^(K-2) carrying such a collision (D odd <= 31) have density
     0.0156, ..., 0.1199 among odd a for K = 9, ..., 20 (a rigorous lower bound for the density of switching a).
  6. Plateaus (FINITE-EXACT): sigma(M_a), a <= AMAX, takes very few values (104 for a <= 6000, 58 for a <= 1200).
  7. Debt resolution grows with size (NUMERICAL, seeded): a random first-reset-2 source n of B bits merges at
     equal odd-step time with (n-1)/2 with frequency rising from about 0.19 (B = 16) to about 0.8 (B = 2048).
Run: python3 sevens_20261006_mersenne_plateaus.py [AMAX]   (AMAX = 1200: about 10 s)
"""
import sys, random, math
from collections import defaultdict

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg, flush=True)
    OK &= bool(cond)


def U2(x):
    y = 3 * x + 1
    e = (y & -y).bit_length() - 1
    return y >> e, e


def U(x):
    return U2(x)[0]


def oddpart(x):
    return x >> ((x & -x).bit_length() - 1)


def v2(x):
    return (x & -x).bit_length() - 1


def sigma(n):
    s = 0
    while n != 1:
        n = U(n)
        s += 1
    return s


def orbit(n):
    xs = [n]
    while xs[-1] != 1:
        xs.append(U(xs[-1]))
    return xs


AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1200

print("1-2. ternary repunits and the reset pairing")
good = True
for a in range(2, 200):
    x = 2 ** a - 1
    for _ in range(a):
        x = U(x)
    good &= x == oddpart(3 ** a - 1)
check(good, "U^a(2^a - 1) = oddpart(3^a - 1) for 2 <= a < 200")
R = lambda a: (3 ** a - 1) // 2
check(all(R(a) == 3 * R(a - 1) + 1 for a in range(1, 300)), "R_a = 3 R_(a-1) + 1 (R_0 = 0)")
check(all(U(R(2 * k - 1)) == oddpart(R(2 * k)) for k in range(2, 300)), "U(R_(2k-1)) = oddpart(R_(2k)) for 2 <= k < 300 (R_(2k-1) is odd)")
check(all(R(a) == 3 * 2 ** (1 + v2(a - 1)) * oddpart(R(a - 1)) + 1 for a in range(3, 300, 2)),
      "odd a: R_a = 3 * 2^v * oddpart(R_(a-1)) + 1 with v = 1 + v_2(a-1) (a debt state (1, v))")
sig = {a: sigma(2 ** a - 1) for a in range(1, AMAX + 1)}
check(all(sig[2 * k] == sig[2 * k - 1] for k in range(2, AMAX // 2 + 1)), f"sigma(M_(2k)) = sigma(M_(2k-1)) for 2 <= k <= {AMAX // 2}")
oeis = [0, 2, 5, 5, 39, 39, 15, 15, 20, 20, 56, 56, 56, 56, 44, 44, 80, 80, 61, 61, 109, 109, 174, 174, 162, 162, 138, 138, 162, 162,
        162, 162, 191, 191, 191, 191, 191, 191, 191, 191, 191, 191, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 309, 309, 210, 210, 309]
check([sig[a] for a in range(1, 58)] == oeis, "agrees with the 57 data terms of OEIS A390816")

print("3. odd shift distances; partners are multiples of 3")
orb = {a: orbit(2 ** a - 1) for a in range(1, min(AMAX, 400) + 1)}


def eq_meet(a, b):
    xa, xb = orb[a], orb[b]
    for j in range(1, min(len(xa), len(xb))):
        if xa[j] == xb[j] and xa[j] != 1:
            return j
    return None


least, lemma = {}, True
for a in range(3, min(AMAX, 400) + 1, 2):
    for b in range(a - 1, 0, -1):
        j = eq_meet(a, b)
        if j is not None:
            least.setdefault(a, a - b)
            if b % 2 == 1 and b + 1 < a:
                lemma &= eq_meet(a, b + 1) is not None
    if a in least:
        b = a - least[a]
        lemma &= (2 ** b - 1) % 3 == 0 and U((4 ** (b // 2) - 1) // 3) == 1
odd = list(range(3, min(AMAX, 400) + 1, 2))
check(lemma and all(D % 2 == 1 for D in least.values()),
      f"odd a <= {min(AMAX, 400)}: every odd partner b also has partner b + 1; least D always odd; the least partner 2^b - 1 = 3 (4^(b/2) - 1)/3")
print(f"   odd a in [3, {min(AMAX, 400)}] with an equal-time Mersenne partner: {len(least)} of {len(odd)}; in [3, 121]: "
      f"{sum(1 for a in least if a <= 121)} of 60")

print("4. 2-adic periodicity of the reduced words and propagation of switches")


def reduced(a, L):
    # the 2-adic word (no stop at 1: U(1) = 1 with exponent 2)
    w, x = [], 2 ** a - 1
    while len(w) < a - 1 + L:
        x, e = U2(x)
        w.append(e)
    return w[a - 1:]


def upto(w, K):
    s, out = 0, []
    for c in w:
        if s + c > K:
            break
        s += c
        out.append(c)
    return out


good = all(upto(reduced(a, 40), K) == upto(reduced(a + 2 ** (K - 2), 40), K) for K in range(3, 12) for a in range(3, 250))
check(good, "reduced (2-adic) word of 2^a - 1 up to total K depends only on a mod 2^(K-2) (3 <= K <= 11, 3 <= a < 250; words continued past 1 with U(1) = 1)")


def eq_meet_big(a, b):
    x, y, j = 2 ** a - 1, 2 ** b - 1, 0
    while x != 1 and y != 1:
        x, y, j = U(x), U(y), j + 1
        if x == y:
            return j if x != 1 else None
    return None


prop = [(95, 1, 3, 9), (159, 1, 3, 10), (31, 1, 3, 11)]   # (a0, D, |u|, K) found by the class search below
good = all(eq_meet_big(a0 + t * 2 ** (K - 2), a0 + t * 2 ** (K - 2) - D) == a0 + t * 2 ** (K - 2) - 1 + p
           for a0, D, p, K in prop for t in range(0, 4))
check(good, "collisions (a0, D, K) = (95, 1, 9), (159, 1, 10), (31, 1, 11) propagate along a0 + t 2^(K-2), t = 0..3, merging at step (a - 1) + |u|")

print("5. certified density of switching exponents (classes mod 2^(K-2), D odd <= 31)")


def prefix_table(c, K):
    M = 1 << (K + 1)
    x = (2 * pow(3, (c - 1) % (1 << (K - 2)), M) - 1) % M
    prec, S, out, N, E, L = K + 1, 0, {}, None, None, 0
    while True:
        y = (3 * x + 1) % (1 << prec)
        if y == 0:
            break
        e = (y & -y).bit_length() - 1
        if S + e > K:
            break
        S += e
        L += 1
        N, E = ((-1, e - 1) if N is None else (3 * N + (1 << E), E + e))
        out[S] = (L, N)
        x, prec = y >> e, prec - e
    return out


dens = {}
for K in range(9, 21):
    P = 1 << (K - 2)
    tab = [prefix_table(c, K) for c in range(P)]
    good = 0
    for c in range(1, P, 2):
        if any(T in tab[(c - D) % P] and tab[(c - D) % P][T] == (L + D, N)
               for D in range(1, 32, 2) for T, (L, N) in tab[c].items()):
            good += 1
    dens[K] = good / (P // 2)
check(abs(dens[9] - 1 / 64) < 1e-12 and abs(dens[20] - 0.1199) < 1e-4 and all(dens[K] <= dens[K + 1] for K in range(9, 20)),
      "certified densities " + ", ".join(f"K={K}: {d:.4f}" for K, d in dens.items()))

print("6. plateaus of sigma on the Mersenne line")
for A in (100, 400, 1200):
    if A <= AMAX:
        print(f"   a <= {A}: {len({sig[a] for a in range(1, A + 1)})} distinct values of sigma(2^a - 1)")
lev = defaultdict(list)
for a in range(1, AMAX + 1):
    lev[sig[a]].append(a)
check(all(len(v) % 2 == 0 for s, v in lev.items() if min(v) >= 3), "every level set {a : sigma(M_a) = s} (a >= 3) is a union of reset pairs {2k-1, 2k}")

print("7. debt resolution grows with size (NUMERICAL, seeded)")
rng = random.Random(2026)
rates = []
for B in (16, 64, 256):
    tot = hit = 0
    while tot < 400:
        n = rng.getrandbits(B) | (1 << (B - 1)) | 1
        r = v2(n + 1) - 1
        if r < 1 or 1 + v2(3 ** r * ((n + 1) >> (r + 1)) - 1) < 3:
            continue
        tot += 1
        x, y = n, (n - 1) // 2
        while x != 1 and y != 1:
            x, y = U(x), U(y)
            if x == y:
                hit += x != 1
                break
    rates.append(hit / tot)
check(rates[0] < rates[1] < rates[2], "equal-time merge rate of reset-2 sources with (n-1)/2: " + ", ".join(f"B={B}: {p:.3f}" for B, p in zip((16, 64, 256), rates)))
print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
