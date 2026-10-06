#!/usr/bin/env python3
"""
collatz_descent_set_vs_rising_cones_20261005.py

OPEN question (S6 shadows note, S7 correction, 2026-10-05): the descent set
    D = { odd n : some odd m < n has T^k(m) = n }
contains the union of the rising cones { n : n = x_w mod 3^l } over rising valuation words w
(2^A < 3^l).  Is it equal to it?  Equivalently: does every n in D have SOME smaller odd ancestor m
whose segment word m -> n is net-rising (2^A < 3^l)?  A smaller ancestor through a NON-rising word
(2^A > 3^l) exists only when m < x_w = c_w/(2^A - 3^l), the positive rational cycle point of w.

Method (exact, forward sweep).  T = Syracuse map on odd integers: n -> (3n+1)/2^v.  For every odd
m < X follow the orbit of m down to 1, maintaining (l, A) = (odd steps, halvings) since m.  Every odd
orbit value n with m < n < X is recorded in D with the flag rising := 2^A < 3^l for that segment.
Because every odd m < n with n in orbit(m) is swept, D is computed exactly below X, together with
    R = { n in D : at least one smaller ancestor word is rising }   (subset of the cone union)
and the candidates  D \ R  = elements of D all of whose smaller-ancestor words are non-rising.
Each candidate is then re-tested directly against the cone criterion: n = x_w mod 3^l for some
rising word w with l <= LMAX (so a candidate surviving this is outside every cone of depth <= LMAX).

Also recorded: for n = 2 mod 3 the depth-1 cone (w = (1), x_w = -1) covers everything, so the
question lives in n = 1 mod 3 (n = 0 mod 3 has no odd ancestor... it has: 3 | n is a leaf of the
inverse Syracuse tree only as a TARGET of nothing? multiples of 3 have no odd preimage, so they are
never in D; checked).

Reproduce: python3 collatz_descent_set_vs_rising_cones_20261005.py [X] [LMAX]
"""
import sys, time
from collections import Counter

X = int(sys.argv[1]) if len(sys.argv) > 1 else 1 << 18
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 24
CHECKS = 0
def check(c, msg):
    global CHECKS
    CHECKS += 1
    if not c:
        print("CHECK FAILED:", msg); sys.exit(1)

def v2(x):
    return (x & -x).bit_length() - 1

t0 = time.time()
inD = bytearray(X)        # inD[n] = 1 if n in D
rising_ok = bytearray(X)  # 1 if some smaller ancestor word is rising
best_nonrising = {}       # n -> (m, l, A) for a non-rising smaller ancestor (first found)
for m in range(1, X, 2):
    n = m; l = 0; A = 0
    while n != 1:
        y = 3 * n + 1
        a = v2(y)
        n = y >> a
        l += 1; A += a
        if m < n < X:
            inD[n] = 1
            if (1 << A) < 3 ** l:
                rising_ok[n] = 1
            elif n not in best_nonrising:
                best_nonrising[n] = (m, l, A)
        if n < m:
            # the orbit is below m: every later odd value n' > m it visits is also in orbit(n'' ) for the
            # smaller ancestor n'' < m that is itself swept; we may stop here.
            break
print(f"sweep to X = {X} done in {time.time()-t0:.1f}s")

D = [n for n in range(1, X, 2) if inD[n]]
check(all(n % 3 != 0 for n in D), "multiples of 3 are never in D")
check(all(inD[n] for n in range(5, X, 6)), "every n = 2 mod 3 (odd: 5 mod 6) is in D")
check(all(rising_ok[n] for n in range(5, X, 6)), "n = 2 mod 3: the one-step word (1) is rising")
cand = [n for n in D if not rising_ok[n]]
odd_count = (X - 1) // 2
units = sum(1 for n in range(1, X, 2) if n % 3)
print(f"odd n < X: {odd_count}; prime to 3: {units}; |D| = {len(D)} (density among odd {len(D)/odd_count:.6f}, among units {len(D)/units:.6f})")
print(f"n in D with a rising smaller-ancestor word: {sum(1 for n in D if rising_ok[n])}")
print(f"CANDIDATES (in D, every smaller-ancestor word non-rising): {len(cand)}")
for n in cand[:40]:
    m, l, A = best_nonrising[n]
    print(f"   n = {n}  (n mod 3 = {n % 3}); witness ancestor m = {m}, word length l = {l}, halvings A = {A}, 2^A/3^l = {2**A/3**l:.4f}")

# direct cone test for the candidates: n = x_w mod 3^l for some rising word w, l <= LMAX
# enumerate rising words by (l, A) with 2^A < 3^l; x_w mod 3^l = c_w * inv(2^A) mod 3^l where
# c_w = sum_{i=1}^{l} 3^{l-i} 2^{A_1 + ... + A_{i-1}}  (A_i = valuation of step i).  We enumerate the
# residues reached by DFS over words, which is exponential; restrict to l <= LMAX with pruning by the
# rising condition on the FULL word only (prefixes need not rise).  Count distinct cone classes.
def cone_classes(lmax):
    classes = {}   # l -> set of residues mod 3^l
    import itertools
    # DFS over valuation sequences a_1..a_l >= 1 with sum < l*log2(3)
    import math
    for l in range(1, lmax + 1):
        mod = 3 ** l
        Amax = math.floor(l * math.log2(3))   # 2^A < 3^l  <=>  A <= Amax (A integer)
        res = set()
        # iterate compositions of A <= Amax into l parts >= 1: A >= l needed
        def rec(i, Asum, c, pow2):
            # i = number of steps done, c = partial carry, pow2 = 2^(A so far)
            if i == l:
                # x_w = c / (2^A - 3^l); modulo 3^l: 2^A invertible; x_w = -c * inv(2^A)?? careful:
                # x_w (2^A - 3^l) = c  =>  x_w * 2^A = c (mod 3^l)  =>  x_w = c * inv(2^A) mod 3^l
                res.add((c * pow(pow2, -1, mod)) % mod)
                return
            remaining = l - i - 1
            for a in range(1, Amax - Asum - remaining + 1):
                # carry: c_w = sum 3^{l-i} 2^{A_1+...+A_{i-1}}; step i (1-based) contributes 3^{l-i} * 2^{A_<i}
                rec(i + 1, Asum + a, c + 3 ** (l - i - 1) * pow2, pow2 << a)
        rec(0, 0, 0, 1)
        classes[l] = res
    return classes

cc = cone_classes(min(LMAX, 20))
newc = {}
covered = set()
for l in sorted(cc):
    mod = 3 ** l
    fresh = 0
    for r in cc[l]:
        # a class mod 3^l is new if not contained in a coarser covered class
        if not any((r % (3 ** k)) in cc[k] for k in range(1, l)):
            fresh += 1
    newc[l] = fresh
print("new primitive cone classes per depth l (should match S6/S7: 1,1,0,1,0,2,8,0,28,0,124,602,0,2498,0,12319):")
print("  ", [newc[l] for l in sorted(newc)])

def in_cone(n):
    for l in sorted(cc):
        if (n % (3 ** l)) in cc[l]:
            return l
    return None
still = []
for n in cand:
    l = in_cone(n)
    if l is None:
        still.append(n)
print(f"candidates not in any rising cone of depth <= {min(LMAX,20)}: {len(still)}")
for n in still[:40]:
    m, l, A = best_nonrising[n]
    print(f"   n = {n}: witness non-rising ancestor m = {m} (l = {l}, A = {A}); n mod 27 = {n % 27}")
# Sanity: every element of D with a rising word must be in a cone of depth = that word's length; test a sample
import random
random.seed(1)
sample = random.sample([n for n in D if rising_ok[n] and n % 3 == 1], min(2000, len(D)))
miss = [n for n in sample if in_cone(n) is None]
print(f"sanity: of {len(sample)} sampled D-elements (1 mod 3) with a rising word, {len(miss)} are outside all cones of depth <= {min(LMAX,20)} (expected: only those whose rising words are longer than {min(LMAX,20)})")
print(f"\nALL {CHECKS} CHECKS PASSED")
