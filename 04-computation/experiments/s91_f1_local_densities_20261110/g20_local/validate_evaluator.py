#!/usr/bin/env python3
"""S91 helper (g20_local): validation of the two-apex evaluator g20_core.py.

Checks (all exact integers):
  1. evaluator = brute-force Hamiltonian-path count (vertex-subset DP) for ALL words of length <= 6 (5,461 words) and
     for 300 random words of lengths 7..12;
  2. evaluator = Session 194's 36-state transfer (session194/two-apex/verify.py, imported, independent code) on all
     words of length <= 7;
  3. one-apex reductions: b a sink (letters in {0,1}) and a a source (letters in {1,3}): H = x(h) = 2 psi(h) + 1 where
     the surviving apex word, stripped of its trailing 1s (sinks), is bits(h) followed by 0 (H = 1 if nothing is left);
     all words of length <= 12 in each alphabet;
  4. reversal symmetry (reverse all arcs and the chain): H(w) = H(reverse of w with letters 0 <-> 3 swapped), all words
     of length <= 7;
  5. the eight "first part seen by one apex, second part by the other" sub-languages A1* A2*:
       {0,1}*{1,3}*  (a over C1, C1 -> b;  b over C2, a -> C2):  ordinal sum, H = x1 x2;
       {0,2}*{2,3}*  (b over C1, C1 -> a;  a over C2, b -> C2):  F1, H = x1 x2 + 2^s y1 y2 - [s = 0] for the word
                      bits(h1) 0 2^s 3 revcomp-bits(h2)  (letters: b-bit of the first part, a-bit of the second);
     exhaustively for h1, h2 < 2^5, s <= 5 (F1) and all pairs of words of length <= 5 (product).
"""
import importlib.util
import itertools
import os
import random
import sys

import numpy as np

import g20_core as g

ok = True


def check(name, cond):
    global ok
    print(("PASS  " if cond else "FAIL  ") + name)
    ok &= bool(cond)


def Hn(word):
    return int(g.nb_H(np.array(word, dtype=np.int64), g.T_TO, g.T_CNT, g.ACC, g.INIT))


# 1. brute force
bad = 0
cnt = 0
for n in range(0, 7):
    for w in itertools.product(range(4), repeat=n):
        arr = np.array(w, dtype=np.int64)
        cnt += 1
        if Hn(w) != int(g.nb_bruteforce(arr)) or Hn(w) != g.H_word(w):
            bad += 1
check(f"1a. evaluator = brute force on all {cnt} words of length <= 6", bad == 0 and cnt == 5461)
rng = random.Random(91)
bad = 0
for _ in range(300):
    n = rng.randint(7, 12)
    w = [rng.randrange(4) for _ in range(n)]
    if Hn(w) != int(g.nb_bruteforce(np.array(w, dtype=np.int64))):
        bad += 1
check("1b. evaluator = brute force on 300 random words of lengths 7..12", bad == 0)

# 2. Session 194 transfer
here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", "..", "..", ".."))
p194 = os.path.join(root, "session194", "two-apex", "verify.py")
if os.path.exists(p194):
    spec = importlib.util.spec_from_file_location("s194verify", p194)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    bad = 0
    cnt = 0
    for n in range(0, 8):
        for w in itertools.product(range(4), repeat=n):
            cnt += 1
            if Hn(w) != mod.value(w):
                bad += 1
    check(f"2. evaluator = Session 194 transfer on all {cnt} words of length <= 7", bad == 0)
else:
    print("SKIP  2. session194/two-apex/verify.py not found at", p194)


# 3. one-apex reductions
def one_apex_value(bits):
    bits = list(bits)
    while bits and bits[-1] == 1:
        bits.pop()
    if not bits:
        return 1
    h = 0
    for b in bits[:-1]:
        h = 2 * h + b
    return g.x_of(h)


bad = 0
cnt = 0
for n in range(0, 13):
    for bits in itertools.product(range(2), repeat=n):
        cnt += 1
        want = one_apex_value(bits)
        if Hn(bits) != want:                        # letters 0/1: a-bit = bit, b a sink
            bad += 1
        if Hn([1 + 2 * b for b in bits]) != want:   # letters 1/3: a a source, b-bit = bit
            bad += 1
check(f"3. one-apex formula 2 psi(h) + 1 (b a sink; a a source) on all {cnt} binary words of length <= 12", bad == 0)

# 4. reversal symmetry
bad = 0
swap = {0: 3, 3: 0, 1: 1, 2: 2}
for n in range(0, 8):
    for w in itertools.product(range(4), repeat=n):
        if Hn(w) != Hn([swap[c] for c in reversed(w)]):
            bad += 1
check("4. H(w) = H(reverse(w) with 0 <-> 3) on all words of length <= 7", bad == 0)


# 5. sub-languages
def bits_of(h, k):
    return [(h >> (k - 1 - i)) & 1 for i in range(k)]


bad = 0
cnt = 0
for k1 in range(0, 6):
    for k2 in range(0, 6):
        for h1 in range(1 << k1):
            for h2 in range(1 << k2):
                for s in range(0, 6):
                    # first part: b-bits of h1 (letter 2 = bit 1, letter 0 = bit 0), then the core's final 0
                    w1 = [2 * b for b in bits_of(h1, k1)] + [0]
                    # second part: a over C2; reversal duality: reading from the right, letter 2 = bit 1, letter 3 = bit 0
                    w2 = [3] + [2 if b else 3 for b in reversed(bits_of(h2, k2))]
                    w = w1 + [2] * s + w2
                    want = g.x_of(h1) * g.x_of(h2) + (1 << s) * g.y_of(h1) * g.y_of(h2) - (1 if s == 0 else 0)
                    cnt += 1
                    if Hn(w) != want:
                        bad += 1
check(f"5a. F1 = {{0,2}}*{{2,3}}*: H(bits(h1) 0 2^s 3 revcomp(h2)) = x1 x2 + 2^s y1 y2 - [s=0] on {cnt} triples "
      f"(h1, h2 < 2^5 with all paddings of leading zeros, s <= 5)", bad == 0)
# degenerate ends: a leading 0 (source) and a trailing 3 (sink) never change H; 2^s alone = 0 2^s 3
bad = 0
for n in range(0, 7):
    for w in itertools.product(range(4), repeat=n):
        if Hn([0] + list(w)) != Hn(w) or Hn(list(w) + [3]) != Hn(w):
            bad += 1
check("5b. a leading letter 0 (source) and a trailing letter 3 (sink) leave H unchanged, all words of length <= 6", bad == 0)
bad = 0
for n1 in range(0, 6):
    for n2 in range(0, 6):
        for u in itertools.product(range(2), repeat=n1):
            for v in itertools.product(range(2), repeat=n2):
                w = list(u) + [1 + 2 * b for b in v]
                if Hn(w) != one_apex_value(u) * one_apex_value(v):
                    bad += 1
check("5c. {0,1}*{1,3}* is the ordinal sum: H = H(a over C1) H(b over C2), all pairs of binary words of length <= 5", bad == 0)

print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
sys.exit(0 if ok else 1)
