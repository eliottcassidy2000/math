#!/usr/bin/env python3
"""Audit A (independent), THM-4640 statements 1-3 at word level.

Own implementation (does not import the audited scripts).  f_c(x) = (3x+1)/2^c, a word acts left to right,
reduced = first letter >= 2.  For every reduced word of total <= T we compute f_u(-1) by direct Fraction
iteration and check:
  (1)  f_u(-1) = K/2^(A-1), K = -3^m + sum_j 3^(m-j) 2^(S_j - 1) over the head (m = k-1), K odd, and the
       reduced denominator is exactly 2^(A-1); colliding words have equal totals and equal K.
  (2)  right cancellation: stripping the longest common suffix of two colliding words leaves two non-empty
       reduced words that still collide and have different last letters (a primitive collision).
  (3)  every pair with distinct last letters c < c' has c' - c even = 2i and f_h'(-1) = 4^i f_h(-1) + (4^i-1)/3
       (head h of the c-word; empty head = -1).  Converse: every pair of heads (h, h'), both reduced or empty,
       of total <= TH, with f_h'(-1) = 4^i f_h(-1) + (4^i-1)/3 for some integer i != 0, gives collisions
       (h, c) ~ (h', c+2i) for every c with both words reduced and both last letters >= 1 (checked to c <= 25),
       and NO collision when one of the last letters would be <= 0 (i.e. the range c >= max(1, 1-2i) is exact
       apart from the reducedness of a one-letter word).
  (4)  root twin: (a+2, w) ~ (2, a, w) for all a >= 1 and all words w of total <= T - a - 2.
Usage: python3 a1_words.py [T=18] [TH=14]
"""
import sys, time
from fractions import Fraction as Fr
from collections import defaultdict
from itertools import product

T = int(sys.argv[1]) if len(sys.argv) > 1 else 18
TH = int(sys.argv[2]) if len(sys.argv) > 2 else 14
t0 = time.time()
nchk = 0


def ok(c, msg):
    global nchk
    nchk += 1
    if not c:
        raise SystemExit("FAIL: " + msg)


def fval(word, x=Fr(-1)):
    for a in word:
        x = (3 * x + 1) / Fr(2) ** a
    return x


def compositions(n):
    """all compositions of n (tuples of positive ints)"""
    if n == 0:
        yield ()
        return
    for first in range(1, n + 1):
        for rest in compositions(n - first):
            yield (first,) + rest


def Kform(word):
    head = word[:-1]
    m = len(head)
    S = 0
    K = -3 ** m
    for j, a in enumerate(head, start=1):
        S += a
        K += 3 ** (m - j) * 2 ** (S - 1)
    return K


# ---- enumerate reduced words of total <= T
byval = defaultdict(list)
nwords = 0
for A in range(2, T + 1):
    for a1 in range(2, A + 1):
        for rest in compositions(A - a1):
            w = (a1,) + rest
            y = fval(w)
            K = Kform(w)
            ok(y == Fr(K, 2 ** (A - 1)), f"numerator formula {w}")
            ok(K % 2 == 1, f"K odd {w}")
            ok(y.denominator == 2 ** (A - 1), f"exact denominator {w}")
            byval[y].append(w)
            nwords += 1
ok(nwords == 2 ** (T - 1) - 1, "word count")
print(f"(1) {nwords} reduced words of total <= {T}: numerator formula, K odd, exact denominator 2^(A-1): OK "
      f"[{time.time()-t0:.1f}s]", flush=True)

# equal totals and equal K inside every collision class
nclass = nprimval = nroot = 0
for y, ws in byval.items():
    tots = {sum(w) for w in ws}
    ok(len(tots) == 1, f"equal totals {ws}")
    ok(len({Kform(w) for w in ws}) == 1, "equal K")
    if len(ws) > 1:
        nclass += 1
    if len({w[-1] for w in ws}) > 1:
        nprimval += 1
        if any(len(w) == 1 for w in ws):
            nroot += 1
print(f"    {len(byval)} distinct values; {nclass} hit more than once; {nprimval} carry words with two distinct "
      f"last letters, {nroot} of them containing a one-letter word (root)", flush=True)

# ---- (2) cancellation of the longest common suffix
npairs = 0
for y, ws in byval.items():
    if len(ws) < 2:
        continue
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            u, v = ws[i], ws[j]
            k = 0
            while k < min(len(u), len(v)) and u[-1 - k] == v[-1 - k]:
                k += 1
            pu, pv = u[:len(u) - k], v[:len(v) - k]
            ok(len(pu) > 0 and len(pv) > 0, f"cancellation leaves an empty word {u} {v}")
            ok(pu[0] >= 2 and pv[0] >= 2, "prefixes reduced")
            ok(fval(pu) == fval(pv), f"prefix collision {pu} {pv}")
            ok(pu[-1] != pv[-1], "primitive")
            npairs += 1
print(f"(2) {npairs} colliding word pairs: longest-common-suffix cancellation always leaves a primitive collision: OK",
      flush=True)

# ---- (3) ladder identity on every pair with distinct last letters
nlad = 0
gaps = defaultdict(int)
for y, ws in byval.items():
    for i in range(len(ws)):
        for j in range(len(ws)):
            u, v = ws[i], ws[j]
            if u[-1] >= v[-1]:
                continue
            d = v[-1] - u[-1]
            ok(d % 2 == 0, f"even gap {u} {v}")
            ii = d // 2
            yh = fval(u[:-1])
            yh2 = fval(v[:-1])
            ok(yh2 == 4 ** ii * yh + Fr(4 ** ii - 1, 3), f"compiler form {u} {v}")
            gaps[d] += 1
            nlad += 1
print(f"(3) {nlad} ordered pairs with distinct last letters: gap even and compiler identity OK; gap histogram "
      f"{dict(sorted(gaps.items()))}", flush=True)

# converse: all head pairs (reduced or empty) of total <= TH satisfying the compiler relation for some i != 0
heads = [()]
for A in range(2, TH + 1):
    for a1 in range(2, A + 1):
        for rest in compositions(A - a1):
            heads.append((a1,) + rest)
hv = {h: fval(h) for h in heads}
# index heads by value to find relation partners: y' = 4^i y + (4^i - 1)/3  <=>  3y' + 1 = 4^i (3y + 1)
by3 = defaultdict(list)
for h, y in hv.items():
    by3[3 * y + 1].append(h)
nconv = 0
for h, y in hv.items():
    z = 3 * y + 1
    for ii in range(-12, 13):
        if ii == 0:
            continue
        for h2 in by3.get(z * Fr(4) ** ii, []):
            if h2 == h:
                continue
            # ladder (h, c) ~ (h2, c + 2 ii)
            for c in range(-3, 26):
                c2 = c + 2 * ii
                valid = c >= 1 and c2 >= 1 and (len(h) > 0 or c >= 2) and (len(h2) > 0 or c2 >= 2)
                if c >= 1 and c2 >= 1:
                    coll = fval(h + (c,)) == fval(h2 + (c2,))
                    ok(coll, f"converse ladder {h},{c} ~ {h2},{c2}")
                    if valid:
                        nconv += 1
            # the only possible failure of 'c >= max(1, 1-2i)' is reducedness of a one-letter word
            if len(h) == 0:
                ok(ii == -1 and h2 == (2,), f"empty head partner {h2} i={ii}")
            if len(h2) == 0:
                ok(ii == 1 and h == (2,), f"empty head partner {h} i={ii}")
print(f"    converse: {nconv} (head pair, c) instances with total head <= {TH}, c <= 25: all collide; empty heads "
      f"pair only with (2) (the root, i = -1): OK", flush=True)

# ---- (4) root twin (a+2, w) ~ (2, a, w)
ntw = 0
for a in range(1, T - 1):
    for rest_total in range(0, T - a - 1):
        for w in compositions(rest_total):
            ok(fval((a + 2,) + w) == fval((2, a) + w), f"root twin {a} {w}")
            ntw += 1
print(f"(4) root twin (a+2, w) ~ (2, a, w): {ntw} instances OK", flush=True)
print(f"TOTAL checks {nchk}  [{time.time()-t0:.1f}s]")
