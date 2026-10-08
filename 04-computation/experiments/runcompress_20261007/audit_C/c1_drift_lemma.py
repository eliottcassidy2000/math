#!/usr/bin/env python3
"""audit_C item 1: THM-4603 (1), the drift lemma, re-tested with exact rational arithmetic (independent code).

Claim: w, w' U-words, c = c_w, c' = c_w' (cycle points at which w, w' are read), Lam = lcm(Sum w, Sum w').
If u - c = 3^k (v - c') and v2(v - c') > j*Lam then for jLam Terras steps u reads w^(jLam/Sum w), v reads w'^(jLam/Sum w'),
T^(jLam)(u) - c = 3^k' (T^(jLam)(v) - c') with k' = k + jLam(|w|/Sum w - |w'|/Sum w'), and the runs end together.

Tests (all exact, Fractions in Q cap Z_(2)):
  A. random word pairs (length <= 4, letters <= 6), k in [-9, 9] (negative k included), j in {1,2,3},
     precision N = v2(v - c') in {jLam+1 (boundary), jLam+random}; u = c + 3^k (v - c').
     Checks: U-words, exact relation after jLam steps, debt formula (sum of par(u)-par(v)), both orbits leave their cycles
     at exactly Terras time N (runs end together).
  B. boundary N = jLam exactly: relation and drift still hold after jLam Terras steps, but the last U-letter is longer
     (the U-word statement needs '>' exactly; the Terras-step statement only needs '>=').
  C. integer instances: integral cycle points {-1, 1, -5, -17} (k >= 0 and k < 0), u, v actual integers.
  D. the halving point 0 (not a U-word cycle point): Terras form of the lemma with c' = 0 (density 0), as used for the
     D <= 2 collapse (state (J, 1): u - 1 = 3^J (v - 0)).
  E. 'words not read from their cycle point': with c a different odd point of the same cycle, u reads the rotation, not w.
"""
import random
from fractions import Fraction as Fr
from math import lcm

rnd = random.Random(20261007)


def par(q):
    return q.numerator & 1


def T(q):
    return (3 * q + 1) / 2 if par(q) else q / 2


def v2q(q):
    """2-adic valuation of a rational with odd denominator (q != 0)"""
    n = q.numerator
    assert q.denominator % 2 == 1
    return (n & -n).bit_length() - 1


def cycle_point(w):
    m, S, B, Q = len(w), sum(w), 0, 1
    for a in w:
        B = 3 * B + Q
        Q *= 2 ** a
    return Fr(B, 2 ** S - 3 ** m)


def uword(x, nsteps):
    """U-letters read by x (odd) in exactly nsteps Terras steps; returns (letters, exact) where exact means the value at
    Terras time nsteps is odd (so the last letter is complete and exact)."""
    letters = []
    t = 0
    while t < nsteps:
        assert par(x) == 1
        x = T(x)
        t += 1
        a = 1
        while t < nsteps and par(x) == 0:
            x = T(x)
            t += 1
            a += 1
        letters.append(a)
    return letters, par(x) == 1, x


def first_exit(x, c, maxsteps):
    """first Terras time at which the parity of x's orbit differs from that of c's orbit"""
    for s in range(maxsteps):
        if par(x) != par(c):
            return s
        x, c = T(x), T(c)
    return None


def random_word():
    L = rnd.randint(1, 4)
    while True:
        w = tuple(rnd.randint(1, 6) for _ in range(L))
        if sum(w) <= 12:
            return w


counts = dict(A=0, B=0, C=0, D=0, E=0)

# ---------------- A ----------------
for trial in range(1500):
    w, wp = random_word(), random_word()
    c, cp = cycle_point(w), cycle_point(wp)
    S, Sp = sum(w), sum(wp)
    Lam = lcm(S, Sp)
    if Lam > 60:
        continue
    j = rnd.choice((1, 2, 3))
    k = rnd.randint(-9, 9)
    N = j * Lam + (1 if trial % 2 == 0 else rnd.randint(1, 12))
    r = Fr(2 * rnd.getrandbits(30) + 1, 2 * rnd.randint(0, 40) + 1)       # 2-adic unit
    v = cp + 2 ** N * r
    u = c + Fr(3) ** k * (v - cp)
    assert v2q(v - cp) == N and v2q(u - c) == N
    n = j * Lam
    lu, exu, Tu = uword(u, n)
    lv, exv, Tv = uword(v, n)
    assert lu == list(w) * (n // S) and exu, (w, lu)
    assert lv == list(wp) * (n // Sp) and exv, (wp, lv)
    rho = Fr(len(w), S) - Fr(len(wp), Sp)
    kp = k + n * rho
    assert kp.denominator == 1
    kp = int(kp)
    assert Tu - c == Fr(3) ** kp * (Tv - cp)
    # debt bookkeeping step by step
    kk, x, y = k, u, v
    for s in range(n):
        kk += par(x) - par(y)
        x, y = T(x), T(y)
    assert kk == kp and x == Tu and y == Tv
    # runs end together: both leave their cycles exactly at Terras time N
    assert first_exit(u, c, N + 5) == N and first_exit(v, cp, N + 5) == N
    counts['A'] += 1

# ---------------- B (boundary N = jLam) ----------------
for trial in range(600):
    w, wp = random_word(), random_word()
    c, cp = cycle_point(w), cycle_point(wp)
    S, Sp = sum(w), sum(wp)
    Lam = lcm(S, Sp)
    if Lam > 60:
        continue
    j = rnd.choice((1, 2))
    k = rnd.randint(-9, 9)
    n = j * Lam
    r = Fr(2 * rnd.getrandbits(30) + 1, 2 * rnd.randint(0, 40) + 1)
    v = cp + 2 ** n * r
    u = c + Fr(3) ** k * (v - cp)
    lu, exu, Tu = uword(u, n)
    lv, exv, Tv = uword(v, n)
    # Terras parities agree for the first n steps, so the relation and the drift hold after n steps ...
    kp = int(k + n * (Fr(len(w), S) - Fr(len(wp), Sp)))
    assert Tu - c == Fr(3) ** kp * (Tv - cp)
    # ... but the value at time n is even, so the last U-letter is not exact (it is longer)
    assert (not exu) and (not exv)
    counts['B'] += 1

# ---------------- C (actual integers) ----------------
int_words = {(1,): -1, (2,): 1, (1, 2): -5, (1, 1, 1, 2, 1, 1, 4): -17}
for w, val in int_words.items():
    assert cycle_point(w) == val, (w, cycle_point(w))
for trial in range(800):
    w, wp = rnd.choice(list(int_words)), rnd.choice(list(int_words))
    c, cp = cycle_point(w), cycle_point(wp)
    S, Sp = sum(w), sum(wp)
    Lam = lcm(S, Sp)
    j = rnd.choice((1, 2))
    k = rnd.randint(-7, 7)
    n = j * Lam
    N = n + rnd.randint(1, 8)
    # v - c' = 2^N * 3^max(0,-k) * odd  -> u integer
    r = (2 * rnd.getrandbits(40) + 1) * 3 ** max(0, -k)
    v = cp + 2 ** N * r
    u = c + Fr(3) ** k * (v - cp)
    assert u.denominator == 1 and v.denominator == 1
    lu, exu, Tu = uword(u, n)
    lv, exv, Tv = uword(v, n)
    assert lu == list(w) * (n // S) and lv == list(wp) * (n // Sp) and exu and exv
    kp = int(k + n * (Fr(len(w), S) - Fr(len(wp), Sp)))
    assert Tu - c == Fr(3) ** kp * (Tv - cp)
    assert first_exit(u, c, N + 3) == N and first_exit(v, cp, N + 3) == N
    counts['C'] += 1

# ---------------- D (halving point 0) ----------------
for trial in range(300):
    w = rnd.choice([(2,), (1,), (1, 2), (3,), (1, 3)])
    c = cycle_point(w)
    S = sum(w)
    k = rnd.randint(-6, 8)
    n = S * rnd.randint(1, 4)
    N = n + rnd.randint(0, 6)
    r = Fr(2 * rnd.getrandbits(30) + 1, 2 * rnd.randint(0, 20) + 1)
    v = 2 ** N * r                      # near 0 (even)
    u = c + Fr(3) ** k * v
    kk, x, y = k, u, v
    for s in range(n):
        kk += par(x) - par(y)
        x, y = T(x), T(y)
    assert kk == k + n * Fr(len(w), S)          # drift = dens(c) - 0
    assert x - c == Fr(3) ** kk * (y - 0)
    counts['D'] += 1
# the D <= 2 collapse: state (J, 1) is u - 1 = 3^J (v - 0): drift +1/2
J = 5
u, v, k = Fr(1) + 3 ** J * 2 ** 40 * 7, Fr(2 ** 40 * 7), J
for s in range(30):
    k += par(u) - par(v)
    u, v = T(u), T(v)
assert k == J + 15

# ---------------- E (rotation) ----------------
w = (1, 2, 3)
c = cycle_point(w)
# odd points of the cycle of c and the words read there
pts = []
x = c
for s in range(sum(w)):
    if par(x):
        pts.append(x)
    x = T(x)
for p in pts[1:]:
    u = p + 2 ** 30 * 7
    lu, exu, _ = uword(u, sum(w))
    assert lu != list(w) and sorted(lu) == sorted(w)   # a rotation of w, not w
    counts['E'] += 1

print("drift lemma exact tests:", counts)
print("A: random word pairs, k in [-9,9], N = v2(v-c') > jLam (incl. jLam+1): U-words, relation, debt formula, common exit time N")
print("B: N = jLam exactly: relation and drift still hold after jLam Terras steps, last U-letter not exact")
print("C: integral cycle points {-1,1,-5,-17}, actual integers, k in [-7,7]")
print("D: Terras form with the halving point 0 (density 0), incl. the (J,1) collapse drift +1/2")
print("E: a different odd point of the cycle reads a rotation of w")
print("ALL CHECKS PASSED")
