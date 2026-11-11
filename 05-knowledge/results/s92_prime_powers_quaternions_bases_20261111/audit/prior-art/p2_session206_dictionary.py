#!/usr/bin/env python3
"""Auditor prior-art (S92).  Is Session206's weighted-capacity source family the fixed-smaller-state family of F1?

Session206 (session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md): a physical source s = (h, R, e, d),
R > d > 0, has tail values H_s(m) = h + R m - d Q(m) - e (m mod 2), Q(0) = 0, Q(2m) = Q(m),
Q(2m+1) = Q(m) + m + (m mod 2).  Grammar (THM-4888 setting; coordinates (D, E, U, d) = (h, e, R, d)):
    O(D,E,U,d) = (D, 0, 2U, d),  Z = (D-E+U, d, 2U-d, d),  J = (D-E+U, D-E, D+U, D),  root t = (1, 0, 2, 1).
F1 (THM-5128 / THM-5134): N = x1 x2 + 2^s y1 y2 - [s = 0], (x, y)(h) = (2 psi(h) + 1, 2h + 1).

Checks (all exact, own code, nothing imported from the session):
 C1  Q = 2q on m < 2^16  (q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2));  x(m) = 1 + 2m - Q(m) (root tail).
 C2  the states of the root O/Z words that end in Z are exactly (D, E, U, d) = (x(h1), 1, y(h1), 1), h1 >= 1,
     each h1 once among words of a given length range (a bijection word <-> h1).
 C3  for such a word W (and for the root itself, h1 = 0, s >= 1) the state of  W O^s J O  is
     (h, e, R, d) = (x1 + 2^s y1 - [s = 0], 0, 2 (x1 + 2^s y1), x1):   a PHYSICAL Session206 source.
 C4  its tail H_s(m) equals the F1 value N(h1, m, s) for all m < 2^9, both by the closed formula and by running
     the grammar on the binary digits of m (least significant digit first).
 C5  gcd(R, d, p) = gcd(x1, y1, p) and h = -[s = 0] (mod gcd) for odd p: Session206's coset is THM-5134 (A)'s.
 C6  Session206's capacity c_0 + c_1 = ln(R/(R-d))/d equals the F1 second-state mass (1/x1) ln((x1 + 2^s y1)/(x1/2 + 2^s y1)).
 C7  a direct count: for three sources and the moduli 9, 15, 45, 49 the classes of H_s(m) <= X are confined to the
     coset and close to uniform on it (illustration only, X = 2^22).
"""
import math
import sys
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def psi(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return psi(h // 2) + h // 2
    k = h // 2
    return psi(k) + k // 2 + 1


def qf(h):
    return h - psi(h)


def X(h):
    return 2 * psi(h) + 1


def Y(h):
    return 2 * h + 1


def Q(m):
    if m == 0:
        return 0
    if m % 2 == 0:
        return Q(m // 2)
    k = m // 2
    return Q(k) + k + (k % 2)


def O(s):
    D, E, U, d = s
    return (D, 0, 2 * U, d)


def Z(s):
    D, E, U, d = s
    return (D - E + U, d, 2 * U - d, d)


def J(s):
    D, E, U, d = s
    return (D - E + U, D - E, D + U, D)


ROOT = (1, 0, 2, 1)
ok = True


def check(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print("%-4s %s  %s" % (name, "PASS" if cond else "FAIL", detail))


# C1
bad = [m for m in range(1 << 16) if Q(m) != 2 * qf(m)]
bad2 = [m for m in range(1 << 16) if X(m) != 1 + 2 * m - Q(m)]
check("C1", not bad and not bad2, "Q = 2q and x(m) = 1 + 2m - Q(m) for m < 2^16")

# C2: enumerate O/Z words from the root
LMAX = 11
states_endZ = {}
layer = {(): ROOT}
allwords = {(): ROOT}
for L in range(1, LMAX + 1):
    new = {}
    for w, s in layer.items():
        new[w + ("O",)] = O(s)
        new[w + ("Z",)] = Z(s)
    layer = new
    allwords.update(new)
h1_of_word = {}
c2 = True
seen_h1 = {}
for w, s in allwords.items():
    if w and w[-1] == "Z":
        D, E, U, d = s
        if not (E == 1 and d == 1 and U % 2 == 1):
            c2 = False
            continue
        h1 = (U - 1) // 2
        if D != X(h1):
            c2 = False
        h1_of_word[w] = h1
        seen_h1.setdefault(h1, []).append(w)
# every h1 in [1, 2^LMAX) should occur, exactly once
missing = [h for h in range(1, 1 << LMAX) if h not in seen_h1]
multi = [h for h, ws in seen_h1.items() if len(ws) > 1]
check("C2", c2 and not missing and not multi,
      "words ending in Z, length <= %d: %d words <-> h1 in [1, %d); state = (x(h1), 1, y(h1), 1); missing %d, repeated %d"
      % (LMAX, len(h1_of_word), 1 << LMAX, len(missing), len(multi)))
# example dictionary
for w in [("Z",), ("Z", "Z"), ("O", "Z"), ("Z", "O", "Z"), ("O", "O", "Z")]:
    print("       word %-6s -> state %s, h1 = %d, (x, y) = (%d, %d)" % ("".join(w), allwords[w], h1_of_word[w], X(h1_of_word[w]), Y(h1_of_word[w])))

# C3, C4, C5, C6
c3 = c4 = c5 = c6 = True
ncase = 0
sources = []
for w, h1 in list(h1_of_word.items()) + [((), 0)]:
    if len(w) > 7:
        continue
    for s in range(0, 5):
        if h1 == 0 and s == 0:
            # root itself is (x1, 2 y1, 0, 1): the case s = 0 is the letter O (N = x2 + y2 - 1), not a J-source
            st = O(ROOT)
            src = st
        else:
            st = allwords[w]
            if h1 == 0:
                # root = (1, 0, 2, 1) = (x1, 0, 2^1 y1, 1): already one O "inside"
                for _ in range(s - 1):
                    st = O(st)
            else:
                for _ in range(s):
                    st = O(st)
            src = O(J(st))
        x1, y1 = X(h1), Y(h1)
        want = (x1 + (1 << s) * y1 - (1 if s == 0 else 0), 0, 2 * (x1 + (1 << s) * y1), x1)
        if src != want:
            c3 = False
            print("   C3 mismatch", w, s, src, want)
        D, E, U, d = src
        hh, R, e = D, U, E
        if not (R > d > 0):
            c3 = False
        sources.append((h1, s, hh, R, e, d))
        for m in range(1 << 9):
            N = x1 * X(m) + (1 << s) * y1 * Y(m) - (1 if s == 0 else 0)
            Hs = hh + R * m - d * Q(m) - e * (m % 2)
            if N != Hs:
                c4 = False
            # grammar on the digits of m, least significant first
            st2 = src
            mm = m
            while mm:
                st2 = Z(st2) if mm & 1 else O(st2)
                mm >>= 1
            if st2[0] != N:
                c4 = False
        for p in (3, 5, 7, 9, 15, 21, 25, 27, 45, 49, 63, 75, 81, 105, 125):
            g = gcd(gcd(R, d), p)
            if g != gcd(gcd(x1, y1), p):
                c5 = False
            if (hh + (1 if s == 0 else 0)) % g != 0:
                c5 = False
        cap = math.log(R / (R - d)) / d
        f1 = math.log((x1 + (1 << s) * y1) / (x1 / 2 + (1 << s) * y1)) / x1
        if abs(cap - f1) > 1e-12:
            c6 = False
        ncase += 1
check("C3", c3, "state of W O^s J O = (x1 + 2^s y1 - [s=0], 0, 2(x1 + 2^s y1), x1), R > d > 0; %d (h1, s) cases" % ncase)
check("C4", c4, "H_s(m) = N(h1, m, s) for m < 2^9, formula and grammar run")
check("C5", c5, "gcd(R, d, p) = gcd(x1, y1, p), h = -[s=0] mod gcd, 15 odd moduli")
check("C6", c6, "capacity ln(R/(R-d))/d = (1/x1) ln((x1 + 2^s y1)/(x1/2 + 2^s y1))")

# C7 illustration
print("C7   illustration (X = 2^22): classes of H_s(m) <= X; max relative deviation from g/p on the coset, mass off the coset")
XB = 1 << 22
for (h1, s) in ((1, 1), (4, 0), (2, 2), (38, 1)):
    x1, y1 = X(h1), Y(h1)
    hh, R, e, d = x1 + (1 << s) * y1 - (1 if s == 0 else 0), 2 * (x1 + (1 << s) * y1), 0, x1
    vals = []
    m = 0
    while True:
        v = hh + R * m - d * Q(m)
        if hh + (R - d) * m > XB:
            break
        if v <= XB:
            vals.append(v)
        m += 1
    for p in (9, 15, 45, 49):
        g = gcd(gcd(R, d), p)
        cnt = [0] * p
        for v in vals:
            cnt[v % p] += 1
        off = sum(c for r, c in enumerate(cnt) if (r - hh) % g)
        on = [c for r, c in enumerate(cnt) if (r - hh) % g == 0]
        exp = len(vals) * g / p
        dev = max(abs(c - exp) / exp for c in on)
        print("       h1 = %3d s = %d (x1, y1) = (%4d, %4d) g(x1,y1) = %3d   p = %2d  g = %2d  #m = %7d  off-coset %d  max rel dev %.4f"
              % (h1, s, x1, y1, gcd(x1, y1), p, g, len(vals), off, dev))
print("ALL PASS" if ok else "SOME FAIL")
