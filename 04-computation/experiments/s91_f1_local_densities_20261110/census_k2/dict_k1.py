#!/usr/bin/env python3
"""S91 helper (census_k2), step 1: the dictionary between the K = 1 words of THM-4916 and the F1 triples.

Canon model (THM-4916, rc_lf.c, audit/census_mechanism.c): states (D, E, U, d), root t = (1, 0, 2, 1), letters
    O: (D, 0, 2U, d)        Z: (D - E + U, d, 2U - d, d)        J: (D - E + U, D - E, D + U, D).
A counted word ends in Z or J; its value is the first coordinate D of its final state (= F = D - E + U of the state
before the last letter).  M_K(N) = number of counted words with at most K letters J and value N.
Tail code m of an O/Z word ending in Z: the letters are the binary digits of m from the LOWEST to the top one,
Z = 1, O = 0 (m = 0 is the empty tail).  Q(1) = 0, Q(2m) = Q(m), Q(2m + 1) = Q(m) + m + (m mod 2) (Q = 2q).

F1 (CONTEXT_S91): N = x1 x2 + 2^s y1 y2 - [s = 0],  x(h) = 1 + 2h - Q(h),  y(h) = 2h + 1, ordered triples (h1, h2, s).

Statements checked here (all exact, universe: every counted word of value < X = 2^14, K = 1):
 D1  every J-free state is (D, E, U, d) = (x(g), [a = 0 and g >= 1], 2^a y(g), 1) for a unique pair (g, a), and
     A -> (g, a) is a bijection from the O/Z words onto {g >= 1, a >= 0} u {g = 0, a >= 1}   (a = trailing O's of A,
     plus one for the empty core: the root t carries one padding step).
 D2  tail m -> (c, h2):  m even: c = 0, h2 = m/2;  m odd: c = v2(m + 1) >= 1, h2 = m >> (c + 1); this is a bijection
     from m >= 0 onto {c >= 0, h2 >= 0}, and value(A J tail m) = D x(h2) + 2^c U y(h2) - [c = 0] E      (THM-4916 (1b)).
 D3  hence value = x1 x2 + 2^s y1 y2 - [s = 0] with (h1, h2, s) = (g, h2, a + c): odd tails lie on the sheet s >= 1;
     an even tail lies on the sheet s = 0 exactly when A ends in Z (a = 0), and on s >= 1 otherwise.
 D4  multiplicity: the one-J words over a triple (h1, h2, s) number s + 1 - [h1 = 0]      (CONTEXT: s + 1; the
     boundary h1 = 0 loses one because the root already carries the padding a = 1).
 D5  M_1(N) = #{h >= 1: x(h) = N} + sum over triples with value N of (s + 1 - [h1 = 0]),  for every odd N < X,
     against a letter-by-letter enumeration that uses none of the formulas above.
"""
import sys
from collections import Counter

X = 1 << 14


def letter_O(s):
    D, E, U, d = s
    return (D, 0, 2 * U, d)


def letter_Z(s):
    D, E, U, d = s
    return (D - E + U, d, 2 * U - d, d)


def letter_J(s):
    D, E, U, d = s
    return (D - E + U, D - E, D + U, D)


ROOT = (1, 0, 2, 1)

Qtab = [0, 0]
for m in range(2, 1 << 18):
    half = m >> 1
    Qtab.append(Qtab[half] + (half + (half & 1) if m & 1 else 0))


def xval(h):
    return 1 + 2 * h - Qtab[h]


def yval(h):
    return 2 * h + 1


def brute_M(K):
    """letter-by-letter: prefixes w with F(t.w) < X; counted words wZ and (if #J(w) < K) wJ, both of value F(t.w)."""
    M = Counter()
    stack = [(ROOT, 0)]
    nwords = 0
    while stack:
        s, nj = stack.pop()
        D, E, U, d = s
        F = D - E + U
        if F >= X:
            continue
        M[F] += 1
        nwords += 1
        if nj < K:
            M[F] += 1
            nwords += 1
            stack.append((letter_J(s), nj + 1))
        stack.append((letter_O(s), nj))
        stack.append((letter_Z(s), nj))
    return M, nwords


def apply_tail(s, m):
    """apply the O/Z word with code m (lowest digit first) to the state s."""
    while m:
        s = letter_Z(s) if m & 1 else letter_O(s)
        m >>= 1
    return s


def main():
    ok = True
    # ---- D1: J-free states -------------------------------------------------------------------------------------
    seen = {}
    frontier = [((), ROOT)]
    nstates = 0
    d1 = True
    maxlen = 16
    for length in range(maxlen + 1):
        nxt = []
        for word, s in frontier:
            D, E, U, d = s
            a = (U & -U).bit_length() - 1
            g = ((U >> a) - 1) // 2
            trailing = 0
            for letter in reversed(word):
                if letter == 0:
                    trailing += 1
                else:
                    break
            core_empty = (trailing == len(word))
            d1 &= (d == 1) and (D == xval(g)) and (U == (yval(g) << a))
            d1 &= (E == (1 if (a == 0 and g >= 1) else 0))
            d1 &= (a == trailing + (1 if core_empty else 0)) and ((g == 0) == core_empty)
            d1 &= (g, a) not in seen
            seen[(g, a)] = word
            nstates += 1
            if length < maxlen:
                nxt.append((word + (0,), letter_O(s)))
                nxt.append((word + (1,), letter_Z(s)))
        frontier = nxt
    # surjectivity on a box: every (g, a) with bitlength(g) + a <= maxlen (g >= 1) and (0, a), 1 <= a <= maxlen + 1
    for g in range(0, 1 << 10):
        for a in range(0, 6):
            if g == 0 and a == 0:
                continue
            d1 &= (g, a) in seen
    print(f"D1 J-free states (x(g), [a=0, g>=1], 2^a y(g), 1), bijection A -> (g, a): {d1}   ({nstates} words of length <= {maxlen})")
    ok &= d1

    # ---- D2, D3, D4: one-J words by (A, tail m) ------------------------------------------------------------------
    d2 = True
    d3 = True
    tri = Counter()          # multiset of triples (h1, h2, s) produced by one-J words
    M_formula = Counter()    # values of the one-J words
    sheet_by_parity = Counter()
    seen_tail = set()
    for (g, a), word in seen.items():
        s = ROOT
        for letter in word:
            s = letter_Z(s) if letter else letter_O(s)
        D, E, U, d = s
        if D - E + U >= X:
            continue
        post = letter_J(s)
        h, e, R, dd = post[0], post[1], post[2], post[3]
        m = 0
        while h + (R - dd) * m < X:
            if m & 1:
                c = ((m + 1) & -(m + 1)).bit_length() - 1
                h2 = m >> (c + 1)
                d2 &= (xval(h2) == m - Qtab[m]) and ((yval(h2) << c) == m + 1)
            else:
                c = 0
                h2 = m >> 1
                d2 &= (xval(h2) == m + 1 - Qtab[m]) and (yval(h2) == m + 1)
            seen_tail.add((c, h2, m))
            val_direct = apply_tail(post, m)[0]                     # letter by letter
            val_seed = h + R * m - dd * Qtab[m] - e * (m & 1)       # rc_lf / decomp.c seed formula
            val_1b = D * xval(h2) + (U << c) * yval(h2) - (E if c == 0 else 0)
            sh = a + c
            val_f1 = xval(g) * xval(h2) + (yval(g) * yval(h2) << sh) - (1 if sh == 0 else 0)
            d2 &= (val_direct == val_seed == val_1b)
            d3 &= (val_1b == val_f1)
            if val_direct < X:
                tri[(g, h2, sh)] += 1
                M_formula[val_direct] += 1
                sheet_by_parity[(m & 1, 0 if sh == 0 else 1, 1 if a == 0 else 0)] += 1
            m += 1
    # tail map injective: distinct m give distinct (c, h2)
    d2 &= len({(c, h2) for (c, h2, m) in seen_tail}) == len({m for (c, h2, m) in seen_tail})
    print(f"D2 tail m <-> (c, h2) and value = D x(h2) + 2^c U y(h2) - [c = 0] E (three ways): {d2}")
    print(f"D3 value = x1 x2 + 2^s y1 y2 - [s = 0], (h1, h2, s) = (g, h2, a + c): {d3}")
    print("   one-J words of value < 2^14 by (tail parity, sheet, A ends in Z):")
    for key in sorted(sheet_by_parity):
        print(f"      tail {'odd ' if key[0] else 'even'}  sheet {'s = 0 ' if key[1] == 0 else 's >= 1'}  "
              f"A ends in {'Z' if key[2] else 'O (or A empty)'}: {sheet_by_parity[key]}")
    d3 &= all(not (par == 1 and sheet == 0) for (par, sheet, az) in sheet_by_parity)          # odd tails: s >= 1
    d3 &= all((sheet == 0) == (az == 1) for (par, sheet, az) in sheet_by_parity if par == 0)   # even: s = 0 iff A ends Z
    ok &= d2 and d3

    # D4: multiplicity of each triple with value < X
    d4 = True
    ntri = 0
    wsum = 0
    plain = 0
    M_tri = Counter()
    s = 0
    while (1 << s) + 1 - (1 if s == 0 else 0) < X:
        h1 = 0
        while True:
            x1, y1 = xval(h1), yval(h1)
            # x is not monotone; x(h1) >= h1 + 1 gives a safe stopping rule
            if (h1 + 1) + (y1 << s) - 1 >= X:
                break
            h2 = 0
            while True:
                N = x1 * xval(h2) + ((y1 * yval(h2)) << s) - (1 if s == 0 else 0)
                if N >= X:
                    # the value is not monotone in h2; x(h2) >= h2 + 1 gives a safe stopping rule
                    if x1 * (h2 + 1) + ((y1 * (2 * h2 + 1)) << s) - 1 >= X:
                        break
                else:
                    w = s + 1 - (1 if h1 == 0 else 0)
                    d4 &= tri.get((h1, h2, s), 0) == w
                    ntri += 1
                    wsum += w
                    plain += s + 1
                    M_tri[N] += w
                h2 += 1
            h1 += 1
        s += 1
    d4 &= wsum == sum(tri.values())
    print(f"D4 every triple of value < 2^14 carries exactly s + 1 - [h1 = 0] one-J words: {d4}   "
          f"({ntri} triples, {wsum} one-J words; the plain rule s + 1 would give {plain}, i.e. {plain - wsum} too many)")
    # the plain s + 1 rule fails exactly on h1 = 0; exhibit the smallest counterexample
    zero_s0 = [(0, h2, 0) for h2 in range(0, 4)]
    print("   boundary: triples (0, h2, 0) are carried by", [tri.get(k, 0) for k in zero_s0], "words (s + 1 = 1 would predict 1);"
          " e.g. (h1, h2, s) = (0, 0, 1), value 3: words", tri.get((0, 0, 1), 0), "(s + 1 = 2)")
    ok &= d4 and tri.get((0, 0, 1), 0) == 1 and all(tri.get(k, 0) == 0 for k in zero_s0)

    # ---- D5: M_1 against the letter-by-letter enumeration ------------------------------------------------------
    M1, nwords = brute_M(1)
    M0, _ = brute_M(0)
    free = Counter()
    h = 1
    while h + 1 < X:
        if xval(h) < X:
            free[xval(h)] += 1
        h += 1
    d5 = True
    for N in range(3, X, 2):
        d5 &= M0[N] == free[N]
        d5 &= M1[N] == free[N] + M_tri[N]
        d5 &= M1[N] == free[N] + M_formula[N]
    d5 &= all(N & 1 for N in M1)
    print(f"D5 M_1(N) = #J-free + sum of (s + 1 - [h1 = 0]) over triples, every odd N in [3, 2^14): {d5}   "
          f"({nwords} counted words; J-free {sum(free.values())}, one-J {sum(M_formula.values())})")
    ok &= d5

    # ---- D6: the J-free words fill the boundary exactly: M_1(N) = sum over ALL triples of (s + 1) --------------
    # the J-free word A Z with A <-> (g, a) has value x(g) + 2^a y(g) - [a = 0] = value of the triple (0, g, a); the
    # triple (0, 0, 0) (value 1) is the empty word.
    M_plain = Counter()
    for (h1, h2, sh), w in tri.items():
        M_plain[xval(h1) * xval(h2) + ((yval(h1) * yval(h2)) << sh) - (1 if sh == 0 else 0)] += sh + 1
    # triples with h1 = 0 and s = 0 carry no one-J word, so they are not keys of tri: add them
    h2 = 0
    while (h2 + 1) + (2 * h2 + 1) - 1 < X:
        N = xval(h2) + yval(h2) - 1
        if N < X:
            M_plain[N] += 1
        h2 += 1
    d6 = True
    for N in range(3, X, 2):
        d6 &= M1[N] == M_plain[N]
    d6 &= M_plain[1] == 1
    boundary = Counter()
    for (g, a) in seen:
        N = xval(g) + (yval(g) << a) - (1 if a == 0 else 0)
        if N < X:
            boundary[N] += 1
    d6 &= all(boundary[N] == free[N] for N in range(3, X, 2))
    print(f"D6 M_1(N) = sum over ALL triples of (s + 1), every odd N in [3, 2^14) (J-free word A Z <-> triple (0, g, a); "
          f"the empty word <-> (0, 0, 0)): {d6}")
    ok &= d6
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
