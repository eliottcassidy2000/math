#!/usr/bin/env python3
"""Checks for the psi-coordinates of the one-apex count and for the three-letter rooted system (S, B, F).

All statements here are restatements or re-derivations of repository results (Session 16 weighted inversions,
THM-4876 (6), THM-4885/4887, the O/Z tail law of Sessions 195-205 / THM-4888, the O/Z/J grammar of Sessions 197-200);
the script exists so that the note's formulas are executable.

(1) One-apex count.  Chain v_1 -> ... -> v_m (v_i -> v_j for i < j), apex x, word u (u_i = 1 iff x -> v_i).  For a
    reduced word (u_1 = 1, u_m = 0) with binary value n = 2h (u_1 the leading bit):
        H = 2 psi(h) + 1,   psi(0) = 0, psi(2h) = psi(h) + h, psi(2h+1) = psi(h) + floor(h/2) + 1,
    and psi(h) = h - q(h), q(h) = sum over pairs of one-bits at distance d of w(d), w(1) = w(2) = 1, w(d) = 2^(d-2).
(2) G(n) := sum over zero bits p of n of floor(n / 2^(p+1)) satisfies G(2n) = G(n) + n, G(2n+1) = G(n),
    G(2h) + rho(h) = 2 psi(h) (rho = number of blocks of ones), and
        G(n) = #{ m < n : phi(m) > phi(n) },     phi = van der Corput radical inverse in base 2.
(3) Splitting law: for t < 2^j,
        psi(m 2^j + t) = psi(m) + psi(t) + kappa_j(t) m - t_(j-1) ceil(m/2),
        kappa_j(t) = 2^j - 1 - sum over one-bits p <= j - 2 of t of 2^(j-p-2).
(4) Rooted system.  State (h, R, e, u) of a pointed tournament (A, a): h = H(A),
    R = sum over X not containing a of H(A[X]) s_a(A - X), e = number of paths ending at a, u = H(A - a).
        S: new sink                         -> (h, 2R, 0, u)
        B: new vertex beating only a        -> (h + R - e, 2R - u, u, u)        (root stays)
        F: new vertex beating only a        -> (h + R - e, h + R, h - e, h)     (root moves to the new vertex)
    verified against brute force; the odd numbers below 6000 that are not an h-coordinate of the orbit of (1,1,1,1)
    are exactly the 29 values of THM-4884 / Session 197.
Run:  python3 psi_sbf.py > psi_sbf.out
"""
import itertools
import random
import sys
from fractions import Fraction


def ham_tables(n, beats):
    size = 1 << n
    cnt = [[0] * n for _ in range(size)]
    for v in range(n):
        cnt[1 << v][v] = 1
    Hs = [0] * size
    Hs[0] = 1
    for mask in range(1, size):
        row = cnt[mask]
        tot = 0
        for v in range(n):
            c = row[v]
            if c:
                tot += c
                nxt = beats[v] & ~mask
                while nxt:
                    w = nxt & -nxt
                    cnt[mask | w][w.bit_length() - 1] += c
                    nxt ^= w
        Hs[mask] = tot
    return Hs, cnt


def psi(h, memo={0: 0}):
    if h in memo:
        return memo[h]
    v = psi(h >> 1) + ((h >> 1) if h % 2 == 0 else ((h >> 1) // 2 + 1))
    memo[h] = v
    return v


def q(h):
    ones = [i for i in range(h.bit_length()) if h >> i & 1]
    s = 0
    for a in range(len(ones)):
        for b in range(a + 1, len(ones)):
            d = ones[b] - ones[a]
            s += 1 if d <= 2 else 1 << (d - 2)
    return s


def G(n):
    g = 0
    k = n
    while k:
        if k & 1 == 0:
            g += k >> 1
        k >>= 1
    return g


def rho(n):
    return len([blk for blk in bin(n)[2:].split("0") if blk])


def phi(m):
    x = Fraction(0)
    d = Fraction(1, 2)
    while m:
        if m & 1:
            x += d
        d /= 2
        m >>= 1
    return x


def one_apex_H(bits):
    m = len(bits)
    beats = [0] * (m + 1)
    for i in range(m):
        for j in range(i + 1, m):
            beats[i] |= 1 << j
        if bits[i]:
            beats[m] |= 1 << i
        else:
            beats[i] |= 1 << m
    return ham_tables(m + 1, beats)[0][-1]


def rooted(n, beats, a):
    Hs, cnt = ham_tables(n, beats)
    rb = [0] * n
    for i in range(n):
        for j in range(n):
            if beats[i] >> j & 1:
                rb[j] |= 1 << i
    Hr, cntr = ham_tables(n, rb)
    full = (1 << n) - 1
    R = 0
    for X in range(1 << n):
        if X >> a & 1:
            continue
        R += Hs[X] * cntr[full ^ X][a]
    return Hs[full], R, cnt[full][a], Hs[full ^ (1 << a)]


def build(word):
    beats = [0]
    root = 0
    n = 1
    for ch in word:
        z = n
        n += 1
        beats.append(0)
        for v in range(z):
            if ch in "BF" and v == root:
                beats[z] |= 1 << v
            else:
                beats[v] |= 1 << z
        if ch == "F":
            root = z
    return n, beats, root


def step(state, ch):
    h, R, e, u = state
    if ch == "S":
        return (h, 2 * R, 0, u)
    if ch == "B":
        return (h + R - e, 2 * R - u, u, u)
    return (h + R - e, h + R, h - e, h)


def main():
    ok = True
    # (1)
    bad = 0
    for m in range(2, 11):
        for bits in itertools.product([0, 1], repeat=m):
            if bits[0] != 1 or bits[-1] != 0:
                continue
            n = int("".join(map(str, bits)), 2)
            if one_apex_H(bits) != 2 * psi(n // 2) + 1:
                bad += 1
    good = (bad == 0) and all(psi(h) == h - q(h) for h in range(5000))
    ok &= good
    print("(1) H = 2 psi(h) + 1 on all 511 reduced words of length <= 10; psi(h) = h - q(h) for h < 5000:", "OK" if good else "FAIL")
    # (2)
    good = all(G(2 * n) == G(n) + n and G(2 * n + 1) == G(n) for n in range(5000))
    good &= all(G(2 * h) + rho(h) == 2 * psi(h) for h in range(1, 5000))
    good &= all(G(n) == sum(1 for mm in range(n) if phi(mm) > phi(n)) for n in range(700))
    ok &= good
    print("(2) G recurrences (n < 5000), G(2h) + rho(h) = 2 psi(h) (h < 5000), G(n) = van der Corput inversion count (n < 700):",
          "OK" if good else "FAIL")
    # (3)
    bad = 0
    for j in range(1, 9):
        for t in range(1 << j):
            kappa = (1 << j) - 1 - sum(1 << (j - p - 2) for p in range(j - 1) if t >> p & 1)
            top = (t >> (j - 1)) & 1
            for m in range(0, 200):
                if psi(m * (1 << j) + t) != psi(m) + psi(t) + kappa * m - top * ((m + 1) // 2):
                    bad += 1
    good = (bad == 0)
    ok &= good
    print("(3) splitting law for j <= 8, all t < 2^j, m < 200:", "OK" if good else "FAIL")
    # (4)
    random.seed(20261110)
    bad = 0
    for trial in range(500):
        L = random.randint(1, 9)
        word = "".join(random.choice("SBF") for _ in range(L))
        st = (1, 1, 1, 1)
        for ch in word:
            st = step(st, ch)
        n, beats, root = build(word)
        if rooted(n, beats, root) != st:
            bad += 1
    good = (bad == 0)
    ok &= good
    print("(4) S/B/F transfer rules against brute force on 500 random words of length <= 9:", "OK" if good else "FAIL")
    X = 6000
    seen = {(1, 1, 1, 1)}
    frontier = [(1, 1, 1, 1)]
    while frontier:
        nxt = []
        for st in frontier:
            for ch in "SBF":
                t = step(st, ch)
                # h never decreases and R >= h along the system; an orbit point with R > 4X cannot lead to h <= X
                # except through states already generated from a smaller R (S only doubles R), so the cap is safe
                if t[0] <= X and t[1] <= 4 * X and t not in seen:
                    seen.add(t)
                    nxt.append(t)
        frontier = nxt
    hit = {s[0] for s in seen}
    miss = [N for N in range(1, X + 1, 2) if N not in hit]
    s197 = [7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583, 637,
            821, 1067, 1339, 1861, 2647]
    good = (miss == s197)
    ok &= good
    print(f"(4) orbit of (1,1,1,1) under S, B, F: {len(seen)} states with h <= {X}, R <= {4*X}; odd h-values missed below {X}: {miss}")
    print("    equal to the 29 values of THM-4884 / Session 197:", "OK" if good else "FAIL")
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
