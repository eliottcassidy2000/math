#!/usr/bin/env python3
"""Word reversal on the integer cycles of x -> (3x + k)/2^v (k odd, x odd), k <= 41 (S20 addendum).
For a cycle with word w (length d, cost A) the cycle point is y_0 = k C_w/(2^A - 3^d); the reversed word gives
y_0' = k C_rev(w)/(2^A - 3^d), an integer cycle of the same map iff (2^A - 3^d) | k C_rev(w).
Cycles found from all odd starts |x| <= 5 10^4 (orbits capped at 10^15 and 1000 steps; the start -k/3 is skipped).
Run: python 04-computation/experiments/collatz_five_mirrors_reversal_20260929.py
"""
from __future__ import annotations

from fractions import Fraction


def carry(word):
    d = len(word)
    s, pref = 0, 0
    for j in range(1, d + 1):
        s += 3 ** (d - j) * 2 ** pref
        pref += word[j - 1]
    return s


def step(x, k):
    t = 3 * x + k
    v = (t & -t).bit_length() - 1
    return t >> v, v


if __name__ == "__main__":
    N = 5 * 10 ** 4
    for k in range(1, 42, 2):
        cycles = {}
        for x0 in range(-N, N + 1, 2):
            if 3 * x0 + k == 0:
                continue
            x = x0
            seen = {}
            for i in range(1000):
                if 3 * x + k == 0:
                    break
                if x in seen:
                    # cycle found: extract it
                    cyc = []
                    y = x
                    while True:
                        cyc.append(y)
                        y, _ = step(y, k)
                        if y == x:
                            break
                    key = min(cyc)
                    if key not in cycles:
                        cycles[key] = cyc
                    break
                if abs(x) > 10 ** 15:
                    break
                seen[x] = i
                x, _ = step(x, k)
        rows = []
        n_sym = 0
        for key, cyc in sorted(cycles.items()):
            # word from the minimum element
            i0 = cyc.index(key)
            cyc = cyc[i0:] + cyc[:i0]
            word = []
            for y in cyc:
                _, v = step(y, k)
                word.append(v)
            d, A = len(word), sum(word)
            den = 2 ** A - 3 ** d
            y0 = Fraction(k * carry(word), den)
            assert y0 == key, (k, word, y0, key)
            y0r = Fraction(k * carry(word[::-1]), den)
            rot = any(word[i:] + word[:i] == word[::-1] for i in range(d))
            integral = y0r.denominator == 1
            verified = False
            partner = None
            if integral:
                # the affine fixed point is an integer; verify that its actual orbit has the reversed word exactly
                x = int(y0r)
                ok = x % 2 == 1
                cyc2 = [x]
                for b in word[::-1]:
                    t = 3 * x + k
                    v = (t & -t).bit_length() - 1 if t else -1
                    ok &= v == b
                    x = t >> b if b <= v and t else 0
                    cyc2.append(x)
                verified = ok and cyc2[-1] == int(y0r)
                if verified:
                    partner = min(cyc2[:-1])
            n_sym += verified
            tag = ('R' if rot else '') + ('I' if integral else '') + ('V' if verified else '')
            rows.append(f"{key}(d={d},A={A}){tag}->{y0r if not integral else int(y0r)}" + (f"[cycle {partner}]" if verified and not rot else ""))
        pairs = sorted(set(tuple(sorted((int(r.split('(')[0]), int(r.split('[cycle ')[1].rstrip(']'))))) for r in rows if '[cycle ' in r))
        print(f"k={k:2d}: {len(cycles)} cycles; reversal gives an integer cycle for {n_sym} (verified orbits); non-rotation dual pairs {pairs}: " + "; ".join(rows))
    print("   (R = reversed word is a rotation of the word; I = affine fixed point of the reversed word is an integer; V = its orbit is a cycle with exactly the reversed word)")
    print("DONE")
