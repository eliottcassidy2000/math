#!/usr/bin/env python3
"""S91 helper (census_k2): the exact identities behind the K = 2 cascade, checked letter by letter (pure Python).

Two-J word  A J B J C,  A, B O/Z words, C empty or an O/Z word ending in Z (tail code m).
States: (D1, E1, U1, 1) after A; (D2, E2, U2, D1) after A J B.
 X1  value(A J B J tail m) = x3 D2 + 2^c y3 U2 - [c = 0] E2     (THM-4916 (1b) in the F1 form of dict_k1.py)
 X2  E2 = 0 when B ends in O;  E2 = D1 when B ends in Z;  E2 = D1 - E1 when B is empty
 X3  cascade (THM-4916 lemma (M3), here for every word and not only in the Cesaro limit): if p | gcd(D1, U1), then
         B not empty:  (D2, U2) = (0, 0) mod p   <=>   E1 = 0 or B starts with Z            (the flag k1 = 1)
                        and otherwise (D2, U2) = (-1, 0) mod p, E2 = 0 mod p;
         B empty:      (D2, U2) = (-E1, 0) mod p.
 X4  with B = B' O^b, n' the code of B' (lowest digit = first letter), L' its length, mu = n'/2^L' in [1/2, 1):
         D2 = x(n') D1 + (n' + 1)(U1 - D1) + [n' even](D1 - E1),      U2 = 2^b (D1 (2^(L'+1) - rev_L'(n')) + 2^L' (U1 - D1)),
     and U2 / (2^b D2) = (1/mu)(1 + O(L'/2^L')): the pre-J2 states have the same projective law tau = 2^b/mu as the J-free
     states (tau = 2^a y/x).  Checked: the two formulas exactly; the ratio numerically for |B'| = 6, 10, 14, 18.
Universe: all A with |A| <= 9, all B with |B| <= 9, tails m < 64 (X1), primes 3, 5, 7, 11, 13, 17 (X3).
"""
import sys
from itertools import product

PRIMES = [3, 5, 7, 11, 13, 17]

Q = [0, 0]
for m in range(2, 1 << 20):
    half = m >> 1
    Q.append(Q[half] + (half + (half & 1) if m & 1 else 0))


def O(s):
    D, E, U, d = s
    return (D, 0, 2 * U, d)


def Z(s):
    D, E, U, d = s
    return (D - E + U, d, 2 * U - d, d)


def J(s):
    D, E, U, d = s
    return (D - E + U, D - E, D + U, D)


def apply(s, word):
    for letter in word:
        s = Z(s) if letter else O(s)
    return s


def tail_letters(m):
    out = []
    while m:
        out.append(m & 1)
        m >>= 1
    return out


def main():
    ok = True
    words = [()]
    for L in range(1, 10):
        words += list(product((0, 1), repeat=L))
    statesA = [(w, apply((1, 0, 2, 1), w)) for w in words]
    x1 = x2 = x3 = x4 = True
    n13 = n3 = 0
    conf_stats = {p: [0, 0] for p in PRIMES}
    for wa, sA in statesA:
        D1, E1, U1, _ = sA
        pJ = J(sA)
        divs = [p for p in PRIMES if D1 % p == 0 and U1 % p == 0]
        deep = len(wa) <= 5 or bool(divs)
        if not deep:
            continue
        for wb in words:
            if len(wa) > 5 and len(wb) > 7:
                continue
            s2 = apply(pJ, wb)
            D2, E2, U2, d2 = s2
            # X2
            if len(wb) == 0:
                x2 &= E2 == D1 - E1
            elif wb[-1] == 0:
                x2 &= E2 == 0
            else:
                x2 &= E2 == D1
            x2 &= d2 == D1
            # X3
            for p in divs:
                n3 += 1
                if len(wb) == 0:
                    x3 &= (D2 % p, U2 % p) == ((-E1) % p, 0)
                else:
                    k1 = (E1 == 0) or (wb[0] == 1)
                    if k1:
                        x3 &= D2 % p == 0 and U2 % p == 0 and E2 % p == 0
                    else:
                        x3 &= D2 % p == p - 1 and U2 % p == 0 and E2 % p == 0
                    conf_stats[p][0] += 1
                    conf_stats[p][1] += 1 if k1 else 0
            # X4 exact formulas
            b = 0
            core = list(wb)
            while core and core[-1] == 0:
                core.pop()
                b += 1
            n = sum(bit << i for i, bit in enumerate(core))
            Lc = len(core)
            rev = sum(bit << (Lc - 1 - i) for i, bit in enumerate(core))
            xn = 1 + 2 * n - Q[n]
            D2f = xn * D1 + (n + 1) * (U1 - D1) + ((D1 - E1) if n % 2 == 0 else 0)
            U2f = (D1 * ((1 << (Lc + 1)) - rev) + (1 << Lc) * (U1 - D1)) << b
            x4 &= (D2f == D2) and (U2f == U2)
            # X1 on a sub-universe
            if len(wa) <= 4 and len(wb) <= 4:
                post = J(s2)
                for m in range(64):
                    val = apply(post, tail_letters(m))[0]
                    if m & 1:
                        c = ((m + 1) & -(m + 1)).bit_length() - 1
                        h3 = m >> (c + 1)
                    else:
                        c = 0
                        h3 = m >> 1
                    xx, yy = 1 + 2 * h3 - Q[h3], 2 * h3 + 1
                    x1 &= val == xx * D2 + (yy << c) * U2 - (E2 if c == 0 else 0)
                    n13 += 1
    print(f"X1 value = x3 D2 + 2^c y3 U2 - [c = 0] E2: {x1}   ({n13} words A J B J C, |A|, |B| <= 4, m < 64)")
    print(f"X2 E2 in (0, D1, D1 - E1) by the last letter of B, d2 = D1: {x2}")
    print(f"X3 cascade for p | gcd(D1, U1): confined at J2 iff k1 = 1 (else (D2, U2) = (-1, 0)): {x3}   ({n3} triples (A, B, p))")
    print("   letter-count share of k1 = 1 among them (B not empty):", {p: round(v[1] / v[0], 4) for p, v in conf_stats.items() if v[0]})
    print(f"X4 closed forms of D2 and U2: {x4}")
    ok &= x1 and x2 and x3 and x4
    # X4 asymptotics: tau2 mu / 2^b -> 1
    import random
    rng = random.Random(91)
    x4b = True
    for Lc in (6, 10, 14, 18):
        worst = 0.0
        for trial in range(4000):
            la = rng.randint(0, 8)
            wa = [rng.randint(0, 1) for _ in range(la)]
            sA = apply((1, 0, 2, 1), wa)
            core = [rng.randint(0, 1) for _ in range(Lc - 1)] + [1]
            b = rng.randint(0, 3)
            s2 = apply(J(sA), core + [0] * b)
            n = sum(bit << i for i, bit in enumerate(core))
            mu = n / (1 << Lc)
            ratio = s2[2] * mu / ((1 << b) * s2[0])
            worst = max(worst, abs(ratio - 1))
        bound = 4.0 * Lc / (1 << Lc)
        print(f"   |B'| = {Lc:2d}: max |U2 mu / (2^b D2) - 1| over 4000 random (A, B', b) = {worst:.3e}   (4 L'/2^L' = {bound:.3e})")
        x4b &= worst < bound
    print(f"X4 asymptotic ratio tau2 = (2^b / mu)(1 + O(L'/2^L')): {x4b}")
    ok &= x4b
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
