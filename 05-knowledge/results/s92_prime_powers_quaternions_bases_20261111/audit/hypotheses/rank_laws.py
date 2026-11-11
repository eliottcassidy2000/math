#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9666 (A), (B). Exact class laws of the hyperbolic form sum_{i<=r} a_i b_i for uniform
# variables modulo q = p^k (convolution of the multiplication table law), against the stated factors
#   (1 - p^-r) sigma_(r-1)(p^v) / p^((r-1) v)   (v < k),
# and against the sum of 2r squares (r = 3: six squares; the two laws agree iff p = 1 mod 4).
# Also: class-0 factor modulo p for r = 4: 1 + p^-3 - p^-4 (bias 2/81 at 3, 4/625 at 5), and the share 17/32 of the
# class 0 modulo 2 for r = 4.
from fractions import Fraction
import numpy as np


def conv(a, b, q):
    out = np.zeros(q, dtype=object)
    for i in range(q):
        if a[i]:
            for j in range(q):
                if b[j]:
                    out[(i + j) % q] += a[i] * b[j]
    return out


def table_law(q):
    t = np.zeros(q, dtype=object)
    for a in range(q):
        for b in range(q):
            t[(a * b) % q] += 1
    return t


def square_law(q):
    t = np.zeros(q, dtype=object)
    for a in range(q):
        t[(a * a) % q] += 1
    return t


def power(law, r, q):
    out = np.zeros(q, dtype=object)
    out[0] = 1
    for _ in range(r):
        out = conv(out, law, q)
    return out


def vp(c, p, k):
    if c == 0:
        return k
    v = 0
    while c % p == 0:
        c //= p
        v += 1
    return v


def main():
    print("# rank_laws.py (auditor hypotheses; own code)")
    for p, k in ((3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1), (7, 2), (13, 1), (2, 1), (2, 3)):
        q = p ** k
        tab = table_law(q)
        sq = square_law(q)
        for r in (2, 3, 4):
            law = power(tab, r, q)
            tot = Fraction(q) ** (2 * r)
            ok = True
            facs = {}
            for c in range(q):
                v = vp(c, p, k)
                f = Fraction(int(law[c]) * q) / tot
                facs.setdefault(v, set()).add(f)
                if v < k:
                    pred = (1 - Fraction(1, p ** r)) * sum(Fraction(1, p ** ((r - 1) * i)) for i in range(v + 1))
                    ok &= f == pred
            sqlaw = power(sq, 2 * r, q)
            same = all(int(law[c]) == int(sqlaw[c]) for c in range(q))
            exp = (r % 2 == 0) or (p % 4 == 1)
            note = ""
            if p == 2:
                note = " (p = 2: not covered by the statement about squares)"
            print("q = %3d r = %d: factor formula for v < k: %s; factors by valuation: %s; equal to the sum of %d squares: %s (expected %s)%s" % (
                q, r, ok, {v: sorted(str(x) for x in s) for v, s in sorted(facs.items())}, 2 * r, same, exp if p != 2 else "-", note))
    # the numbers of HYP-9666 (A)
    for p in (3, 5):
        f0 = 1 + Fraction(1, p ** 3) - Fraction(1, p ** 4)
        print("rank 4, class 0 modulo %d: factor %s = 1 + %s" % (p, f0, f0 - 1))
    tab2 = table_law(2)
    l4 = power(tab2, 4, 2)
    print("rank 4 modulo 2: share of the class 0 = %s (a Hamiltonian path count is odd: share 0)" % Fraction(int(l4[0]), 2 ** 8))
    # THM-4916 two root changes, class 0: 1 + (3p-2)(p-1)/(2p^3)
    for p in (3, 5, 7):
        print("THM-4916 K = 2, class 0 modulo %d: %s; uniform rank-4 value %s" % (
            p, 1 + Fraction((3 * p - 2) * (p - 1), 2 * p ** 3), 1 + Fraction(1, p ** 3) - Fraction(1, p ** 4)))


if __name__ == "__main__":
    main()
