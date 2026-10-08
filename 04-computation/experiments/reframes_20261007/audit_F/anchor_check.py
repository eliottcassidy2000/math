#!/usr/bin/env python3
"""audit_F (E5): independent anchor census (section 1.5 of the results note), integer arithmetic only.

Anchors: the odd points c (letter a = v_2(3c+1) >= 2, i.e. c = 1 mod 4) of the Syracuse cycles of primitive U-words
(letters >= 1) with sum <= 9 and length <= 4.  Cycle point of a word w: fixed point of F_w = composition of
x -> (3x+1)/2^a; computed here from the affine composition (alpha, beta), not with the author's helper.
For each anchor c and 1 <= D <= 40 the pair (u, v) = (c, (c+1)/3^D - 1) with debt k = D (u = 3^D v + 3^D - 1) is iterated
under the Terras map with a common odd denominator Q (numerators follow x -> x/2, (3x+Q)/2), k += par(u) - par(v),
until absorption (u = v and k = 0) or a repeat of (u, v) (then the future is periodic):
  ABSORB; ANCHOR (repeat with u = v, k != 0: a numeric meeting with unpaid debt); SHIFT (repeat, u != v, zero drift);
  DRIFT (repeat, nonzero drift); UNRESOLVED (no repeat within MAXIT).
Also: the three-step merge for c = 5 mod 8, checked symbolically on residues mod 2^12 and on all anchors.
"""
import sys, itertools
from fractions import Fraction as Fr
from collections import Counter, defaultdict

MAXIT = 20000


def v2(n):
    n = abs(n); v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    return v


def word_point(w):
    alpha, beta = Fr(1), Fr(0)
    for a in w:                       # x -> (3x + 1)/2^a applied in order
        alpha = alpha * 3 / 2 ** a
        beta = (3 * beta + 1) / 2 ** a
    return beta / (1 - alpha)


def necklaces(maxsum, maxlen):
    out = []
    for L in range(1, maxlen + 1):
        for w in itertools.product(range(1, maxsum + 1), repeat=L):
            if sum(w) > maxsum:
                continue
            if any(L % p == 0 and w == w[:p] * (L // p) for p in range(1, L)):
                continue
            rots = [w[i:] + w[:i] for i in range(L)]
            if w != min(rots):
                continue
            out.append(w)
    return out


def syr(x):
    """Syracuse step on a rational with odd denominator, odd numerator"""
    y = 3 * x + 1
    return y / 2 ** v2(y.numerator), v2(y.numerator)


def run_pair(c, D):
    Q = c.denominator * 3 ** D
    v = (c + 1) / Fr(3 ** D) - 1
    U = int(c * Q); V = int(v * Q)
    assert Fr(U, Q) == c and Fr(V, Q) == v
    k = D
    seen = {}
    for s in range(MAXIT):
        if U == V and k == 0:
            return 'ABSORB', s, k
        key = (U, V)
        if key in seen:
            s0, k0 = seen[key]
            drift = k - k0
            if drift != 0:
                return 'DRIFT', s0, k0
            return ('ANCHOR' if U == V else 'SHIFT'), s0, k0
        seen[key] = (s, k)
        pu, pv = U & 1, V & 1
        k += pu - pv
        U = (3 * U + Q) >> 1 if pu else U >> 1
        V = (3 * V + Q) >> 1 if pv else V >> 1
    return 'UNRESOLVED', MAXIT, k


if __name__ == '__main__':
    neck = necklaces(9, 4)
    anchors = []          # (necklace, rotation word, point)
    for w in neck:
        L = len(w)
        for i in range(L):
            r = w[i:] + w[:i]
            c = word_point(r)
            # sanity: the Syracuse orbit of c follows the rotation r and returns to c
            x = c
            for a in r:
                x, aa = syr(x)
                assert aa == a, (r, c)
            assert x == c
            if r[0] >= 2:
                assert (c.numerator * pow(c.denominator, -1, 4)) % 4 == 1
                anchors.append((w, r, c))
    print(f"cycles (primitive necklaces, sum <= 9, length <= 4): {len(neck)}; legal anchors (first letter >= 2): {len(anchors)}; "
          f"first letter 2: {sum(1 for w, r, c in anchors if r[0] == 2)}; >= 3: {sum(1 for w, r, c in anchors if r[0] >= 3)}")
    outc = Counter(); absorbD = defaultdict(list)
    for w, r, c in anchors:
        for D in range(1, 41):
            kind, t, k = run_pair(c, D)
            outc[kind] += 1
            if kind == 'ABSORB':
                absorbD[c].append((D, t))
    print(f"outcomes over (anchor, D <= 40): {dict(outc)}  (total {sum(outc.values())})")
    fl3 = [(r, c) for w, r, c in anchors if r[0] >= 3]
    fl2 = [(r, c) for w, r, c in anchors if r[0] == 2]
    ok3 = [c for r, c in fl3 if (1, 3) in absorbD[c]]
    print(f"first letter >= 3: {len(fl3)} anchors; absorbed at D = 1 at time 3: {len(ok3)}")
    bar2 = [(r, c) for r, c in fl2 if not absorbD[c]]
    print(f"first letter 2: {len(fl2)} anchors; never absorbed for D <= 40: {len(bar2)}; absorbed for some D: {len(fl2) - len(bar2)}")
    print(f"   barriers: {[str(c) for r, c in bar2]}")
    for name in ('7/23', '7/503', '1'):
        c = Fr(name)
        print(f"   anchor {name}: absorbing (D, time) = {absorbD.get(c, [])}")
    # three-step merge, symbolic on residues: for c = 5 mod 8 the pair (c, (c-2)/3) merges at Terras time 3 (not earlier)
    def T(x): return (3 * x + 1) // 2 if x & 1 else x // 2
    MOD = 2 ** 14
    bad = 0; pat = Counter()
    for c0 in range(5, 3 * MOD, 8):
        if (c0 - 2) % 3:
            continue
        u, v = c0, (c0 - 2) // 3
        par = []
        for s in range(3):
            par.append((u & 1, v & 1)); u, v = T(u), T(v)
        pat[tuple(par)] += 1
        if u != v or T(T(c0)) == T(T((c0 - 2) // 3)) or T(c0) == T((c0 - 2) // 3):
            bad += 1
    print(f"three-step merge c = 5 mod 8 (integers c < {3*MOD} with 3 | c-2): failures {bad}; parity pattern (u,v) per step: {dict(pat)}")
