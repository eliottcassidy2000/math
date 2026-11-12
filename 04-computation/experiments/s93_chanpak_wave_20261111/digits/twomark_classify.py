#!/usr/bin/env python3
"""Lane "digits": exact affine classes of the two-mark one-vertex system (17-dimensional integer representation
of twomark.py; readout = first coordinate = H).

(1) Check of the representation on fresh data: random marked tournaments of orders 6 and 7 and random words.
(2) Integer eigenvalues of every word of length <= LEN.
(3) Single words: class spaces C(word; a), census states in them, offsets.
(4) Pairs of words: categories (i) C = 0, (ii) the two maps commute on all of C, (iii) kappa != 0 somewhere on C;
    for (iii): census states in C (actual non-commuting exact pairs).
Census: the states of word(base) for every base of order <= 5 and every word of length <= CEN.

Usage: twomark_classify.py REPFILE LEN CEN
"""
import itertools
import json
import random
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

from affine import class_space, dot, int_eigs, matmul, matvec
from tour import hcount, iso_classes
from twomark import OPS, H_of


def main():
    repfile = sys.argv[1] if len(sys.argv) > 1 else "twomark_rep_5_3.json"
    maxlen = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    cen = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    t0 = time.time()
    rep = json.load(open(repfile))
    d = rep["d"]
    mats = [[[int(Fraction(x)) for x in row] for row in M] for M in rep["mats"]]
    assert all(Fraction(x).denominator == 1 for M in rep["mats"] for row in M for x in row)
    swords = [tuple(w) for w in rep["state_words"]]
    L = [1] + [0] * (d - 1)
    # (1) fresh check
    rng = random.Random(93)
    bad = 0
    ntest = 0
    for t in range(300):
        n = rng.choice((6, 7))
        adj = np.zeros(n, np.int64)
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < 0.5:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
        z1, z2 = rng.sample(range(n), 2)
        s = [H_of(n, adj, z1, z2, w) for w in swords]
        w = [rng.randrange(12) for _ in range(rng.randint(1, 5))]
        x = s
        for o in w:
            x = matvec(mats[o], x)
        ntest += 1
        if x[0] != H_of(n, adj, z1, z2, w):
            bad += 1
    print(f"(1) representation (d = {d}) against direct counts on {ntest} random marked tournaments of orders 6, 7 with random words of length <= 5: {bad} failures")
    # census
    classes = iso_classes(5)
    base_states = set()
    for n in range(2, 6):
        for key, adj in classes[n].items():
            for z1 in range(n):
                for z2 in range(n):
                    if z1 != z2:
                        base_states.add(tuple(H_of(n, adj, z1, z2, w) for w in swords))
    cur = set(base_states)
    census = set(cur)
    for ln in range(cen):
        nxt = set()
        for s in cur:
            for M in mats:
                nxt.add(tuple(matvec(M, list(s))))
        cur = nxt - census
        census |= nxt
    S = np.array(sorted(census), dtype=object)
    Sl = [tuple(int(x) for x in s) for s in S]
    print(f"census: {len(base_states)} base states (orders 2..5), closed under words of length <= {cen}: {len(Sl)} distinct states  [{time.time() - t0:.1f} s]")
    # (2) words
    words = {}
    cur = {(): [[1 if i == j else 0 for j in range(d)] for i in range(d)]}
    ec = Counter()
    for ln in range(1, maxlen + 1):
        nxt = {}
        for w, M in cur.items():
            for o in range(12):
                nxt[w + (o,)] = matmul(mats[o], M)
        cur = nxt
        for w, M in cur.items():
            e = [a for a in int_eigs(M) if a >= 1]
            for a in e:
                ec[(ln, a)] += 1
            if e:
                words[w] = (M, e)
    print(f"(2) words of length <= {maxlen} with an integer eigenvalue >= 1: {len(words)}; (length, eigenvalue): count = {dict(sorted(ec.items()))}  [{time.time() - t0:.1f} s]")

    def wn(w):
        return ".".join(f"({OPS[o][0]:02b},{OPS[o][1]})" for o in w)

    Sarr = np.array(Sl, dtype=np.int64)
    assert np.abs(Sarr).max() < 2 ** 40

    def in_space(V):
        """census states annihilated by every row of V (rows scaled to integers; int64 exact: entries checked)."""
        rows = V.rows
        if not rows:
            return list(Sl)
        mask = np.ones(len(Sl), dtype=bool)
        for r in rows:
            den = 1
            for x in r:
                den = den * x.denominator // int(np.gcd(den, x.denominator))
            ri = [int(x * den) for x in r]
            if max(abs(x) for x in ri) < 2 ** 20:
                mask &= (Sarr @ np.array(ri, dtype=np.int64)) == 0
            else:
                idx = np.nonzero(mask)[0]
                for i in idx:
                    if dot(ri, Sl[i]) != 0:
                        mask[i] = False
            if not mask.any():
                return []
        return [Sl[i] for i in np.nonzero(mask)[0]]

    # (3) singles
    print("(3) SINGLE WORDS with an exact class and moving census states (op = (S, mu): S = set of marks beaten, mu = 0 stay, 1: z1 := v, 2: z2 := v)")
    keep = []
    summ = defaultdict(int)
    for w, (M, e) in sorted(words.items(), key=lambda t: (len(t[0]), t[0])):
        for a in e:
            V, lam = class_space(L, [M], [a])
            dimC = d - V.dim()
            if dimC == 0:
                summ[(len(w), a, "C = 0")] += 1
                continue
            keep.append((w, a))
            inC = in_space(V)
            mov = [s for s in inC if a * s[0] + dot(lam[0], s) != s[0]]
            bs = Counter(dot(lam[0], s) for s in mov)
            kind = "no moving census state" if not mov else ("b = 0" if all(b == 0 for b in bs) else "AFFINE b != 0")
            summ[(len(w), a, kind)] += 1
            if mov and kind == "AFFINE b != 0" and len(w) <= 2:
                bl = sorted(bs.items())
                bshow = bl if len(bl) <= 8 else bl[:4] + ["..."] + bl[-3:]
                print(f"   {wn(w)}: a = {a}, dim C = {dimC}, moving census states {len(mov)}, b: {bshow}")
    print("   summary (length, multiplier, kind): count")
    for k2, v in sorted(summ.items()):
        print(f"      {k2}: {v}")
    print(f"   [{time.time() - t0:.1f} s]", flush=True)
    # (4) pairs
    cat = Counter()
    free = []
    for x in range(len(keep)):
        w1, a1 = keep[x]
        for y in range(x + 1, len(keep)):
            w2, a2 = keep[y]
            if w1 == w2:
                continue
            V, lam = class_space(L, [words[w1][0], words[w2][0]], [a1, a2])
            dimC = d - V.dim()
            if dimC == 0:
                cat["(i) C = 0"] += 1
                continue
            kappa = [(a2 - 1) * p - (a1 - 1) * q for p, q in zip(lam[0], lam[1])]
            if V.contains(kappa):
                cat["(ii) commute on all of C"] += 1
                continue
            cat["(iii) kappa != 0 somewhere on C"] += 1
            inC = in_space(V)
            good = [(s[0], dot(lam[0], s), dot(lam[1], s)) for s in inC if dot(kappa, s) != 0]
            free.append((w1, a1, w2, a2, dimC, len(inC), good))
    print(f"(4) PAIRS OF WORDS (length <= {maxlen}) with nonzero single class spaces: {dict(cat)}")
    nst = sum(1 for f in free if f[6])
    print(f"    category (iii): {len(free)} pairs; with a census state in C and kappa != 0 (an actual non-commuting exact pair): {nst}")
    byt = Counter((a1, a2, dimC) for (w1, a1, w2, a2, dimC, ninC, good) in free)
    print("    category (iii) by (a1, a2, dim C):", dict(sorted(byt.items())))
    for (w1, a1, w2, a2, dimC, ninC, good) in [f for f in free if f[6]][:40]:
        print(f"    ACTUAL: {wn(w1)} (a={a1}) & {wn(w2)} (a={a2}): dim C = {dimC}, census states in C {ninC}; (H, b1, b2): {good[:6]}")
    for (w1, a1, w2, a2, dimC, ninC, good) in free[:25]:
        print(f"    virtual: {wn(w1)} (a={a1}) & {wn(w2)} (a={a2}): dim C = {dimC}, census states in C: {ninC}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
