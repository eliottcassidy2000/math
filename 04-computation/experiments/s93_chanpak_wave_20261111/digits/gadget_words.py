#!/usr/bin/env python3
"""Lane "digits", part (2): words in the gadget system of THM-5130 (A) with the sink operation.

Alphabet: sigma = "add a sink to the base, keep the selected set"  (r, M, h, E) -> (r, 2M, h, 0)
          and the gadget matrices G(X, P), X of order <= NMAX.
Readouts: L_h = (0, 0, 1, 0)  (H of the base)  and  L_1 = (1, 2, 0, 0)  (H of the closed tournament B => pt => v).
For every word of length <= LEN: integer eigenvalues >= 1 (tested exactly for 1..64 and through float roots).
For every word with such an eigenvalue a: the class space C(word; a) for each readout, the census states in it
and their offsets.  For every pair of such words: categories (i) C = 0, (ii) commuting on all of C, (iii) not.

Usage: gadget_words.py NMAX LEN
"""
import itertools
import json
import sys
import time
from collections import Counter, defaultdict

from affine import class_space, dot, int_eigs, matmul
from gadgets import pair_state
from tour import iso_classes

SIGMA = [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]]


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    maxlen = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    t0 = time.time()
    gad = json.load(open(f"gadget_matrices_{nmax}.json"))
    letters = [("s", SIGMA)] + [(f"g{i}", g["G"]) for i, g in enumerate(gad)]
    classes = iso_classes(6)
    states = {(1, 0, 1, 0): (0, [], 0)}
    for n in range(1, 7):
        for key, adj in classes[n].items():
            for O in range(1, 1 << n):
                s = pair_state(n, adj, O)
                states.setdefault(s, (n, [int(x) for x in adj], O))
    S = sorted(states)
    print(f"alphabet: sigma + {len(gad)} gadget matrices (X of order <= {nmax}); words of length <= {maxlen}; {len(S)} census states")
    words = {}
    cur = {(): [[1 if i == j else 0 for j in range(4)] for i in range(4)]}
    eigcount = Counter()
    nwords = 0
    for ln in range(1, maxlen + 1):
        nxt = {}
        for w, M in cur.items():
            for name, G in letters:
                if ln > 1 and w[-1] == "s" and name == "s" and False:
                    continue
                M2 = matmul(G, M)          # apply w first, then the new letter
                w2 = w + (name,)
                nxt[w2] = M2
        cur = nxt
        for w, M in cur.items():
            nwords += 1
            if len(w) > 1 and all(x == "s" for x in w):
                continue          # sigma^k, k >= 2: powers of the letter sigma, which is kept
            e = [a for a in int_eigs(M) if a >= 1]
            if e:
                words[w] = (M, e)
                for a in e:
                    eigcount[(ln, a)] += 1
        print(f"length {ln}: {len(cur)} words; words (sigma^k, k >= 2, omitted) with an integer eigenvalue >= 1 so far: {len(words)}  [{time.time() - t0:.1f} s]", flush=True)
    print("integer eigenvalues >= 1 by (word length, eigenvalue): count")
    for k2, v in sorted(eigcount.items()):
        print(f"   {k2}: {v}")

    def wname(w):
        out = []
        for x in w:
            if x == "s":
                out.append("s")
            else:
                g = gad[int(x[1:])]
                out.append(f"G(n={g['n']},X={g['adj']},P={g['P']:0{g['n']}b})")
        return " . ".join(out)

    for Lname, L in (("L_h", [0, 0, 1, 0]), ("L_1", [1, 2, 0, 0])):
        print()
        print(f"=== readout {Lname} = {L}")
        summ = defaultdict(int)
        shown = 0
        keep = []
        for w, (M, e) in words.items():
            for a in e:
                V, lam = class_space(L, [M], [a])
                dimC = 4 - V.dim()
                if dimC == 0:
                    summ[(len(w), a, "C = 0")] += 1
                    continue
                inC = [s for s in S if all(dot(r, s) == 0 for r in V.rows)]
                mov = [s for s in inC if a * dot(L, s) + dot(lam[0], s) != dot(L, s)]
                bs = Counter(dot(lam[0], s) for s in mov)
                kind = "no moving census state" if not mov else ("b = 0" if all(b == 0 for b in bs) else "AFFINE b != 0")
                summ[(len(w), a, kind)] += 1
                keep.append((w, a))
                if mov and kind != "b = 0" and shown < 60:
                    shown += 1
                    bl = sorted(bs.items())
                    bshow = bl if len(bl) <= 8 else bl[:4] + ["..."] + bl[-3:]
                    print(f"  word [{wname(w)}]: a = {a}, dim C = {dimC}, moving census states {len(mov)}, b: {bshow}, e.g. {mov[:2]}")
        print(f"  summary (length, multiplier, kind): count")
        for k2, v in sorted(summ.items()):
            print(f"     {k2}: {v}")
        cat = Counter()
        free = []
        for x in range(len(keep)):
            w1, a1 = keep[x]
            for y in range(x + 1, len(keep)):
                w2, a2 = keep[y]
                if w1 == w2:
                    continue
                V, lam = class_space(L, [words[w1][0], words[w2][0]], [a1, a2])
                dimC = 4 - V.dim()
                if dimC == 0:
                    cat["(i) C = 0"] += 1
                    continue
                kappa = [(a2 - 1) * p - (a1 - 1) * q for p, q in zip(lam[0], lam[1])]
                if V.contains(kappa):
                    cat["(ii) commute on all of C"] += 1
                    continue
                cat["(iii) kappa != 0 somewhere on C"] += 1
                inC = [s for s in S if all(dot(r, s) == 0 for r in V.rows)]
                good = [(s, dot(lam[0], s), dot(lam[1], s)) for s in inC if dot(kappa, s) != 0]
                free.append((w1, a1, w2, a2, dimC, len(inC), good, V.annihilator()))
        print(f"  pairs of words with class spaces: {dict(cat)}")
        for (w1, a1, w2, a2, dimC, ninC, good, basis) in free[:40]:
            print(f"   (iii): [{wname(w1)}] a={a1} & [{wname(w2)}] a={a2}: dim C = {dimC}, basis {basis}, census in C: {ninC}, non-commuting states: {good[:3]}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
