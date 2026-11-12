#!/usr/bin/env python3
"""Lane "digits": the two-apex family G(2,0) - minimal representation and exact affine classes.

G(2,0): transitive spine s_1 -> ... -> s_n, apexes u -> v, one letter (a_i, b_i) per spine vertex (a_i = 1 iff
u -> s_i, b_i = 1 iff v -> s_i).  Appending a letter = adding a new last spine vertex, i.e. a one-vertex extension
"v is beaten by everything except possibly the two marks u, v".  H by the 36-state transfer of the seed probe
(re-derived here and checked against a direct count for n <= 6).

Output: dimension of the minimal representation, integer letter matrices, integer eigenvalues of all words of
length <= LEN, class spaces of single words and of pairs of words (categories (i), (ii), (iii) as in
onemark_classify.py), census = states of all words of length <= CENSUS.

Usage: twoapex_classify.py LEN CENSUS
"""
import itertools
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

from affine import class_space, dot, int_eigs_q, matmul, matvec
from minrep import minimal_representation
from tour import hcount


def step_matrix(a, b):
    """36 x 36 matrix of one letter (state index arr*18 + r1*6 + r2*2 + r3, as in seed/two_apex_probe.py)."""
    M = [[0] * 36 for _ in range(36)]
    for arr in range(2):
        bx = a if arr == 0 else b
        by = b if arr == 0 else a
        for r1 in range(3):
            for r2 in range(3):
                for r3 in range(2):
                    src = arr * 18 + r1 * 6 + r2 * 2 + r3
                    n1 = 1 if bx == 0 else 2
                    M[arr * 18 + n1 * 6 + r2 * 2 + r3][src] += 1
                    if r2 != 0 or bx == 1:
                        n2 = 1 if by == 0 else 2
                        M[arr * 18 + r1 * 6 + n2 * 2 + r3][src] += 1
                    if r3 != 0 or by == 1:
                        M[arr * 18 + r1 * 6 + r2 * 2 + 1][src] += 1
    return M


def final_row():
    phi = [0] * 36
    for arr in range(2):
        for r1 in range(3):
            if r1 == 2:
                continue
            for r2 in range(3):
                if r2 == 2:
                    continue
                if r2 == 0 and arr == 1:
                    continue
                for r3 in range(2):
                    phi[arr * 18 + r1 * 6 + r2 * 2 + r3] += 1
    return phi


def direct(word):
    """direct count: vertices spine 0..n-1, u = n, v = n + 1."""
    n = len(word)
    V = n + 2
    u, v = n, n + 1
    adj = np.zeros(V, np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
        a, b = word[i]
        if a:
            adj[u] |= 1 << i
        else:
            adj[i] |= 1 << u
        if b:
            adj[v] |= 1 << i
        else:
            adj[i] |= 1 << v
    adj[u] |= 1 << v
    return int(hcount(V, adj))


def main():
    maxlen = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    census_len = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    t0 = time.time()
    letters = [(0, 0), (1, 0), (0, 1), (1, 1)]
    T = [step_matrix(a, b) for (a, b) in letters]
    v0 = [0] * 36
    v0[0] = 1
    v0[18] = 1
    phi = final_row()
    # check the transfer against the direct count
    ok = True
    for n in range(0, 7):
        for w in itertools.product(range(4), repeat=n):
            vec = v0
            for l in w:
                vec = matvec(T[l], vec)
            if dot(phi, vec) != direct([letters[l] for l in w]):
                ok = False
    print("36-state transfer = direct count for all words of length <= 6:", ok)
    rep = minimal_representation(T, v0, phi)
    d = rep["d"]
    print(f"reachable dimension {rep['dR']}, minimal dimension d = {d}")
    C = rep["C"]
    dens = set(x.denominator for M in C for row in M for x in row) | set(Fraction(x).denominator for x in rep["x0"])
    print("denominators in the minimal representation:", sorted(dens))
    Ci = [[[Fraction(x) for x in row] for row in M] for M in C]      # rational matrices (NOT truncated)
    x0 = [Fraction(x) for x in rep["x0"]]
    nonint = sum(1 for M in Ci for row in M for x in row if x.denominator != 1)
    print(f"non-integer entries in the letter matrices: {nonint} (the coordinates are integer on actual states)")
    print("coordinates: x_k = H(word . u_k) for the suffix words u_k =", ["".join("0123"[l] for l in w) or "empty" for w in rep["obs_words"]])
    print("x0 =", x0)
    for l, M in zip(letters, Ci):
        print(f"letter {l}: integer eigenvalues {int_eigs_q(M, -8, 64)}")
    L = [Fraction(1)] + [Fraction(0)] * (d - 1)
    # census of actual states
    states = [tuple(x0)]
    cur = [tuple(x0)]
    for n in range(census_len):
        nxt = []
        for s in cur:
            for M in Ci:
                nxt.append(tuple(matvec(M, list(s))))
        cur = nxt
        states.extend(cur)
    S = sorted(set(states))
    assert all(x.denominator == 1 for s in S for x in s)
    print(f"census: {len(states)} words of length <= {census_len}, {len(S)} distinct states")
    # words
    words = {}
    cur = {(): [[Fraction(1 if i == j else 0) for j in range(d)] for i in range(d)]}
    for ln in range(1, maxlen + 1):
        nxt = {}
        for w, M in cur.items():
            for li in range(4):
                nxt[w + (li,)] = matmul(Ci[li], M)
        cur = nxt
        for w, M in cur.items():
            e = [a for a in int_eigs_q(M) if a >= 1]
            if e:
                words[w] = (M, e)
    print(f"words of length <= {maxlen} with an integer eigenvalue >= 1: {len(words)}")
    ec = Counter()
    for w, (M, e) in words.items():
        for a in e:
            ec[(len(w), a)] += 1
    print("  (length, eigenvalue): count", dict(sorted(ec.items())))

    def wn(w):
        return "".join("0123"[l] for l in w)

    print()
    print("SINGLE WORDS: exact classes h -> a h + b with moving census states (letters: 0 = (0,0) sink, 1 = (1,0), 2 = (0,1), 3 = (1,1))")
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
            inC = [s for s in S if all(dot(r, s) == 0 for r in V.rows)]
            mov = [s for s in inC if a * s[0] + dot(lam[0], s) != s[0]]
            bs = Counter(dot(lam[0], s) for s in mov)
            kind = "no moving census state" if not mov else ("b = 0" if all(b == 0 for b in bs) else "AFFINE b != 0")
            summ[(len(w), a, kind)] += 1
            if mov and len(w) <= 2:
                bl = sorted(bs.items())
                bshow = bl if len(bl) <= 10 else bl[:5] + ["..."] + bl[-4:]
                print(f"  word {wn(w)}: a = {a}, dim C = {dimC}, moving census states {len(mov)}, b: {bshow} [{kind}]")
    print("  summary (length, multiplier, kind): count")
    for k2, v in sorted(summ.items()):
        print(f"     {k2}: {v}")
    print()
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
            inC = [s for s in S if all(dot(r, s) == 0 for r in V.rows)]
            good = [(s, dot(lam[0], s), dot(lam[1], s)) for s in inC if dot(kappa, s) != 0]
            free.append((w1, a1, w2, a2, dimC, len(inC), good))
    print(f"PAIRS OF WORDS (length <= {maxlen}): {dict(cat)}")
    nst = sum(1 for f in free if f[6])
    print(f"category (iii): {len(free)} pairs; with a census state in C and kappa != 0 (a non-commuting exact pair): {nst}")
    for (w1, a1, w2, a2, dimC, ninC, good) in free[:60]:
        print(f"   {wn(w1)} (a={a1}) & {wn(w2)} (a={a2}): dim C = {dimC}, census states in C: {ninC}, non-commuting: {[(s[0], b1, b2) for (s, b1, b2) in good[:5]]}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
