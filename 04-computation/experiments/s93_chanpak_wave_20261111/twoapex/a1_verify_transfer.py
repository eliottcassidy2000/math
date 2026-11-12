#!/usr/bin/env python3
"""Lane twoapex, deliverable (a), part 1: the 36-state transfer count equals the number of Hamiltonian paths.
Three independent counts: subset DP on the tournament (all words, n <= 6), literal enumeration of vertex orders
(all words, n <= 4), literal enumeration of the 3-colourings of the run decomposition (all words, n <= 5), and
random words of length 7..10 against the subset DP.  Also the elementary identities used later:
H(0 w) = H(w) (a source), H(w 3) = H(w) (a sink), and the converse symmetry.
"""
import random
import sys

from ta_core import H36, H_direct, brute, brute_perm, letters_ab, words


def main():
    ok = True
    total = 0
    for n in range(0, 7):
        bad = 0
        for w in words(n):
            a, b = letters_ab(w)
            h = H36(w)
            if h != brute(n, a, b):
                bad += 1
            total += 1
        print(f"n = {n}: all {4 ** n} words, transfer != subset DP on {bad} words")
        ok &= bad == 0
    for n in range(0, 5):
        bad = sum(1 for w in words(n) if H36(w) != brute_perm(n, *letters_ab(w)))
        print(f"n = {n}: all {4 ** n} words, transfer != literal vertex orders on {bad} words")
        ok &= bad == 0
    for n in range(0, 6):
        bad = sum(1 for w in words(n) if H36(w) != H_direct(w))
        print(f"n = {n}: all {4 ** n} words, transfer != literal 3-colourings on {bad} words")
        ok &= bad == 0
    rng = random.Random(20261111)
    bad = 0
    cnt = 0
    for n in (7, 8, 9, 10):
        for _ in range(60 if n < 10 else 25):
            w = [rng.randrange(4) for _ in range(n)]
            a, b = letters_ab(w)
            cnt += 1
            if H36(w) != brute(n, a, b):
                bad += 1
    print(f"random words of length 7..10: {cnt} tested, mismatches {bad}")
    ok &= bad == 0
    # identities
    bad = 0
    for n in range(0, 7):
        for w in words(n):
            h = H36(w)
            if H36((0,) + w) != h or H36(w + (3,)) != h:
                bad += 1
            # converse: reverse all arcs and the spine order; the apexes swap roles (v -> u becomes the marked arc
            # u' = v, v' = u): a'_i = 1 - b_(n+1-i), b'_i = 1 - a_(n+1-i)
            a, b = letters_ab(w)
            w2 = tuple((1 - b[n - 1 - i]) + 2 * (1 - a[n - 1 - i]) for i in range(n))
            if H36(w2) != h:
                bad += 1
    print(f"H(0 w) = H(w), H(w 3) = H(w), converse symmetry w -> w*: violations {bad} (all words, n <= 6)")
    ok &= bad == 0
    print("ALL CHECKS PASS" if ok else "FAILURE")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
