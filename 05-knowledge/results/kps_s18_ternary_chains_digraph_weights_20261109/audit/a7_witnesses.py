#!/usr/bin/env python3
"""AUDIT: recheck the explicit one-defect witnesses of out/one_defect_witness.out and the order-4 H = 7 matrix of the
README: pair types (T / B / Z), loops, and H by brute-force permutation count."""
import itertools
import re
import sys


def H(A):
    n = len(A)
    return sum(1 for p in itertools.permutations(range(n)) if all(A[p[i]][p[i + 1]] for i in range(n - 1)))


def types(A):
    n = len(A)
    t = {"T": 0, "B": 0, "Z": 0}
    for i in range(n):
        for j in range(i + 1, n):
            s = A[i][j] + A[j][i]
            t["Z" if s == 0 else ("T" if s == 1 else "B")] += 1
    return t, sum(A[i][i] for i in range(n))


for line in open(sys.argv[1]):
    m = re.search(r"rows ([01 ]+)$", line.strip())
    if not m:
        continue
    rows = m.group(1).split()
    A = [[int(c) for c in r] for r in rows]
    t, loops = types(A)
    print(f"{line.strip()[:60]}... -> H = {H(A)}, pair types {t}, loops {loops}")
A = [[int(c) for c in r] for r in "00001 10000 01000 11100 01110".split()]
print("README order-5 H=7 one-missing-pair example:", H(A), types(A))
