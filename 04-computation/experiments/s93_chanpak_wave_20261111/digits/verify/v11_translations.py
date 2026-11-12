#!/usr/bin/env python3
"""Verifier (lane digits), check V11: exact translations H -> H + b (multiplier 1).

Commuting systems give no covering ONLY when some multiplier differs from 1; a single exact translation with b != 0
would give an arithmetic progression of path counts.  The lane excludes translations on census states.  Here,
census-free: for every operation / word M with the eigenvalue 1, is the offset lambda = L M - L identically 0 on the
class space C(M; 1)?  (If yes: H is FIXED on every state of the class, actual or not.)
Systems: gadget words (sigma + gadgets, both readouts), two-mark words of length <= 2 (contains the two-apex family).
"""
import json, os, sys, time
from collections import Counter
from fractions import Fraction
from vlin import charpoly_int, class_space, matmul
from vroots import int_roots

here = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
SIGMA = [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]]
for nmax, maxlen in ((4, 2), (3, 3)):
    gad = json.load(open(f"{here}/../gadget_matrices_{nmax}.json"))
    letters = [SIGMA] + [g["G"] for g in gad]
    cur = [[[1 if i == j else 0 for j in range(4)] for i in range(4)]]
    words = []
    for ln in range(maxlen):
        cur = [matmul(G, M) for M in cur for G in letters]
        words += cur
    for Lname, L in (("L_h", [0, 0, 1, 0]), ("L_1", [1, 2, 0, 0])):
        res = Counter()
        for M in words:
            if 1 in int_roots(charpoly_int(M)):
                V, lam = class_space(L, [M], [1])
                if V.dim() == 4:
                    res["C = 0"] += 1
                elif V.contains(lam[0]):
                    res["offset = 0 on all of C (H fixed)"] += 1
                else:
                    res["offset != 0 somewhere on C (virtual translation)"] += 1
        print(f"gadget words, alphabet sigma + {len(gad)}, length <= {maxlen}, readout {Lname}: words with eigenvalue 1: {dict(res)}   [{time.time() - t0:.1f} s]")
rep = json.load(open(here + "/../twomark_rep_5_3.json"))
d = rep["d"]
A = [[[int(Fraction(x)) for x in row] for row in M] for M in rep["mats"]]
names = [(o,) for o in range(12)] + [(o1, o2) for o1 in range(12) for o2 in range(12)]
L = [1] + [0] * (d - 1)
res = Counter()
res_stay = Counter()
for w in names:
    M = A[w[0]] if len(w) == 1 else matmul(A[w[1]], A[w[0]])
    if 1 in int_roots(charpoly_int(M)):
        V, lam = class_space(L, [M], [1])
        key = "C = 0" if V.dim() == d else ("offset = 0 on all of C (H fixed)" if V.contains(lam[0]) else "offset != 0 somewhere on C (virtual translation)")
        res[key] += 1
        if all(o % 3 == 0 for o in w):
            res_stay[key] += 1
print(f"two-mark words of length <= 2 with eigenvalue 1 (exact rational linear algebra): {dict(res)}; among them words in the mu = 0 letters (two-apex): {dict(res_stay)}   [{time.time() - t0:.1f} s]")
