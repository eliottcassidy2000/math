#!/usr/bin/env python3
"""Verifier: the gadget words with eigenvalue 1 whose offset is not identically 0 on C(M; 1): which word, which class space."""
import json, os
from vlin import charpoly_int, class_space, matmul, rowmat
from vroots import int_roots
from v2_onemark_pairs import ann_basis
here = os.path.dirname(os.path.abspath(__file__))
SIGMA = [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]]
for nmax, maxlen in ((4, 2), (3, 3)):
    gad = json.load(open(f"{here}/../gadget_matrices_{nmax}.json"))
    letters = [("s", SIGMA)] + [(f"G(n={g['n']},X={g['adj']},P={g['P']:b})", g["G"]) for g in gad]
    cur = [((), [[1 if i == j else 0 for j in range(4)] for i in range(4)])]
    words = []
    for ln in range(maxlen):
        cur = [(w + (nm,), matmul(G, M)) for (w, M) in cur for (nm, G) in letters]
        words += cur
    for Lname, L in (("L_h", [0, 0, 1, 0]), ("L_1", [1, 2, 0, 0])):
        for w, M in words:
            if 1 in int_roots(charpoly_int(M)):
                V, lam = class_space(L, [M], [1])
                if V.dim() < 4 and not V.contains(lam[0]):
                    print(f"alphabet {nmax}, readout {Lname}: word {' . '.join(w)}: matrix {M}, eigenvalues {int_roots(charpoly_int(M))}, lambda {lam[0]}, C basis (r, M, h, E): {ann_basis(V, 4)}")
