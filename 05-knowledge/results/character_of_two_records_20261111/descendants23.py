#!/usr/bin/env python3
"""After the audit refuted the class-wide form of HYP-9667 at N = 24: the doubly regular form on ALL descendants.
For each of the 13 classes of drt23_classes.json and the double of P_11: the skew Hadamard matrix of (source) + T,
its 24 descendants (all doubly regular), their odd cycle parities by length, H mod 256, and the number of
isomorphism classes among them; and, for contrast, the 23 single-vertex switchings of (source) + T.
Usage: descendants23.py CYC LAYERED JSON"""
import sys, json
from drt_zoo import descendant, double, paley_matrix, is_drt, is_skew_hadamard, classes, layered_H
from cycle_parities import cycles, code_to_out, odd_vector
from switch_class import plus, switch
def matrix_of(out):
    """skew Hadamard matrix of (source) + T: M = [[1, 1], [-1, I + S_T]]."""
    n = len(out); N = n + 1
    M = [[1] * N] + [[-1] + [1 if a == b else (1 if (out[a] >> b) & 1 else -1) for b in range(n)] for a in range(n)]
    assert is_skew_hadamard(M); return M
if __name__ == "__main__":
    cyc, lay, js = sys.argv[1], sys.argv[2], sys.argv[3]
    d = json.load(open(js))
    rows = [(rec['label'], code_to_out(rec['code'], 23)) for rec in d] + [("double of P_11", descendant(double(paley_matrix(11)), 0))]
    tot = bad_par = bad8 = 0; allreps = []
    for name, T in rows:
        M = matrix_of(T); desc = [descendant(M, c) for c in range(24)]; assert all(is_drt(o) for o in desc)
        vec = {}; h8 = {}
        for o in desc:
            v = odd_vector(cycles(cyc, o), 23); vec[v] = vec.get(v, 0) + 1
            h = layered_H(lay, o); h8[h % 8] = h8.get(h % 8, 0) + 1
            tot += 1; bad_par += (v != '0' * 11); bad8 += (h % 8 != 5)
        ncl = len(classes(desc)); classes(desc, allreps)
        U = plus(T); sv = {}
        for v in range(1, 24):
            w = odd_vector(cycles(cyc, switch(U, 1 << v)), 24); sv[w] = sv.get(w, 0) + 1
        print("%-15s 24 descendants in %d isomorphism classes: odd cycle parities %s; H mod 8 %s | 23 single-vertex switchings of (source) + T: %s" % (
            name, ncl, sorted(vec.items()), sorted(h8.items()), sorted(sv.items()))); sys.stdout.flush()
    print("descendants: %d doubly regular tournaments of order 23 in %d isomorphism classes; odd cycle parities not all even: %d; H != 5 (mod 8): %d" % (tot, len(allreps), bad_par, bad8))
