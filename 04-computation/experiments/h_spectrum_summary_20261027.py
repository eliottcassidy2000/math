#!/usr/bin/env python3
"""h_spectrum_summary_20261027.py -- opus-2026-10-27-S63
Summaries over the exhaustive spectrum outputs h_spectrum_sweep_20261027_n{3..8}.out (labelled counts, class counts,
strong-class counts per Hamiltonian-path count H): the sets A_n and S_n, first-instance orders n(v) and n_s(v), the
complete strong fibres F_s(v) (number of strong classes with H = v over all orders) for v <= 74 (complete because
f(9) = 75 > 74 by augmentation / Busch), the values with a unique strong source, the largest fibres, and the odd values
below 610 that no strong tournament of order <= 8 attains."""
import os, sys, re
D = os.path.dirname(os.path.abspath(__file__))
tab = {}   # n -> {v: (labelled, classes, strong_classes, strong_labelled)}
for n in range(3, 9):
    t = {}
    for line in open(os.path.join(D, f"h_spectrum_sweep_20261027_n{n}.out"), encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) == 5 and p[0].isdigit(): t[int(p[0])] = tuple(int(x) for x in p[1:])
    tab[n] = t
A = {n: sorted(tab[n]) for n in tab}; S = {n: sorted(v for v in tab[n] if tab[n][v][2] > 0) for n in tab}
print("order n: |A_n| |S_n| f(n) max(n) classes strong_classes")
for n in tab:
    cl = sum(x[1] for x in tab[n].values()); st = sum(x[2] for x in tab[n].values())
    print(f"{n}: {len(A[n])} {len(S[n])} {S[n][0]} {max(A[n])} {cl} {st}")
first = {}; firsts = {}
for n in tab:
    for v in A[n]: first.setdefault(v, n)
    for v in S[n]: firsts.setdefault(v, n)
print("\nfirst-instance orders n(v) (any tournament) and n_s(v) (strong), odd v <= 189:")
print(" ".join(f"{v}:{first.get(v,'-')}/{firsts.get(v,'-')}" for v in range(1, 190, 2)))
print("\nodd v < 610 missing from every S_n, n <= 8:", [v for v in range(1, 610, 2) if v not in firsts])
print("odd v < 610 missing from every A_n, n <= 8:", [v for v in range(1, 610, 2) if v not in first])
print("values first attained at order 8 (not in A_7):", [v for v in A[8] if v not in first or first[v] == 8][:60])
Fs = {}
for n in tab:
    for v, x in tab[n].items():
        if x[2]: Fs[v] = Fs.get(v, 0) + x[2]
print("\ncomplete strong fibres F_s(v) = number of strong classes with H = v, all orders (complete for v <= 74):")
print(" ".join(f"{v}:{Fs.get(v, 0)}" for v in range(1, 75, 2)))
print("unique strong source (F_s = 1), v <= 74:", [v for v in range(3, 75, 2) if Fs.get(v, 0) == 1])
print("empty strong fibre, v <= 74:", [v for v in range(1, 75, 2) if Fs.get(v, 0) == 0])
big = sorted(((Fs[v], v) for v in Fs if v <= 74), reverse=True)[:5]
print("largest strong fibres, v <= 74:", [(v, c) for c, v in big])
for n in (7, 8):
    cl = tab[n]; top = sorted(((x[1], v) for v, x in cl.items()), reverse=True)[:5]
    print(f"order {n}: most populous values (classes):", [(v, c) for c, v in top], "; values with a single class:", sum(1 for x in cl.values() if x[1] == 1), "of", len(cl))
