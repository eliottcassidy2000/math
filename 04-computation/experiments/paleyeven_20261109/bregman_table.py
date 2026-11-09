#!/usr/bin/env python3
"""bregman_table.py - Bregman upper bound for H of doubles D(B) of doubly regular B (order q) versus the mean
E[H(m)] = m!/2^(m-1) at m = 2q - 1 (for H(D)) and m = 2q - 2 (for every deletion; worst case deletes a b-vertex).
H(D) <= (n!)^(1/n) prod_v ((d+(v)+1)!)^(1/(d+(v)+1)); out-degrees (3q-3)/2 on A and (q+1)/2 on B'.
Prints q, log10 ratio for D, log10 ratio for the worst deletion; and the first q (= 3 mod 4) from which both stay < 0."""
from math import lgamma, log
def lf(k): return lgamma(k + 1)
def breg_log(outdegs):
    n = len(outdegs); s = lf(n) / n
    for d in outdegs: s += lf(d + 1) / (d + 1)
    return s
def logE(m): return lf(m) - (m - 1) * log(2)
rows = []
for q in range(7, 20001, 4):
    da = [(3 * q - 3) // 2] * q; db = [(q + 1) // 2] * q
    rows.append((q, (breg_log(da + db) - logE(2 * q - 1)) / log(10), (breg_log(da + db[1:]) - logE(2 * q - 2)) / log(10)))
first = next(q for q, a, b in rows if a < 0 and b < 0)
print("first q with both ratios < 0:", first)
print("all q >= first (to 20000) negative:", all(a < 0 and b < 0 for q, a, b in rows if q >= first))
print("deletion ratio strictly decreasing from first on:", all(rows[i + 1][2] < rows[i][2] for i in range(len(rows) - 1) if rows[i][0] >= first))
for q, a, b in rows:
    if q <= 103 or q in (199, 499, 999, 1999, 4999, 9999, 19999):
        print(q, round(a, 3), round(b, 3))
