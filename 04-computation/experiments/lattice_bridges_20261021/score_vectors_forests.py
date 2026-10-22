# Score vectors of labelled n-tournaments = integer points of the permutohedron conv{sigma(0,1,...,n-1)} (Landau's theorem),
# and their number equals the number of labelled forests on n vertices (Stanley).  The permutohedron is (up to scaling) the
# Voronoi cell of the dual root lattice A_(n-1)^*; reversing one arc u -> v moves the score vector by the root e_v - e_u.
# Check: count distinct score vectors by Landau's criterion for n <= 9 and compare with forests (A001858).
# mac-mini-2026-10-21-fanopsp.
from itertools import product
from math import comb
def landau_count(n):
    cnt = 0
    for s in product(range(n), repeat=n):
        if sum(s) != comb(n, 2): continue
        t = sorted(s); ok = True; acc = 0
        for k in range(1, n + 1):
            acc += t[k - 1]
            if acc < comb(k, 2): ok = False; break
        if ok: cnt += 1
    return cnt
forests = [1, 1, 2, 7, 38, 291, 2932, 36961, 561948, 10026505]   # A001858, n = 0..9
for n in range(1, 9):
    c = landau_count(n)
    print(f"n={n}: score vectors (Landau) = {c}, labelled forests = {forests[n]}, equal: {c == forests[n]}")
s25 = sum(k * k for k in range(25))
print("TT_25 score vector (0..24): sum of squares =", s25, "= 70^2:", s25 == 70 ** 2)
