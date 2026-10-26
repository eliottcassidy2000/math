# Character sums of graphs on 4 and 5 points, the building blocks of the subtournament census of Paley structures.
# For a graph G on vertices 0..k-1 (edges i<j), T_G(p) = sum over distinct v in F_p^k with v_0 = 0, v_1 = 1 of
# prod_(ij in G) chi(v_j - v_i)  (translation and dilation normalised; chi^|E| = 1 needs |E| even, else use chi(-1)).
# (4) K4: T = sum_{x,y} chi(x)chi(y)chi(x-1)chi(y-1)chi(y-x) -- Greene's 3F2(1) / the weight-3 CM form eta(4z)^6,
#     which vanishes for p = 3 mod 4: the 4-census of Paley TOURNAMENTS is polynomial, of Paley GRAPHS is not.
# (5) the dense 5-point graphs K5, K5 - e, K5 - 2K2 whose sums carry the 5-census arithmetic.
import sys
import numpy as np
from sympy import primerange

def chi_table(p):
    t = np.full(p, -1, dtype=np.int64); t[0] = 0
    for x in range(1, p): t[(x * x) % p] = 1
    return t

def T_K4(p):
    c = chi_table(p); xs = np.arange(p)
    tot = 0
    for x in range(2, p):
        y = xs[(xs != 0) & (xs != 1) & (xs != x)]
        tot += int(c[x] * c[x - 1] * np.sum(c[y] * c[(y - 1) % p] * c[(y - x) % p]))
    return tot

def T5(p, edges):
    """vertices 0..4 with v0 = 0, v1 = 1, v2 = x, v3 = y, v4 = z"""
    c = chi_table(p); xs = np.arange(p)
    tot = 0
    Y, Z = np.meshgrid(xs, xs, indexing='ij')
    base = (Y != Z) & (Y != 0) & (Y != 1) & (Z != 0) & (Z != 1)
    for x in range(2, p):
        V = [0, 1, x, Y, Z]
        mask = base & (Y != x) & (Z != x)
        prod = np.ones_like(Y)
        for i, j in edges:
            prod = prod * c[(np.asarray(V[j]) - np.asarray(V[i])) % p] if (i >= 2 or j >= 2) else prod * c[(V[j] - V[i]) % p]
        tot += int(np.sum(prod[mask]))
    return tot

K5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
K5e = [e for e in K5 if e != (3, 4)]                  # K5 minus the edge {3,4} (both ends free)
K5mm = [e for e in K5 if e not in ((1, 2), (3, 4))]   # K5 minus a matching {1,2},{3,4}
if __name__ == '__main__':
    Pmax = int(sys.argv[1]) if len(sys.argv) > 1 else 110
    for p in primerange(5, Pmax):
        t4 = T_K4(p)
        t5 = T5(p, K5); t5e = T5(p, K5e); t5m = T5(p, K5mm)
        print(f"p={p} mod4={p % 4} K4={t4} K5={t5} K5-e={t5e} K5-2K2={t5m}", flush=True)
