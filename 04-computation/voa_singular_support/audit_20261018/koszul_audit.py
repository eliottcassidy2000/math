#!/usr/bin/env python3
"""Audit of section 1.3 of the results note (derived arc space of the fat point C[x]/(x^N)).

Checks, for given N and weight bound W:
  (a) prod_{n>=0}(1-q^{2N+n}) / prod_{j>=0}(1-q^{j+2})  ==  1/prod_{m=2}^{2N-1}(1-q^m)   (series);
  (b) HS(C[x]/J), J = (x_0^N)_d, computed by linear algebra, == prod_{n != 0,+-1 mod 2N+1}(1-q^n)^{-1};
  (c) Koszul homology H_i (i=0..3) of the sequence (d^n f)_{n>=0} on C[x_0,x_1,...] by weight,
      sum_i (-1)^i HS(H_i) == chi, and the first weight where H_1 != 0 (claimed 2N+3).
Exact ranks over Q (flint.fmpq_mat).  Independent of the author's code.
Usage: python koszul_audit.py N W
"""
import sys
from itertools import combinations
import flint
from alin import rank
from virmod import partitions_min2, poly_mul, poly_d

Q = flint.fmpq
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2
W = int(sys.argv[2]) if len(sys.argv) > 2 else 18


def series_mul(a, b, M):
    out = [0] * (M + 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b[:M + 1 - i]):
            out[i + j] += x * y
    return out


def inv_one_minus(e, M):
    s = [0] * (M + 1)
    for i in range(0, M + 1, e):
        s[i] = 1
    return s


M = W
# (a)
num = [1] + [0] * M
for n in range(0, M + 1):
    e = 2 * N + n
    if e > M:
        break
    t = [0] * (M + 1)
    t[0] = 1
    t[e] = -1
    num = series_mul(num, t, M)
den_inv = [1] + [0] * M
for j in range(0, M + 1):
    e = j + 2
    if e > M:
        break
    den_inv = series_mul(den_inv, inv_one_minus(e, M), M)
chi = series_mul(num, den_inv, M)
fin = [1] + [0] * M
for m in range(2, 2 * N):
    fin = series_mul(fin, inv_one_minus(m, M), M)
print(f"N={N}: chi        =", chi)
print(f"N={N}: finite prod=", fin, " equal:", chi == fin)

# (b) Andrews-Gordon product
ag = [1] + [0] * M
for n in range(1, M + 1):
    r = n % (2 * N + 1)
    if r in (0, 1, 2 * N):
        continue
    ag = series_mul(ag, inv_one_minus(n, M), M)

f = {(2,) * N: Q(1)}
dfs = [f]
for n in range(1, M + 1):
    dfs.append(poly_d(dfs[-1]))
bases = {w: partitions_min2(w) for w in range(M + 1)}


def Jspan(w):
    gens = []
    for n in range(0, w - 2 * N + 1):
        for m in bases[w - 2 * N - n]:
            gens.append(poly_mul({m: Q(1)}, dfs[n]))
    return gens


hs0 = []
rk1 = {}
for w in range(M + 1):
    g = Jspan(w)
    r = rank(g, bases[w]) if g else 0
    rk1[w] = r
    hs0.append(len(bases[w]) - r)
print(f"N={N}: HS(arc)    =", hs0)
print(f"N={N}: AG product =", ag, " equal:", hs0 == ag)


# (c) Koszul complex ranks
def K_basis(i, w):
    """basis of K_i in weight w: (tuple of distinct n's (increasing), monomial)."""
    out = []
    # minimal weight of e_{n1}..e_{ni}: sum (2N+n)
    def rec(start, cur, wt):
        if len(cur) == i:
            for m in bases[w - wt]:
                out.append((tuple(cur), m))
            return
        n = start
        while wt + 2 * N + n + (i - len(cur) - 1) * (2 * N + n + 1) <= w:
            cur.append(n)
            rec(n + 1, cur, wt + 2 * N + n)
            cur.pop()
            n += 1
    if i == 0:
        return [((), m) for m in bases[w]]
    rec(0, [], 0)
    return out


def boundary(i, w):
    """images of K_i basis (weight w) in K_{i-1}, as dicts keyed by (ns, monomial)."""
    imgs = []
    for ns, m in K_basis(i, w):
        img = {}
        for pos, n in enumerate(ns):
            sign = (-1) ** pos
            rest = ns[:pos] + ns[pos + 1:]
            for mm, cc in poly_mul({m: Q(1)}, dfs[n]).items():
                key = (rest, mm)
                v = img.get(key, 0) + sign * cc
                if v == 0:
                    img.pop(key, None)
                else:
                    img[key] = v
        imgs.append(img)
    return imgs


H = {i: [] for i in range(4)}
for w in range(M + 1):
    dims = [len(K_basis(i, w)) for i in range(5)]
    rks = [0]
    for i in range(1, 5):
        if dims[i] == 0:
            rks.append(0)
            continue
        imgs = boundary(i, w)
        tgt = K_basis(i - 1, w)
        rks.append(rank(imgs, tgt))
    for i in range(4):
        H[i].append(dims[i] - rks[i] - rks[i + 1])
for i in range(4):
    print(f"N={N}: dim H_{i} by weight =", H[i])
alt = [H[0][w] - H[1][w] + H[2][w] - H[3][w] for w in range(M + 1)]
print(f"N={N}: sum(-1)^i HS(H_i) (i<=3) == chi :", alt == chi,
      " (valid while K_4 = 0, i.e. w <", 4 * 2 * N + 6, ")")
first = next((w for w in range(M + 1) if H[1][w] != 0), None)
print(f"N={N}: first weight with H_1 != 0: {first} (claimed 2N+3 = {2 * N + 3}), dim there = "
      f"{H[1][first] if first is not None else None}")
print(f"N={N}: HS(arc) - chi =", [hs0[w] - chi[w] for w in range(M + 1)])
