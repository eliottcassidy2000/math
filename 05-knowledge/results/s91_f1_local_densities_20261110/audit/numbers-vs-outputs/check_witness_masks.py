#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: the three out-mask witnesses printed in THM-5129 (levels 7, 8, 9), recounted with an
own pure-Python DP: tournament check, H(B), U_B(y), E_B(y), U_even/U_odd, and H of the closed tournaments
B => TT_t => z (z -> y) for t = 0 (head), 1, 2, 3."""
from functools import lru_cache

W = {7: ((444, 273, 258, 6, 268, 286, 319, 382, 8), 8, 129, 144),
     8: ((82, 36, 129, 135, 142, 157, 190, 3), 7, 257, 288),
     9: ((402, 360, 3, 261, 270, 285, 317, 382, 4), 8, 513, 576)}


def ham_counts(out, n):
    """st[S][v] = number of Hamiltonian paths of T[S] starting at v"""
    full = (1 << n) - 1
    st = [[0] * n for _ in range(full + 1)]
    h = [0] * (full + 1)
    h[0] = 1
    for S in range(1, full + 1):
        for v in range(n):
            if (S >> v) & 1:
                R = S ^ (1 << v)
                if R == 0:
                    st[S][v] = 1
                else:
                    m = out[v] & R
                    c = 0
                    while m:
                        w = (m & -m).bit_length() - 1
                        m &= m - 1
                        c += st[R][w]
                    st[S][v] = c
        h[S] = sum(st[S])
    return st, h


def H(out, n):
    return ham_counts(out, n)[1][(1 << n) - 1]


for s, (masks, y, D, U) in W.items():
    n = len(masks)
    out = list(masks)
    ok = all(((out[i] >> j) & 1) + ((out[j] >> i) & 1) == 1 for i in range(n) for j in range(i + 1, n)) and all(not (out[i] >> i) & 1 for i in range(n)) and all(out[i] < (1 << n) for i in range(n))
    st, h = ham_counts(out, n)
    full = (1 << n) - 1
    u_even = sum(st[S][y] * h[full ^ S] for S in range(1, full + 1) if (S >> y) & 1 and bin(S).count("1") % 2 == 0)
    u_odd = sum(st[S][y] * h[full ^ S] for S in range(1, full + 1) if (S >> y) & 1 and bin(S).count("1") % 2 == 1)
    conv = [sum(1 << j for j in range(n) if j != i and not (out[i] >> j) & 1) for i in range(n)]
    e = ham_counts(conv, n)[0][full][y]          # paths of B ending at y = paths of the converse starting at y
    closed = []
    for t in range(0, 4):
        N = n + t + 1
        z = n + t
        o = [0] * N
        for v in range(n):
            o[v] = out[v]
            for i in range(t):
                o[v] |= 1 << (n + i)
            if v != y:
                o[v] |= 1 << z
        for i in range(t):
            for k in range(i + 1, t):
                o[n + i] |= 1 << (n + k)
            o[n + i] |= 1 << z
        o[z] = 1 << y
        closed.append(H(o, N))
    print(f"s = {s}: tournament {ok}; order {n}; H(B) = {h[full]} (claimed {D}); U = {u_even + u_odd} (claimed {U}); "
          f"U_even = {u_even}, U_odd = {u_odd}, E = {e} (U_odd - U_even = {u_odd - u_even}); outdeg(y) = {bin(out[y]).count('1')}")
    print(f"       closed counts t = 0..3: {closed};  D - E + U = {h[full] - e + u_even + u_odd}, D + 2 U_even = {h[full] + 2 * u_even}, D + 2^t U = {[h[full] + (1 << t) * (u_even + u_odd) for t in (1, 2, 3)]}")
