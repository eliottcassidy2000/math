#!/usr/bin/env python3
"""Exact bias of the parity of the number of cyclic triangles: E (-1)^{c_3(T)} over labelled tournaments of order n.
c_3 = C(n,3) - sum_v C(s_v,2) and C(s,2) mod 2 is the degree-2 elementary symmetric polynomial of the out-arc
indicators at v, so c_3 mod 2 is a quadratic polynomial over F_2 in the C(n,2) arc variables.  Its character sum is
computed exactly by reducing the quadratic form (and by brute force for n <= 6 as a control)."""
from itertools import combinations
from fractions import Fraction
import sys
def quad_c3(n):
    prs = list(combinations(range(n), 2)); idx = {p: k for k, p in enumerate(prs)}; N = len(prs)
    # polynomial: const c, linear lin[k], quadratic Q[k] (bitmask of k' > k)
    c = (n * (n - 1) * (n - 2) // 6) & 1
    lin = [0] * N; Q = [0] * N
    def addprod(a, b):   # a, b = (var index, const) meaning x_var + const ; add (x_a + ca)(x_b + cb)
        nonlocal c
        (va, ca), (vb, cb) = a, b
        # x_a x_b + cb x_a + ca x_b + ca cb
        lo, hi = min(va, vb), max(va, vb)
        Q[lo] ^= 1 << hi
        if cb: lin[va] ^= 1
        if ca: lin[vb] ^= 1
        if ca and cb: c ^= 1
    for v in range(n):
        terms = []
        for w in range(n):
            if w == v: continue
            if v < w: terms.append((idx[(v, w)], 0))     # out-arc indicator x_vw
            else: terms.append((idx[(w, v)], 1))         # out-arc indicator 1 + x_wv
        for a, b in combinations(terms, 2): addprod(a, b)
    return N, c, lin, Q
def gauss(N, c, lin, Q):
    """sum over x in F_2^N of (-1)^(c + lin.x + x^T Q x) by elimination; returns exact integer"""
    # symmetric polar matrix
    B = [0] * N
    for i in range(N):
        for j in range(i + 1, N):
            if Q[i] >> j & 1: B[i] |= 1 << j; B[j] |= 1 << i
    lin = lin[:]; alive = list(range(N)); total_sign = -1 if c else 1; scale = 1
    # iterative reduction: pick i, j with B[i][j] = 1; summing over x_i, x_j gives factor 2 * (-1)^(l_i l_j) after
    # substituting; standard: f = x_i x_j + x_i A + x_j Bf + rest, A, Bf linear in others (incl. constants)
    alive_set = set(alive)
    while True:
        pair = None
        for i in list(alive_set):
            m = B[i]
            if m:
                j = (m & -m).bit_length() - 1
                # ensure j alive (B only among alive vars)
                pair = (i, j); break
        if pair is None: break
        i, j = pair
        Ai = B[i] & ~(1 << j); Aj = B[j] & ~(1 << i)   # linear forms (other vars) multiplying x_i and x_j
        li, lj = lin[i], lin[j]
        # sum_{x_i,x_j} (-1)^(x_i x_j + x_i(A+li) + x_j(Bf+lj)) = 2 * (-1)^((A+li)(Bf+lj))
        scale *= 2
        # add product (Ai + li)(Aj + lj) to the remaining form
        for a in range(N):
            if Ai >> a & 1:
                # x_a * (Aj + lj)
                B[a] ^= Aj & ~(1 << a)
                if Aj >> a & 1: lin[a] ^= 1          # x_a^2 = x_a
                if lj: lin[a] ^= 1
        if li:
            for a in range(N):
                if Aj >> a & 1: lin[a] ^= 1
            if lj: total_sign = -total_sign
        # remove i, j
        for a in range(N):
            B[a] &= ~((1 << i) | (1 << j))
        B[i] = B[j] = 0; alive_set.discard(i); alive_set.discard(j)
        # symmetrise is preserved because we added Ai x Aj + Aj x Ai entries: fix symmetry explicitly
        for a in range(N):
            if Aj >> a & 1:
                B[a] ^= Ai & ~(1 << a)
        # the two loops together add the symmetric matrix Ai Aj^T + Aj Ai^T (diagonal handled through lin)
    # remaining variables: linear form only
    for a in alive_set:
        if lin[a]: return 0
        scale *= 2
    return total_sign * scale
def brute(n):
    N, c, lin, Q = quad_c3(n); tot = 0
    for x in range(1 << N):
        v = c
        for k in range(N):
            if x >> k & 1:
                v ^= lin[k] ^ (bin(Q[k] & x).count('1') & 1)
        tot += -1 if v else 1
    return tot
for n in range(3, 7):
    N, c, lin, Q = quad_c3(n); g = gauss(N, c, lin, Q); b = brute(n)
    assert g == b, (n, g, b)
print("control: elimination = brute force for n = 3..6")
for n in range(3, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 15):
    N, c, lin, Q = quad_c3(n); g = gauss(N, c, lin, Q)
    print("n=%2d  sum (-1)^c3 = %d  ;  E(-1)^c3 = %s" % (n, g, Fraction(g, 1 << N)))
