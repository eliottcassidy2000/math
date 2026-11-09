#!/usr/bin/env python3
"""paley_build.py - Paley doubles, chains and back-arc variants (paleyeven lane, mac-mini-2026-11-09).

Output lines "n code tag" in the tprop.c / mdefect.c format: code is the upper triangle in row-major order
(pairs (i, j), i < j), '1' iff i -> j.

Objects (vertex a_x = x, b_x = p + x, chain block c: c*p + x):
  double  p          T_p: two copies of QR_p, all arcs a_x -> b_y (x != y), back arcs b_x -> a_x.
  double  p sigma    back arcs b_x -> a_sigma(x) for a permutation sigma (affine x -> u*x + t).
  kreg    p E        back arcs b_x -> a_y for (x, y) in E (any bipartite pattern).
  chain   p r        r copies of QR_p in a row; block c sends all arcs forward to blocks c' > c, except the matched
                     back arcs from block c+1 to block c (x -> x); non-adjacent blocks have no back arcs.
usage: python3 paley_build.py   (prints the standard list used in the lane)
"""
import itertools, sys


def qr(p):
    return {(x * x) % p for x in range(1, p)}


def paley_arc(p, Q, x, y):
    return ((y - x) % p) in Q


def code_of(n, arc):
    return ''.join('1' if arc(i, j) else '0' for i in range(n) for j in range(i + 1, n))


def double(p, back):
    """back: set of (x, y) meaning b_x -> a_y."""
    Q = qr(p)
    n = 2 * p

    def arc(i, j):  # i < j
        if j < p:                    # both in A
            return paley_arc(p, Q, i, j)
        if i >= p:                   # both in B
            return paley_arc(p, Q, i - p, j - p)
        x, y = i, j - p              # a_x vs b_y
        return (y, x) not in back    # a_x -> b_y unless b_y -> a_x is a back arc
    return n, code_of(n, arc)


def chain(p, r, perm=None):
    Q = qr(p)
    n = r * p

    def arc(i, j):
        ci, xi = divmod(i, p)
        cj, xj = divmod(j, p)
        if ci == cj:
            return paley_arc(p, Q, xi, xj)
        # ci < cj (since i < j): forward unless matched back arc from block ci+1 to block ci
        if cj == ci + 1 and xj == xi:
            return False
        return True
    return n, code_of(n, arc)


def fano_incidence():
    # points 0..6, lines {x, x+1, x+3} mod 7 (difference set {0, 1, 3})
    return {(L, (L + d) % 7) for L in range(7) for d in (0, 1, 3)}


def main():
    out = []
    for p in (3, 7, 11):
        n, c = double(p, {(x, x) for x in range(p)})
        out.append((n, c, f'double_p{p}_identity'))
    # affine matchings for p = 7: sigma(x) = u x + t
    for u in (1, 2, 3, 6):
        for t in (0, 1):
            n, c = double(7, {(x, (u * x + t) % 7) for x in range(7)})
            out.append((n, c, f'double_p7_affine_u{u}_t{t}'))
    # 2-regular back-arc patterns for p = 7: b_x -> a_x, a_{x+d}
    for d in (1, 3):
        n, c = double(7, {(x, x) for x in range(7)} | {(x, (x + d) % 7) for x in range(7)})
        out.append((n, c, f'double_p7_2reg_d{d}'))
    # 3-regular: Fano incidence (b_line -> a_point)
    n, c = double(7, fano_incidence())
    out.append((n, c, 'double_p7_fano'))
    # chains
    for r in (3,):
        n, c = chain(7, r)
        out.append((n, c, f'chain_p7_r{r}'))
    for n, c, tag in out:
        print(n, c, tag)


if __name__ == '__main__':
    main()
