#!/usr/bin/env python3
"""chain_build.py - block chains of circulant tournaments glued by difference-set back arcs (paleyeven lane).

A chain  C(m, S; r; D)  has r blocks B_0..B_{r-1}, each a copy of the circulant tournament on Z_m with connection set S
(x -> x + s for s in S; S cap -S = empty, S cup -S = Z_m \ {0}; S = squares gives Paley QR_p).  Between blocks c < c':
every arc goes forward B_c -> B_c', except the back arcs  (c+1, z) -> (c, z + d)  for d in D  (adjacent blocks only).
With D = {0} and S = QR_p this is the Paley chain of chain_build/paley_build; r = 2 is the Paley double T_p.

usage: python3 chain_build.py m "s1,s2,.." r "d1,d2,.." [tag]   -> prints "n code tag"
       python3 chain_build.py --family                           -> prints the lane's standard list
"""
import sys


def build(m, S, r, D):
    S = {s % m for s in S}
    D = {d % m for d in D}
    n = m * r

    def arc(i, j):  # i < j: True iff i -> j
        ci, xi = divmod(i, m)
        cj, xj = divmod(j, m)
        if ci == cj:
            return ((xj - xi) % m) in S
        if cj == ci + 1 and ((xi - xj) % m) in D:   # back arc (cj, xj) -> (ci, xj + d) with xi = xj + d
            return False
        return True
    code = ''.join('1' if arc(i, j) else '0' for i in range(n) for j in range(i + 1, n))
    return n, code


def qr(p):
    return sorted({(x * x) % p for x in range(1, p)})


FAMILY = [
    # (m, S, r, D, tag)
    (5, [1, 2], 3, [0, 2], 'chain_R5_r3_D02'),
    (5, [1, 2], 3, [0, 1], 'chain_R5_r3_D01'),
    (5, [1, 2], 3, [0, 1, 3], 'chain_R5_r3_D013'),
    (5, [1, 2], 4, [0, 2], 'chain_R5_r4_D02'),
    (3, [1], 5, [0, 1], 'chain_C3_r5_D01'),
    (3, [1], 6, [0, 1], 'chain_C3_r6_D01'),
    (7, qr(7), 2, [0], 'double_P7'),
    (7, qr(7), 3, [0], 'chain_P7_r3_D0'),
    (7, qr(7), 2, [0, 1], 'double_P7_D01'),
    (7, [1, 2, 3], 2, [0], 'double_R7_D0'),
    (7, [1, 2, 3], 3, [0], 'chain_R7_r3_D0'),
    (9, [1, 2, 3, 4], 2, [0], 'double_R9_D0'),
    (5, [1, 2], 2, [0, 2], 'double_R5_D02'),
    (5, [1, 2], 2, [0, 1, 3], 'double_R5_D013'),
]


def main():
    if len(sys.argv) > 1 and sys.argv[1] != '--family':
        m = int(sys.argv[1]); S = [int(t) for t in sys.argv[2].split(',')]
        r = int(sys.argv[3]); D = [int(t) for t in sys.argv[4].split(',')]
        tag = sys.argv[5] if len(sys.argv) > 5 else f'chain_m{m}_r{r}'
        n, code = build(m, S, r, D)
        print(n, code, tag)
        return
    for m, S, r, D, tag in FAMILY:
        n, code = build(m, S, r, D)
        print(n, code, tag)


if __name__ == '__main__':
    main()
