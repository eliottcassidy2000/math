#!/usr/bin/env python3
"""Audit P1: independent builders for THM-4995 (5) (negative examples) and the chains; output "name n out_0..".
  double with back pattern P: a_x = x, b_x = q + x; within halves B; a_x -> b_y unless (x, y) in P means b_y -> a_x.
  chain(B, r, D): blocks c = 0..r-1 (vertex c*q + x); within block B; between blocks c < c' all arcs forward,
                  except back arcs (c+1, z) -> (c, z + d) for d in D (adjacent blocks only)."""
from build_tournaments import qr_tour, circ_tour, is_tournament, masks

def double_pattern(B, back):            # back: set of (y, x): b_y -> a_x
    q = len(B); n = 2 * q
    A = [[False] * n for _ in range(n)]
    for x in range(q):
        for y in range(q):
            if B[x][y]:
                A[x][y] = True; A[q + x][q + y] = True
    for x in range(q):
        for y in range(q):
            if (y, x) in back: A[q + y][x] = True
            else: A[x][q + y] = True
    return A

def chain(B, r, D):
    q = len(B); n = q * r
    A = [[False] * n for _ in range(n)]
    for c in range(r):
        for x in range(q):
            for y in range(q):
                if B[x][y]: A[c * q + x][c * q + y] = True
    for c in range(r):
        for c2 in range(c + 1, r):
            for x in range(q):
                for z in range(q):
                    if c2 == c + 1 and any((x - z - d) % q == 0 for d in D):
                        A[c2 * q + z][c * q + x] = True     # back arc (c+1, z) -> (c, z + d)
                    else:
                        A[c * q + x][c2 * q + z] = True
    return A

Q7 = qr_tour(7)
objs = [
    ('T7_2reg_d1', double_pattern(Q7, {(x, x) for x in range(7)} | {(x, (x + 1) % 7) for x in range(7)})),
    ('T7_2reg_d3', double_pattern(Q7, {(x, x) for x in range(7)} | {(x, (x + 3) % 7) for x in range(7)})),
    ('T7_fano', double_pattern(Q7, {(L, (L + d) % 7) for L in range(7) for d in (0, 1, 3)})),
    ('chain_QR7_r3', chain(Q7, 3, [0])),
    ('chain_R5_r3_D02', chain(circ_tour(5, [1, 2]), 3, [0, 2])),
    ('D_C9_1234', chain(circ_tour(9, [1, 2, 3, 4]), 2, [0])),
    ('T7_check_as_chain', chain(Q7, 2, [0])),
]
for name, A in objs:
    assert is_tournament(A), name
    print(name, len(A), *masks(A))
