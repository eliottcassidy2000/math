#!/usr/bin/env python3
"""Audit P1 addendum (independent): the order-16 example D8_144 = D(B8) of THM-4995 (4b).

usage: python3 addendum_d8.py EXPDIR > addendum_d8_structure.out
  - parses d8_144.txt and d8_144_block.txt (code = upper triangle row-major, '1' iff i -> j);
  - rebuilds D(B8) with the audit's own double() (build_tournaments.py) and compares arc by arc;
  - B8 scores and twin distance (brute force over pairs);
  - writes d8_144.txt in the audit's "name n out_0 .. out_15" format to stdout-side file d8_144_masks.txt
    for tstruct / ham_dp / ham_ie;
  - checks the 'unequal double' structure of every order-15 tournament in d87_list.txt and writes them in mask format
    (d87_masks.txt) for the optional negative check.
"""
import sys
import os
from build_tournaments import double, is_tournament, masks


def parse(n, code):
    assert len(code) == n * (n - 1) // 2
    A = [[False] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[k] == '1':
                A[i][j] = True
            elif code[k] == '0':
                A[j][i] = True
            else:
                raise SystemExit('bad code char')
            k += 1
    return A


def twin_distance(B):
    q = len(B)
    best = None
    for u in range(q):
        for v in range(u + 1, q):
            dist = sum(1 for z in range(q) if z not in (u, v) and B[z][u] != B[z][v])
            best = dist if best is None else min(best, dist)
    return best


def main():
    ex = sys.argv[1]
    n16, code16, tag = open(os.path.join(ex, 'd8_144.txt')).read().split()
    n8, code8, tag8 = open(os.path.join(ex, 'd8_144_block.txt')).read().split()
    T = parse(int(n16), code16)
    B8 = parse(int(n8), code8)
    print(f"block code given in the task: 1010100100111001101001010101 ; file: {code8} ; equal: "
          f"{code8 == '1010100100111001101001010101'}")
    print("T is a tournament:", is_tournament(T), "; B8 is a tournament:", is_tournament(B8))
    D = double(B8)
    same = all(T[i][j] == D[i][j] for i in range(16) for j in range(16))
    print("d8_144 == D(B8) (A = 0..7, B' = 8..15, back arcs b_x -> a_x, other cross arcs A -> B'):", same)
    print("B8 scores (out-degrees):", [sum(r) for r in B8])
    print("B8 twin distance:", twin_distance(B8))
    with open('d8_144_masks.txt', 'w') as f:
        f.write('D8_144 16 ' + ' '.join(map(str, masks(T))) + '\n')
    # optional: d87 unequal doubles
    Q = {1, 2, 4}
    bad = 0
    cnt = 0
    with open('d87_masks.txt', 'w') as f:
        for line in open(os.path.join(ex, 'd87_list.txt')):
            p = line.split()
            if len(p) < 3:
                continue
            n, code, t = int(p[0]), p[1], p[2]
            U = parse(n, code)
            cnt += 1
            ok = (n == 15) and is_tournament(U)
            sc = [sum(U[x][y] for y in range(8)) for x in range(8)]
            ok &= sorted(sc) == [3] * 4 + [4] * 4
            ok &= all(U[8 + x][8 + y] == (((y - x) % 7) in Q) for x in range(7) for y in range(7) if x != y)
            back = [(y, a) for y in range(7) for a in range(8) if U[8 + y][a]]
            f_map = {}
            for y, a in back:
                f_map.setdefault(y, []).append(a)
            ok &= len(back) == 7 and all(len(v) == 1 for v in f_map.values()) and len(f_map) == 7
            if ok:
                img = {v[0] for v in f_map.values()}
                ok &= len(img) == 7 and {x for x in range(8) if sc[x] == 4} <= img
            if not ok:
                bad += 1
            f.write(f'{t} 15 ' + ' '.join(map(str, masks(U))) + '\n')
    print(f"d87_list.txt: {cnt} tournaments; structure violations (A near-regular, B' = QR_7, 7 back arcs b_y -> a_f(y),"
          f" f injective covering the 4 in-degree-3 vertices of A): {bad}")


if __name__ == '__main__':
    main()
