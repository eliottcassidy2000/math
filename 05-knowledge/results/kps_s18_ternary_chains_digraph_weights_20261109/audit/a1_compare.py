#!/usr/bin/env python3
"""AUDIT: compare my black-box matrix sets (a1_blocks.py) with the packet's chain_blocks.py output, and rebuild
every generator in out/generators_base3_15.txt and out/generators_base2_13.txt from its block representative.

For each generator line "matrix | k rows w u'" the block is checked to be a tournament, glued after the test
tournaments of a1_blocks.py, and its transfer matrix recomputed from scratch; it must equal the printed matrix,
and the block order must be <= 4.
Usage: a1_compare.py MYDIR PACKETDIR
"""
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
from a1_blocks import TESTS, glue, state_polys, ev, solve3  # noqa: E402


def load_mine(path):
    d = {}
    for line in open(path):
        kmin, mat, rep = line.split("|")
        d[tuple(map(int, mat.split()))] = int(kmin)
    return d


def load_theirs(path):
    d = {}
    for line in open(path):
        mat, rep = line.split("|")
        d[tuple(map(int, mat.split()))] = int(rep.split()[0])
    return d


def matrix_of(B, w, up, x):
    old = [[ev(p, x) for p in state_polys(len(T), T, u)] for (T, u) in TESTS]
    new = []
    for (T, u) in TESTS:
        A, port = glue(T, u, B, w, up)
        new.append([ev(p, x) for p in state_polys(len(A), A, port)])
    S = [[old[c][r] for c in range(3)] for r in range(3)]
    Sp = [[new[c][r] for c in range(3)] for r in range(3)]
    N = solve3(S, Sp)
    for t in range(3, len(TESTS)):
        assert [sum(N[r][c] * old[t][c] for c in range(3)) for r in range(3)] == new[t]
    assert all(a.denominator == 1 for row in N for a in row)
    return tuple(int(a) for row in N for a in row)


def main():
    mydir, pk = sys.argv[1], sys.argv[2]
    for x in (2, 3):
        mine = load_mine(f"{mydir}/mats_x{x}.txt")
        theirs = load_theirs(f"{mydir}/theirs_x{x}_k5.txt.reps")
        for kk in (4, 5):
            a = {m for m, k in mine.items() if k <= kk}
            b = {m for m, k in theirs.items() if k <= kk}
            print(f"x={x} k<={kk}: mine {len(a)}, packet {len(b)}, identical sets: {a == b}")
    for x, fn in ((3, "generators_base3_15.txt"), (2, "generators_base2_13.txt")):
        mine = load_mine(f"{mydir}/mats_x{x}.txt")
        n = ok = 0
        for line in open(f"{pk}/out/{fn}"):
            if line.startswith("#"):
                continue
            mat, rep = line.split("|")
            M = tuple(map(int, mat.split()))
            k, rows, w, up = rep.split()
            k, w, up = int(k), int(w), int(up)
            B = [[int(rows[i * k + j]) for j in range(k)] for i in range(k)]
            is_t = all(B[i][i] == 0 for i in range(k)) and all(B[i][j] + B[j][i] == 1
                                                               for i in range(k) for j in range(k) if i != j)
            N = matrix_of(B, w, up, x)
            good = is_t and k <= 4 and N == M and M in mine and mine[M] <= 4
            n += 1
            ok += good
            if not good:
                print("GENERATOR MISMATCH", fn, line.strip(), N, is_t)
        print(f"{fn}: {n} generators, {ok} rebuilt from their block (tournament, order <= 4, matrix equal)")


if __name__ == "__main__":
    main()
