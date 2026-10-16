"""Referee check (A): the double chains D(n, j) and the failure of D_in > 0.

D(n, j) = transitive on 0..n-1 (higher beats lower) with the arcs reversed so that 0 -> j and j -> n-1.
Part 1: brute force over ALL n! permutations (no DFS, no DP) for D(11, 5), D(12, 6) and also D(11, 6), D(12, 7):
        H, every p(u before 0), E[N] - d/2 at w = 0, best in-neighbour, plus vertex-lemma / majority-source /
        F strong / M transitive for the whole tournament.
Part 2: independent explicit DFS for every j and 6 <= n <= 18: min vertex-lemma value (and ALL j attaining it),
        min majority-source value, min D_in (and j), M transitive (triangle test), F strong (Warshall closure),
        strongness of D(n, j) (Landau score test).  Also D(n,1) vs T_n (H = 2^(n-2) + 1).
"""
import sys, time
from fractions import Fraction as Fr
import numpy as np
from pa_core import (transitive_plus, brute_force_perms, enumerate_paths_dfs, strong_landau, strong_warshall,
                     transitive_by_triangles)


def stats(n, A, H, C):
    H = int(H)
    vl = None; vl_w = None; ms = None; dmin = None; dmin_w = None
    for w in range(n):
        ins = [u for u in range(n) if u != w and A[u, w]]
        best = max(Fr(int(C[u, w]), H) for u in ins)
        anyb = max(Fr(int(C[x, w]), H) for x in range(n) if x != w)
        D = sum(Fr(2 * int(C[u, w]) - H, H) for u in ins)
        if vl is None or best < vl:
            vl, vl_w = best, w
        if ms is None or anyb < ms:
            ms = anyb
        if dmin is None or D < dmin:
            dmin, dmin_w = D, w
    M = np.zeros((n, n), dtype=np.int8)
    Fm = np.zeros((n, n), dtype=np.int8)
    for a in range(n):
        for b in range(n):
            if a != b and 2 * int(C[a, b]) > H:
                M[a, b] = 1
                if A[a, b]:
                    Fm[a, b] = 1
    return vl, vl_w, ms, dmin, dmin_w, transitive_by_triangles(n, M), strong_warshall(n, Fm)


def part1(n, j):
    A = transitive_plus(n, [(0, j), (j, n - 1)])
    t0 = time.time()
    H, C, Nd, L = brute_force_perms(n, A)
    H = int(H)
    w = 0
    U = [u for u in range(n) if u != w and A[u, w]]
    outs = [x for x in range(n) if A[w, x]]
    pu = {u: Fr(int(C[u, w]), H) for u in U}
    EN = sum(pu.values())
    EN2 = Fr(sum(k * int(Nd[w, k]) for k in range(n + 1)), H)
    d = len(U)
    print(f'D({n},{j}) [brute force over {n}! permutations, {time.time() - t0:.1f} s]: H = {H}; out-neighbours of 0: {outs}; d = {d}; '
          f'strong (Landau) = {strong_landau(n, A)}')
    print('   p(u before 0): ' + ', '.join(f'{u}: {pu[u]}' for u in U))
    print(f'   E[N] = {EN} (from N-distribution: {EN2}); E[N] - d/2 = {EN - Fr(d, 2)}; D_in(0) = {2 * EN - d}; '
          f'best in-neighbour p = {max(pu.values())}')
    vl, vl_w, ms, dmin, dmin_w, mtrans, fstrong = stats(n, A, H, C)
    print(f'   whole tournament: vertex-lemma value {vl} (at w = {vl_w}); majority-source value {ms}; min D_in {dmin} (at w = {dmin_w}); '
          f'M transitive: {mtrans}; F strong: {fstrong}')
    return H, C


def part2(nmax):
    for n in range(6, nmax + 1):
        t0 = time.time()
        rows = []
        for j in range(1, n - 1):
            A = transitive_plus(n, [(0, j), (j, n - 1)])
            C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
            H = enumerate_paths_dfs(n, A, C, Nd, L)
            # sanity: Nd row sums equal H, and E[N] from Nd equals sum of p over in-neighbours
            for w in range(n):
                assert int(Nd[w].sum()) == H
                assert sum(k * int(Nd[w, k]) for k in range(n + 1)) == sum(int(C[u, w]) for u in range(n) if u != w and A[u, w])
            vl, vl_w, ms, dmin, dmin_w, mtrans, fstrong = stats(n, A, H, C)
            rows.append((j, int(H), vl, ms, dmin, dmin_w, mtrans, fstrong, strong_landau(n, A)))
        vlmin = min(r[2] for r in rows)
        js_vl = [r[0] for r in rows if r[2] == vlmin]
        msmin = min(r[3] for r in rows)
        dm = min(r[4] for r in rows)
        js_d = [r[0] for r in rows if r[4] == dm]
        Tn = Fr(2 ** (n - 3) + 1, 2 ** (n - 2) + 1)
        print(f'n = {n}: H(D(n,1)) = {rows[0][1]} (2^(n-2)+1 = {2 ** (n - 2) + 1}); all strong: {all(r[8] for r in rows)}; '
              f'min vertex-lemma value {vlmin} = T_n value {Tn}: {vlmin == Tn}, attained at j = {js_vl}; '
              f'min majority-source value {msmin} ({float(msmin):.6f}); min D_in {dm} = {float(dm):.4f} at j = {js_d} (w = {[r[5] for r in rows if r[4] == dm]}); '
              f'#j with D_in <= 0 somewhere: {sum(1 for r in rows if r[4] <= 0)}; M transitive: {sum(r[6] for r in rows)}; '
              f'F not strong: {sum(not r[7] for r in rows)}  [{time.time() - t0:.1f} s]', flush=True)


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    print('Part 1: brute force over all permutations')
    part1(11, 5)
    part1(12, 6)
    part1(11, 6)
    part1(12, 7)
    print('Part 2: explicit DFS, every j')
    part2(nmax)
