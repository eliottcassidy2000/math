#!/usr/bin/env python3
"""Independent check (explicit enumeration of Hamiltonian paths, exact fractions) that the averaged form D_in(w) > 0
of HYP-9280 fails, while the vertex lemma, 'no majority source' and HYP-9280 itself hold (session opus-2026-10-14-S36, post-audit).

The 'double chain' D(n, j): transitive order on positions 0..n-1 (higher beats lower) with the two arcs 0 -> j and
j -> n-1 reversed.  It is strong for 1 <= j <= n-2.  At w = 0 (out-degree 1, out-neighbour j) the top vertex n-1 is
an in-neighbour that precedes w only rarely, which drags E[N] below d/2.
Part 1: D(11, 5) and D(12, 6) by pure-Python DFS: all p(u before 0), E[N] - d/2, best in-neighbour.
Part 2: for 6 <= n <= NMAX and every j, a numba DFS computes the full precedence matrix and reports, over the family,
        the minimum vertex-lemma value, the minimum majority-source value, the minimum D_in, M(T) transitive?, F strong?
"""
import sys, time
from fractions import Fraction
import numpy as np
import numba as nb


def build(n, arcs):
    beats = [[False] * n for _ in range(n)]
    for a in range(n):
        for b in range(a):
            beats[a][b] = True
    for lo, hi in arcs:
        beats[hi][lo] = False
        beats[lo][hi] = True
    return beats


def ham_paths(n, beats):
    paths = []

    def rec(path, used):
        if len(path) == n:
            paths.append(tuple(path))
            return
        v = path[-1]
        for x in range(n):
            if not used[x] and beats[v][x]:
                used[x] = True
                path.append(x)
                rec(path, used)
                path.pop()
                used[x] = False
    for s in range(n):
        used = [False] * n
        used[s] = True
        rec([s], used)
    return paths


def part1(n, j):
    beats = build(n, [(0, j), (j, n - 1)])
    P = ham_paths(n, beats)
    H = len(P)
    w = 0
    U = [u for u in range(n) if u != w and beats[u][w]]
    pos = [{v: k for k, v in enumerate(p)} for p in P]
    pu = {u: Fraction(sum(1 for q in pos if q[u] < q[w]), H) for u in U}
    EN = sum(pu.values())
    print(f'D({n}, {j}): H = {H}; w = 0 has out-neighbour {j} only, in-degree d = {len(U)}')
    print('   p(u before 0): ' + ', '.join(f'{u}: {pu[u]}' for u in U))
    print(f'   E[N] = {EN}, d/2 = {Fraction(len(U), 2)}, E[N] - d/2 = {EN - Fraction(len(U), 2)} (< 0: {EN < Fraction(len(U), 2)}); '
          f'best in-neighbour p = {max(pu.values())} (> 1/2: {max(pu.values()) > Fraction(1, 2)})')


@nb.njit(cache=True)
def enum_counts(n, adj):
    C = np.zeros((n, n), dtype=np.int64)
    H = 0
    path = np.zeros(n, dtype=np.int64)
    nxt = np.zeros(n, dtype=np.int64)
    for s in range(n):
        depth = 0
        path[0] = s
        used = 1 << s
        nxt[0] = 0
        while depth >= 0:
            if depth == n - 1:
                H += 1
                for a in range(n):
                    for b in range(a + 1, n):
                        C[path[a], path[b]] += 1
                used &= ~(1 << path[depth])
                depth -= 1
                continue
            v = path[depth]
            x = nxt[depth]
            while x < n and not (not (used >> x) & 1 and (adj[v] >> x) & 1):
                x += 1
            if x < n:
                nxt[depth] = x + 1
                depth += 1
                path[depth] = x
                used |= 1 << x
                nxt[depth] = 0
            else:
                if depth == 0:
                    break
                used &= ~(1 << path[depth])
                depth -= 1
    return H, C


def strong(n, M):
    def reach(fwd):
        seen = {0}
        st = [0]
        while st:
            v = st.pop()
            for x in range(n):
                if x not in seen and (M[v][x] if fwd else M[x][v]):
                    seen.add(x)
                    st.append(x)
        return len(seen) == n
    return reach(True) and reach(False)


def part2(n):
    res = []
    for j in range(1, n - 1):
        beats = build(n, [(0, j), (j, n - 1)])
        adj = np.array([sum(1 << b for b in range(n) if beats[a][b]) for a in range(n)], dtype=np.int64)
        H, C = enum_counts(n, adj)
        H = int(H)
        M = [[a != b and 2 * C[a, b] > H for b in range(n)] for a in range(n)]
        trans = sorted(sum(r) for r in M) == list(range(n))
        Fs = strong(n, [[beats[a][b] and M[a][b] for b in range(n)] for a in range(n)])
        vl = min(max(Fraction(int(C[u, w]), H) for u in range(n) if u != w and beats[u][w]) for w in range(n))
        ms = min(max(Fraction(int(C[x, w]), H) for x in range(n) if x != w) for w in range(n))
        dmin = min(sum(Fraction(2 * int(C[u, w]) - H, H) for u in range(n) if u != w and beats[u][w]) for w in range(n))
        res.append((j, H, vl, ms, dmin, trans, Fs))
    return res


if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    part1(11, 5)
    part1(12, 6)
    for n in range(6, NMAX + 1):
        t0 = time.time()
        res = part2(n)
        vl = min(r[2] for r in res); ms = min(r[3] for r in res); dm = min(r[4] for r in res)
        jd = min(res, key=lambda r: r[4])[0]
        print(f'n = {n}: double chains j = 1..{n - 2}: min vertex-lemma value {vl} ({float(vl):.6f}); min majority-source value '
              f'{float(ms):.6f}; min D_in {float(dm):.4f} (at j = {jd}); M(T) transitive: {sum(r[5] for r in res)}; '
              f'F(T) not strong: {sum(not r[6] for r in res)}  [{time.time() - t0:.1f} s]', flush=True)
