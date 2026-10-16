"""Referee check (E): (SS) at k = 1, Pr[N >= d-1] >= Pr[N <= 1], at every vertex w of in-degree d >= 3.

N = number of in-neighbours of w before w in a uniform random Hamiltonian path.  Counts by explicit DFS (pa_core),
not the session's subset DP.  Parts:
  algebra   the claim '(SS) at k = 0 and k = 1 give E[N] >= 3/2 for d = 3 and E[N] >= 2 for d = 4': exact LP over all
            distributions on {0..d} (vertex enumeration with Fractions), also d = 5 (expected NOT to reach 5/2).
  random    uniformly random strong tournaments, n = 8..13 (unplanted), and tournaments with a planted vertex of
            out-degree 1 or 2 (where all known (SS) failures live), n = 8..11; every vertex with d >= 3.
  climb     greedy local search (single-arc flips keeping strongness) minimising the smallest relative k = 1 slack
            [#(N >= d-1) - #(N <= 1)] / H over vertices with d >= 3, n = 8, 9, 10.
For every (T, w) also: (SS) at k = 0, and the consequences E[N] >= 3/2 (d = 3), E[N] >= 2 (d = 4) checked directly.
"""
import sys, time, itertools
from fractions import Fraction as Fr
import numpy as np
import numba as nb
from pa_core import enumerate_paths_dfs, strong_landau, random_tournament


def lp_min_EN(d, ks):
    """min E[N] over probability vectors P on {0..d} with Pr[N >= d-k] >= Pr[N <= k] for k in ks.
    Exact: enumerate all vertices of the polytope (choose d+1 tight constraints among the inequalities), Fractions."""
    m = d + 1
    rows = []  # inequality rows a.P >= 0
    for k in ks:
        a = [Fr(0)] * m
        for j in range(d - k, d + 1):
            a[j] += 1
        for j in range(0, k + 1):
            a[j] -= 1
        rows.append(a)
    for j in range(m):
        a = [Fr(0)] * m; a[j] = Fr(1); rows.append(a)
    best = None; arg = None
    for S in itertools.combinations(range(len(rows)), m - 1):
        M = [rows[i][:] + [Fr(0)] for i in S] + [[Fr(1)] * m + [Fr(1)]]
        # Gaussian elimination
        ok = True
        for c in range(m):
            piv = next((r for r in range(c, m) if M[r][c] != 0), None)
            if piv is None:
                ok = False; break
            M[c], M[piv] = M[piv], M[c]
            for r in range(m):
                if r != c and M[r][c] != 0:
                    f = M[r][c] / M[c][c]
                    M[r] = [x - f * y for x, y in zip(M[r], M[c])]
        if not ok:
            continue
        P = [M[i][m] / M[i][i] for i in range(m)]
        if all(sum(a[j] * P[j] for j in range(m)) >= 0 for a in rows):
            EN = sum(j * P[j] for j in range(m))
            if best is None or EN < best:
                best, arg = EN, P
    return best, arg


@nb.njit(cache=False)
def ss_eval(n, A, Nd, H, res, mins):
    """res: [pairs d>=3, k1 fail, k0 fail, d=3 pairs, d=3 E[N] < 3/2, d=4 pairs, d=4 E[N] < 2, d=3 E[N] < 1 + Pr[N>=2]];
    mins[0] = min k1 slack / H;  mins[1] = min over d=3 of E[N] - 3/2;  mins[2] = min over d=4 of E[N] - 2."""
    for w in range(n):
        d = 0
        for u in range(n):
            if u != w and A[u, w]:
                d += 1
        if d < 3:
            continue
        res[0] += 1
        hi = Nd[w, d] + Nd[w, d - 1]
        lo = Nd[w, 0] + Nd[w, 1]
        if hi < lo:
            res[1] += 1
        if Nd[w, d] < Nd[w, 0]:
            res[2] += 1
        v = (hi - lo) / H
        if v < mins[0]:
            mins[0] = v
        s1 = 0
        for j in range(d + 1):
            s1 += j * Nd[w, j]
        if d == 3:
            res[3] += 1
            if 2 * s1 < 3 * H:
                res[4] += 1
            ge2 = Nd[w, 2] + Nd[w, 3]
            if s1 < H + ge2:
                res[7] += 1
            v = s1 / H - 1.5
            if v < mins[1]:
                mins[1] = v
        if d == 4:
            res[5] += 1
            if s1 < 2 * H:
                res[6] += 1
            v = s1 / H - 2.0
            if v < mins[2]:
                mins[2] = v


def show(tag, r, mins, t0):
    print(f'  {tag}: (T, w) with d >= 3: {int(r[0])}; (SS) k = 1 failures {int(r[1])}; k = 0 failures {int(r[2])}; '
          f'min k=1 slack/H = {mins[0]:.6f}; d = 3: {int(r[3])} pairs, E[N] < 3/2 in {int(r[4])}, E[N] < 1 + Pr[N>=2] in {int(r[7])}, '
          f'min E[N]-3/2 = {mins[1]:.4f}; d = 4: {int(r[5])} pairs, E[N] < 2 in {int(r[6])}, min E[N]-2 = {mins[2]:.4f}  '
          f'[{time.time() - t0:.1f} s]', flush=True)


def plant(A, n, w, outs):
    for x in range(n):
        if x != w:
            A[w, x] = 0; A[x, w] = 1
    for o in outs:
        A[w, o] = 1; A[o, w] = 0


def part_algebra():
    print('PART algebra: exact LP min of E[N] under (SS) constraints', flush=True)
    for d, ks in [(3, [0, 1]), (4, [0, 1]), (5, [0, 1]), (5, [0, 1, 2]), (6, [0, 1, 2]), (3, [0]), (4, [0])]:
        best, arg = lp_min_EN(d, ks)
        print(f'  d = {d}, (SS) at k in {ks}: min E[N] = {best} (d/2 = {Fr(d, 2)}; >= d/2: {best >= Fr(d, 2)}); a minimiser P = {[str(x) for x in arg]}', flush=True)


def part_random(rng):
    print('PART random: unplanted random strong tournaments', flush=True)
    for n, trials in [(8, 400), (9, 400), (10, 300), (11, 120), (12, 40), (13, 12)]:
        t0 = time.time(); r = np.zeros(8, dtype=np.int64); mins = np.full(3, 9.0); got = 0
        while got < trials:
            A = random_tournament(n, rng)
            if not strong_landau(n, A):
                continue
            got += 1
            C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
            H = enumerate_paths_dfs(n, A, C, Nd, L)
            ss_eval(n, A, Nd, H, r, mins)
        show(f'n = {n}, {trials} random strong', r, mins, t0)
    print('PART random: planted vertex of out-degree 1 or 2 (rest random), strong only', flush=True)
    for n, trials in [(8, 600), (9, 500), (10, 300), (11, 100)]:
        for od in (1, 2):
            t0 = time.time(); r = np.zeros(8, dtype=np.int64); mins = np.full(3, 9.0); got = 0
            while got < trials:
                A = random_tournament(n, rng)
                w = int(rng.integers(n))
                others = [x for x in range(n) if x != w]
                outs = list(rng.choice(others, size=od, replace=False))
                plant(A, n, w, outs)
                if not strong_landau(n, A):
                    continue
                got += 1
                C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
                H = enumerate_paths_dfs(n, A, C, Nd, L)
                ss_eval(n, A, Nd, H, r, mins)
            show(f'n = {n}, {trials} with a planted out-degree-{od} vertex', r, mins, t0)


def objective(n, A):
    C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
    H = enumerate_paths_dfs(n, A, C, Nd, L)
    r = np.zeros(8, dtype=np.int64); mins = np.full(3, 9.0)
    ss_eval(n, A, Nd, H, r, mins)
    return mins[0], int(r[1]), int(H)


def part_climb(rng):
    print('PART climb: greedy single-arc-flip search minimising the relative (SS) k = 1 slack', flush=True)
    for n, restarts, steps in [(8, 6, 1500), (9, 5, 1200), (10, 4, 700)]:
        t0 = time.time(); overall = 9.0; fails = 0; bestA = None; bestH = None
        for rs in range(restarts):
            while True:
                A = random_tournament(n, rng)
                if strong_landau(n, A):
                    break
            cur, f, H = objective(n, A)
            for st in range(steps):
                i, j = rng.choice(n, size=2, replace=False)
                B = A.copy(); B[i, j], B[j, i] = B[j, i], B[i, j]
                if not strong_landau(n, B):
                    continue
                val, f, H = objective(n, B)
                fails += f
                if val <= cur:
                    A, cur = B, val
                    if val < overall:
                        overall, bestA, bestH = val, A.copy(), H
        print(f'  n = {n}: {restarts} restarts x {steps} steps: min relative k=1 slack found {overall:.6f} (H = {bestH}; '
              f'1/H(T_n) = {1 / (2 ** (n - 2) + 1):.6f}); (SS) k = 1 failures met during search: {fails}; '
              f'best tournament out-sets {[sorted(int(x) for x in range(n) if bestA[v, x]) for v in range(n)]}  [{time.time() - t0:.1f} s]', flush=True)


def part_k2deg():
    """In-degree d of the vertices where (SS) fails at k = 2 in near-transitive sweeps (does the k <= 2 route to the
    vertex lemma for d <= 6 survive?)."""
    print('PART k2deg: in-degree of (SS) k = 2 failures in near-transitive sweeps (transitive + k reversed arcs)', flush=True)
    for n, k in [(10, 2), (11, 2), (12, 2), (10, 3), (11, 3), (10, 4)]:
        t0 = time.time()
        cand = [(lo, hi) for lo in range(n) for hi in range(lo + 2, n)]
        degs = {}
        first = None
        for arcs in itertools.combinations(cand, k):
            A = np.zeros((n, n), dtype=np.int8)
            for a in range(n):
                for b in range(a):
                    A[a, b] = 1
            for lo, hi in arcs:
                A[hi, lo] = 0; A[lo, hi] = 1
            if not strong_landau(n, A):
                continue
            C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
            H = int(enumerate_paths_dfs(n, A, C, Nd, L))
            for w in range(n):
                d = int(sum(A[u, w] for u in range(n) if u != w))
                if d < 5:
                    continue
                hi_ = int(Nd[w, d] + Nd[w, d - 1] + Nd[w, d - 2]); lo_ = int(Nd[w, 0] + Nd[w, 1] + Nd[w, 2])
                if hi_ < lo_:
                    degs[d] = degs.get(d, 0) + 1
                    if first is None:
                        first = (arcs, w, d, Fr(hi_ - lo_, H))
        print(f'  n = {n}, k = {k}: (SS) k = 2 failures by in-degree d: {dict(sorted(degs.items()))}; first: {first}  [{time.time() - t0:.1f} s]', flush=True)


if __name__ == '__main__':
    part = sys.argv[1]
    if part == 'k2deg':
        part_k2deg()
    if part == 'algebra':
        part_algebra()
    elif part == 'random':
        part_random(np.random.default_rng(20261015))
    elif part == 'climb':
        part_climb(np.random.default_rng(15102026))
