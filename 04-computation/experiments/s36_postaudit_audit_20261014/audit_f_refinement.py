"""Referee check (F): the per-in-neighbour refinement of (SS) at k = 1.

For a vertex w of in-degree d >= 3 and an in-neighbour u:
  a_u = # Hamiltonian paths in which u is the ONLY in-neighbour of w before w,
  b_u = # Hamiltonian paths in which u is the ONLY in-neighbour of w after w.
Claimed (session output): a_u > b_u in 720 of 2280 triples (T, w, u) at n = 5, 82320 of 222720 at n = 6, 11,734,800 of
28,413,840 at n = 7 (labelled strong tournaments).  Method: explicit DFS over Hamiltonian paths with per-path position
bookkeeping (no subset DP), all 2^C(n,2) labelled tournaments, strongness by the Landau score test.  Sanity: for every
(T, w), sum_u a_u = #(N = 1) and sum_u b_u = #(N = d - 1) (N-distribution from pa_core's DFS).
Usage: python audit_f_refinement.py NMAX
"""
import sys, time
import numpy as np
import numba as nb
from pa_core import strong_landau, enumerate_paths_dfs


@nb.njit(cache=False)
def ab_counts(n, A, a, b):
    """a[w, u], b[w, u] as above (only for w with in-degree >= 3); returns H."""
    outm = np.zeros(n, dtype=np.int64)
    indeg = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if A[i, j]:
                outm[i] |= 1 << j
                indeg[j] += 1
    for i in range(n):
        for j in range(n):
            a[i, j] = 0
            b[i, j] = 0
    path = np.zeros(n, dtype=np.int64)
    cand = np.zeros(n, dtype=np.int64)
    pos = np.zeros(n, dtype=np.int64)
    H = 0
    for s in range(n):
        path[0] = s
        used = 1 << s
        depth = 0
        cand[0] = outm[s] & ~used
        while depth >= 0:
            if depth == n - 1:
                H += 1
                for t in range(n):
                    pos[path[t]] = t
                for w in range(n):
                    d = indeg[w]
                    if d < 3:
                        continue
                    nb_ = 0; ub = -1; ua = -1; na = 0
                    for u in range(n):
                        if u != w and A[u, w]:
                            if pos[u] < pos[w]:
                                nb_ += 1
                                ub = u
                            else:
                                na += 1
                                ua = u
                    if nb_ == 1:
                        a[w, ub] += 1
                    if na == 1:
                        b[w, ua] += 1
                used &= ~(1 << path[depth])
                depth -= 1
                continue
            c = cand[depth]
            if c == 0:
                used &= ~(1 << path[depth])
                depth -= 1
                continue
            lb = c & (-c)
            cand[depth] = c & ~lb
            v = 0
            while (lb >> v) != 1:
                v += 1
            depth += 1
            path[depth] = v
            used |= lb
            cand[depth] = outm[v] & ~used
    return H


@nb.njit(parallel=True, cache=False)
def sweep(n, pu, pv, nchunk, sanity):
    E = pu.shape[0]
    total = 1 << E
    res = np.zeros((nchunk, 6), dtype=np.int64)   # triples, a > b, a == b, a < b, strong count, sanity failures
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        A = np.zeros((n, n), dtype=np.int8)
        a = np.zeros((n, n), dtype=np.int64); b = np.zeros((n, n), dtype=np.int64)
        C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for i in range(n):
                for j in range(n):
                    A[i, j] = 0
            for e in range(E):
                if (code >> e) & 1:
                    A[pu[e], pv[e]] = 1
                else:
                    A[pv[e], pu[e]] = 1
            if not strong_landau(n, A):
                continue
            res[c, 4] += 1
            H = ab_counts(n, A, a, b)
            if sanity:
                H2 = enumerate_paths_dfs(n, A, C, Nd, L)
                bad = H2 != H
            for w in range(n):
                d = 0
                for u in range(n):
                    if u != w and A[u, w]:
                        d += 1
                if d < 3:
                    continue
                sa = 0; sb = 0
                for u in range(n):
                    if u != w and A[u, w]:
                        res[c, 0] += 1
                        sa += a[w, u]; sb += b[w, u]
                        if a[w, u] > b[w, u]:
                            res[c, 1] += 1
                        elif a[w, u] == b[w, u]:
                            res[c, 2] += 1
                        else:
                            res[c, 3] += 1
                if sanity:
                    if sa != Nd[w, 1] or sb != Nd[w, d - 1]:
                        bad = True
            if sanity and bad:
                res[c, 5] += 1
    tot = np.zeros(6, dtype=np.int64)
    for c in range(nchunk):
        for i in range(6):
            tot[i] += res[c, i]
    return tot


if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    claimed = {5: (2280, 720), 6: (222720, 82320), 7: (28413840, 11734800)}
    for n in range(4, NMAX + 1):
        t0 = time.time()
        prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in prs], dtype=np.int64); pv = np.array([p[1] for p in prs], dtype=np.int64)
        r = sweep(n, pu, pv, 512, True)
        cl = claimed.get(n)
        print(f'n = {n}: {1 << len(prs)} labelled tournaments, {int(r[4])} strong; triples (T, w, u) with d >= 3: {int(r[0])}; '
              f'a_u > b_u: {int(r[1])} ({100 * r[1] / max(1, r[0]):.2f}%), a_u = b_u: {int(r[2])}, a_u < b_u: {int(r[3])}; '
              f'sanity failures (sum a = #(N=1), sum b = #(N=d-1), H): {int(r[5])}; '
              f'session claim {cl}: {"MATCH" if cl and (int(r[0]), int(r[1])) == cl else ("-" if not cl else "MISMATCH")}  [{time.time() - t0:.1f} s]', flush=True)
