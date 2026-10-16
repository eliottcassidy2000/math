"""Shared core for the S36 post-audit referee checks (independent of the session's subset DP).

Conventions: a tournament on n vertices is an n x n int8 matrix A with A[i, j] = 1 iff i beats j (i -> j).
'transitive + reversed arcs': positions 0..n-1, higher position beats lower; an arc (lo, hi) in the list is reversed,
so that lo -> hi.

Methods used here (all different from the session's subset DP):
  * enumerate_paths_dfs: explicit-stack depth-first enumeration of every Hamiltonian path, accumulating per path
      C[x, y]  (# paths with x before y), Nd[w, j] (# paths with exactly j in-neighbours of w before w),
      L[w, u]  (# paths with u immediately before w).
  * brute_force_perms: lexicographic enumeration of ALL n! permutations, each tested arc by arc (no pruning
      beyond the first failing arc), parallel over the first vertex.
  * strong_landau: score-sequence criterion (sorted scores s_1 <= ... <= s_n; strong iff sum_{i<=k} s_i > C(k,2)
      for all k < n).   strong_warshall: boolean transitive closure (Floyd-Warshall) of an arbitrary digraph.
  * transitive_by_triangles: a tournament is transitive iff it has no directed 3-cycle.
No numba on-disk cache is used (cache=False), so no cache files are written.
"""
import numpy as np
import numba as nb


def transitive_plus(n, arcs):
    A = np.zeros((n, n), dtype=np.int8)
    for a in range(n):
        for b in range(a):
            A[a, b] = 1
    for lo, hi in arcs:
        A[hi, lo] = 0
        A[lo, hi] = 1
    return A


@nb.njit(cache=False)
def strong_landau(n, A):
    s = np.zeros(n, dtype=np.int64)
    for i in range(n):
        t = 0
        for j in range(n):
            t += A[i, j]
        s[i] = t
    s = np.sort(s)
    acc = 0
    for k in range(1, n):
        acc += s[k - 1]
        if acc == k * (k - 1) // 2:
            return False
    return True


@nb.njit(cache=False)
def strong_warshall(n, M):
    R = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(n):
            R[i, j] = 1 if (M[i, j] != 0 or i == j) else 0
    for k in range(n):
        for i in range(n):
            if R[i, k]:
                for j in range(n):
                    if R[k, j]:
                        R[i, j] = 1
    for i in range(n):
        for j in range(n):
            if R[i, j] == 0:
                return False
    return True


@nb.njit(cache=False)
def transitive_by_triangles(n, M):
    for a in range(n):
        for b in range(n):
            if b == a or M[a, b] == 0:
                continue
            for c in range(n):
                if c == a or c == b:
                    continue
                if M[b, c] and M[c, a]:
                    return False
    return True


@nb.njit(cache=False)
def _popc(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@nb.njit(cache=False)
def enumerate_paths_dfs(n, A, C, Nd, L):
    """Explicit DFS over all Hamiltonian paths.  Fills C (n x n), Nd (n x (n+1)), L (n x n); returns H."""
    outm = np.zeros(n, dtype=np.int64)
    inm = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if A[i, j]:
                outm[i] |= 1 << j
                inm[j] |= 1 << i
    for i in range(n):
        for j in range(n):
            C[i, j] = 0
            L[i, j] = 0
        for j in range(n + 1):
            Nd[i, j] = 0
    path = np.zeros(n, dtype=np.int64)
    cand = np.zeros(n, dtype=np.int64)
    H = 0
    for s in range(n):
        path[0] = s
        used = 1 << s
        depth = 0
        cand[0] = outm[s] & ~used
        while depth >= 0:
            if depth == n - 1:
                H += 1
                pm = 0
                for t in range(n):
                    y = path[t]
                    for t2 in range(t):
                        C[path[t2], y] += 1
                    Nd[y, _popc(inm[y] & pm)] += 1
                    if t >= 1:
                        L[y, path[t - 1]] += 1
                    pm |= 1 << y
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


@nb.njit(cache=False)
def enumerate_paths_dfs_first(n, A, C, Nd, L, first):
    """As enumerate_paths_dfs, and first[v] = # Hamiltonian paths starting at v."""
    for v in range(n):
        first[v] = 0
    outm = np.zeros(n, dtype=np.int64)
    inm = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if A[i, j]:
                outm[i] |= 1 << j
                inm[j] |= 1 << i
    for i in range(n):
        for j in range(n):
            C[i, j] = 0
            L[i, j] = 0
        for j in range(n + 1):
            Nd[i, j] = 0
    path = np.zeros(n, dtype=np.int64)
    cand = np.zeros(n, dtype=np.int64)
    H = 0
    for s in range(n):
        path[0] = s
        used = 1 << s
        depth = 0
        cand[0] = outm[s] & ~used
        while depth >= 0:
            if depth == n - 1:
                H += 1
                first[s] += 1
                pm = 0
                for t in range(n):
                    y = path[t]
                    for t2 in range(t):
                        C[path[t2], y] += 1
                    Nd[y, _popc(inm[y] & pm)] += 1
                    if t >= 1:
                        L[y, path[t - 1]] += 1
                    pm |= 1 << y
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


@nb.njit(cache=False)
def count_paths_dfs(n, A):
    """Number of Hamiltonian paths only (DFS)."""
    if n == 0:
        return 1
    outm = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if A[i, j]:
                outm[i] |= 1 << j
    path = np.zeros(n, dtype=np.int64)
    cand = np.zeros(n, dtype=np.int64)
    H = 0
    for s in range(n):
        path[0] = s
        used = 1 << s
        depth = 0
        cand[0] = outm[s] & ~used
        while depth >= 0:
            if depth == n - 1:
                H += 1
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


@nb.njit(cache=False)
def _next_perm(a, lo):
    """Lexicographic next permutation of a[lo:]; returns False after the last one."""
    n = a.shape[0]
    i = n - 2
    while i >= lo and a[i] >= a[i + 1]:
        i -= 1
    if i < lo:
        return False
    j = n - 1
    while a[j] <= a[i]:
        j -= 1
    t = a[i]; a[i] = a[j]; a[j] = t
    l = i + 1; r = n - 1
    while l < r:
        t = a[l]; a[l] = a[r]; a[r] = t
        l += 1; r -= 1
    return True


@nb.njit(parallel=True, cache=False)
def brute_force_perms(n, A):
    """All n! permutations, parallel over the first vertex.  Returns H, C (x before y), Nd, L (as in DFS)."""
    Cs = np.zeros((n, n, n), dtype=np.int64)
    Nds = np.zeros((n, n, n + 1), dtype=np.int64)
    Ls = np.zeros((n, n, n), dtype=np.int64)
    Hs = np.zeros(n, dtype=np.int64)
    for f in nb.prange(n):
        a = np.zeros(n, dtype=np.int64)
        a[0] = f
        k = 1
        for v in range(n):
            if v != f:
                a[k] = v
                k += 1
        pos = np.zeros(n, dtype=np.int64)
        while True:
            ok = True
            for t in range(n - 1):
                if A[a[t], a[t + 1]] == 0:
                    ok = False
                    break
            if ok:
                Hs[f] += 1
                for t in range(n):
                    pos[a[t]] = t
                for x in range(n):
                    for y in range(n):
                        if x != y and pos[x] < pos[y]:
                            Cs[f, x, y] += 1
                for w in range(n):
                    cnt = 0
                    for u in range(n):
                        if u != w and A[u, w] and pos[u] < pos[w]:
                            cnt += 1
                    Nds[f, w, cnt] += 1
                for t in range(1, n):
                    Ls[f, a[t], a[t - 1]] += 1
            if not _next_perm(a, 1):
                break
    H = 0
    C = np.zeros((n, n), dtype=np.int64)
    Nd = np.zeros((n, n + 1), dtype=np.int64)
    L = np.zeros((n, n), dtype=np.int64)
    for f in range(n):
        H += Hs[f]
        C += Cs[f]
        Nd += Nds[f]
        L += Ls[f]
    return H, C, Nd, L


def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A
