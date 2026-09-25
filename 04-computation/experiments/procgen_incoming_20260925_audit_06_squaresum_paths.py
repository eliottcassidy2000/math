#!/usr/bin/env python3
"""Independent audit of decoder_prime_square / square-sum correction:
Hamiltonian path counts of Q_N (up to reversal) for N=2..27 by plain DFS over
all start vertices (each undirected path counted twice, then halved), with a
simple connectivity prune; check every Q_25 path uses 11-25-24; check the Q_24
endpoint-aware lemma: in Q_24 no Hamiltonian path exists (0), and the claimed
forced 11-cycle 1-8-17-19-6-10-15-21-4-12-24-1 lies in Q_24.
"""
import sys
from math import isqrt


def adjacency(N):
    adj = {v: [] for v in range(1, N + 1)}
    for x in range(1, N + 1):
        for y in range(x + 1, N + 1):
            s = x + y
            if isqrt(s) ** 2 == s:
                adj[x].append(y)
                adj[y].append(x)
    return adj


def count_paths(N, record=False):
    adj = adjacency(N)
    total = 0
    uses_ear = True
    found = []
    visited = [False] * (N + 1)
    path = []

    def reachable_ok(cur):
        # all unvisited vertices must be reachable from cur through unvisited vertices
        seen = {cur}
        stack = [cur]
        cnt = 0
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if not visited[w] and w not in seen:
                    seen.add(w)
                    cnt += 1
                    stack.append(w)
        return cnt == N - len(path)

    def dfs(v):
        nonlocal total, uses_ear
        if len(path) == N:
            total += 1
            if record and len(found) < 3:
                found.append(list(path))
            if N == 25:
                i = path.index(25)
                nb = set()
                if i > 0:
                    nb.add(path[i - 1])
                if i < N - 1:
                    nb.add(path[i + 1])
                if nb != {11, 24}:
                    uses_ear = False
            return
        if not reachable_ok(v):
            return
        for w in adj[v]:
            if not visited[w]:
                visited[w] = True
                path.append(w)
                dfs(w)
                path.pop()
                visited[w] = False

    for s in range(1, N + 1):
        visited[s] = True
        path.append(s)
        dfs(s)
        path.pop()
        visited[s] = False
    if N == 1:
        return 1, uses_ear, found
    assert total % 2 == 0
    return total // 2, uses_ear, found


def main():
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 27
    res = []
    for N in range(2, Nmax + 1):
        c, ear, f = count_paths(N, record=(N in (15, 23, 25)))
        res.append((N, c))
        extra = ""
        if N == 25:
            extra = " every path uses 11-25-24: %s" % ear
        print("N=%d Hamiltonian paths up to reversal: %d%s" % (N, c, extra), flush=True)
    adj = adjacency(24)
    cyc = [1, 8, 17, 19, 6, 10, 15, 21, 4, 12, 24, 1]
    print("forced 11-cycle edges in Q24:", all(b in adj[a] for a, b in zip(cyc, cyc[1:])))
    print("degrees in Q24 of 18,9,11,22,2,20,5,4:", {v: sorted(adj[v]) for v in (18, 9, 11, 22, 2, 20, 5, 4, 8, 17, 19, 10, 21, 24)})


if __name__ == "__main__":
    main()
