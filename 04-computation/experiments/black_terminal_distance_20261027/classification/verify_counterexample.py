#!/usr/bin/env python3
"""Independent brute-force verifier (pure Python, no shared code with the C++ checkers).

Usage: python3 verify_counterexample.py n CODE [CODE ...]
CODE = upper triangle row by row; character k is the pair (i,j), i<j, lexicographic; '1' means i->j.

For each tournament it prints: scores, strong?, regular?, lambda, the list of minimum out-cuts, the number of
Hamiltonian paths, the number of descent candidates (start in some minimum out-cut X, end outside X), how many
of them are black, and whether the 'stronger form' certificate exists.
"""
import sys
from itertools import combinations


def parse(n, code):
    assert len(code) == n * (n - 1) // 2, (len(code), n)
    adj = [[False] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[k] == '1':
                adj[i][j] = True
            elif code[k] == '0':
                adj[j][i] = True
            else:
                raise ValueError(code[k])
            k += 1
    return adj


def is_strong(n, adj):
    def reach(src, fwd):
        seen = {src}
        stack = [src]
        while stack:
            x = stack.pop()
            for y in range(n):
                if y not in seen and (adj[x][y] if fwd else adj[y][x]):
                    seen.add(y)
                    stack.append(y)
        return seen
    return len(reach(0, True)) == n and len(reach(0, False)) == n


def ham_paths(n, adj):
    paths = []
    path = []
    used = [False] * n

    def dfs():
        if len(path) == n:
            paths.append(tuple(path))
            return
        last = path[-1]
        for y in range(n):
            if not used[y] and adj[last][y]:
                used[y] = True
                path.append(y)
                dfs()
                path.pop()
                used[y] = False

    for s in range(n):
        used[s] = True
        path.append(s)
        dfs()
        path.pop()
        used[s] = False
    return paths


def is_blue(n, adj, Q):
    r = [None] * n
    for i in range(n):
        r[Q[i]] = Q[n - 1 - i]
    for u in range(n):
        for v in range(n):
            if u != v and adj[u][v] and not adj[r[v]][r[u]]:
                return False
    return True


def main():
    n = int(sys.argv[1])
    for code in sys.argv[2:]:
        adj = parse(n, code)
        s = [sum(adj[v]) for v in range(n)]
        strong = is_strong(n, adj)
        regular = all(2 * x == n - 1 for x in s)
        verts = list(range(n))
        best = None
        cuts = []
        for size in range(1, n):
            for X in combinations(verts, size):
                Xs = set(X)
                f = sum(1 for x in X for y in verts if y not in Xs and adj[x][y])
                if best is None or f < best:
                    best = f
                    cuts = [Xs]
                elif f == best:
                    cuts.append(Xs)
        P = ham_paths(n, adj)
        cand = [Q for Q in P if any(Q[0] in X and Q[-1] not in X for X in cuts)]
        black = [Q for Q in cand if not is_blue(n, adj, Q)]
        closed_cert = [Q for Q in cand if adj[Q[-1]][Q[0]] and s[Q[0]] + s[Q[-1]] != n - 1]
        print(f"code {code}")
        print(f"  scores {s}  strong={strong} regular={regular}  lambda={best}")
        print(f"  minimum out-cuts ({len(cuts)}): {[sorted(X) for X in cuts]}")
        print(f"  Hamiltonian paths: {len(P)}; descent candidates: {len(cand)}; black candidates: {len(black)}; "
              f"closed score-certified candidates: {len(closed_cert)}")
        for Q in cand[:12]:
            print(f"    candidate {Q}  blue={is_blue(n, adj, Q)}  closed={adj[Q[-1]][Q[0]]}")
        verdict = "KEY LEMMA HOLDS" if black else "KEY LEMMA FAILS (no black descent candidate)"
        print(f"  verdict: {verdict}")


if __name__ == '__main__':
    main()
