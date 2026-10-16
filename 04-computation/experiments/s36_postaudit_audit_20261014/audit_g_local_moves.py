"""Referee check (G): '(SS) at k = 1 cannot be realised by local moves' (session script
tournament_ss_k1_local_move_matching_20261014.py; claim: in a seeded random sample, a matching from X = {N <= 1} into
Y = {N >= d-1} along single-vertex relocations is impossible for 717 of 740 (T, w) pairs at n = 8; 504 with 2-blocks).

Independent code: Hamiltonian paths by brute force over all n! permutations (itertools), moves generated afresh,
maximum bipartite matching by Hopcroft-Karp (own implementation), and a CERTIFICATE for every deficient case: a Hall
violator S subset X with |N(S)| < |S|, obtained from the alternating-reachability set of the maximum matching (Koenig),
then re-checked directly from the edge lists (N(S) recomputed, |S| - |N(S)| = deficiency).
Part 1: 5 explicit (T, w) cases at n = 6 and 7 (numpy seed 20261015), printed in full.
Part 2: re-run of ONLY the n = 6 part of the session's random sample (random.seed(1), same generation order), counting
deficient (T, w) with this code, to compare with the session's 1203 / 654 / 99.  (n = 7, 8 not rerun, as instructed.)
"""
import itertools, random, time
from collections import deque
import numpy as np


def ham_paths_bruteforce(n, A):
    return [p for p in itertools.permutations(range(n)) if all(A[p[i]][p[i + 1]] for i in range(n - 1))]


def strong(n, A):
    R = [[bool(A[i][j]) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]:
                        R[i][j] = True
    return all(all(r) for r in R)


def relocations(seq, maxblock):
    out = set()
    n = len(seq)
    for L in range(1, maxblock + 1):
        for i in range(n - L + 1):
            blk = seq[i:i + L]
            rest = seq[:i] + seq[i + L:]
            for j in range(len(rest) + 1):
                s = rest[:j] + blk + rest[j:]
                if s != seq:
                    out.add(s)
    return out


def hopcroft_karp(X, adj):
    """X: list of left vertices; adj: dict left -> list of right vertices.  Returns (matchL, matchR)."""
    INF = float('inf')
    matchL = {x: None for x in X}
    matchR = {}
    while True:
        dist = {}
        q = deque()
        for x in X:
            if matchL[x] is None:
                dist[x] = 0
                q.append(x)
            else:
                dist[x] = INF
        found = False
        while q:
            x = q.popleft()
            for y in adj[x]:
                x2 = matchR.get(y)
                if x2 is None:
                    found = True
                elif dist[x2] == INF:
                    dist[x2] = dist[x] + 1
                    q.append(x2)
        if not found:
            break

        def dfs(x):
            for y in adj[x]:
                x2 = matchR.get(y)
                if x2 is None or (dist[x2] == dist[x] + 1 and dfs(x2)):
                    matchL[x] = y
                    matchR[y] = x
                    return True
            dist[x] = INF
            return False

        for x in X:
            if matchL[x] is None:
                dfs(x)
    return matchL, matchR


def hall_certificate(X, adj, matchL, matchR):
    """Alternating reachability from unmatched left vertices; returns (S, NS)."""
    S = set(x for x in X if matchL[x] is None)
    NS = set()
    q = deque(S)
    while q:
        x = q.popleft()
        for y in adj[x]:
            if y not in NS:
                NS.add(y)
                x2 = matchR.get(y)
                if x2 is not None and x2 not in S:
                    S.add(x2)
                    q.append(x2)
    return S, NS


def analyse(n, A, w, paths=None, verbose=False):
    if paths is None:
        paths = ham_paths_bruteforce(n, A)
    I = [u for u in range(n) if u != w and A[u][w]]
    d = len(I)
    Iset = set(I)
    Nof = {}
    for p in paths:
        k = p.index(w)
        Nof[p] = sum(1 for u in p[:k] if u in Iset)
    X = [p for p in paths if Nof[p] <= 1]
    Y = set(p for p in paths if Nof[p] >= d - 1)
    out = {'d': d, 'H': len(paths), 'X': len(X), 'Y': len(Y)}
    for blk in (1, 2):
        adj = {x: sorted(y for y in relocations(x, blk) if y in Y) for x in X}
        matchL, matchR = hopcroft_karp(X, adj)
        # validity of the matching
        used = [y for y in matchL.values() if y is not None]
        assert len(used) == len(set(used))
        assert all(matchL[x] is None or matchL[x] in adj[x] for x in X)
        nu = len(used)
        iso = sum(1 for x in X if not adj[x])
        S, NS = hall_certificate(X, adj, matchL, matchR)
        # direct re-check of the certificate from the edge lists
        NS_direct = set(y for x in S for y in adj[x])
        cert_ok = (NS_direct == NS) and (len(S) - len(NS) == len(X) - nu)
        out[blk] = {'nu': nu, 'deficient': nu < len(X), 'isolated': iso, 'S': len(S), 'NS': len(NS), 'cert_ok': cert_ok,
                    'S_N0': sum(1 for x in S if Nof[x] == 0), 'S_N1': sum(1 for x in S if Nof[x] == 1),
                    'example_isolated': next((x for x in X if not adj[x]), None)}
    return out


def outsets(n, A):
    return [[x for x in range(n) if A[v][x]] for v in range(n)]


def part1():
    print('PART 1: explicit (T, w) cases (numpy default_rng(20261015) random strong tournaments, first deficient ones)', flush=True)
    rng = np.random.default_rng(20261015)
    shown = 0
    shown2 = 0
    t0 = time.time()
    for n, want in ((6, 2), (7, 3)):
        got = 0
        while got < want:
            A = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    if rng.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            if not strong(n, A):
                continue
            paths = ham_paths_bruteforce(n, A)
            for w in range(n):
                if got >= want:
                    break
                d = sum(1 for u in range(n) if u != w and A[u][w])
                if d < 3:
                    continue
                r = analyse(n, A, w, paths)
                need_block_case = (n == 7 and got == want - 1 and shown2 == 0)
                if not r[1]['deficient']:
                    continue
                if need_block_case and not r[2]['deficient']:
                    continue
                got += 1
                shown += 1
                if r[2]['deficient']:
                    shown2 += 1
                print(f'  case {shown}: n = {n}, w = {w}, out-sets {outsets(n, A)}; d = {r["d"]}, H = {r["H"]}, |X| = |N<=1| = {r["X"]}, '
                      f'|Y| = |N>=d-1| = {r["Y"]}', flush=True)
                for blk in (1, 2):
                    q = r[blk]
                    name = 'single-vertex relocations' if blk == 1 else 'vertex or 2-block relocations'
                    print(f'     {name}: maximum matching {q["nu"]} < |X| = {r["X"]}: {q["deficient"]}; X-paths with no move into Y: {q["isolated"]}'
                          f'{" (e.g. " + str(q["example_isolated"]) + ")" if q["isolated"] else ""}; Hall violator |S| = {q["S"]} '
                          f'(N = 0: {q["S_N0"]}, N = 1: {q["S_N1"]}) with |N(S)| = {q["NS"]} (certificate re-checked: {q["cert_ok"]})', flush=True)
    print(f'  [{time.time() - t0:.1f} s]', flush=True)


def part2():
    print('PART 2: the session sample at n = 6 only (random.seed(1), same generation order), recounted with this code', flush=True)
    t0 = time.time()
    random.seed(1)
    n, trials = 6, 400
    pairs = 0; f1 = 0; f2 = 0; cert_bad = 0; iso1 = 0
    done = 0
    while done < trials:
        adjb = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5:
                    adjb[i] |= 1 << j
                else:
                    adjb[j] |= 1 << i
        A = [[(adjb[i] >> j) & 1 for j in range(n)] for i in range(n)]
        if not strong(n, A):
            continue
        paths = ham_paths_bruteforce(n, A)
        for w in range(n):
            d = sum(1 for u in range(n) if u != w and A[u][w])
            if d < 3:
                continue
            r = analyse(n, A, w, paths)
            pairs += 1
            f1 += r[1]['deficient']
            f2 += r[2]['deficient']
            iso1 += r[1]['isolated'] > 0 and r[1]['deficient']
            cert_bad += (not r[1]['cert_ok']) + (not r[2]['cert_ok'])
        done += 1
    print(f'  n = 6: (T, w) pairs with d >= 3: {pairs} (session 1203); single-vertex matching deficient: {f1} (session 654); '
          f'vertex-or-2-block deficient: {f2} (session 99); deficient single-vertex cases with an isolated X-path: {iso1}; '
          f'certificate failures: {cert_bad}  [{time.time() - t0:.1f} s]', flush=True)


if __name__ == '__main__':
    part1()
    part2()
