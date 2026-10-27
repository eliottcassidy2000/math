"""Constructive B8: from the transitive tournament, repeatedly move along a Hamiltonian path from a maximum-score to a
minimum-score vertex (THM-4815 B7); after exactly (n^2-1)/8 moves the tournament is regular.  Randomised choices;
reports which regular classes are reached (identified by the number H of Hamiltonian paths) and checks R drops by one
per move.  Usage: python3 greedy_regular.py n trials seed"""
import sys, random
n, trials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); random.seed(seed)
def ham_paths_between(out, a, b):
    # one random Hamiltonian path from a to b (DP over subsets with random tie-breaking), or None
    full = (1 << n) - 1
    reach = {}  # (S, v) reachable from a covering S ending at v
    import functools
    @functools.lru_cache(maxsize=None)
    def ok(S, v):   # exists path from v covering complement-of-S ... we do forward DP instead
        return None
    # forward DP: dp[S] = set of ends
    dp = [0] * (1 << n); dp[1 << a] = 1 << a
    for S in range(1 << n):
        e = dp[S]
        if not e: continue
        for v in range(n):
            if e >> v & 1:
                nx = out[v] & ~S
                while nx:
                    w = (nx & -nx).bit_length() - 1; nx &= nx - 1
                    dp[S | 1 << w] |= 1 << w
    if not (dp[full] >> b & 1): return None
    # backtrack a random path ending at b
    path = [b]; S = full; v = b
    while S != (1 << a):
        S2 = S & ~(1 << v)
        cands = [u for u in range(n) if (dp[S2] >> u & 1) and (out[u] >> v & 1)]
        u = random.choice(cands); path.append(u); S = S2; v = u
    return path[::-1]
def move(out, Q):
    full = (1 << n) - 1
    r = [(~out[u]) & full & ~(1 << u) for u in range(n)]
    for i in range(n - 1):
        u, v = Q[i], Q[i + 1]; r[u] |= 1 << v; r[v] &= ~(1 << u)
    return r
def count_paths(out):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c: continue
            nx = out[v] & ~S
            while nx:
                w = (nx & -nx).bit_length() - 1; nx &= nx - 1
                dp[S | 1 << w][w] += c
    return sum(dp[full])
m = (n - 1) // 2
seen = {}
for t in range(trials):
    out = [sum(1 << j for j in range(i + 1, n)) for i in range(n)]   # transitive
    steps = 0
    while True:
        s = [bin(o).count("1") for o in out]
        R = sum(abs(x - m) for x in s) // 2
        if R == 0: break
        mx, mn = max(s), min(s)
        pairs = [(a, b) for a in range(n) for b in range(n) if s[a] == mx and s[b] == mn]
        random.shuffle(pairs)
        Q = None
        for a, b in pairs:
            Q = ham_paths_between(out, a, b)
            if Q: break
        assert Q, "B7 violated"
        out2 = move(out, Q); s2 = [bin(o).count("1") for o in out2]
        assert sum(abs(x - m) for x in s2) // 2 == R - 1, "R did not drop"
        out = out2; steps += 1
    H = count_paths(out)
    seen[H] = seen.get(H, 0) + 1
    assert steps == (n * n - 1) // 8
print(f"n={n}: {trials} greedy descents from the transitive tournament, each exactly (n^2-1)/8 = {(n*n-1)//8} moves; regular classes reached, by Hamiltonian-path count H: {dict(sorted(seen.items()))}")
