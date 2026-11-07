import itertools, random, sys
sys.path.insert(0, '.')
from sheaf_check import all_tournaments, random_tournament, ham_paths

def gate(adj):
    n = len(adj)
    a, b = n, n + 1
    G = [[False] * (n + 2) for _ in range(n + 2)]
    for i in range(n):
        for j in range(n):
            G[i][j] = adj[i][j]
        G[b][i] = True      # b -> T
        G[i][a] = True      # T -> a
    G[a][b] = True          # a -> b
    return G

rng = random.Random(5)
cases = [adj for n in range(1, 6) for adj in all_tournaments(n)] + [random_tournament(n, rng) for n in (6, 7) for _ in range(150)]
ok = bad = 0
for adj in cases:
    n = len(adj); V = list(range(n))
    lhs = ham_paths(gate(adj), list(range(n + 2)))
    s = 0
    for k in range(1, n + 1, 2):
        for S in itertools.combinations(V, k):
            rest = [v for v in V if v not in S]
            s += ham_paths(adj, S) * ham_paths(adj, rest)
    rhs = ham_paths(adj, V) + 2 * s
    if lhs == rhs: ok += 1
    else: bad += 1
print('H(Gamma T) = H(T) + 2 sum_{|S| odd} H(T[S]) H(T[V-S]):', ok, 'ok,', bad, 'bad')
