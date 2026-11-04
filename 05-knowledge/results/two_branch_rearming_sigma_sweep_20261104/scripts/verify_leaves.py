import json
d = json.load(open('direct_children_receipt.json'))
def adj_from_rows(rows):
    n = len(rows); return [[1 if rows[i] >> j & 1 else 0 for j in range(n)] for i in range(n)]
def hcount(adj, verts, start=None, end=None):
    n = len(verts); idx = {v: i for i, v in enumerate(verts)}
    f = [[0] * n for _ in range(1 << n)]
    for i, v in enumerate(verts):
        if start is None or v == start: f[1 << i][i] = 1
    for S in range(1, 1 << n):
        row = f[S]
        for i in range(n):
            c = row[i]
            if not c: continue
            vi = verts[i]
            for j in range(n):
                if adj[vi][verts[j]] and not (S >> j & 1): f[S | 1 << j][j] += c
    full = (1 << n) - 1
    return sum(f[full][i] for i, v in enumerate(verts) if end is None or v == end)
def strong(adj):
    n = len(adj)
    def reach(A):
        seen = {0}; st = [0]
        while st:
            u = st.pop()
            for v in range(n):
                if A[u][v] and v not in seen: seen.add(v); st.append(v)
        return len(seen) == n
    return reach(adj) and reach([[adj[j][i] for j in range(n)] for i in range(n)])
root = d['root_index']
for leaf in d['four_leaves']:
    adj = adj_from_rows(leaf['rows']); n = len(adj); V = list(range(n))
    H = hcount(adj, V)
    ready = []
    for a in range(n):
        if a == root: continue
        outs = [v for v in range(n) if adj[a][v]]
        if outs == [root]:
            e = hcount(adj, [v for v in V if v != a], end=root)
            ready.append((a, e))
    print('leaf H %d (receipt %d): strong %s; sole-exit auxiliaries into root %d with e_root(T-a): %s' % (H, leaf['H'], strong(adj), root, ready))
