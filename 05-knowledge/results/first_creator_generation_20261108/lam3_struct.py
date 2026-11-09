"""Structure of the 3-arc-strong 11-tournaments with H in S(10) (non-first-creators at order 11)."""
import itertools
S9 = set(int(x) for x in open('S9.txt').read().split()); S10 = set(int(x) for x in open('S10_klein.txt').read().split())
def decode(code, n):
    adj = [[0]*n for _ in range(n)]; b = 0
    for i in range(n):
        for j in range(i+1, n):
            if code >> b & 1: adj[i][j] = 1
            else: adj[j][i] = 1
            b += 1
    return adj
def ham(adj, verts):
    verts = list(verts); m = len(verts)
    f = [[0]*m for _ in range(1 << m)]
    for i in range(m): f[1 << i][i] = 1
    for S in range(1, 1 << m):
        for i in range(m):
            x = f[S][i]
            if not x: continue
            for j in range(m):
                if not S >> j & 1 and adj[verts[i]][verts[j]]: f[S | 1 << j][j] += x
    return sum(f[(1 << m) - 1])
def modules(adj, n):
    full = (1 << n) - 1; res = set()
    for a in range(n):
        for b in range(a+1, n):
            M = (1 << a) | (1 << b); ch = True
            while ch:
                ch = False
                for z in range(n):
                    if not M >> z & 1:
                        outs = sum(adj[z][w] for w in range(n) if M >> w & 1); k = bin(M).count('1')
                        if 0 < outs < k: M |= 1 << z; ch = True
            if M != full: res.add(M)
    return res
seen = set()
for line in open('lambda3_witnesses_order11.txt'):
    _, code, x, h, st = line.split(); code = int(code); x = int(x); h = int(h)
    if st != 'inS10' or h in seen: continue
    seen.add(h)
    n = 10; T = decode(code, n)
    adj = [row + [0] for row in T] + [[0]*11]
    for a in range(10):
        if x >> a & 1: adj[a][10] = 1
        else: adj[10][a] = 1
    Hp = ham(T, range(10)); Ht = ham(adj, range(11))
    dels = [ham(adj, [u for u in range(11) if u != v]) for v in range(11)]
    fcdel = [d for d in dels if d not in S9]
    mods = modules(adj, 11)
    sc = sorted(sum(r) for r in adj)
    print('H=%d (check %d) parent H=%d parentFC=%s scores=%s FC-deletions=%d/11 nontrivial modules=%d deletion H=%s' % (h, Ht, Hp, Hp not in S9, sc, len(fcdel), len(mods), sorted(dels)))
