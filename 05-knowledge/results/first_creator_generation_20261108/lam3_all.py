"""Every 3-arc-strong 11-tournament with H in S(10) (the complete list from lam11_all): is it EXT (a first-creator
deletion) or SUB (a nontrivial module)?  Decides the generation grammar at order 11 for lambda >= 3."""
exec(open('lam3_struct.py').read().split('seen = set()')[0])
tot = ext = sub = both = neither = 0; byH = {}
for line in open('lambda3_witnesses_order11.txt'):
    _, code, x, h, st = line.split(); code = int(code); x = int(x); h = int(h)
    if st != 'inS10': continue
    T = decode(code, 10)
    adj = [row + [0] for row in T] + [[0]*11]
    for a in range(10):
        if x >> a & 1: adj[a][10] = 1
        else: adj[10][a] = 1
    dels = [ham(adj, [u for u in range(11) if u != v]) for v in range(11)]
    e = any(d not in S9 for d in dels); s = len(modules(adj, 11)) > 0
    tot += 1; ext += e; sub += s; both += e and s; neither += (not e and not s)
    byH.setdefault(h, set()).add(('EXT' if e else '') + ('SUB' if s else ''))
print('pairs with H in S(10): %d; EXT %d, SUB %d, both %d, neither %d' % (tot, ext, sub, both, neither))
print(byH)
