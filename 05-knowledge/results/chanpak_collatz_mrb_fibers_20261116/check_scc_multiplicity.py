from collections import Counter
from itertools import combinations, permutations, product

def fg(n):
    perms = list(permutations(range(n)))
    pairs = list(combinations(range(n), 2))
    seen = set()
    f, a, g, c = [Counter() for _ in range(4)]
    for bits in product([0,1],repeat=len(pairs)):
        if bits in seen: continue
        E = {((i,j) if b else (j,i)) for (i,j),b in zip(pairs,bits)}
        h = sum(all((i,j) in E for i,j in zip(p,p[1:])) for p in perms)
        orbit = {tuple(int((p[i],p[j]) in E) for i,j in pairs) for p in perms}
        seen.update(orbit)
        f[h] += 1
        self_converse = tuple(1-b for b in bits) in orbit
        if self_converse: g[h] += 1
        strong = n > 0
        for v in range(n):
            reaches = {v}
            for _ in range(n):
                reaches |= {j for i,j in E if i in reaches}
            strong &= len(reaches) == n
        if strong:
            a[h] += 1
            if self_converse: c[h] += 1
    return f,a,g,c

rows = [fg(n) for n in range(6)]
F = []
for n, (f,a,g,c) in enumerate(rows):
    ff = Counter({1:1}) if n==0 else Counter()
    for i in range(1,n+1):
        for d,q in rows[i][1].items():
            for e,r in F[n-i].items():
                ff[d*e] += q*r
    gg = Counter()
    for i in range(n//2+1):
        middle = rows[n-2*i][3] if n != 2*i else Counter({1:1})
        for d,q in (F+[ff])[i].items():
            for e,r in middle.items():
                gg[d*d*e] += q*r
    assert f == ff and g == gg, (n, f, ff, g, gg)
    F.append(ff)
    print(n, sum(f.values()), sum(g.values()), sorted(f.items()))

# Literal certificates for the finite multiplicative atoms.
for h, masks in {
    49: [42,96,67,118,103,68,1],
    63: [188,189,240,244,160,64,147,32],
    343: [206,192,170,162,15,147,188,16],
}.items():
    n=len(masks)
    assert all(m < (1<<n) and not (m>>i&1) for i,m in enumerate(masks))
    assert all((masks[i]>>j&1)+(masks[j]>>i&1)==1 for i,j in combinations(range(n),2))
    assert sum(all(masks[i]>>j&1 for i,j in zip(p,p[1:])) for p in permutations(range(n)))==h
    for v in range(n):
        reaches={v}
        for _ in range(n):
            reaches |= {j for i in reaches for j in range(n) if masks[i]>>j&1}
        assert len(reaches)==n
    print("finite atom", h, "at order", n)
