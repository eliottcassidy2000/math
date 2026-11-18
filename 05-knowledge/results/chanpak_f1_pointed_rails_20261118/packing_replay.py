from collections import Counter
from itertools import combinations

def cycles(clones, high):
    names = '0'+clones+'BC456'
    edges = set(combinations(names,2))
    changes = [('0','B'),('0','C')] + ([('0','4'),('4','6')] if high else [])
    for i,j in changes:
        edges.remove((i,j))
        edges.add((j,i))
    out = {x:{y for i,y in edges if i==x} for x in names}
    found = []
    def dfs(word, free):
        if len(word)>=3 and len(word)%2 and word[0] in out[word[-1]]:
            found.append(word)
        for y in free & out[word[-1]]:
            if y>word[0]: dfs(word+y,free-{y})
    for x in names: dfs(x,set(names)-{x})
    return found

for a in range(1,8):
    names='abcdefg'[:a]
    p,q=cycles(names,False), cycles(names,True)
    d=len(p)
    assert d==3*(1<<(a-1))-1 and len(q)==5*d
    for found, target in [(p,0),(q,d)]:
        masks=list(map(set,found))
        assert sum(not (x&y) for x,y in combinations(masks,2))==target
        disjoint = [(x,y) for x,y in combinations(masks,2) if not x&y]
        assert not any(not (z & (x|y)) for x,y in disjoint for z in masks)
    fiber=lambda pack: Counter(''.join(x for x in names[1:] if x in r) for r in pack)
    low,high=fiber(p),fiber(q)
    assert high==Counter({u:5*c for u,c in low.items()})
    print(a, (1,d,0), (1,5*d,d))
