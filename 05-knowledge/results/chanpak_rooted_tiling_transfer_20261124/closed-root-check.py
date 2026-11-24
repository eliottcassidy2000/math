"""Exact replay for the closed-root obstruction."""
from pathlib import Path
from runpy import run_path

lib=run_path(str(Path(__file__).with_name("rooted-transfer-check.py")))
tiles, rows, enc, paths, ref, classify, adjacent=(
    lib[k] for k in ("tiles","rows","encode","all_paths","reflect","classify","adjacent"))

def closed(n,x):
    return not(x>>tiles(n).index((0,n-1))&1)

def closed_parts(n,fibre):
    unseen={x for x in fibre if closed(n,x)}
    parts=[]
    while unseen:
        a=unseen.pop()
        todo=[a]
        for x in todo:
            for y,_ in adjacent(n,x):
                if y in unseen:
                    unseen.remove(y)
                    todo.append(y)
        parts.append(set(todo))
    return parts

for n,expected in ((3,(0,0)),(4,(0,0)),(5,(1,0)),(6,(3,2))):
    unseen=set(range(1<<len(tiles(n))))
    split=rigid=0
    while unseen:
        fibre,g=classify(n,min(unseen))
        unseen-=fibre
        if len(closed_parts(n,fibre))>1:
            split+=1
            rigid+=len(g)==1
    assert (split,rigid)==expected

for a,expected in (
 ((18,12,9,17,6),({18},{41})),
 ((34,4,9,19,39,14),({49,771,407,672,757,978},{740,470,834,117,32,263}))):
    n=len(a)
    pc=paths(a)
    m=enc(a,pc[0])
    fibre,g=classify(n,m)
    assert sorted(map(sorted,closed_parts(n,fibre)))==sorted(map(sorted,expected))
    if n==5:
        assert fibre=={12,18,41} and len(pc)==15 and len(g)==5
        p=(18,12,41)
    else:
        assert len(pc)==29 and len(g)==1
        assert [enc(a,tuple(map(int,w))) for w in ("340512","345120","534120")]==[757,250,470]
        p=(757,250,470)
    for x,y in zip(p,p[1:]):
        assert y in {t for t,_ in adjacent(n,x)}
    assert all(ref(n,x)==x and closed(n,x) for x in p[::2])
    assert not closed(n,p[1]) and (ref(n,p[1])==p[1])==(n==5)

def ok(c,edge):
    return all(edge(x,y) for x,y in zip(c,c[1:]+c[:1]))

def first_erasure(c,edge):
    n=len(c)
    return min(d-1 for d in range(2,n-1) if any(
        edge(c[i],c[(i+d)%n]) for i in range(n)))

for n in range(5,62):
    m=n//2
    if n%2:
        k=m-1
        edge=lambda a,b: (b-a)%n==1 or n-m<=(b-a)%n<=n-2
        cut=set(range(k+2,n))
        p=([0,3,1,4,2] if k==1 else
           list(range(k+1))+[2*k+2]+list(range(k+2,2*k+2))+[k+1])
        second=[(-2*i)%n for i in range(n)]
        assert ok(second,edge)
        assert any(edge(a,b)!=edge((-2*a)%n,(-2*b)%n) for a in range(n) for b in range(n))
    else:
        edge=lambda a,b: ((b-a)%n==1 or m<(b-a)%n<=n-2 or
                          ((b-a)%n==m and 1<=a<=m))
        cut=set(range(m+1,n))
        p=([0,1,4,2,5,3] if m==3 else
           list(range(m))+[n-1]+list(range(m+1,n-1))+[m])
    c=list(range(n))
    assert set(p)==set(c) and len(p)==n
    base=[x for x in c if x not in cut]
    assert base==[x for x in p if x not in cut]
    assert ok(c,edge) and ok(base,edge) and ok(p,edge)
    assert len(cut)==first_erasure(c,edge)==(n-1)//2
    assert first_erasure(p,edge)==1

from itertools import permutations
for bits in range(8):
    def e6(a,b):
        d=(b-a)%6
        if d==3:
            low=min(a,b)
            return ((bits>>low)&1)==(a==low)
        return d in (1,4)
    base=list(range(6))
    assert first_erasure(base,e6)==2
    alternatives=[list(p) for p in permutations(base[1:])]
    assert any(ok(p,e6) and first_erasure(p,e6)==1 and any(
        ok([x for x in base if x not in (a,b)],e6)
        and [x for x in p if x not in (a,b)]==[x for x in base if x not in (a,b)]
        for a in base for b in base if a<b) for p in ([0]+x for x in alternatives))

print("Closed midpoint bridges, census, and sharp deletion examples checked.")
