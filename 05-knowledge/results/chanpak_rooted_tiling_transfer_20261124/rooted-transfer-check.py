"""Rooted tiling classifier and small exact replay. No third-party packages."""
from collections import Counter
from functools import lru_cache
from itertools import permutations
from math import comb
from pathlib import Path
from runpy import run_path

@lru_cache(None)
def tiles(n):
    return tuple((i,j) for i in range(n) for j in range(i+2,n))

def rows(n,m):
    a=[0]*n
    for i in range(n-1):
        a[i]|=1<<(i+1)
    for k,(i,j) in enumerate(tiles(n)):
        if m>>k&1:
            a[i]|=1<<j
        else:
            a[j]|=1<<i
    return tuple(a)

def encode(a,path):
    return sum(((a[path[i]]>>path[j])&1)<<k for k,(i,j) in enumerate(tiles(len(a))))

def valid(a,p):
    return all(a[i]>>j&1 for i,j in zip(p,p[1:]))

def all_paths(a):
    return [p for p in permutations(range(len(a))) if valid(a,p)]

def reflect(n,m):
    lookup={pair:i for i,pair in enumerate(tiles(n))}
    return sum((m>>lookup[n-1-j,n-1-i]&1)<<k for k,(i,j) in enumerate(tiles(n)))

def adjacent(n,m):
    a=rows(n,m)
    for v in range(n):
        p0=tuple(i for i in range(n) if i!=v)
        if not valid(a,p0):
            continue
        for j in range(n):
            if j!=v:
                p=p0[:j]+(v,)+p0[j:]
                if valid(a,p):
                    yield encode(a,p),p

def comp(p,q):
    return tuple(p[i] for i in q)

def inv(p):
    r=list(p)
    for a,b in enumerate(p):
        r[b]=a
    return tuple(r)

def classify(n,m,voltage=True):
    """Return the complete ordinary fibre and its automorphism group."""
    e=tuple(range(n))
    tree={m:e}
    queue=[m]
    loops=set()
    for x in queue:
        for y,q in adjacent(n,x):
            p=comp(tree[x],q)
            if y not in tree:
                queue.append(y)
                tree[y]=p
            elif voltage:
                loops.add(comp(p,inv(tree[y])))
    group={e}
    todo=[e]
    for a in todo:
        for b in loops:
            y=comp(a,b)
            if y not in group:
                group.add(y)
                todo.append(y)
    return set(tree),group

def count_cycles(a):
    return sum(p[0]==0 and a[p[-1]]&1==1 for p in all_paths(a))

expected_classes={3:2,4:4,5:12,6:56}
for n,expected in expected_classes.items():
    unseen=set(range(1<<len(tiles(n))))
    found=[]
    while unseen:
        m=min(unseen)
        masks,group=classify(n,m)
        paths=all_paths(rows(n,m))
        assert {encode(rows(n,m),p) for p in paths}==masks
        assert len(paths)==len(masks)*len(group)
        assert {p for p in paths if encode(rows(n,m),p)==m}==group
        assert len(masks)%2==1
        if len(group)==1:
            for x in masks:
                if reflect(n,x)==x:
                    assert all(reflect(n,y)!=y for y,_ in adjacent(n,x))
        found.append((len(paths),masks))
        unseen-=masks
    assert len(found)==expected
    eleven=[m for h,m in found if h==11]
    if n<=4:
        assert not eleven
    if n==5:
        assert len(eleven)==1
        E=eleven[0]

assert E=={4,6,9,13,20,22,33,40,43,44,57}
assert {x for x in E if reflect(5,x)==x}=={4,22,33}
assert {x for x in E if not(x>>tiles(5).index((0,4))&1)}=={9,33,40,43,57}
assert {(x,63^x) for x in E if x<63^x and (63^x) in E}=={(6,57),(20,43)}
assert reflect(5,6)==20 and reflect(5,57)==43
for x,y,h,doubled in ((9,40,5,True),(13,44,13,False)):
    choices=[classify(5,63^m)[0] for m in (x,y)]
    assert {reflect(5,m) for m in choices[0]}==choices[1]
    assert (choices[0]!=choices[1])==doubled
    assert all(len(all_paths(rows(5,min(s))))==h for s in choices)
assert rows(5,19)==(14,20,8,18,5)
assert (35,(0,2,3,1,4)) in set(adjacent(5,19))
assert [len(all_paths(rows(5,63^x))) for x in (19,35)]==[11,13]
assert not (4&8) and (9&8)

for n in (3,4):
    old=(1<<len(tiles(n)))-1
    new=(1<<len(tiles(n+2)))-1
    rank={p:i for i,p in enumerate(tiles(n+2))}
    def restriction(m):
        return sum(((m>>rank[i+1,j+1])&1)<<k for k,(i,j) in enumerate(tiles(n)))
    def defects(s,m):
        return (m^reflect(s,m)).bit_count()//2
    def line(m,x):
        return min(m,m^x)
    def orbit(s,m,x):
        return min(line(m,x),line(reflect(s,m),x))
    raw=Counter()
    projected=Counter()
    side=Counter()
    closebit=1<<rank[0,n+1]
    for m in range((new+1)//2):
        root=restriction(m)
        parent=line(root,old)
        k=defects(n+2,m)-defects(n,root)
        raw[parent,k]+=1
        side[parent,k,restriction(m if m&closebit==0 else m^new)]+=1
        if orbit(n+2,m,new)==m:
            projected[orbit(n,root,old),k]+=1
    for p in range((old+1)//2):
        for k in range(n):
            target=2**n*comb(n-1,k)
            assert raw[p,k]==target
            assert side[p,k,p]==side[p,k,p^old]==target//2
            if orbit(n,p,old)==p:
                fixed=reflect(n,p)==p
                assert projected[p,k]==target//(2 if fixed and k else 1)

for a,h in (((18,36,9,19,38,77,31),7),((114,68,9,19,38,78,156,63),21)):
    assert count_cycles(a)==h
for n in (3,4,5):
    for m in range(1<<len(tiles(n))):
        a=rows(n,m)
        b=n
        c=n+1
        wrap=tuple(x|1<<b for x in a)+(1<<c,(1<<n)-1)
        assert count_cycles(wrap)==len(all_paths(a))

folder=Path(__file__).resolve().parent.parent/"chanpak_lrc_two_edit_ports_20261123"
bank=run_path(str(folder/"third-slot-bank.py"))["data"]
data={str(i)+tag:(d,q,r[0]) for i,(d,q,a,b) in enumerate(bank)
      for tag,r in zip("ab",(a,b))}
words={frozenset(data)}
for p in (3,5,7,13):
    words={s & frozenset(k for k,t in data.items() if t[j]%p)
           for s in words for j,c in enumerate((11,42,20)) if c%p}
maxset={s for s in words if len(s)==max(map(len,words))}
assert maxset=={frozenset(t.split()) for t in
               ("0a 0b 1a 2a 2b 3a 3b 6a 6b 7a 7b 8a 8b 11a 11b",
                "0a 1a 2a 2b 3a 3b 6a 6b 7a 7b 8a 8b 11a 13a 13b")}
for support in maxset:
    for p in (3,5,7,13,11):
        assert sum(c%p!=0 and all(data[k][j]%p for k in support)
                   for j,c in enumerate((11,42,20)))==(p!=11)
print("Root class and voltage, blue/black shells, H=11 ports, cycle lift, and wheel checked.")
