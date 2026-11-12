"""Small exact engines; no imported historical numeric claims."""
from itertools import permutations,product
from collections import defaultdict

def hp(rows):
    n=len(rows)
    if not n:return 1
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,a in enumerate(D[mask]):
            if not a:continue
            todo=rows[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1
                D[mask|bit][j]+=a
    return sum(D[-1])

def hp_permutations(rows):
    return sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(rows)-1)) for p in permutations(range(len(rows))))

def switched_tt(n,mask):
    return tuple(sum(1<<j for j in range(n) if j!=i and ((i<j)^bool((mask>>i&1)^(mask>>j&1)))) for i in range(n))

def canonical(rows):
    cells=defaultdict(list)
    for i,r in enumerate(rows):cells[r.bit_count()].append(i)
    groups=[cells[k] for k in sorted(cells)]
    best=None
    for pieces in product(*(permutations(g) for g in groups)):
        p=sum((tuple(g) for g in pieces),())
        r=tuple(sum(1<<j for j,v in enumerate(p) if rows[u]>>v&1) for u in p)
        if best is None or r<best:best=r
    return best

def locally_transitive(rows):
    full=(1<<len(rows))-1
    for i,r in enumerate(rows):
        for hood in (r,full^r^(1<<i)):
            vs=[j for j in range(len(rows)) if hood>>j&1]
            scores=sorted((rows[j]&hood).bit_count() for j in vs)
            if scores!=list(range(len(vs))):return False
    return True

def strong(rows):
    if not rows:return False
    for reverse in (False,True):
        reached=1;todo=1
        while todo:
            bit=todo&-todo;todo-=bit;i=bit.bit_length()-1
            r=(sum(1<<j for j in range(len(rows)) if rows[j]>>i&1) if reverse else rows[i])
            new=r&~reached;reached|=new;todo|=new
        if reached!=(1<<len(rows))-1:return False
    return True

def cyclic_blocks(a,b,c):
    groups=[0]*a+[1]*b+[2]*c
    return tuple(sum(1<<j for j in range(a+b+c) if (groups[i]==groups[j] and i<j) or (groups[j]-groups[i])%3==1) for i in range(a+b+c))
