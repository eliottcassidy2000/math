from itertools import product, permutations
from math import factorial

def per(M):
    return sum(all(M[i][j] for i,j in enumerate(p)) for p in permutations(range(len(M))))
def bit(d,h):return d>>(1-h)&1
def connector_columns(pq,x,y):
    if pq is None:return (1, int((x,y)==(0,1)),1)
    p,q=pq
    return (1-bit(q,x),bit(p,x)*(1-bit(q,y)),bit(p,y))
def count(types):
    t=len(types);types=tuple(types)+(None,)*(3-t)
    ps=[per(tuple(zip(*(connector_columns(pq,x,y) for pq in types)))) for x,y in ((0,1),(1,0))]
    return sum(ps)//factorial(3-t)
if __name__=='__main__':
    for t in range(4):
        for d in range(t+1):
            # d singleton runs, all remaining runs at least two vertices
            vals={count(p) for p in product(*([(tuple((u,u) for u in range(4)))]*d+[tuple(product(range(4),repeat=2))]*(t-d)))}
            print(t,d,sorted(vals))

def fibre_load(types):
    "Indexed by cut bits of three components and y; x fixed at 0."
    assert len(types)==3
    weights=[0]*16
    for x,y in ((0,1),(1,0)):
        def internal(id,h,into):
            p,q=types[id]
            return bit(p,h) if into else 1-bit(q,h)
        for a,b,c in permutations(range(3)):
            forward=(internal(a,x,False),internal(b,x,True),internal(b,y,False),internal(c,y,True))
            ids=(a,3+x,b,3+y,c)
            bits={ids[0]:0}
            for u,v,good in zip(ids,ids[1:],forward):
                bits[v]=bits[u]^(1-good)
            index=sum((bits[i]^bits[3])<<i for i in (0,1,2)) + 8*(bits[4]^bits[3])
            weights[index]+=1
    assert weights[0]==count(types) and sum(weights)==12
    return tuple(weights)

if __name__=='__main__':
    dist={}
    for t in product(tuple(product(range(4),repeat=2)),repeat=3):
        w=fibre_load(t);s=sum(z>0 for z in w)
        dist.setdefault(s,0);dist[s]+=1
    print('nonempty switch charts for one fixed 3-piece forest:',sorted(dist.items()))
