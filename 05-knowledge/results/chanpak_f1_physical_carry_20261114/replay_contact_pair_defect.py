"""Check marked-pair contractions, minimum frames, and exact gap-two rays."""
from itertools import product, combinations
from replay_ordinal_defect import paths, profile, trans, join

def graphs(n):
    pairs=list(combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        rows=[0]*n
        for k,(i,j) in enumerate(pairs):
            if bits>>k&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
        yield tuple(rows)

def data(rows,u,v):
    W=[i for i in range(len(rows)) if i not in (u,v)]
    Q=tuple(sum(((rows[i]>>j)&1)<<k for k,j in enumerate(W)) for i in W)
    n=len(W)
    c=[((rows[i]>>u)&1)+((rows[i]>>v)&1) for i in W]
    J=0
    H,end=paths(Q)
    _,start=paths(tuple(((1<<n)-1)^(1<<i)^r for i,r in enumerate(Q)))
    for X in range(1,(1<<n)-1):
        Y=((1<<n)-1)^X
        for i in range(n):
            C=((rows[W[i]]>>u)&1) | (((rows[W[i]]>>v)&1)<<1)
            if C not in (1,2): continue
            for j in range(n):
                B=((rows[W[j]]>>u)&1) | (((rows[W[j]]>>v)&1)<<1)
                if C==B: J+=end[X][i]*start[Y][j]
    mask=1<<n
    def contracted(orient):
        return tuple(r|(mask if flag else 0) for r,flag in zip(Q,orient))+(sum((not a)<<i for i,a in enumerate(orient)),)
    A=paths(contracted([q==2 for q in c]))[0][-1]
    O=paths(contracted([q!=0 for q in c]))[0][-1]
    order=sorted(range(n),key=lambda i:-Q[i].bit_count())
    typ=[((rows[W[i]]>>u)&1) | (((rows[W[i]]>>v)&1)<<1) for i in order]
    minimum=(H[-1]==1 and [c[i] for i in order]==sorted(c,reverse=True) and typ.count(1)<=1 and typ.count(2)<=1)
    suffix=bool(typ and typ[-1]==0)
    return A,O,J,minimum,suffix

if __name__=='__main__':
    checked=minima=0
    for n in range(2,6):
        for rows in graphs(n):
            for u,v in combinations(range(n),2):
                h,p,e,m,s,t=profile(rows,(1<<u)|(1<<v))
                A,O,J,b,z=data(rows,u,v)
                assert e==A+O+J and e>=2 and e%2==0
                assert b==(e==2)
                if b: assert (s==0)==z
                checked+=1
                minima+=b
    assert (checked,minima)==(10650,1922)
    for rows,h,p in [((2,0),1,3),((0,5,1),1,5),((2,4,1),3,6),((8,5,9,2),5,13),((10,4,9,2),5,14)]:
        for s in range(1,6):
            H,P,E,M,S,T=profile(join(rows,trans('0'*s)[0]),3)
            assert (H,P,E,S)==(h,p<<s,2,0)
    print('PASS:',checked,'marked pairs;',minima,'minimum cases; 25 literal ray witnesses')
