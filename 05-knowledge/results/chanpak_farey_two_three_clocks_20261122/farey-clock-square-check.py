"""Literal arc-flip squares and ordered-slot reads on small asymmetric cores."""
from itertools import product
from random import Random

def ordinal_pad(rows,m):
    n=len(rows)
    return tuple(((1<<(n+m))-(1<<(i+1))) for i in range(m))+tuple(r<<m for r in rows)

def wrapper(rows, closed):
    n=len(rows)
    return (((1<<(n+2))-2)-(int(closed)<<(n+1)),)+tuple(
        (r<<1)|(1<<(n+1)) for r in rows)+(int(closed),)

def subpaths(rows):
    n=len(rows); M=1<<n
    end=[[0]*n for _ in range(M)]
    ans=[1]+[0]*(M-1)
    for mask in range(1,M):
        for b in range(n):
            if not (mask>>b)&1:continue
            rest=mask^(1<<b)
            end[mask][b]=1 if not rest else sum(end[rest][a] for a in range(n) if (rows[a]>>b)&1)
        ans[mask]=sum(end[mask])
    return ans

def physical(rows):
    reads=subpaths(rows)
    ret=[]
    for slots in (1,2,3):
        total=0
        for seq in product(range(slots),repeat=len(rows)):
            ms=[0]*slots
            for i,j in enumerate(seq):ms[j]|=1<<i
            weight=1
            for m in ms:weight*=reads[m]
            total+=weight
        ret.append(total)
    return ret

if __name__=="__main__":
    rng=Random(6123)
    tests=0
    for n in range(1,6):
        for _ in range(5):
            rows=[0]*n
            for i in range(n):
                for j in range(i):
                    if rng.randrange(2):rows[i]|=1<<j
                    else:rows[j]|=1<<i
            h,Q,R=physical(rows)
            for m in range(3):
                base=ordinal_pad(rows,m)
                X,Y=(1<<m)*Q,3**m*R
                w=[subpaths(wrapper(wrapper(base,a),b))[-1] for a,b in ((0,0),(1,0),(0,1),(1,1))]
                assert w==[h,h+X,h+4*X,-h+5*X+2*Y],(n,m,rows,w)
                tests+=1
    print(tests,"literal flip-square controls PASS")
