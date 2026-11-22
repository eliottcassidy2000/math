from itertools import combinations, product
from functools import lru_cache
@lru_cache(None)
def tiles(n): return [(a,b) for a in range(1,n+1) for b in range(1,a-1)]
@lru_cache(None)
def H(n,t):
    w=[0]*n
    for i in range(n-1): w[i]|=1<<(i+1)
    for k,(a,b) in enumerate(tiles(n)):
        if t>>k&1: w[a-1]|=1<<(b-1)
        else: w[b-1]|=1<<(a-1)
    a=[[0]*n for _ in range(1<<n)]
    for i in range(n): a[1<<i][i]=1
    for s,row in enumerate(a):
        for i,c in enumerate(row):
            if not c: continue
            z=w[i]&~s
            while z:
                b=z&-z;z^=b;a[s|b][b.bit_length()-1]+=c
    return sum(a[-1])
@lru_cache(None)
def restricted_index(n,x,y,z):
    return [tiles(n).index((a+y+z,b+z)) for a,b in tiles(n-x-y-z)]
fq=[abc for abc in product(range(2),repeat=3) if any(abc)]
def faces(n,t):
    return tuple(H(n-sum(s),sum(((t>>i)&1)<<j for j,i in enumerate(restricted_index(n,*s)))) for s in fq)
for n in range(4,8):
    mem={}
    rec=0
    apex=1<<tiles(n).index((n,1))
    for t in range(1<<len(tiles(n))):
        if t&apex: continue
        p=(faces(n,t),faces(n,t^apex))
        hh=H(n,t),H(n,t^apex)
        if p in mem:
            u,v=mem[p]
            if abs(hh[1]-hh[0])==2 and abs(v[1]-v[0])!=2:
                print('n',n,'old',u,v,'new',t,hh,'face vectors',p,'order',fq,flush=True);break
            if abs(hh[1]-hh[0])!=2 and abs(v[1]-v[0])==2:
                print('n',n,'old',u,v,'new',t,hh,'face vectors',p,'order',fq,flush=True);break
        mem[p]=(t,hh)
    else: print('none',n,len(mem),flush=True)
