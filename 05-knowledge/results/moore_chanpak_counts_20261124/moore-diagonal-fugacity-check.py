from itertools import combinations

ORDERS=["02341","02314","02431","02134","01432","01342","01243","01423"]

def rows(order, mask):
    order=tuple(map(int,order));m=len(order);n=2*m; r=[0]*n
    for i,j in combinations(range(n),2):
        if j<m:
            forward=True
        elif i>=m:
            forward=False
        else:
            j0=j-m
            forward=bool(mask>>i&1) if i==j0 else (order[i]-order[j0])%m in [1,m-1]
        if forward:r[i]|=1<<j
        else:r[j]|=1<<i
    return r

def count(r):
    n=len(r);full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,full):
        for i,x in enumerate(dp[s]):
            if x:
                nxt=r[i]&~s
                while nxt:
                    b=nxt&-nxt;nxt-=b
                    dp[s|b][b.bit_length()-1]+=x
    return sum(dp[full])

def starts(r,start):
    n=len(r);full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    dp[1<<start][start]=1
    for s in range(1,full):
        if not s>>start&1:continue
        for u,v in enumerate(dp[s]):
            if not v:continue
            nxt=r[u]&~s
            while nxt:
                b=nxt&-nxt;nxt-=b
                dp[s|b][b.bit_length()-1]+=v
    return sum(dp[-1])

def quotient(o):
    a=list(map(int,o));n=len(a);r=[0]*n
    for i,j in combinations(range(n),2):
        if (a[i]-a[j])%n in [1,n-1]:r[i]|=1<<j
        else:r[j]|=1<<i
    return r

def full_cycles(r):
    n=len(r);full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    dp[1][0]=1
    for mask in range(1,full,2):
        for i,ways in enumerate(dp[mask]):
            if not ways:continue
            nbs=r[i]&~mask
            while nbs:
                b=nbs&-nbs;nbs-=b
                dp[mask|b][b.bit_length()-1]+=ways
    return sum(ways for i,ways in enumerate(dp[full]) if r[i]&1)

def mobius(v):
    f=v[:]
    for j in range((len(f)-1).bit_length()):
        for i in range(len(f)):
            if i>>j&1:f[i]-=f[i^(1<<j)]
    return f

def cycles(r):
    n=len(r);out=[]
    for start in range(n):
        def visit(path,S):
            last=path[-1]
            if len(path)>=3 and len(path)%2 and r[last]>>start&1:
                out.append(tuple(path))
            cand=(r[last]&~S)&~((1<<start)-1)
            while cand:
                b=cand&-cand;cand-=b
                visit(path+[b.bit_length()-1],S|b)
        visit([start],1<<start)
    return out

def orbit_parity(r):
    n=len(r);m=n//2;rtn=0
    for c in cycles(r):
        d=tuple((i+m)%n for i in reversed(c))
        j=d.index(min(d));d=d[j:]+d[:j]
        assert d!=c
        if c<d and set(c)&set(d):rtn^=1
    return rtn

if __name__=='__main__':
    best=None; steps4=[];joint=set()
    for o in ORDERS:
        hs=[count(rows(o,s)) for s in range(32)]
        assert all(h%4==1 for h in hs)
        for s,h in enumerate(hs):
            if o==ORDERS[0]:
                assert h%8==1+4*orbit_parity(rows(o,s))
            for j in range(5):
                if not(s>>j&1):
                    t=s^(1<<j);d=abs(hs[t]-h)
                    if d==4:steps4.append((o,s,j,h,hs[t]))
                    if d and (best is None or d<best[0]):best=(d,o,s,j,h,hs[t])
        assert all(v%4==0 for v in mobius(hs)[1:])
        cs=[full_cycles(rows(o,s)) for s in range(32)]
        co=mobius(cs)
        assert all(c%2==0 for s,c in enumerate(co) if s.bit_count()>1)
        q=quotient(o)
        slope=[starts(q,j)%2 for j in range(5)]
        assert [co[1<<j]%2 for j in range(5)]==slope
        assert sum(slope)%2 and sum(c%2 for c in cs)==16
        print(o,hs[-1],sorted(set(h%8 for h in hs)),'C parity slopes',slope,flush=True)
    assert len(steps4)==2 and best==(4,'01342',1,2,4225,4221)
    assert steps4[-1]==('01342',18,2,8369,8365)
    assert all(8369%i for i in range(2,92))
    print('nonzero minimum' ,best,'steps of size 4',len(steps4),steps4[:5])
