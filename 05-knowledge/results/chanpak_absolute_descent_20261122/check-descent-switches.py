from itertools import product

def hops(rows):
    n=len(rows)
    f=[[0]*n for _ in range(1<<n)]
    for j in range(n):f[1<<j][j]=1
    for s in range(1,1<<n):
        for i,v in enumerate(f[s]):
            avail=rows[i]&~s
            while v and avail:
                e=avail&-avail
                f[s|e][e.bit_length()-1]+=v
                avail-=e
    return sum(f[-1])

def H(x):
    return 1+sum((2 if j-i==1 else 2**(j-i-1))*x[i]*(1-x[j])
                 for i in range(len(x)) for j in range(i+1,len(x)))

def rows(x):
    n=len(x)
    z=[sum(1<<j for j in range(i+1,n))|((1-x[i])<<n) for i in range(n)]
    return z+[sum(x[i]<<i for i in range(n))]

def kind(d):
    t=max((i for i,v in enumerate(d) if v),default=-1)
    if t==0:return 1
    if t>0 and d[0] and d[t]==-d[0] and all(v==d[0] for v in d[:t]):
        return -1
    return 0

for n in range(1,7):
    for x in product(range(2),repeat=n):
        assert H(x)==hops(rows(x))
for n in range(2,9):
    for p in range(1,n):
        for d in product((-1,0,1),repeat=p):
            if not any(d):continue
            for e in product((-1,0,1),repeat=n-p):
                if not any(e):continue
                u=tuple(int(v<0) for v in d);v=tuple(int(w<0) for w in e)
                up=tuple(a+z for a,z in zip(u,d));vp=tuple(a+z for a,z in zip(v,e))
                D=H(up+vp)-H(up+v)-H(u+vp)+H(u+v)
                A=sum(z<<i for i,z in enumerate(d[::-1]))
                B=sum(z<<i for i,z in enumerate(e))
                assert D==-A*B-d[-1]*e[0]
                assert (D==0)==(kind(d[::-1])*kind(e)==-1)
print("exact spine cuts and adjacent carry classification PASS")
