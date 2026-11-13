A=(2,28,25,81,97,79,7)
P=tuple(x|128 for x in A)+(0,)
B=(2,80,3,7,173,143,189,15)
def counts(R):
    n=len(R);N=1<<n
    d=[[0]*n for _ in range(N)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,N):
        for i,a in enumerate(d[s]):
            if not a:continue
            t=R[i]&~s
            while t:
                b=t&-t;t^=b;d[s|b][b.bit_length()-1]+=a
    return [1]+[sum(x) for x in d[1:]]
def extend(R,v,a,b=None):
    n=len(R);F=((1<<n)-1)^(1<<v)
    x=(F|((1-a)<<v))|((1<<(n+1)) if b is not None else 0)
    res=[r|((a and j==v)<<n) for j,r in enumerate(R)]+[x]
    if b is not None:
        for j in range(n):
            if j!=v or not b:res[j]|=1<<(n+1)
        res.append(b<<v)
    return res
def require(condition):
    if not condition: raise ValueError('exact marked-port check failed')
for R,one,single,expected in (
    (P,[297,325,405,405,409,741,445,95],
       [411,245,239,239,197,165,255,1371],[752,416,608,608,416,752,672,1276]),
    (B,[145,191,215,299,609,577,673,413],
       [901,219,633,451,193,239,159,327],[868,128,868,868,612,868,612,868])):
    old=counts(R)
    for v in range(len(R)):
        require(counts(extend(R,v,1))[ -1 ]==single[v])
        X=tuple(r|((i!=v)<<8) for i,r in enumerate(R))+(1<<v,)
        require(counts(X)[-1]==one[v])
        h={(a,b):counts(extend(R,v,a,b))[-1] for a in (0,1) for b in (0,1)}
        F=255^(1<<v)
        q=sum(old[s]*old[F^s] for s in range(256) if s&~F==0)
        require(h[1,1]-h[1,0]-h[0,1]+h[0,0]==q==expected[v])
print('PASS one- and two-exterior selective port detectors, 16 marked cores')
