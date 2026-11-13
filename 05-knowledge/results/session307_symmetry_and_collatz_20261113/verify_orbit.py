from math import factorial
from random import Random

def covers(a):
    n=len(a); m=(1<<n)-1
    paths=[[0]*n for _ in range(m+1)]
    for i in range(n): paths[1<<i][i]=1
    for s in range(1,m+1):
        for i in range(n):
            if paths[s][i]:
                for j in range(n):
                    if not s>>j&1 and a[i][j]:
                        paths[s|(1<<j)][j]+=paths[s][i]*a[i][j]
    h=list(map(sum,paths)); f=[[0]*(n+1) for _ in range(m+1)]; f[0][0]=1
    for s in range(1,m+1):
        t=s
        while t:
            if t&(s&-s):
                for k in range((s^t).bit_count()+1):
                    f[s][k+1]+=h[t]*f[s^t][k]
            t=(t-1)&s
    return f[m]

rng=Random(307)
for p,m,f,count in ((3,1,3,20),(3,2,1,20),(3,3,1,35),(5,1,2,10),(5,2,0,10)):
    for _ in range(count):
        n=p*m+f; a=[[0]*n for _ in range(n)]; d=[[0]*m for _ in range(m)]
        def arc(x,y,b): a[x][y]=b; a[y][x]=1-b
        for i in range(m):
            for x in range(p):
                for y in range(x+1,p):
                    arc(i*p+x,i*p+y,int((y-x)%p<=(p-1)//2))
            for j in range(i+1,m):
                b=[rng.randrange(2) for _ in range(p)]
                d[i][j]=sum(b); d[j][i]=p-d[i][j]
                for x in range(p):
                    for y in range(p): arc(i*p+x,j*p+y,b[(y-x)%p])
            for z in range(p*m,n):
                b=rng.randrange(2)
                for x in range(p): arc(i*p+x,z,b)
        for x in range(p*m,n):
            for y in range(x+1,n): arc(x,y,rng.randrange(2))
        v=covers(a); r=covers([row[p*m:] for row in a[p*m:]]); w=covers(d)
        z=[0]*(n+1)
        for i,x in enumerate(r):
            for j,y in enumerate(w): z[i+p*j]+=x*y
        if [(x-y)%p for x,y in zip(v,z)]!=[0]*(n+1): raise ValueError('cover fold')
        if any(w[j]%p and (j-m)%2 for j in range(m+1)): raise ValueError('skew parity')
print('PASS: weighted invariant-cover fold and skew support on 95 exact controls')

for fixed in (0,1):
    for _ in range(3):
        p=3; sizes=[9,3]+([1] if fixed else [])
        offsets=[0,9,12]; n=sum(sizes); a=[[0]*n for _ in range(n)]
        for i,k in enumerate(sizes):
            for x in range(k):
                for y in range(x+1,k):
                    u,v=offsets[i]+x,offsets[i]+y
                    a[u][v]=int(y-x<=k//2);a[v][u]=1-a[u][v]
            for j,t in enumerate(sizes[i+1:],i+1):
                from math import gcd
                g=gcd(k,t); mask=[rng.randrange(2) for _ in range(g)]
                for x in range(k):
                    for y in range(t):
                        u,v=offsets[i]+x,offsets[j]+y
                        a[u][v]=mask[(y-x)%g];a[v][u]=1-a[u][v]
        c=covers(a)
        if any((v-(i==12+fixed))%p for i,v in enumerate(c)):
            raise ValueError('prime-power tier fold')
print('PASS: independent prime-power tiers on six mixed 9+3(+1) controls')
