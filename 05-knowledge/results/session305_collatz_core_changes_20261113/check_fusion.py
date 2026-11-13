from fractions import Fraction as Q
from random import Random
from collections import Counter

def paths(ar):
    n=len(ar)
    dp=[Counter() for _ in range(1<<n)]
    for i in range(n): dp[1<<i][(i,i)]=1
    for m in range(1,1<<n):
        for (s,t),z in dp[m].items():
            for i in range(n):
                if not m>>i&1 and ar[t][i]:
                    dp[m|(1<<i)][s,i]+=z
    return dp

def stats(ar, out):
    dp=paths(ar)
    full=len(dp)-1
    P=S=E=F=A=J=D=0
    for (s,t),z in dp[full].items():
        P+=z; S+=z*out[s]; E+=z*(1-out[t])
    for m in range(1,full):
        for (s,t),z in dp[m].items():
            for (u,v),w in dp[full^m].items():
                F+=z*w; A+=z*w*(1-out[t]); J+=z*w*out[u]; D+=z*w*(1-out[t])*out[u]
    return P,*(Q(x,P) for x in (S,E,F,A,J,D))

def check(blocks, ports):
    from math import prod
    t=[stats(ar,out) for ar,out in zip(blocks,ports)]
    p,s,e,f,a,j,d=zip(*t)
    r=len(t)
    hh=s[0]+e[-1]+sum(d)+sum(e[i]*s[i+1] for i in range(r-1))
    hh+=sum((s[i]+j[i])*(e[k]+a[k])*prod(2+f[z] for z in range(i+1,k))
            for i in range(r) for k in range(i+1,r))
    ids=[(i,u) for i,b in enumerate(blocks) for u in range(len(b))]
    ar=[[int(i<k) if i!=k else blocks[i][u][v] for k,v in ids]+[1-ports[i][u]] for i,u in ids]
    ar.append([ports[i][u] for i,u in ids]+[0])
    dp=paths(ar); h=sum(dp[-1].values())
    if h!=prod(p)*hh: raise ValueError((blocks,ports,h,prod(p)*hh))
    if r>=2:
        cy=sum(z for (u,v),z in dp[(1<<len(ids))-1].items() if ar[-1][u] and ar[v][-1])
        if cy!=prod(p)*s[0]*e[-1]: raise ValueError('cycle')
    if all(len(set(out))==1 for out in ports):
        w=[out[0] for out in ports]
        hh2=1+sum((k==i+1)+prod(2+f[z] for z in range(i+1,k))*(1+f[i])*(1+f[k])
                  for i in range(r) for k in range(i+1,r) if w[i] and not w[k])
        if hh2!=hh: raise ValueError('uniform')

rng=Random(305)
for trial in range(500):
    sizes=[rng.randrange(1,5) for _ in range(rng.randrange(1,5))]
    if sum(sizes)>10: sizes=[1]+sizes[:2]
    bs=[]; ps=[]
    for n in sizes:
        ar=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                ar[i][j]=rng.randrange(2); ar[j][i]=1-ar[i][j]
        bs.append(ar)
        ps.append([rng.randrange(2)]*n if trial%2 else [rng.randrange(2) for _ in range(n)])
    check(bs,ps)
print('500 fixed-seed shared-apex path/cycle cases passed')

def block(n, code):
    ar=[[0]*n for _ in range(n)]
    k=0
    for i in range(n):
        for j in range(i+1,n):
            ar[i][j]=(code>>k)&1
            ar[j][i]=1-ar[i][j]
            k+=1
    return ar

def cyclic_blowup(t):
    sizes=list(map(len,t))
    ids=[(i,j) for i,s in enumerate(sizes) for j in range(s)]
    return [[t[i][j][k] if i==u else int((u-i)%3==1)
             for u,k in ids] for i,j in ids]

for n,c,P,F,G,q in ((7,245600,29,522,3402,20),(7,168491,33,594,3714,20),
                    (7,1104621,37,666,4026,20),(3,0,1,6,6,8),(4,2,3,18,36,8),
                    (5,41,13,78,198,8),(5,76,15,90,210,8)):
    ar=block(n,c)
    z=[sum(a.values()) for a in paths(ar)]
    m=len(z)-1
    f=sum(z[a]*z[m^a] for a in range(1,m))
    g=sum(z[a]*z[b]*z[m^a^b] for a in range(1,m)
          for b in range(1,m^a) if not (a&b))
    if (z[-1],f,g,Q(2*z[-1]+f,z[-1]))!=(P,F,G,q):
        raise ValueError('clock bank')
    t2=[[0,1],[0,0]]
    got=sum(paths(cyclic_blowup([ar,t2,t2]))[-1].values())
    if got!=11*P+17*F+4*G:
        raise ValueError('three-visit contrast')
print('Seven equal-clock witnesses and three-visit contrasts passed')
