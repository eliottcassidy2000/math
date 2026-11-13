# check cap(Q_T) = s(T) + d(T) for every tournament of order <= 9:
#   s = number of singleton strong components, d = number of even strong components with a detachable partner
import collections
from trade import decode
DET={(6,(15,186,828,1752,1800,720)),(8,(95,2362,18756,73464,162000,205920,141120,40320))}
prof={}
for n in range(3,10):
    for line in open('prof_%d.txt'%n):
        p=line.split(); prof[(n,int(p[1]))]=tuple(map(int,p[7:]))
def sccs(out):
    n=len(out)
    reach=[0]*n
    for v in range(n):
        seen=1<<v; fr=1<<v
        while fr:
            nf=0
            for w in range(n):
                if fr>>w&1: nf|=out[w]
            nf&=~seen; seen|=nf; fr=nf
        reach[v]=seen
    comps=[]; done=0
    for v in range(n):
        if done>>v&1: continue
        c=sum(1<<w for w in range(n) if reach[v]>>w&1 and reach[w]>>v&1)
        comps.append(c); done|=c
    return comps
def sub(out,mask):
    idx=[v for v in range(len(out)) if mask>>v&1]; pos={v:i for i,v in enumerate(idx)}
    return [sum(1<<pos[w] for w in idx if out[v]>>w&1) for v in idx]
# profile of a strong component: compute via small DP (orders <= 8)
def profile(out):
    n=len(out); N=1<<n; end=[[0]*n for _ in range(N)]; Hs=[0]*N
    for S in range(1,N):
        tot=0
        for v in range(n):
            if not S>>v&1: continue
            P=S&~(1<<v)
            a=1 if P==0 else sum(end[P][w] for w in range(n) if P>>w&1 and out[w]>>v&1)
            end[S][v]=a; tot+=a
        Hs[S]=tot
    g=Hs[:]; F=[g[N-1]]
    for j in range(2,n+1):
        ng=[0]*N
        for S in range(N):
            s=0; B=S
            while B:
                A=S&~B
                if A: s+=g[A]*Hs[B]
                B=(B-1)&S
            ng[S]=s
        g=ng; F.append(g[N-1])
    return tuple(F)
bad=0; groups=collections.defaultdict(set); checked=0
cache={}
for (n,code),F in prof.items():
    out=decode(code,n); comps=sccs(out)
    s=sum(1 for c in comps if bin(c).count('1')==1)
    d=0
    for c in comps:
        k=bin(c).count('1')
        if k in (6,8):
            Y=sub(out,c); key=(k,tuple(Y))
            if key not in cache: cache[key]=profile(Y)
            if (k,cache[key]) in DET: d+=1
    groups[(n,F)].add((s,s+d))
    checked+=1
viol=[(k,v) for k,v in groups.items() if len(set(x[1] for x in v))>1 or max(x[0] for x in v)!=next(iter(v))[1]]
print('classes checked',checked,'profile groups',len(groups),'violations',len(viol))
for k,v in viol[:5]: print(k[0],k[1][:3],sorted(v))
