import numpy as np, sys, math
from collections import defaultdict, Counter
def run(c,k,phi,L):
    N=1<<L
    p=np.zeros(N,dtype=np.int64)
    for u in range(1,N):
        w=u>>1; p[u]=p[w] if u&1 else p[w]+w+(w&1)
    U=np.arange(N,dtype=np.int64); G=c*U+k*p+phi*(U&1)
    order=np.argsort(G,kind='stable'); Gs=G[order]
    groups=np.split(order,np.nonzero(np.diff(Gs))[0]+1)
    def rev(x,a):
        r=0
        for i in range(a): r=(r<<1)|((x>>i)&1)
        return r
    lines=defaultdict(int); tot=0
    for g in groups:
        if len(g)<2: continue
        g=sorted(int(x) for x in g)
        for i in range(len(g)):
            for j in range(i+1,len(g)):
                u,v=g[i],g[j]; l=(u^v).bit_length()
                b=u&((1<<l)-1); bp=v&((1<<l)-1)
                if b<bp: b,bp=bp,b
                dl=b^bp; s=(dl&-dl).bit_length()-1
                lines[(l,s,b>>s,bp>>s)]+=1; tot+=1
    cls=Counter()
    for (l,s,H0,H1),m in lines.items():
        a=l-s; h=H0-H1; d=rev(H0,a)-rev(H1,a); g=math.gcd(h,d)
        if h==d: key='palindromic (h=d)'
        elif m>=2 and (d*1.0/h).is_integer() or (h*1.0/d).is_integer() and m>=2: key='other rational, multi'
        elif m>=2: key='other multi'
        else: key='single-point line'
        cls[key]+=m
    print('row %dU+%dp+%d[odd], L=%d: unordered collisions %d (%.4f per 2^L)'%(c,k,phi,L,tot,tot/N))
    for kk,v in sorted(cls.items(),key=lambda x:-x[1]): print('   %-24s %7d  %.1f%%'%(kk,v,100*v/max(tot,1)))
for args in [(18,3,2,18),(2,1,0,18),(189,11,0,18)]: run(*args)
