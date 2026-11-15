"""Independent literal replay of the signed binary triangle-seam decoder."""
from functools import cache
from random import Random

def check(x):
    if not x:
        raise AssertionError("triangle decoder failed")

@cache
def reads(rows):
    n=len(rows); f=(1<<n)-1
    h=[1]+[0]*f
    g=[[0]*n for _ in range(f+1)]
    for s in range(1,f+1):
        for i in range(n):
            if s>>i&1:
                t=s^(1<<i)
                g[s][i]=sum(g[t][j] for j in range(n) if rows[i]&t&(1<<j)) if t else 1
        h[s]=sum(g[s])
    return h[f], sum(h[s]*h[s^f] for s in range(f+1)), [sum(g[s][i]*h[f^s] for s in range(1,f+1)) for i in range(n)]

def inv(rows):
    return tuple(sum(1<<i for i,v in enumerate(rows) if v>>j&1) for j in range(len(rows)))

def glue(a,b):
    return tuple(v|(((1<<len(b))-1)<<len(a)) for v in a)+tuple(v<<len(a) for v in b)

def apex(b,u):
    return tuple(v|((1<<len(b)) if not u>>i&1 else 0) for i,v in enumerate(b))+(u,)

def replay(deck):
    B=(); w=mark=0; m=r=0
    seams={}
    for triangle,j in deck:
        mark|=((1<<j)-1)<<len(B)
        B=glue(B,(2,4,1) if triangle else (0,))
        if triangle and j in (1,2):
            seams[m]=3-2*j
        width=2 if triangle else 1
        w=(w<<width)|j
        m+=width
        r+=triangle
    h,Q,W=reads(B)
    U=reads(inv(B))[2]
    s=sum(x for i,x in enumerate(W) if mark>>i&1)
    t=sum(x for i,x in enumerate(U) if not mark>>i&1)
    rc=sum((1-(w>>i&1))<<(m-1-i) for i in range(m))
    check((h,Q,s)==(3**r,(3**r)<<m,(3**r)*w) and t%h==0)
    d=t//h-rc
    pos=0
    found={}
    while d:
        z=0 if not d%2 else (1 if d%4==1 else -1)
        if z: found[pos]=z
        d=(d-z)//2
        pos+=1
    check(found==seams)
    C=(); c=0; p=0
    while p<m:
        if p in found:
            j=(3-found[p])//2
            width=2
            V=(2,4,1)
        else:
            j=(w>>(m-p-1))&1
            width=1
            V=(0,)
        c|=((1<<j)-1)<<len(C)
        C=glue(C,V)
        p+=width
    check(reads(apex(B,mark))[0]==3**(r-len(found))*reads(apex(C,c))[0])

G=Random(513)
for n in range(1,10):
    for _ in range(15):
        deck=[]; k=n
        while k:
            triangle=k>=3 and G.randrange(2)
            deck.append((triangle,G.randrange(4 if triangle else 2)))
            k-=3 if triangle else 1
        replay(deck)
print("135 literal apex-DP, reverse-complement and signed-seam decoder controls.")
