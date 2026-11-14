"""Independent literal checks of the finite-set support decoder and minimal switch."""
from random import Random
from math import prod

def paths(rows):
    n=len(rows);M=1<<n
    start=[[0]*n for _ in range(M)]
    H=[1]+[0]*(M-1)
    for s in range(1,M):
        for i,r in enumerate(rows):
            if not (s>>i&1):continue
            t=s^(1<<i)
            start[s][i]=sum(start[t][j] for j in range(n) if r>>j&1) if t else 1
        H[s]=sum(start[s])
    return H,start

def sig(rows,U):
    n=len(rows);M=(1<<n)-1
    H,A=paths(rows)
    _,B=paths(tuple((M^a)^(1<<i) for i,a in enumerate(rows)))
    P=sum(sum(A[s][i] for i in range(n) if U>>i&1)*H[M^s] for s in range(1,M+1))
    Q=sum(sum(B[s][i] for i in range(n) if not (U>>i&1))*H[M^s] for s in range(1,M+1))
    S=sum(B[M][i] for i in range(n) if U>>i&1)
    T=sum(A[M][i] for i in range(n) if not (U>>i&1))
    return H[M],sum(H[s]*H[M^s] for s in range(M+1)),P,Q,S,T

def ordinal(blocks):
    n=sum(len(rows) for rows,_ in blocks)
    R=[];U=0
    for rows,a in blocks:
        j=len(R);U|=a<<j
        R += [(r<<j)|((1<<n)-(1<<(j+len(rows)))) for r in rows]
    return R,U

def augmented(rows,U):
    n=len(rows)
    return [r|((0 if U>>i&1 else 1)<<n) for i,r in enumerate(rows)]+[U]

def tt(l):
    return tuple((1<<l)-(1<<(i+1)) for i in range(l))

def rev(x,n):
    a=0
    for _ in range(n):
        a=2*a+x%2;x//=2
    return a

def recover(bits_shape, P):
    # For a singleton run, bits_shape has its width. For a fixed
    # nonsingleton mask the precomputed (h,f,P) entry replaces the run.
    hs=[item[0] if isinstance(item,tuple) else 1 for item in bits_shape]
    h=prod(hs)
    out=[]
    for item in reversed(bits_shape):
        if isinstance(item,int):
            Q=1<<item
            d=P*pow(h,-1,Q)%Q
            out.append(rev(d,item))
            d*=h
        else:
            h0,Q,p0=item
            h//=h0
            d=h*p0
        P,r=divmod(P-d,Q)
        if r or P<0:return None
    if P!=0:return None
    return list(reversed(out))

def naf(x):
    S={}
    while x:
        j=(abs(x)&-abs(x)).bit_length()-1
        s=1 if (x>>j)%4==1 else -1
        S[j]=s
        x-=s*(1<<j)
    return S

rng=Random(590320122)
C3=(2,4,1);C5=(8,1,3,6)
for step in range(180):
    cores=[C3]*rng.randrange(4)
    if step%5==0:cores=[C5]+cores[:1]
    widths=[rng.randrange(2) for _ in range(len(cores)+1)]
    blocks=[];shape=[];word='';marked={};expected=[]
    for j,w in enumerate(widths):
        a=rng.randrange(1<<w)
        blocks.append((tt(w),a));shape.append(w);expected.append(a)
        if w:word+=format(rev(a,w),f'0{w}b')
        if j<len(cores):
            C=cores[j];v=rng.randrange(1<<len(C));blocks.append((C,v))
            shape.append(sig(C,v)[:3])
            k=v.bit_count()
            if k in (1,2):marked[len(word)]=1 if k==1 else -1
            word+=format(k,'02b')
    T,V=ordinal(blocks)
    h,f,P,Q,_,_=sig(T,V)
    assert recover(shape,P)==expected,(shape,P,expected,blocks,T,V,recover(shape,P))
    if C5 in cores:continue
    x=int(word or '0',2)
    assert (f,P)==(h<<len(word),h*x)
    assert naf(Q//h-((1<<len(word))-1-rev(x,len(word))))==marked
    # Retain only marked mixed triangles; replace uniform triangles
    # by individual slots in the compressed binary support.
    filtered=[];j=0
    while j<len(word):
        if j in marked:
            filtered.append((C3,(1<<int(word[j:j+2],2))-1));j+=2
        else:
            filtered.append(((0,),int(word[j])));j+=1
    R,U=ordinal(filtered)
    K=paths(augmented(T,V))[0][-1]
    Knew=paths(augmented(R,U))[0][-1]
    assert K*3**len(marked)==Knew*h

A=(510,1008,994,486,1004,448,384,768,512,105)
B=(510,1020,504,496,896,272,816,864,512,45)
for G,K in [(A,579),(B,581)]:
    assert sig(tuple(x&511 for x in G[:-1]),G[-1])==(5,1024,716,810,0,0)
    assert paths(G)[0][-1]==K

def wordblock(x,m,rc=False):
    if rc:x=((1<<m)-1)^rev(x,m)
    return tt(m),rev(x,m)

for x,y in [(1,1),(1,2),(2,1),(2,2),(4,1),(1,4),(3,5),(5,3)]:
    yl=8*y+4; xr=8*x+3
    L,u=ordinal([wordblock(x,x.bit_length()),(C5,4),wordblock(yl,yl.bit_length(),True)])
    R,v=ordinal([wordblock(xr,xr.bit_length()),(C5,2),wordblock(y,y.bit_length(),True)])
    assert sig(L,u)==sig(R,v)
    assert paths(augmented(R,v))[0][-1]-paths(augmented(L,u))[0][-1]==2*(x%2+y%2-1)
print("PASS: 180 support decodings, literal minimum switch, growing parity branches.")
