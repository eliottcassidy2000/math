"""Exact bit ceilings and a full-size tilted-interval codec replay."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, ceil, factorial
from random import Random

N=295
D=55

def bound(c):
    t=c//5
    r=(2065+t)//15
    M=8395-c
    U=M+2*(675-c)-2*(135-t)
    Q=F((295-c)*49+268*t,225)
    L=ceil(F(11*r*r,100)-4*Q-F(268*t+U,225)+F(175+72*t,15))
    return r,57*r-L

def caps():
    out=[]
    for c in (25,100,175):
        r,R=bound(c)
        nrm=[]; lab=[]
        for s in range(26):
            k=r-s
            a,b=divmod(s,5)
            L=11*(k*(539-k)+s*(49-s))//2
            top=56**L*(5+a)**((5-b)*(270-k))*(6+a)**(b*(270-k))
            den=55**(L-R)*5**(5*(270-k))
            # A dyadic point fits any half-open interval of length >= its mesh.
            lab.append(((top-1)//den).bit_length())
            nrm.append(((top-1)//(den*11**294)).bit_length())
        out.append((c,r,R,max(nrm),max(lab)))
    print('bit ceilings:',out)
    assert max(a[-2] for a in out)+17 <= comb(290,2)
    assert max(a[-2] for a in out)+297 <= comb(291,2)
    assert max(a[-1] for a in out)+297 <= comb(295,2)
    assert max(a[-2] for a in out)+17+42==comb(290,2)
    assert max(a[-2] for a in out)+297+52==comb(291,2)
    lo,hi,scale=27182818,27182819,10**7
    partial=sum(F(1,factorial(i)) for i in range(15))
    assert F(lo,scale)<partial and partial+F(16,15*factorial(15))<F(hi,scale)
    assert all(factorial(j)*hi**j<=((j+5)*scale)**j for j in range(1,292))
    for n,bits in [(290,1687),(291,1694)]:
        a=(n*(n+3)+10*(n+1))*scale
        b=2*(n+1)*lo
        assert a**(n+1)<(b**(n+1))<<bits
    return out

SS={k:list(combinations(range(11),k)) for k in (1,2,3)}
def norm(S):
    return min(tuple(sorted((x+a)%11 for x in S)) for a in range(11))
NN={k:sorted({norm(s) for s in ss}) for k,ss in SS.items()}

def factors(S):
    """For each shape factor yield choices: (weight, pair, multiplicity)."""
    w=[comb(11,m)*D**(3-m) for m in range(4)]
    S=set(S)
    for lo,hi in [(0,270),(270,N)]:
        for i,j in combinations(range(lo,hi),2):
            if i in S or j in S:
                yield tuple((w[m],(i,j),m) for m in range(4))
    for i in range(270):
        for c in range(5):
            ports=range(270+5*c,275+5*c)
            if i in S:
                yield tuple((1,(i,j),1) for j in ports)
            else:
                yield ((5,None,0),)+tuple((1,(i,j),1) for j in ports if j in S)

def tree(edges):
    p=list(range(N))
    def f(a):
        while a!=p[a]:
            p[a]=p[p[a]]
            a=p[a]
        return a
    T=set()
    for i,j in sorted(edges):
        a,b=f(i),f(j)
        if a!=b:
            p[a]=b
            T.add((i,j))
    return T

def normalize(data):
    T=tree(data)
    assert len(T)==N-1
    adj=[[] for _ in range(N)]
    for i,j in T:
        a=next(a for a in range(11) if tuple(sorted((x+a)%11 for x in data[i,j]))==norm(data[i,j]))
        adj[i].append((j,a));adj[j].append((i,-a))
    shift={0:0}; q=[0]
    for i in q:
        for j,a in adj[i]:
            if j not in shift:
                shift[j]=(shift[i]+a)%11
                q.append(j)
    return {(i,j):tuple(sorted((v+shift[j]-shift[i])%11 for v in ss))
            for (i,j),ss in data.items()}

def encode(S,data,b):
    data=normalize(data)
    a,w,den=0,1,1
    for f in factors(S):
        d=sum(opt[0] for opt in f)
        off=0
        for weight,pair,m in f:
            if (m and len(data.get(pair,()))==m) or (not m and not any(p in data for _,p,_ in f)):
                break
            off+=weight
        else:
            raise ValueError('no shape option')
        a=a*d+w*off
        w*=weight
        den*=d
    phase=0;K=1;T=tree(data)
    for p,ss in sorted(data.items()):
        alphabet=(NN if p in T else SS)[len(ss)]
        phase=phase*len(alphabet)+alphabet.index(ss)
        K*=len(alphabet)
    a=a*K+w*phase
    den*=K
    assert (w<<b)>=den
    return ((a<<b)+den-1)//den

def decode(S,code,b):
    fs=list(factors(S))
    den=1
    ds=[sum(opt[0] for opt in f) for f in fs]
    for d in ds:
        den*=d
    width=den<<b
    point=code*den
    shape={}
    for f,d in zip(fs,ds):
        unit=width//d
        bin_no=point//unit
        off=0
        for weight,pair,m in f:
            if bin_no<off+weight:
                break
            off+=weight
        point-=off*unit
        width=unit*weight
        if m:
            shape[pair]=m
    T=tree(shape)
    items=[(p,(NN if p in T else SS)[m]) for p,m in sorted(shape.items())]
    K=1
    for _,A in items:
        K*=len(A)
    phase=point*K//width
    data={}
    for p,A in reversed(items):
        phase,r=divmod(phase,len(A))
        data[p]=A[r]
    assert not phase
    return data

def fixture():
    # Full-size synthetic touching chart, not a degree-57 solution.
    c=25;r,R=bound(c)
    S=set(range(r-5))|{270+5*i for i in range(5)}
    k=r-5
    rng=Random(540)
    data={(0,j):(0,) for j in range(1,270)}
    for i in range(270):
        for c in range(5):
            js=list(range(270+5*c,275+5*c))
            if i<k:
                j=js[i%5] if i>=5 else js[i]
            elif (i+c)%3==0:
                j=js[0]
            else:
                continue
            data[i,j]=(rng.randrange(11),)
    left=R-sum(map(len,data.values()))
    for (i,j) in combinations(range(270),2):
        if left<=0: break
        if i>=k or (i,j) in data or rng.randrange(4)!=0:continue
        m=min(left,1+(rng.randrange(20)==0)+(rng.randrange(36)==0))
        data[i,j]=tuple(sorted(rng.sample(range(11),m)))
        left-=m
    assert left==0 and any(len(x)==3 for x in data.values())
    b=41846
    q=encode(S,data,b)
    recovered=decode(S,q,b)
    assert recovered==normalize(data)
    print('full-size round trip:',len(data),'pair supports and',R,'shifts')

if __name__=='__main__':
    caps()
    fixture()
