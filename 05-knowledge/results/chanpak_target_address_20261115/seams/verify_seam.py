from itertools import combinations
from random import Random

def tournament(n,code):
    a=[0]*n
    for i,(x,y) in enumerate(combinations(range(n),2)):
        if code>>i&1:a[x]|=1<<y
        else:a[y]|=1<<x
    return a

def tt(n):
    return tournament(n,(1<<(n*(n-1)//2))-1)

def join(*B):
    a=[]; size=sum(map(len,B))
    for b in B:
        left=len(a); right=left+len(b)
        a.extend((row<<left)|(((1<<size)-1)^((1<<right)-1)) for row in b)
    return a

def reverse(a,x,y):
    a=list(a);a[x]^=1<<y;a[y]^=1<<x
    return a

def dual(a):
    f=(1<<len(a))-1
    return [f^x^(1<<i) for i,x in enumerate(a)]

def stats(a):
    n=len(a);f=(1<<n)-1
    D=[[0]*n for _ in range(f+1)]
    for i in range(n):D[1<<i][i]=1
    for m in range(1,f):
        for j,c in enumerate(D[m]):
            if not c:continue
            x=a[j]&~m
            while x:
                p=x&-x;x-=p
                D[m|p][p.bit_length()-1]+=c
    h=[sum(z) for z in D];h[0]=1
    return h,D

def profile(a):
    h,D=stats(a);_,G=stats(dual(a));f=len(h)-1
    E=[sum(G[m][i]*h[f^m] for m in range(1,f+1) if not m.bit_count()%2) for i in range(len(a))]
    O=[sum(G[m][i]*h[f^m] for m in range(1,f+1) if m.bit_count()%2) for i in range(len(a))]
    assert [o-e for o,e in zip(O,E)]==D[-1]
    return h[-1],sum(h[m]*h[f^m] for m in range(f+1)),E,O,D[-1]

def strong(a):
    f=(1<<len(a))-1
    for i in range(len(a)):
        m=1<<i
        while True:
            z=m
            for j,r in enumerate(a):
                if m>>j&1:z|=r
            if m==z:break
            m=z
        if m!=f:return False
    return True

small={i:set() for i in [1,2]}
strong4=set();total=0
for n in range(1,6):
    for c in range(1<<(n*(n-1)//2)):
        a=tournament(n,c)
        h,q,E,O,end=profile(a);total+=n
        for i in range(n):
            if a[i]:
                assert 2*end[i]<=h-1
                assert h<=4*E[i]-1
            else:
                assert E[i]==0 and end[i]==h
            if E[i] in small:
                small[E[i]].add(h)
        if n==4 and strong(a):
            strong4.add((h,tuple(sorted(E))))
        if n==5 and strong(a):
            assert 2 not in E
assert small=={1:{1,3},2:{1,5}}
assert strong4=={(5,(2,2,3,4))}
print(total,'rooted endpoint identities; exact E=1,2 controls',small)
rng=Random(20261115)
for _ in range(120):
    n,m=rng.randint(1,5),rng.randint(1,5)
    a,b=tournament(n,rng.getrandbits(n*(n-1)//2)),tournament(m,rng.getrandbits(m*(m-1)//2))
    i,j=rng.randrange(n),rng.randrange(m)
    ha,_,E,O,d=profile(a);hb,_,F,P,e=profile(dual(b))
    T=reverse(join(a,b),i,n+j)
    actual=stats(T)[0][-1]
    inc=(E[i]+O[i])*(F[j]+P[j])-d[i]*e[j]
    assert actual==ha*hb+inc
    if not (a[i] or dual(b)[j]):assert inc==0
    else:assert 2*inc>ha*hb
    if a[i] and dual(b)[j]:assert 4*inc>3*ha*hb
print('120 literal seam controls')
for n in range(4,12):
    a=reverse(tt(n),0,n-1);h,q,E,O,end=profile(a)
    x=n-2
    assert h==2**(n-2)+1 and E[x]==1+2**(n-4)
    assert E[x]+O[x]==h+1 and end[x]==(h-1)//2
print('sharpness controls through order 11')
for _ in range(60):
    n,m=rng.randint(1,5),rng.randint(1,3)
    d,b=tournament(n,rng.getrandbits(n*(n-1)//2)),tournament(m,rng.getrandbits(m*(m-1)//2))
    h,q,*_=profile(d); hb,p,*_=profile(b)
    a=reverse(join(tt(1),d,tt(2)),0,n+2)
    ha,_,E,O,end=profile(a); assert ha==h+2*q and E[-1]+O[-1]==ha and E[-2]+O[-2]==ha+h
    for root,expected in ((len(a)-1,ha*(hb+p)),(len(a)-2,ha*(hb+p)+h*p)):
        X=reverse(join(a,b,tt(1)),root,len(a)+len(b))
        assert strong(X) and stats(X)[0][-1]==expected
print('120 strong two-arm literal controls')

def decode(N,h,q,b,p):
    k=(N+h*b)//2
    out=[]
    for d in range(1,k+1):
        if k%d:continue
        y=k//d
        if d<=h or y<=b or (d-h)%q or (y-b)%p:continue
        z,w=(d-h)//q,(y-b)//p
        if z&(z-1)==0 and w&(w-1)==0:out.append((z.bit_length()-1,w.bit_length()-1))
    return out

def build(d,b):
    n=len(d)
    a=reverse(join(tt(1),d,tt(2)),0,n+2)
    return reverse(join(a,b,tt(1)),len(a)-2,len(a)+len(b))
for d,b in [(tt(1),tt(2)),(tournament(3,2),tt(1)),(tt(1),tournament(3,2))]:
    h,q,*_=profile(d);v,p,*_=profile(b)
    for i in range(3):
        for j in range(3):
            N=2*(h+q*2**i)*(v+p*2**j)-h*v
            witnessed=decode(N,h,q,v,p)
            assert (i,j) in witnessed
            for r,s in witnessed:
                assert stats(build(join(d,tt(r)),join(b,tt(s))))[0][-1]==N
print('27 factoring-decoder seeds, with every returned witness counted')

from itertools import permutations
def gap_poly(w,a,v,rooted=False):
    if not w:return [1,1]
    gaps=[(1,0 if rooted else a[v]>>w[0]&1)]
    gaps += [(a[i]>>j&1,(a[i]>>v&1)*(a[v]>>j&1)) for i,j in zip(w,w[1:])]
    gaps += [(1,a[w[-1]]>>v&1)]
    p=[1]
    for z,one in gaps:
        q=[0]*(len(p)+1)
        for k,c in enumerate(p):
            q[k]+=z*c;q[k+1]+=one*c
        p=q
    return p

def compressed(a):
    n=len(a);v=n-1;H=[0]*(n+2);S=[[0]*(n+2) for _ in range(v)]
    for w in permutations(range(v)):
        for k,c in enumerate(gap_poly(w,a,v)):H[k]+=c
        for k in range(1,v+1):
            for i,x in enumerate(gap_poly(w[:k],a,v,True)):
                for j,y in enumerate(gap_poly(w[k:],a,v)):
                    if x and y:S[w[0]][i+j]+=x*y
    def trim(w):
        return tuple(w[1:max(k for k,z in enumerate(w) if z)+1])
    return trim(H),[trim(p) for p in S]
def inflate(a,D):
    n=len(a);k=len(D);v=n-1;full=(1<<k)-1
    return [((r&~(1<<v)) | ((r>>v&1)*full<<v)) for r in a[:-1]]+[(a[-1]|r<<v) for r in D]
def three(a):
    h,_=stats(a);n=len(a);f=(1<<n)-1
    Q=sum(h[x]*h[f^x] for x in range(f+1))
    R=0
    for x in range(f+1):
        r=f^x;y=r
        while True:
            R+=h[x]*h[y]*h[r^y]
            if y==0:break
            y=(y-1)&r
    return h[-1],Q-2*h[-1],R-3*Q+3*h[-1]

for n,code,lo,hi,H,L,U in [
        (5,108,3,2,(11,6),(12,6),(13,6)),
        (6,2330,4,3,(43,54,16),(74,66,16),(75,66,16))]:
    a=tournament(n,code)
    assert strong(a)
    hh,ss=compressed(a)
    assert (hh,ss[lo],ss[hi])==(H,L,U)
    for d in [tt(k) for k in range(1,6)]+[tournament(3,2),tournament(4,17),tournament(5,108)]:
        P=three(d);A=inflate(a,d)
        evalp=lambda z:sum(c*p for c,p in zip(z,P))
        hh,_,E,O,_=profile(A)
        assert hh==evalp(H) and E[lo]+O[lo]==evalp(L) and E[hi]+O[hi]==evalp(U)
        for i,s in [(lo,L),(hi,U)]:
            T=reverse(join(A,tt(2)),i,len(A)+1)
            assert strong(T) and stats(T)[0][-1]==hh+2*evalp(s)
print('two compressed visit tables and 32 direct inflated switch controls')
