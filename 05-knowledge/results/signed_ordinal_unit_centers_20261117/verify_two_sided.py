from itertools import product
from functools import cache

def words(n):
    yield ()
    for k in range(n):
        for w in product((0,1),repeat=k):
            yield w+(1,)

def host(l,r,flip):
    L=list(range(len(l)));base=len(l)
    b,c,y,x,a=range(base,base+5); R=list(range(base+5,base+5+len(r)))
    n=len(l)+len(r)+5; A=[0]*n
    def add(i,j): A[i]|=1<<j
    for T in (L,R):
        for i,v in enumerate(T):
            for w in T[i+1:]: add(v,w)
    for i in L:
        for j in [b,c,y]+R: add(i,j)
    for i in [b,c,y]:
        for j in R: add(i,j)
    for i,j in [(b,c),(c,y),(y,b),(x,b),(a,b),(c,x),(a,c),(x,y),(y,a)]:add(i,j)
    add(a,x) if flip else add(x,a)
    for T,W in ((L,l),(R,r)):
        for i,k in zip(T,W):
            for p in (x,a):
                if k==(T is L):add(p,i)
                else:add(i,p)
    assert all((A[i]>>j&1)+(A[j]>>i&1)==1 for i in range(n)for j in range(i))
    return A,L,R,(b,c,y,x,a)

def count(A):
    @cache
    def d(S,v):
        if not S:return 1
        M=S&A[v];out=0
        while M:
            b=M&-M;M-=b
            out+=d(S^b,b.bit_length()-1)
        return out
    return sum(d(((1<<len(A))-1)^(1<<v),v)for v in range(len(A)))

def fillings(w,T,left):
    W=list(zip(w,T))
    if left:W.reverse()
    i,k=[j for j,(x,_)in enumerate(W)if x][-2:]
    options=[(0,1,2)if j<i else (0,)if j==i else (1,2)if j<k else (1,)for j in range(k+1)]
    for ii in product(*options):
        parts=[[]for _ in range(3)]
        for ind,(_,v)in zip(ii,W):parts[ind].append(v)
        if left:parts=[z[::-1]for z in parts[::-1]]
        yield tuple(parts)

cnt=0
for LL,RR in product(words(3),repeat=2):
    l=LL[::-1];r=RR;hs=[]
    for f in (0,1):
        A,L,R,(b,c,y,x,a)=host(l,r,f);h=count(A);hs.append(h)
        def check(z):
            global cnt
            assert sorted(z)==list(range(len(A)))
            assert all(A[i]>>j&1 for i,j in zip(z,z[1:]))
            cnt+=1
        if sum(l)>1 and sum(r)>1:
            paths=set()
            for (l0,l1,l2),(r0,r1,r2) in product(fillings(l,L,True),fillings(r,R,False)):
                z=tuple(l0+[b]+r0+[x]+l1+[c]+r1+[a]+l2+[y]+r2)
                check(z);paths.add(z)
            assert len(paths)<=h
        if sum(l)>1:
            for l0,l1,l2 in fillings(l,L,True):check(l0+[c,x]+l1+[y,a]+l2+[b]+R)
        if sum(r)>1:
            for r0,r1,r2 in fillings(r,R,False):check(L+[c]+r0+[x,y]+r1+[a,b]+r2)
        if sum(l)==1:
            for group in product((0,1),repeat=len(L)-1):
                l0=[i for i,g in zip(L[1:],group)if g==0];l1=[L[0]]+[i for i,g in zip(L[1:],group)if g==1]
                check(l0+[c,x]+l1+[y,a,b]+R)
        if sum(r)==1:
            for group in product((0,1),repeat=len(R)-1):
                r0=[i for i,g in zip(R[:-1],group)if g==0]+[R[-1]];r1=[i for i,g in zip(R[:-1],group)if g==1]
                check(L+[c,x,y]+r0+[a,b]+r1)
    assert hs[1]-hs[0]==2
print(cnt,"literal traces checked on 64 two-sided addresses")
