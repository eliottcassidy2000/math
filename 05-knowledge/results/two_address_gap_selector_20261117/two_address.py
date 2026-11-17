"""Literal oracle for two-port data; subset counts independent of the cover formula."""
from itertools import permutations, product

def profile(T, X, A):
    n=len(T); S=[0]*n; E=[0]*n; M=[[0]*n for _ in T]
    for p in permutations(range(n)):
        bad=[j for j in range(1,n) if not (T[p[j-1]]>>p[j]&1)]
        if not bad:S[p[0]]+=1;E[p[-1]]+=1
        for j in (range(1,n) if not bad else bad if len(bad)==1 else ()):
            M[p[j-1]][p[j]]+=1
    diff=[(X>>j&1)-(A>>j&1) for j in range(n)]
    row=list(map(sum,M)); col=list(map(sum,zip(*M)))
    inner=lambda v: sum(a*b for a,b in zip(diff,v))
    skew=sum(M[i][j]-M[j][i] for i in range(n) for j in range(n) if X>>i&1 and A>>j&1)
    BS,BE=inner(S),inner(E)
    D=skew-BS-BE-inner(row)
    L=-BE-inner(row); R=-BS-inner(col)
    return D,L,R,BS,BE

def attach(T,X,A,left,right):
    n=len(T); p=len(left); q=len(right); z=(1<<(p+n+q))-1
    masks=[]
    for i in range(p):masks.append(z^((1<<(i+1))-1))
    for r in T:masks.append((r<<p) + (z^((1<<(p+n))-1)))
    for j in range(q):masks.append(z^((1<<(p+n+j+1))-1))
    # left letters in birth order (reverse final listing), right in birth order.
    x=(X<<p)+sum((c==0)<<j for j,c in enumerate(reversed(left)))+sum((c==1)<<(p+n+j) for j,c in enumerate(right))
    a=(A<<p)+sum((c==0)<<j for j,c in enumerate(reversed(left)))+sum((c==1)<<(p+n+j) for j,c in enumerate(right))
    return masks,x,a
def host(T,X,A,plus):
    n=len(T)
    a=[r+((X>>i&1)<<n)+((A>>i&1)<<(n+1)) for i,r in enumerate(T)]
    return a+[ ((1<<n)-1)^X | ((not plus)<<(n+1)), ((1<<n)-1)^A | (plus<<n)]
def count(T):
    n=len(T);N=1<<n;dp=[[0]*n for _ in range(N)]
    for j in range(n):dp[1<<j][j]=1
    for m in range(1,N):
        for j,v in enumerate(dp[m]):
            if not v:continue
            w=T[j]&~m
            while w:
                b=w&-w;w-=b;dp[m|b][b.bit_length()-1]+=v
    return sum(dp[-1])

def tri(prefix,suffix):
    t,x,a=attach((2,4,1),2,4,[1]*len(prefix),[0]*len(suffix))
    for i,c in enumerate(prefix+suffix):
        j=i+(3 if i>=len(prefix) else 0)
        x=(x&~(1<<j)) | (int(not c&1)<<j)
        a=(a&~(1<<j)) | (int(not c&2)<<j)
    return t,x,a

if __name__=='__main__':
    from random import Random
    rng=Random(961)
    word=[(),(0,),(1,),(0,1),(1,0),(1,1)]
    cases=0
    for n in range(1,6):
        for _ in range(12):
            T=[0]*n
            for i in range(n):
                for j in range(i+1,n):
                    a,b=(i,j) if rng.randrange(2) else (j,i)
                    T[a]|=1<<b
            X,A=(rng.randrange(1<<n) for _ in range(2))
            D,L,R,BL,BR=profile(T,X,A)
            for l,r in product(word,repeat=2):
                v,w=(sum(b<<i for i,b in enumerate(u)) for u in (l,r))
                N,x,a=attach(T,X,A,l,r)
                low,high=(count(host(N,x,a,+p)) for p in (False,True))
                if high-low!=D+L*v+R*w+BL*(v&1)+BR*(w&1):raise AssertionError((T,l,r,low,high,(D,L,R,BL,BR)))
                cases+=1
    print('independent literal signed-gap checks',cases)

    T,X,A=tri((2,1),())
    assert profile(T,X,A)==(-2,3,-12,3,0)
    for l,r,h,j in [(4,1,3635,3633),(3,1,3443,3441),(8,2,13183,13181),(7,2,20347,20345)]:
        seq=lambda z: tuple(z>>i&1 for i in range(z.bit_length()))
        G,x,a=attach(T,X,A,seq(l),seq(r))
        assert (count(host(G,x,a,False)),count(host(G,x,a,True)))==(h,j)
    T,X,A=(40,1,3,6,15,30),46,23
    assert profile(T,X,A)[0]==2 and profile(T,X,A)[2]==0
    for k in range(3):
        G,x,a = attach(T,X,A,(), (0,)*k+(1,) if k else ())
        h,j=(count(host(G,x,a,b)) for b in (False,True))
        assert h== (2610*2**(k-1)+345 if k else 509) and j==h+2
    print("committed absolute centers and one-sided tile checked")
