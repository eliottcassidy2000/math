from functools import cache
from physical_prefix import row,greedy,prime

def table(R):
    n=len(R)
    @cache
    def ends(mask,v):
        u=mask^(1<<v)
        if not u:return 1
        return sum(ends(u,i) for i,x in enumerate(R) if u>>i&1 and x>>v&1)
    A=[1]+[sum(ends(m,v) for v in range(n) if m>>v&1) for m in range(1,1<<n)]
    return A

def h_f2(R):
    h=table(R)
    return h[-1],sum(x*y for x,y in zip(h,h[::-1]))

def rooted_profile(R,root):
    path=table(R)
    n=len(R)
    full=(1<<n)-1
    @cache
    def starts(S,v):
        X=S^(1<<v)
        if not X:return 1
        return sum(starts(X,j) for j in range(n) if X>>j&1 and R[v]>>j&1)
    Rplus=Rminus=0
    for S in range(1,1<<n):
        if S>>root&1:
            c=starts(S,root)*path[full^S]
            Rplus+=c
            Rminus+=(-1)**S.bit_count()*c
    return path[-1],Rplus,Rminus

def bridge(R,root,r):
    """Add r transitive followers and a last terminal beating only root."""
    assert r>=0 and 0<=root<len(R)
    n=len(R)
    m=n+r+1
    N=((1<<m)-1)^((1<<n)-1)
    A=[u|N for u in R]
    A[root]^=1<<(m-1)
    for j in range(n,m-1):
        A.append(((1<<m)-1)^((1<<(j+1))-1))
    A.append(1<<root)
    return A

def source(C):
    # Prepend one root that dominates every vertex in the arbitrary core.
    return [((1<<(1+len(C)))-1)^1]+[x<<1 for x in C]

def pad_pair(C,k):
    h,L=h_f2(C)
    u,v=row(k); A=u+2*v
    if h*2>A or L>3*A:return None
    r=1
    while (L<<r)<3*A:r+=1
    b=L<<r
    # r=1 can have equality only if L=3A; in fact L is even, A odd.
    assert b<6*A
    p=u*h+v*b
    assert greedy(p)==(k,h,b)
    return p,h,b,r

def test_examples():
    C0=[150,172,472,433,226,69,11,352,115]
    from random import Random
    R=Random(271)
    C=[C0]
    for _ in range(24):
        n=6
        A=[0]*n
        for i in range(n):
            for j in range(i+1,n):
                if R.randrange(2):A[i]|=1<<j
                else:A[j]|=1<<i
        c3=sum(A[i]>>j&1 and A[j]>>k&1 and A[k]>>i&1 or
               A[i]>>k&1 and A[k]>>j&1 and A[j]>>i&1
               for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        if c3>n-2:C.append(A)
    examples=[]
    for k in (1,2,3):
        for c in C:
            z=pad_pair(c,k)
            if z is None:continue
            p,a,b,r=z
            if not prime(p):continue
            s=bridge(source(c),0,r)
            assert table(s)[-1]==a+b and table(s[:-1])[((1<<(len(s)-1))-1)]==a
            assert rooted_profile(s,len(s)-1)==(a+b,a+b,-a)
            examples.append((k,p,a,b,r,c,s))
            break
    return examples

if __name__=='__main__':
    for answer in test_examples():print(answer)

if __name__=='__main__':
    for n in range(1,5):
        pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
        for bcode in range(1<<len(pairs)):
            R=[0]*n
            for d,(i,j) in enumerate(pairs):
                if bcode>>d&1:R[i]|=1<<j
                else:R[j]|=1<<i
            for v in range(n):
                h,a,b=rooted_profile(R,v)
                for r in (0,1,2):
                    W=bridge(R,v,r)
                    D=rooted_profile(W,len(W)-1)
                    want=(h+a+b,h+a,-h-b) if r==0 else (h+(a<<r),h+(a<<r),-h)
                    assert D==want,(R,v,r,D,want)
    C=(2,24,3,53,37,7)
    assert h_f2(C)==(23,296)
    W=bridge(source(C),0,4)
    assert W==[2046,3972,4016,3974,4074,4042,3982,3840,3584,3072,2048,1]
    assert rooted_profile(W,11)==(4759,4759,-23)
    assert greedy(1512281)==(1,23,4736)
    assert 47*(23+2*4736)+56*4*(23+4736)==1512281 and prime(1512281)
    print('all rooted tournaments through order 4, exact positive and zero bridge interfaces')
