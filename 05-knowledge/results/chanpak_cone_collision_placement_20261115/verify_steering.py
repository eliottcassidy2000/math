from math import isqrt

def p(n):
    ans=0
    while n:
        u,b=divmod(n,2)
        if not b: ans+=u+(u&1)
        n=u
    return ans

def rev(v,r):
    z=0
    for _ in range(r):
        z=2*z+(v&1)
        v//=2
    return z

def f(x,y):
    return 100*x*y+89*x+11*(p(x)+p(y)+1+(x&1))

def parameters(n):
    base=isqrt(max(0,(n-11)//100))+1
    q=1<<((base.bit_length()-1)//3)
    assert q>=2
    u=(2*base+q+1+2*q-1)//(2*q)
    s=2*q*u+q-1
    return q,u,s

def chart(n):
    q,u,s=parameters(n)
    K=2*u+2*q+1
    C=44+100*s*s+178*s+44*u*(q-1)
    def smooth(z):
        return C-100*z*z+178*z
    def lo(z):
        return smooth(z)-22*z-44*(K+1)
    def hi(z):
        return smooth(z)+22*z+44*(K+1)
    def literal(z):
        return f((s+z)//2,(s-z)//2)
    J=(s-3)//2
    assert lo(1)>4*n and hi(2*J+1)<4*n
    def first_fall(g):
        a,b=0,J
        while b-a>1:
            mid=(a+b)//2
            if g(2*mid+1)>=4*n: a=mid
            else: b=mid
        return b
    def first_weak(g):
        a,b=0,J
        while b-a>1:
            mid=(a+b)//2
            if g(2*mid+1)>4*n: a=mid
            else: b=mid
        return b
    raised,last=first_weak(lo),first_fall(hi)-1
    hits=[]
    for j in range(raised,last+1):
        z=2*j+1
        if literal(z)==n: hits.append(((s+z)//2,(s-z)//2))
    upper=literal(2*(raised-1)+1)
    lower=literal(2*(last+1)+1)
    assert lower<n<upper
    return lower,upper,hits,max(0,last-raised+1)

def test():
    from random import Random
    rng=Random(937)
    for _ in range(1000):
        r=rng.randrange(1,18); q=1<<r; u=rng.randrange(2,100000)
        j=rng.randrange(-u,u+1); v=rng.randrange(q)
        U,V=u+j,u-j
        x,y=q*U+v,q*V+q-1-v
        s=2*q*u+q-1
        rv=rev(v,r)
        rhs=(p(U)+p(V)+p(v)+p(q-1-v)+u*(q-1)
             +j*(q-1-2*rv)+(U&1))
        assert p(x)+p(y)==rhs
    max_tests=0
    n_min=100*1000**2+1
    for n in range(n_min,n_min+999,18):
        q,u,s=parameters(n)
        lower,upper,hits,t=chart(n)
        assert hits==[(x,s-x) for x in range((s+1)//2,s) if f(x,s-x)==n]
        max_tests=max(max_tests,t)
    for k in range(32,257,8):
        for _ in range(10):
            n=rng.getrandbits(k)|1<<(k-1)|1
            lower,upper,hits,t=chart(n)
            max_tests=max(max_tests,t)
    print("1000 exact split identities, brute fixed-sum inverse controls, and tests to 256 bits; max candidate tests:",max_tests)

if __name__=="__main__": test()
