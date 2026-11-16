from functools import cache

SEEDS=((2,4,25,19,35,15),(6,4,24,19,35,15))

def rev(b,n):
    return int(f'{b:0{n}b}'[::-1],2) if n else 0

def interaction(b):
    sites=[i for i in range(b.bit_length()) if b>>i&1]
    return sum(2**(j-i-1)+(j==i+1) for j in sites for i in sites if i<j)

def score(b,n,Z,I):
    return I+Z*(2**(n+1)-rev(b,n))+9*(2*b-interaction(b))

def parameters(u,d):
    V=2**u
    if u<2 or d%2 or abs(d)>V:
        raise ValueError('need u>=2 and an even gap |d|<=2**u')
    for m in range(6):
        for s in range(2):
            A=62*2**m*V
            K,r=divmod(A-d+2*s,9)
            if not r:
                assert K%2==0 and V<K<8*2**m*V
                return m,s,K
    raise AssertionError('no phase')

def addresses(u,d,U):
    m,s,K=parameters(u,d)
    V=2**u; i=u+m+3; n=i+u+1
    assert 0<=U<V and U%2==0
    X=K-U
    b=(U<<(i+1))|rev(X,i)
    return m,s,n,b,b|(1<<i)

def literal(b,n,m,s):
    a=SEEDS[s]
    for e in [0]*(1+m)+[int(c) for c in f'{b:0{n}b}']+[0]:
        k=len(a)
        a=tuple(x|((int(j!=3 or e))<<k) for j,x in enumerate(a))+((0 if e else 1<<3),)
    return a

def count(a):
    @cache
    def dp(mask,j):
        if mask==1<<j:return 1
        other=mask^(1<<j)
        return sum(dp(other,k) for k in range(len(a)) if (other&a[j])>>k&1)
    F=(1<<len(a))-1
    return sum(dp(F,j) for j in range(len(a)))

if __name__=='__main__':
    for n in range(1,9):
        for b in range(2**n):
            assert score(2*b,n+1,82,-11)-score(b,n,82,-11)==164*2**n+18*b
    total=0
    for u in range(2,9):
        V=2**u
        for d in range(-V,V+1,2):
            prev=None
            for U in range(0,V,2):
                m,s,n,b,c=addresses(u,d,U)
                Z=82*2**m
                h=score(b,n,Z,-11)
                assert score(c,n,Z,-11+2*s)-h==d
                if prev is not None:
                    assert 381*2**m*V <= h-prev <= 1211*2**m*V
                if u==2 and U==0 and m<=1:
                    a=literal(b,n,m,0); A=literal(c,n,m,s)
                    assert count(a)==h and count(A)==h+d
                    assert sum((x^y).bit_count() for x,y in zip(a,A))//2==1+s
                prev=h
                total+=1
    print(total,'physical-rule comparisons passed')
