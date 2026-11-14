from functools import cache
from random import Random

def check(x):
    if not x: raise AssertionError('marked cut state')

def ordinal(A,B):
    return tuple(x | (((1<<len(B))-1)<<len(A)) for x in A)+tuple(
        x<<len(A) for x in B)

def apex(B,U):
    return tuple(x|((1<<len(B)) if not (U>>i)&1 else 0)
        for i,x in enumerate(B))+(U,)

def profile(B,U):
    n=len(B); full=(1<<n)-1
    reverse=[full^x^(1<<j) for j,x in enumerate(B)]
    @cache
    def p(s,j,end):
        z=s^(1<<j)
        if not z: return 1
        return sum(p(z,k,end) for k in range(n)
                   if (z&(reverse if end else B)[j])>>k&1)
    @cache
    def h(s):
        return sum(p(s,j,False) for j in range(n) if s>>j&1) if s else 1
    S=D=T=F=0
    for s in range(1,full+1):
        sign=1 if s.bit_count()%2 else -1
        q=h(full^s)
        for i in range(n):
            if s>>i&1:
                if U>>i&1:
                    x=q*p(s,i,False); S+=x; D+=sign*x
                else:
                    y=q*p(s,i,True); T+=y; F+=sign*y
    Q=sum(h(s)*h(full^s) for s in range(full+1))
    return h(full),Q,S,D,T,F

def total(B,U):
    A=apex(B,U)
    # Independent unmarked DP.
    @cache
    def path(S,v):
        T=S^(1<<v)
        return sum(path(T,w) for w in range(len(A)) if (T&A[v])>>w&1) if T else 1
    return sum(path((1<<len(A))-1,i) for i in range(len(A)))

def random_t(n,rng):
    B=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if rng.randrange(2): B[i]|=1<<j
            else: B[j]|=1<<i
    return tuple(B)

rng=Random(221)
R5=tuple(1<<((i+3)%5)|1<<((i+4)%5) for i in range(5))
for _ in range(72):
    n,m=rng.randrange(1,5),rng.randrange(1,5)
    B,R=random_t(n,rng),random_t(m,rng)
    U,V=rng.randrange(1<<n),rng.randrange(1<<m)
    h,_,S,d,_,_=profile(B,U)
    a,q,s,e,t,f=profile(R,V)
    N=total(B,U)
    big=ordinal(B,R); bits=U|(V<<n)
    ah,_,AS,Ad,_,_=profile(big,bits)
    check((ah,AS,Ad)==(a*h,q*S+h*s,h*e))
    check(total(big,bits)==a*N+(total(R,V)-a)*h+t*S+f*d)
    check(total(ordinal(B,(0,)),U)==N+S+d)
    check(total(ordinal(B,(2,0)),U)==N+3*S+d)
    check(total(ordinal(B,(2,0)),U|(1<<n))==N+2*h+2*S)

for _ in range(20):
    n,m=rng.randrange(1,4),rng.randrange(1,4)
    L,R=random_t(n,rng),random_t(m,rng)
    U,V=rng.randrange(1<<n),rng.randrange(1<<m)
    base=ordinal(ordinal(L,R5),R)
    s=U|(V<<(n+5))
    check(total(base,s|(9<<n))-total(base,s|(17<<n))==
          2*profile(L,U)[0]*profile(R,V)[0])
print('72 literal cut-state and short-readout controls; 20 independent C5 switches in arbitrary painted contexts.')

# Include complete physical losing-block weights, not just cases
# where the base H is invertible modulo three.
for _ in range(60):
    base=(); mask=0; n=0
    bits=()
    for j in range(rng.randrange(1,4)):
        s=rng.randrange(0,4)
        if s:
            base=ordinal(base,random_t(s,rng))
            n+=s
            bits+=(0,)*(1 if s%2 else 2)
        mask|=1<<n
        base=ordinal(base,(0,)); n+=1; bits+=(1,)
    s=rng.randrange(0,4)
    if s:
        base=ordinal(base,random_t(s,rng))
        bits+=(0,)*(1 if s%2 else 2)
    surrogate=()
    for bit in bits:
        surrogate=ordinal(surrogate,(0,))
    bmp=sum(bit<<i for i,bit in enumerate(bits))
    check((total(base,mask)-profile(base,mask)[0]*total(surrogate,bmp))%3==0)
print('60 literal mod-three singleton-contact skeletons.')
