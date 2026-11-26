"""Exact cut bounds for Moore graph state sums; no degree-57 fixture assumed."""
from fractions import Fraction as Q
from itertools import combinations
from math import ceil, sqrt

def refin(n,d,r,s,t):
    """s is cut size; [t,r] is the centered spectral interval."""
    k=Q((d-r)*s,n)
    a=k-int(k)
    bound=(d-r)*Q(s*(n-s),n)+Q(n,r-t)*a*(1-a)
    # The cut and d*s have the same parity.
    raw=ceil(bound)
    return raw+((d*s-raw)&1),bound

def report(name,n,d,r,t,clique=False):
    third=min((refin(n,d,r,s,t)[0],s)
              for s in range(n//3+1,2*n//3+1))
    half=refin(n,d,r,n//2,t)[0]
    plain=ceil(min((d-r)*Q(s*(n-s),n)
              for s in range(n//3+1,2*n//3+1)))
    print(name,n,d,"plain/corrected balanced",plain,third,"half",half,
          "price exponents",round(sqrt(third[0]),6),round(sqrt(half),6))
    return third[0],half

assert report("Petersen",10,3,1,-2)==(5,5)
assert report("Hoffman-Singleton",50,7,2,-3)
full=report("conditional full",3250,57,7,-8)
leaf=report("conditional leaf",3192,56,7,-8)
frozen=report("conditional exterior",3024,53,7,-8)

def make_petersen():
    V=list(combinations(range(5),2))
    return [set(j for j,v in enumerate(V) if not set(v).intersection(u)) for u in V]
def make_hs():
    A=[set() for _ in range(50)]
    def add(u,v): A[u].add(v);A[v].add(u)
    for i in range(5):
        for x in range(5):
            add(5*i+x,5*i+(x+1)%5)
            add(25+5*i+x,25+5*i+(x+2)%5)
            for j in range(5):add(5*i+x,25+5*j+(x+i*j)%5)
    return A

def check(A,S,r,t):
    n,d,s=len(A),len(A[0]),len(S)
    y=Q((d-r)*s,n)
    z=[len(v&S)-r*(i in S)-y for i,v in enumerate(A)]
    b=sum(len(A[i]-S) for i in S)
    assert sum(z)==0
    assert sum(a*a for a in z)==(r-t)*(b-(d-r)*Q(s*(n-s),n))
    assert b>=refin(n,d,r,s,t)[0]
    return b

P=make_petersen()
exact=[min(check(P,set(S),1,-2) for S in combinations(range(10),s))
       for s in range(11)]
print("Petersen exact profile",exact)
# Check the supplied HS fixture satisfies the full Moore equation, then
# verify the residual-square identity on many independently chosen cuts.
H=make_hs()
for i in range(50):
    assert len(H[i])==7
    for j in range(i):
        assert len(H[i]&H[j])+(j in H[i])==1
from random import Random
rng=Random(57)
for s in range(1,50):
    for _ in range(10):check(H,set(rng.sample(range(50),s)),2,-3)
print("Hoffman-Singleton identity: 490 exact controls")

for k in (3,7,57):
    d=k-1
    half=d*(k//2)*(k-k//2)
    low=min(d*j*(k-j) for j in range(k//3+1,2*k//3+1))
    print("branch unions",k,"balanced",low,"half",half)


def branch_refin(k,L,a,r,t):
    # k branches of L sheets; compare all size-a cuts.
    n=k*L; d=k-1
    l,u=divmod(a,k)
    def piece(i):
        rho=Q((a-(r+1)*i)%L,L)
        return Q(r+1,L)*i*i+Q(L,r-t)*rho*(1-rho)
    bound=(d-r)*Q(a*(n-a),n)-Q(r+1,n)*a*a
    bound+=(k-u)*piece(l)+u*piece(l+1)
    raw=ceil(bound)
    return raw+((d*a-raw)&1)
for k,L,r,t in [(7,6,2,-3),(57,56,7,-8),(54,56,7,-8)]:
    n=k*L
    bmin=min((max(branch_refin(k,L,s,r,t),refin(n,k-1,r,s,t)[0]),s)
             for s in range(n//3+1,2*n//3+1))
    bhalf=max(branch_refin(k,L,n//2,r,t),refin(n,k-1,r,n//2,t)[0])
    print("branch-aware",k,L,bmin,bhalf)

def separator(n,d,r,t):
    half=Q(n,2)
    # d-midpoint vs radius of the centered spectral interval.
    mid,rad=Q(r+t,2),Q(r-t,2)
    K=n*rad**2/((d-mid)**2-rad**2)
    assert K<half
    return K,ceil(half**2/(half+K))
for args in [(3250,57,7,-8),(3192,56,7,-8),(3024,53,7,-8)]:
    K,s=separator(*args)
    print("separator",args[0],"anticomplete constant",K,"bag size",s,"treewidth",s-1)

def chart_check(A,k,L,r,t,nfree):
    roots=sorted(A[0])
    branches=[set(A[i])-{0} for i in roots][:nfree]
    assert len(branches)==nfree and all(len(B)==L for B in branches)
    X=set().union(*branches)
    assert len(X)==nfree*L
    for B in branches:
        for v in B:
            assert all(len(A[v]&C)==(B is not C) for C in branches)
    for _ in range(250):
        a=rng.randrange(len(X)+1)
        S=set(rng.sample(sorted(X),a))
        ii=[len(S & B) for B in branches]
        lam=nfree-1
        base=(lam-r)*Q(a*(len(X)-a),len(X))+(r+1)*sum(Q((v-Q(a,nfree))**2,L) for v in ii)
        z=[len(A[v]&S)-r*(v in S)-Q(a-(r+1)*h,L)
           for h,B in zip(ii,branches) for v in B]
        penalty=Q(0)
        for h in ii:
            frac=Q((a-(r+1)*h)%L,L)
            penalty+=Q(L,r-t)*frac*(1-frac)
        cut=sum(len((A[v]&X)-S) for v in S)
        assert (r-t)*(cut-base)>=sum(x*x for x in z)>=(r-t)*penalty
        assert cut>=branch_refin(nfree,L,a,r,t)
chart_check(H,7,6,2,-3,7)
chart_check(H,7,6,2,-3,4)
print("Hoffman-Singleton branch profiles: 250 full + 250 frozen exact controls")
