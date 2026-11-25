"""Exact integer enclosures for the native Moore path-count interval."""
from fractions import Fraction as Q
from math import ceil
from itertools import combinations

S=10**18
terms=18
def atanh2(a,b):
    # floor each rational summand after scaling, and bound the entire tail
    d=sum(2*S*a**(2*i+1)//((2*i+1)*b**(2*i+1)) for i in range(terms))
    # z=a/b <= 1/3, and ignored sum <= 3*(1/3)^(2*terms+1)
    return d,d+terms+ceil(Q(S,3**(2*terms)))
L2,U2=atanh2(1,3)
def iln(k):
    j=k.bit_length()-1; b=1<<j
    lo,hi=atanh2(k-b,k+b)
    return j*L2+lo,j*U2+hi

n=3250
lo=[0]*(n+1);hi=lo.copy()
for i in range(1,n+1):
    a,b=iln(i);lo[i]=lo[i-1]+a;hi[i]=hi[i-1]+b
def cap(c):
    V=Q(hi[n],n)  # the universal apex has row sum n
    for i in range(n):
        x=n-i-c*Q(n-2*i-1,n)
        j=x.numerator//x.denominator
        V+=(j+1-x)*Q(hi[j],j)+(x-j)*Q(hi[j+1],j+1)
    return V/L2
U=cap(Q(65));W=cap(65+Q(n,n-1))
mu=Q(lo[n],U2)-(n-1)
m=ceil(Q(1625*520,538));D=58;r=4*D+3
assert (m,D,r)==(1571,58,235)
lower=Q(lo[m-D]-hi[r-D],U2)
assert lower>12734 and U<28860 and W<28863
assert mu-U>1125 and mu-W>1122
print("balanced matching >= ",m)
print("lower log2 bound",float(lower),"for",m-D,"!/",r-D,"!")
print("sparse log2 upper bound",float(U),"and gap >",float(mu-U))
print("mate log2 upper bound",float(W),"and gap >",float(mu-W))

# Small falsification control: equal triangle parity does not force equal
# Hamilton-path parity modulo 4 even among native Moore complement spines.
p=[0,2,4,1,9,8,7,3,5,6]
def edge(a,b):
    if a//5!=b//5:return a%5==b%5
    d=(a-b)%5
    return d in ((1,4) if a<5 else (2,3))
assert not any(edge(a,b) for a,b in zip(p,p[1:]))
def chart(mate):
    N=len(p);out=[0]*N
    for i,j in combinations(range(N),2):
        rev=edge(p[i],p[j])^(mate and j==i+1) # mate after converse
        if rev: out[j]|=1<<i
        else: out[i]|=1<<j
    return out
def report(a):
    N=len(a)
    P=[[0]*N for _ in range(1<<N)]
    for i in range(N):P[1<<i][i]=1
    for s in range(1,1<<N):
        for j in range(N):
            if s>>j&1:
                u=s^(1<<j)
                if u:P[s][j]=sum(P[u][i] for i in range(N) if u>>i&1 and a[i]>>j&1)
    C=N*(N-1)*(N-2)//6-sum(d*(d-1)//2 for d in (x.bit_count() for x in a))
    return sum(P[-1]),C
A,B=(report(chart(m)) for m in [False,True])
assert A[0]==5625 and B[0]==5627 and B[1]-A[1]==2
print("Petersen native modulo-4 control: (paths,triangles)",A,B)
