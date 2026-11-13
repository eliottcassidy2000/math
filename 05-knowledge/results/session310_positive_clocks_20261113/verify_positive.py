from collections import Counter
from itertools import permutations, combinations

def req(ok):
    if not ok:
        raise AssertionError("positive-clock certificate")

def spectrum(rows):
    n=len(rows)
    # Literal independent Hamilton path counts on every induced subset.
    S=list(range(1<<n));H=[1]*(1<<n)
    for s in S[1:]:
        H[s]=sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:]))
                 for p in permutations(i for i in range(n) if s>>i&1))
    # Literal odd directed cycles, anchored at their least vertex.
    C=[]
    for s in S[1:]:
        if s.bit_count()<3 or s.bit_count()%2==0:
            continue
        i=(s&-s).bit_length()-1
        for p in permutations(j for j in range(i+1,n) if s>>j&1):
            q=(i,)+p+(i,)
            if all(rows[u]>>v&1 for u,v in zip(q,q[1:])):
                C.append(s)
    P=Counter()
    def add(j,covered,k):
        P[n-covered.bit_count()+k]+=2**k
        for i in range(j,len(C)):
            if not covered&C[i]:
                add(i+1,covered|C[i],k+1)
    add(0,0,0)
    def colors(s,t):
        if not t:
            return int(not s)
        val=0;b=s
        while True:
            val+=H[b]*colors(s^b,t-1)
            if not b:
                return val
            b=(b-1)&s
    for t in range(1,5):
        req(sum(c*t**e for e,c in P.items())==colors(S[-1],t))
    return H[-1],colors(S[-1],2)

def decode(n,x):
    rows=[0]*n
    for i,j in combinations(range(n),2):
        if x&1:
            rows[i]|=1<<j
        else:
            rows[j]|=1<<i
        x>>=1
    return rows

tested=0
for n in range(1,6):
    for code in range(1<<(n*(n-1)//2)):
        spectrum(decode(n,code));tested+=1
h,p=spectrum([4,9,82,21,35,143,187,31])
j,q=spectrum([6,12,48,21,131,219,159,15])
e=6*p*j+h*j-6*h*q
Y=36*p*q+6*h*q+18*p*j-h*j
req((h,p,j,q,e,Y)==(181,3424,179,3416,-1,435779329))
req(Y+e==12*p*(3*q+2*j) and ((Y+e)&-(Y+e))==256)
print(tested,"literal positive expansions and determinant fixture passed")
