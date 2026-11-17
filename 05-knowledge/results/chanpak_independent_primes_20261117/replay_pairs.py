from math import isqrt
from slot_profiles import external, ordinal, profile, strong
from slot_profiles import end_profiles

def tt(k):
    return tuple((1<<k)-(2<<i) for i in range(k))

def coerce(A,B):
    n=len(A)-1
    v=[]
    for i in range(n):
        v.append((A[i]&((1<<n)-1))|(((1<<len(B))-1)<<n if A[i]>>n&1 else 0))
    return tuple(v)+tuple((a<<n)|A[n] for a in B)

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

R=(20,57,50,21,32,73,31)
if sum(a>>6&1 for a in R)!=1 or not strong(R) or not strong(external(R,42,1)):
    raise ValueError("unimodular exterior")
for B in (tt(1),tt(2),(2,4,1),(52,1,10,35,14,22)):
    h,P,*_=profile(B,0)
    p,q=44*P-13*h,71*P-21*h
    R1=coerce(R,B)
    if (end_profiles(R1)[0][-1],end_profiles(external(R1,42,1))[0][-1])!=(p,7*q):
        raise ValueError((p,q))
    if (21*p-13*q,71*p-44*q)!=(P,h):
        raise ValueError("inverse")
    print("Unimodular",(h,P),p,q,prime(p),prime(q))

B=(48,89,35,85,4,26,53)
A=ordinal(tt(3),B)
n=len(A)
F=tuple(v|(1<<n) for v in A)+(1<<(n+1),(1<<n)-1)
g,_,W,_,E,_,K=profile(F,5)
D=E-K
p=3*W+2*g-D
q=p+4*D//7
if D%7 or (p,q)!=(731761,732073) or not all(map(prime,(p,q))):
    raise ValueError((D,p,q))
print("Signed",D,p,q,7*q)

p,q=1747,2819
if not all(map(prime,(p,q))) or (21*p-13*q,71*p-44*q)!=(40,1):
    raise ValueError("ghost target")
