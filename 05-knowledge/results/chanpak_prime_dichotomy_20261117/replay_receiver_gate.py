from slot_profiles import profile, end_profiles, ordinal
from itertools import product

def tt(k):
    return tuple((1<<k)-(2<<i) for i in range(k))

def slots(a):
    H,_=end_profiles(a); f=len(H)-1
    h=H[f]; p=sum(H[m]*H[f^m] for m in range(f+1)); q=0
    for m in range(f+1):
        c=f^m; s=c
        while True:
            q+=H[m]*H[s]*H[c^s]
            if not s:break
            s=(s-1)&c
    return h,p,q

def construct(B,C):
    A=ordinal(tt(3),B); b=len(A); c=len(C)
    F=ordinal(A,C)
    return tuple(v|((i>=b)<<(b+c)) for i,v in enumerate(F))+((1<<b)-1,)

def read(B,C):
    F=construct(B,C)
    g,_,w,_,e,_,k=profile(F,5)
    return 3*w+2*g,e-k

K=((9,-5,0),(-46,91,-30),(0,-78,78))
def model(B,C):
    a,b=slots(B),slots(C)
    return sum(a[i]*b[j]*K[i][j] for i,j in product(range(3),repeat=2)),a[0]*(8*b[0]-b[1])

cases=list(product((tt(1),tt(2),tt(3)),repeat=2))
cases+=list(product(((2,4,1),(118,88,26,97,104,6,36)),
                  ((22,12,24,17,2),ordinal((2,4,1),tt(1)))))
for B,C in cases:
    observed=read(B,C)
    if observed!=model(B,C):raise ValueError((B,C,observed))
    print(slots(B),slots(C),observed)

from math import isqrt
V=(118,88,26,97,104,6,36)
for C,p in ((tt(3),13936607),((22,12,24,17,2),265187651)):
    if model(V,C)!=(p,0) or any(p%d==0 for d in range(2,isqrt(p)+1)):
        raise ValueError(p)
    print("exact prime pair",p,7*p)
for B,p in (((58,356,81,198,490,332,129,295,77),440866433),
            ((178,428,97,309,390,400,315,76,133),441222113)):
    if slots(B)[:2]!=(2207,35824) or read(B,tt(3))!=(p,0) or model(B,tt(3))!=(p,0):
        raise ValueError(B)
    print("third slot distinguishes",slots(B),p)
