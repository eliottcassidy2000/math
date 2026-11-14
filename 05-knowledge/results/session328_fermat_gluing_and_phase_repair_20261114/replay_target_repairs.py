from profile_tools import canon,extend,profile,strong,hp_all
from math import gcd,isqrt
from collections import defaultdict
def need(b):
    if not b:raise RuntimeError('target-specific repair replay failed')
def tt(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def join(A,B):return tuple(r|(((1<<len(B))-1)<<len(A)) for r in A)+tuple(r<<len(A) for r in B)
def reachable(T,u):
    x=1<<u
    while True:
        y=x
        for i,r in enumerate(T):
            if (x>>i)&1:y|=r
        if x==y:return x
        x=y
def core(T):
    rr=[reachable(T,u) for u in range(len(T))];z=len(T)-1
    S=[i for i in range(len(T)) if (rr[i]>>z)&1 and (rr[z]>>i)&1]
    for i in set(range(len(T)))-set(S):
        need(all(i==j or not ((rr[i]>>j)&(rr[j]>>i)&1) for j in range(len(T))))
    return tuple(sum(((T[i]>>j)&1)<<a for a,j in enumerate(S)) for i in S)

letters={(1,2)}
A={(0,)}
for n in range(2,7):
    A={canon(extend(T,c)) for T in A for c in range(1<<(n-1))}
    letters.update(profile(T) for T in A if strong(T) and profile(T)[1]<=256)
cl={(1,1)}
for z in range(1,257):
    cl.update((h*a,q*b) for h,q in tuple(cl) if q==z for a,b in letters if q*b<=256)
D=defaultdict(set)
for h,q in cl:
    if q>1 and gcd(h,99)==1:D[q].add(h)
need({q:h for q,h in D.items() if len(h)>1}=={32:{1,5},64:{1,5},128:{1,5},256:{1,5,19}})
R=(2,4,9,3); B=(2,20,9,3,45,15)
need(profile(R)==(5,32) and profile(join(R,tt(3)))==(5,256))
for n,q,old,c,U,I,S in (
 (3169,32,5,1,tt(5),(2,48,3,7,13,29,63),33),
 (25667,256,5,17,B,(8,9,3,100,15,23,55),9),
):
 h,b=profile(U);need(q==b and n==99*q+c*h)
 d=n*pow(old,-1,198)%198;w=(n-d*old)//99
 need((d,w)==((119,26) if n==3169 else (25,258)))
 T=core(extend(join(I,U),S))
 need(strong(T) and hp_all(T)[-1]==n and all(n%d for d in range(2,isqrt(n)+1)))
print('first exact budget-free phase repairs: Q32 at prime3169; nontrivial Q256 at prime25667')
