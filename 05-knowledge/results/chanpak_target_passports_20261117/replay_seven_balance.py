from math import gcd,isqrt
from functools import reduce
from operator import mul
from replay_ordinal_digits import end_profiles,ordinal,strong

L9=(130,28,281,97,361,71,391,318,43)
U7=(94,60,120,48,96,65,10)
U8=(2,8,19,53,195,87,143,47)
def need(x):
    if not x:raise AssertionError("sevenfold balance replay")
def prime(p):
    return p>1 and all(p%d for d in range(2,isqrt(p)+1))
def tt(n):
    return tuple((1<<n)-(1<<(i+1)) for i in range(n))

def profile(a):
    h,_=end_profiles(a)
    F=len(h)-1
    P=[0]*(F+1)
    for S in range(F+1):
        U=S
        while True:
            P[S]+=h[U]*h[S^U]
            if not U:break
            U=(U-1)&S
    return h[F],P[F],sum(h[S]*P[F^S] for S in range(F+1))

def cycle(A,B,C):
    x,y,z=len(A),len(B),len(C)
    return tuple(a|(((1<<y)-1)<<x) for a in A) + tuple((b<<x)|(((1<<z)-1)<<(x+y)) for b in B) + tuple((c<<(x+y))|((1<<x)-1) for c in C)

tiles=[profile(x) for x in (L9,U7,U8,tt(4))]
need(tiles==[(2197,34928,291327),(91,1184,8121),(169,3304,30249),(1,16,81)])
h,P,Q=(a*b for a,b in zip(tiles[0],tiles[3]))
j,V,Z=(a*b for a,b in zip(tiles[1],tiles[2]))
need((h,P,Q)==(2197,558848,23597487) and (j,V,Z)==(15379,3911936,245652129))
need(j==7*h and V==7*P and Z!=7*Q and gcd(h,P)==1)
need(h==13**3 and P==2**8*37*59)
B=ordinal(L9,tt(4));C=ordinal(U7,U8)
for source,coef in [(B,1),(C,7)]:
    witness=cycle(source,(0,),(0,))
    need(strong(witness) and end_profiles(witness)[0][-1]==coef*(P+h))
    witness=cycle(tt(2),source,(0,))
    need(strong(witness) and end_profiles(witness)[0][-1]==coef*((P-h)*4-P+3*h))
for d,p in [(7,71534741),(13,4578085013),(15,18312333461)]:
    need(p==P*(1<<d)+h and prime(p))
    print('padded-module',d,p,7*p)
for s,p in [(8,141950399),(16,36480127679)]:
    need(p==(P-h)*(1<<s)-P+3*h and prime(p))
    need(strong(cycle(tt(s),B,(0,))) and strong(cycle(tt(s),C,(0,))))
    print('cyclic-reader',s,p,7*p)
print('Exact primitive H/P₂ balance and unequal P₃: PASS')
