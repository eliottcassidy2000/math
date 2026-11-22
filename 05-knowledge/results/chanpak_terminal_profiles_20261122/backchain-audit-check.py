"""Exact positive-boundary inverse and exponent certificate (stdlib)."""
from fractions import Fraction as Q
from functools import lru_cache
from collections import defaultdict

A=[(1,0),(2,0),(4,-1)]  # F+U at a post-positive state after j zeros
B=[(1,0),(2,-1),(3,0)]  # scalar F there
def make(j):
    while len(A)<=j:
        A.append(tuple(map(sum,zip(*A[-3:]))))
        B.append(tuple(map(sum,zip(*B[-3:]))))
def v2(t):
    assert t>0
    return (t&-t).bit_length()-1

@lru_cache(None)
def scalar(h):
    if h==1:return ()
    if h<3 or not h&1:return None
    for g in range(1,h//3+1,2):
        w=back(h,g)
        if w is not None:return w
    return None

@lru_cache(None)
def back(h,g):
    if h<3 or not h&1 or g<1 or not g&1 or 3*g>h:return None
    d=h-g
    if d%g==0:
        x=d//g
        if x>=2 and x&(x-1)==0:
            pre=scalar(g)
            if pre is not None:return pre+(x.bit_length()-1,)
    j=1
    while True:
        make(j)
        a,b=A[j];c,e=B[j]
        if 3*c+e>g:break  # g=B_j(h0,g0) >= (3c+e)
        det=a*e-b*c
        for m in range(1,v2(d)+1):
            u=d>>m
            if det:
                top=(u*e-b*g,a*g-u*c)
                if any(x%det for x in top):continue
                cand=[tuple(x//det for x in top)]
            elif a*g!=c*u or b*g!=e*u:
                continue
            elif not e:
                if g%c:continue
                cand=((g//c,t) for t in range(1,(g//c)//3+1,2))
            else:
                # Enumerate bounded scalar progression; testing divisibility
                # as written is equivalent to a direct modular progression.
                lo=(g+c-1)//c
                hi=3*g//(3*c+e)
                cand=((p,(g-c*p)//e) for p in range(lo,hi+1)
                      if (g-c*p)%e==0)
            for h0,g0 in cand:
                if h0>=h or (a*h0+b*g0,c*h0+e*g0)!=(u,g):continue
                pre=back(h0,g0)
                if pre is not None:return pre+(0,)*j+(m,)
        j+=1
    return None

def eval_word(word):
    D,U,V=0,0,1
    for m in word:
        D,U,V=(D+2*U,V,D+U) if not m else ((1<<m)*(D+2*U+V),0,D+U+V)
    return D,U,V

def forward(N):
    roots=[]
    m=1
    while (1<<m)+1<=N:
        roots.append((1<<m,0,1));m+=1
    seen=set()
    F=defaultdict(set)
    while roots:
        D,U,V=roots.pop()
        if (D,U,V) in seen:continue
        seen.add((D,U,V))
        h=D+U+V
        if not U:F[h].add(V)
        z=(D+2*U,V,D+U)
        if sum(z)<=N:roots.append(z)
        m=1
        while h+(1<<m)*(h+U)<=N:
            roots.append(((1<<m)*(h+U),0,h));m+=1
    return seen,F

seen,f=forward(1001)
for h in range(3,1002,2):
    found=set()
    for g in range(1,h//3+1,2):
        w=back(h,g)
        if w is not None:
            assert eval_word(w)==(h-g,0,g)
            found.add(g)
    assert found==f[h],(h,found,f[h])
print('All terminal boundary fibers through 1001 independently matched.')

make(13)
assert A[13]==(3136,-504) and B[13]==(2632,-423)
alpha=[Q(a)+Q(b,3) for a,b in A]
beta=[Q(c)+Q(e,3) for c,e in B]
assert all(alpha[j+1]>=Q(7,4)*alpha[j] for j in (1,2,3))
assert alpha[13]==2968
r=Q(13,25)
g=Q(3,5)
assert r**25*2**24>1 and g**25*Q(7,4)**24>1 and alpha[13]<2**25
K=10**6
S=0
for j in range(1,13):
    for m in range(1,13):
        a=2**m*alpha[j]+beta[j]
        lo,hi=0,K
        rhs=K**25*a.denominator**24
        d=a.numerator**24
        while hi-lo>1:
            k=(lo+hi)//2
            if k**25*d>=rhs:hi=k
            else:lo=k
        S+=hi
assert S==990915
tail_m=2*sum(1/alpha[j] for j in range(1,13))*r**13/(1-r)
tail_j=2/alpha[13]/(1-g)*r/(1-r)
assert Q(S,K)+tail_m+tail_j<Q(199,200)
print('Exact 24/25 primitive-branch weight bound below 199/200 certified.')
