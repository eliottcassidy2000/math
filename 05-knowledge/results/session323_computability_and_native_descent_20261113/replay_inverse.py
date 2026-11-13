from math import isqrt
from functools import cache

def require(ok):
    if not ok:raise RuntimeError('replay mismatch')

def stats(w):
    n=len(w); h=1;q=1<<(n+1)
    for i in range(n):
        if not w[i]:continue
        for j in range(i+1,n):
            if w[j]:continue
            d=j-i
            h+=(1<<(d-1))+(d==1)
            q+=(1<<(n-d))*(3**(d-1)+1)
    return h,q

def candidates(n):
    # Complete for strong one-apex dyadic words. The last full check is exact.
    if n==3:
        yield (1,0,0);yield(1,1,0);return
    if n<5 or n%4!=1:return
    J=(n-1)//4
    def extend(j,a,b,carry):
        if j==J:
            for v in (0,1):
                yield tuple(a+[v]+[1-i for i in b[::-1]])
            return
        alpha=(3**(n-2-2*j)+1)//4
        beta=(3**(n-3-2*j)+1)//2
        good=0
        for e in (0,1):
            for f in (0,1):
                for g in (0,1):
                    for h in (0,1):
                        aa=a+[e,g];bb=b+[f,h]
                        C=sum(aa[i]*bb[2*j-i] for i in range(2*j+1))
                        D=sum(aa[i]*bb[2*j+1-i] for i in range(2*j+2))
                        v=carry+alpha*C+beta*D
                        if v%4==0:
                            good+=1
                            yield from extend(j+1,aa,bb,v//4)
        require(good <= (6 if j%2==0 else 4))
    alpha=(3**(n-2)+1)//4;beta=(3**(n-3)+1)//2
    require((alpha+beta)%4==0)
    for v in (0,1):
        yield from extend(1,[1,v],[1,1-v],(alpha+beta)//4)

def dyadic(n):
    for w in candidates(n):
        h,q=stats(w)
        if not q&(q-1):yield (w,h,q.bit_length()-1)

def divisors(n):
    ds=[]
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            ds.append(d)
            if d*d!=n:ds.append(n//d)
    return sorted(ds)

def decode(N,A,h,extra=None):
    # Q_V(2)=2^m is part of the grammar, not inferred for arbitrary tournaments.
    if min(N,A,h)<=0 or A<=h:return None
    m=1
    while A<<m<=N:m+=1
    t=(A<<m)-N
    if t<=0 or t%h or not ((t//h)&1):return None
    D=t//h
    available={} if extra is None else dict(extra)
    all_d=divisors(D)
    lengths={n for d in all_d if d>1 for n in (d.bit_length(),d.bit_length()+1) if n+2<=m}
    need=set(all_d)
    for n in sorted(lengths):
        # For each n, generate once, not for every possible H divisor.
        if 2*(3**(n-2)+1)>1<<m:continue
        for w,d,cost in dyadic(n):
            if d in need and cost<=m and (d not in available or cost<available[d][0]):
                available[d]=(cost,w)
    @cache
    def best(t):
        if t==1:return (0,())
        got=(m+1,())
        for d,(cost,w) in available.items():
            if d<=1 or t%d:continue
            u,parts=best(t//d)
            if cost+u<got[0]:got=(cost+u,((d,w),)+parts)
        return got
    cost,parts=best(D)
    if cost>m:return None
    return (m-cost,parts)

# Exact exhaustive check of pruning on all small transitive chain cuts.
for n in range(2,15):
    brute=set()
    for x in range(1<<(n-2)):
        w=(1,)+tuple(x>>i&1 for i in range(n-2))+(0,)
        h,q=stats(w)
        if not q&(q-1):brute.add((w,h,q.bit_length()-1))
    require(brute==set(dyadic(n)))

# Nontrivial absolute examples, forced via descending factors 5 and 19.
for e,f,p in [(1,0,3),(0,1,0),(2,1,4),(3,2,0)]:
    D=5**e*19**f;m=5*e+8*f+p; N=15*(1<<m)-D
    ans=decode(N,15,1);require(ans is not None)
    pad,parts=ans
    hs=1;Q=1<<pad
    for d,w in parts:
        h,q=stats(w);require(h==d)
        hs*=h;Q*=q
    require(15*Q-hs==N)
print('all strong one-apex cuts through 14; physical decreasing compiler controls pass')
