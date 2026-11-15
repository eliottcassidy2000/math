from functools import lru_cache
from itertools import combinations
import random
def allrows(n):
    pairs=list(combinations(range(n),2))
    for m in range(1<<len(pairs)):
        r=[0]*n
        for k,(i,j) in enumerate(pairs):
            if m>>k&1:r[i]|=1<<j
            else:r[j]|=1<<i
        yield tuple(r)
def counts(r):
    @lru_cache(None)
    def start(s,i):
        t=s^(1<<i)
        return sum(start(t,j) for j in range(len(r)) if t&r[i]&(1<<j)) if t else 1
    @lru_cache(None)
    def H(s):
        return sum(start(s,i) for i in range(len(r)) if s>>i&1) if s else 1
    f=(1<<len(r))-1
    return H(f),sum(H(s)*H(f^s) for s in range(f+1)),[(H(f^(1<<q)),sum(start(s,q)*H(f^s) for s in range(f+1) if s>>q&1)) for q in range(len(r))]
def glue(r,q,w,u):
    n,t,m=len(r),len(w),len(u)
    mw=((1<<t)-1)<<n
    mu=((1<<m)-1)<<(n+t)
    return tuple(v|mw|(mu if i!=q else 0) for i,v in enumerate(r))+tuple((v<<n)|mu for v in w)+tuple((v<<(n+t))|(1<<q) for v in u)
def sc(v):
    if not v:return 0
    s=-1 if v<0 else 1
    v=abs(v)
    d=2
    while d*d<=v:
        while v%(d*d)==0:v//=d*d
        d+=1
    return s*v
def main():
    rng=random.Random(391)
    pools={n:list(allrows(n)) for n in range(1,5)}
    trials=0
    for _ in range(600):
        r,w,u=(rng.choice(pools[rng.randrange(1,5)]) for i in range(3))
        q=rng.randrange(len(r))
        A,_,marks=counts(r)
        K,B=marks[q]
        hu,qu,_=counts(u)
        hw,qw,_=counts(w)
        actual=counts(glue(r,q,w,u))[0]
        expected=(A+K)*hu*hw+(qu-hu)*(B*qw-K*hw)
        if actual!=expected:raise AssertionError((r,q,w,u,actual,expected))
        trials+=1
    print('600 literal anchor glue checks:',trials)
    atlas=[
        ((4,1,2),0,3,(3,1,24,145),(51,289,1341)),
        ((8,1,3,6),2,2,(5,3,36,-63),(77,431,2003)),
        ((16,1,3,7,14),1,2,(9,5,40,-1015),(89,479,2219)),
    ]
    def tt(n):
        return tuple((1<<n)-(1<<(i+1)) for i in range(n))
    squareclasses=[]
    for raw,q,k,profile,outputs in atlas:
        n=len(raw)
        r=tuple(a|(((1<<k)-1)<<n) for a in raw)+tuple(a<<n for a in tt(k))
        A,_,marks=counts(r)
        K,B=marks[q]
        D=(B-K)**2-4*B*(A+K)
        reachable=1<<q
        for step in range(len(r)):
            for i,a in enumerate(r):
                if reachable>>i&1: reachable|=a
        got=tuple(counts(glue(r,q,tt(t),tt(t)))[0] for t in (1,2,3))
        if (A,K,B,D)!=profile or got!=outputs or reachable!=(1<<len(r))-1:
            raise AssertionError((r,q,profile,got,reachable))
        if r[q]==reachable^(1<<q):
            raise AssertionError('source gate accidentally collapsed to singleton wrap')
        squareclasses.append(sc(D))
        print('Anchor:',r,q,'profile:',profile,'literal counts:',got)
    a,b,c=squareclasses
    if a==1 or b==1 or c==1 or sc(a*b*c)!=1:raise AssertionError(squareclasses)
    print('Odd square-class relation:',squareclasses)
if __name__=='__main__':main()
