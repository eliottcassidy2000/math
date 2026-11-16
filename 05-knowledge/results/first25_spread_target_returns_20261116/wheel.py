from math import isqrt,gcd,lcm
from collections import defaultdict

def primes(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,isqrt(p)+1))]
P=primes(100)
O={2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}
B={2,3,5,7,11,13,17,19,23,31,47}
C={2,3,5,7,11,13}
def wheel(p):
    if p==2: return None
    j,k,h=1,0,None
    while True:
        if (99*j+23)%p==0:h=k
        j=2*j%p;k+=1
        if j==1:return None if h is None else (h,k)

def merge(x,y):
    a,m=x;b,n=y
    d=gcd(m,n)
    if (b-a)%d:return None
    l=m*n//d
    j=(b-a)//d*pow(m//d,-1,n//d)%(n//d)
    return (a+m*j)%l,l

def union(clauses):
    c={}
    for y in clauses:
        out=defaultdict(int,c)
        out[y]+=1
        for x,k in c.items():
            m=merge(x,y)
            if m is not None:out[m]-=k
        c={x:k for x,k in out.items() if k}
    return c

if __name__=='__main__':
    W={p:w for p in P if (w:=wheel(p)) is not None}
    print('Hits',W)
    reduced={p:(a,m) for p,(a,m) in W.items() if not any(q!=p and m%n==0 and (a-b)%n==0 for q,(b,n) in W.items())}
    print('Redundant',sorted(W.keys()-reduced.keys()))
    l=lcm(*(m for _,m in reduced.values()))
    C=union(reduced.values())
    def size(c):return sum(k*l//m for (_,m),k in c.items())
    print('period',l,'accepted',l-size(C),'length of exact IE',len(C))
    for p in P:
        if l%p==0:
            L=l//p
            D=union(list(reduced.values())+[((a+L)%m,m) for (a,m) in reduced.values()])
            print('essential',p,size(D)-size(C))
            assert size(D)>size(C)
    l0=l
    CO=union(w for p,w in W.items() if p in O)
    # Allow same common period for all three wheels
    for name,part in [('Fi22', {2,3,5,7,11,13}),('Baby',B),('Monster',O)]:
        print(name,'accepted',l-size(union(w for p,w in W.items() if p in part)))
    for k in range(1,17):
        n=99*(1<<k)+23
        if any(n%p==0 for p in P):continue
        d=next((d for d in range(2,isqrt(n)+1) if n%d==0),None)
        assert d is None or (k,d)==(16,1301)
        if d is not None:print('first survivor composite',k,n,d,n//d)
    ord1301=next(i for i in range(1,1301) if pow(2,i,1301)==1)
    assert ord1301==1300 and lcm(l,ord1301)==306048600
