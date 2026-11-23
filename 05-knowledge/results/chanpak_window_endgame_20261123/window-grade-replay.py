from collections import Counter
from itertools import combinations,product
from math import comb

def height(R):
    n=len(R); M=1<<n
    E=[[0]*n for _ in range(M)]
    for i in range(n): E[1<<i][i]=1
    for s in range(1,M):
        for i,c in enumerate(E[s]):
            if not c:continue
            v=R[i]&(M-1-s)
            while v:
                b=v&-v;v-=b
                E[s+b][b.bit_length()-1]+=c
    return sum(E[-1])

BANK={3:(4,1,2),5:(8,1,3,6),35:(4,1,10,99,15,23,55),
      49:(4,1,74,35,15,23,59),63:(8,9,3,148,7,31,63,119),
      343:(130,20,193,7,13,31,59,122)}
def words(N):
    if N==1:
        yield ()
    else:
        for d in BANK:
            if N%d==0:
                for r in words(N//d):yield (d,)+r
def a(N,n):
    total=0
    for w in words(N):
        s=sum(len(BANK[d]) for d in w)
        if n>=s:total+=comb(n-s+len(w),len(w))
    return total
def diff(src,r,n):
    return sum((-1)**(r-j)*comb(r,j)*src(n+j) for j in range(r+1))

def apex(R,S):
    n=len(R)
    return tuple(r|(0 if S>>i&1 else 1<<n) for i,r in enumerate(R))+(S,)
def polycoll(weights,a,M=0):
    n=len(a);total=0
    for d in product((-1,0,1),repeat=n):
        I=[i for i in range(n) if d[i]]
        J=[j for j in range(n) if not d[j]]
        L=sum(a[i]*d[i] for i in I)-sum(weights[i,j]*(d[i]+d[j])//2
                                       for i,j in combinations(I,2))
        counts={0:1}
        for j in J:
            w=sum(weights[min(i,j),max(i,j)]*d[i] for i in I)
            acc=Counter(counts)
            for k,p in counts.items():acc[(k+w)%M if M else k+w]+=p
            counts=acc
        total+=counts.get(L%M if M else L,0)
    return total

if __name__=="__main__":
    assert all(height(t)==h for h,t in BANK.items())
    for r in range(1,8):
        for p,want in ((3,0),(5,r)):
            n=100
            assert diff(lambda i:a(p**r,i),r,n)==1
            assert diff(lambda i:a(7*p**r,i),r,n)==want
        assert diff(lambda i:a(7**r,i),r,100)==0
    print("six inherited literal seeds and restricted-bank finite differences PASS")
    for R in [(4,1,2),(8,1,3,6),(16,1,3,7,14)]:
        n=len(R);H=height(R)
        hs=[(height(apex(R,S))-H)//2 for S in range(1<<n)]
        a=[hs[1<<i] for i in range(n)]
        w={(i,j):a[i]+a[j]-hs[(1<<i)|(1<<j)] for i,j in combinations(range(n),2)}
        for m in (0,2,3,5,11):
            b=Counter(h%m if m else h for h in hs)
            assert sum(v*v for v in b.values())==polycoll(w,a,m)
        actual=Counter(2*v+H for v in hs)
        print("core",H,"height counts",dict(sorted(actual.items())),
              "collision energy",sum(v*v for v in actual.values()))
