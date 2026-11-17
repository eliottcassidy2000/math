"""Standalone verification of signed ordinal tiles and an exact target."""
from itertools import permutations, product
from math import prod

def profile(rows, x, a):
    n=len(rows); S=[0]*n; E=[0]*n; M=[[0]*n for _ in rows]
    for p in permutations(range(n)):
        breaks=[k for k in range(1,n) if not rows[p[k-1]]>>p[k]&1]
        if not breaks: S[p[0]]+=1; E[p[-1]]+=1
        for k in (range(1,n) if not breaks else breaks if len(breaks)==1 else ()):
            M[p[k-1]][p[k]]+=1
    def sig(mask):
        vv=(S,E,[sum(r) for r in M],[sum(r) for r in zip(*M)])
        return tuple(sum(v[i] for i in range(n) if mask>>i&1) for v in vv)
    if sig(x)!=sig(a): raise AssertionError((sig(x),sig(a)))
    delta=sum(M[i][j]-M[j][i] for i in range(n) for j in range(n) if x>>i&1 and a>>j&1)
    return sum(S),delta//2,sig(x)

def graph(tiles, orientation):
    N=2+sum(len(t) for t,x,a in tiles)
    row=[0]*N
    def arc(u,v):row[u]|=1<<v
    arc(1,0) if orientation else arc(0,1)
    i=2
    for t,x,a in tiles:
        for r,e in enumerate(t,i):
            row[r]|=e<<i
            for w in range(2,i): arc(w,r)
            for p,m in enumerate((x,a)):
                arc(r,p) if m>>(r-i)&1 else arc(p,r)
        i+=len(t)
    return row

def count(T):
    n=len(T); N=1<<n
    assert all(((T[u]>>v)&1)+((T[v]>>u)&1)==1 for u in range(n) for v in range(u))
    dp=[[0]*n for _ in range(N)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,N):
        for u,h in enumerate(dp[s]):
            if not h:continue
            r=T[u]&~s
            while r:
                b=r&-r;r-=b;dp[s|b][b.bit_length()-1]+=h
    return sum(dp[-1])

tri=((2,4,1),4,2)
t161=((2,12,17,37,75,87,143,63),37,67)
t155=((14,28,88,96,73,23,35),5,18)
t89=((30,60,120,48,64,81,11),68,33)
hit=((126,76,88,112,98,70,0),19,7)
gate=((0,),1,1)

def test(tiles, expected=None):
    P=[profile(*t) for t in tiles]
    expected_gap=2*sum(s*prod(P[j][0] for j in range(len(P)) if j!=i)
                       for i,(h,s,_) in enumerate(P))
    h,j=(count(graph(tiles,b)) for b in (0,1))
    assert j-h==expected_gap,(P,h,j)
    if expected is not None:assert (h,j)==expected,(h,j)
    return h,j

if __name__=='__main__':
    for t in (tri,t161,t155,t89,hit):print('profile', profile(*t))
    print('signed sum',test([tri,t161],(224143,224145)))
    print('reordered',test([t161,tri],(178695,178697)))
    print('bare 155',test([t155],(2071,2073)))
    print('prescribed Q5',test([hit,gate],(1089,1091)))
    for k in range(4):
        for j in range(k+1):
            for w in product(range(2),repeat=k):
                a=[((0,),x,x) for x in w]
                assert test(a[:j]+[t155]+a[j:])[1]-test(a[:j]+[t155]+a[j:])[0]==2
