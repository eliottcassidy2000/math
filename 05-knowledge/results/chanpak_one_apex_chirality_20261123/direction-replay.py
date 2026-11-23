from itertools import combinations,permutations

def product(x,y):
    a,b=x%7,x//7
    c,d=y%7,y//7
    return (a+(2**b)*c)%7+7*((b+d)%3)

D={1,2,3,14,8,9,10,11,12,13}
def core():
    return tuple(sum(1<<product(x,y) for y in D) for x in range(21))

def contrast(R):
    n=len(R)
    return tuple((1<<n)-1-(1<<i)-r for i,r in enumerate(R))

def apex(R,S):
    n=len(R)
    return tuple(r|((1<<n) if not S>>i&1 else 0) for i,r in enumerate(R))+(S,)

def cycles(R,n,need=()):
    lo=len(R); U=(1<<lo)-1
    ans=0
    if n==3:
        mask=0 if not need else sum(1<<v for v in need)
        for S in combinations(range(lo),3):
            if any(v not in S for v in need):continue
            i,j,k=S
            ans+= ((R[i]>>j&1) and (R[j]>>k&1) and (R[k]>>i&1)) or (
                   (R[i]>>k&1) and (R[k]>>j&1) and (R[j]>>i&1))
        return ans
    if need:
        i=need[0]
        for P in permutations([j for j in range(lo) if j!=i],n-1):
            if any(v not in P for v in need[1:]):continue
            if all(R[a]>>b&1 for a,b in zip((i,)+P,P+(i,))):ans+=1
        return ans
    for S in combinations(range(lo),n):
        for P in permutations(S[1:]):
            if all(R[a]>>b&1 for a,b in zip(S[:1]+P,P+S[:1])):ans+=1
    return ans

def h_ends(R):
    n=len(R);M=(1<<n)-1
    I=[sum(1<<i for i,r in enumerate(R) if r>>j&1) for j in range(n)]
    dp=[[0]*n for _ in range(M+1)]
    for u in range(n):dp[1<<u][u]=1
    for s in range(1,M+1):
        if not s&(s-1):continue
        for v in range(n):
            if s>>v&1:
                t=s^(1<<v)
                dp[s][v]=sum(dp[t][u] for u in range(n) if (t&I[v]) >> u & 1)
    return [1]+[sum(x) for x in dp[1:]],dp

def submasks(s):
    x=s
    while x:
        yield x
        x=(x-1)&s
    yield 0

def q_by_subset(h,m):
    M=len(h)-1
    q=[[1]+[0]*M]
    for _ in range(m):
        prev=q[-1]
        q.append([sum(h[t]*prev[s^t] for t in submasks(s)) for s in range(M+1)])
    return q

def ctest():
    R=core()
    assert all(r.bit_count()==10 and not r>>i&1 for i,r in enumerate(R))
    assert all((R[i]>>j&1)+(R[j]>>i&1)==1 for i in range(21) for j in range(i))
    C=contrast(R)
    def cert(T):
        N=[j for j in range(21) if T[0]>>j&1]
        d={x:sum(T[x]>>j&1 for j in N) for x in N}
        a,b=sorted(N,key=d.get)[-2:]
        assert d[b]==7>d[a]
        return b,sorted(d.values()),sorted(d[j] for j in N if T[b]>>j&1)
    left,right=cert(R),cert(C)
    assert left==(1,[3,3,4,4,4,4,5,5,6,7],[3,4,4,4,5,5,6])
    assert right==(6,[3,3,4,4,4,4,5,5,6,7],[3,4,4,4,4,5,6])
    for T,a,b in ((R,116,117),(C,117,116)):
        Q=apex(apex(T,1),2)
        P=apex(apex(T,2),1)
        assert cycles(Q,5,(21,22))==a and cycles(P,5,(21,22))==b
    print("Cayley chirality, two-apex cycle contrasts PASS")
    from random import Random
    rr=Random(210323);tested=0
    for n in range(2,6):
        for _ in range(12):
            R=[0]*n;M=(1<<n)-1
            for i in range(n):
                for j in range(i):
                    if rr.randrange(2):R[i]|=1<<j
                    else:R[j]|=1<<i
            C=contrast(R);h,end=h_ends(R);hc,sta=h_ends(C)
            assert h==hc
            q=q_by_subset(h,n)
            D={m:[sum((sta[t][u]-end[t][u])*q[m-1][M^t]
                 for t in range(1,M+1)) for u in range(n)] for m in range(1,n+2)}
            for s in range(1<<n):
                ha,_=h_ends(apex(R,s)); hb,_=h_ends(apex(C,s))
                aa,bb=q_by_subset(ha,n),q_by_subset(hb,n)
                for m in range(1,n+1):
                    goal=m*sum(D[m][u]+D[m+1][u] for u in range(n) if s>>u&1)
                    assert aa[m][-1]-bb[m][-1]==goal,(n,s,m,R,D,goal)
                    tested+=1
    print(tested,"literal cover-polynomial divergence identities PASS")

if __name__=="__main__":ctest()
