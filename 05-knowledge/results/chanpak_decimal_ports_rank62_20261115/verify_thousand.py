from itertools import combinations
from math import gcd,isqrt

def factors(n):
    out={}
    d=2
    while d*d<=n:
        while n%d==0:
            out[d]=1+out.get(d,0)
            n//=d
        d+=1
    if n>1: out[n]=1
    return out

def sig(n):
    prod=1
    for p,k in factors(n).items(): prod*=sum(p**j for j in range(k+1))
    return prod

def s(n): return sig(n)-n
def sq(n): return isqrt(n)**2==n

def orbit(n):
    row=[n]
    while n:
        n=s(n)
        assert n not in row
        row.append(n)
    return row

def dp(adj):
    n=len(adj)
    h=[[0]*n for _ in range(1<<n)]
    for u in range(n): h[1<<u][u]=1
    for S in range(1,1<<n):
        for u,count in enumerate(h[S]):
            if not count: continue
            choices=adj[u]&~S
            while choices:
                bit=choices&-choices
                choices-=bit
                h[S|bit][bit.bit_length()-1]+=count
    return h

def apex(rows,mask):
    n=len(rows)
    return [row|((1<<n) if not(mask>>i&1) else 0) for i,row in enumerate(rows)]+[mask]

def profile(rows,mask):
    n=len(rows)
    full=(1<<n)-1
    ends=dp(rows)
    starts=dp([sum(1<<j for j,r in enumerate(rows) if r>>i&1) for i in range(n)])
    hs=[1]+[sum(r) for r in ends[1:]]
    def weights(A,accept):
        p=d=0
        for S in range(1,full+1):
            q=sum(A[S][v] for v in range(n) if accept>>v&1)*hs[full^S]
            p+=q
            d+=q if S.bit_count()%2 else -q
        return p,d
    return ((hs[-1],sum(hs[S]*hs[full^S] for S in range(full+1)))
            +weights(starts,mask)+weights(ends,full^mask)+(sum(dp(apex(rows,mask))[-1]),))

def count_transfer(blocks):
    h,N,S,d=1,1,0,0
    for a,q,s,e,t,f,n in blocks:
        h,N,S,d=a*h,a*N+(n-a)*h+t*S+f*d,q*S+h*s,h*e
    return N

def compile_host(x,y,rows,mask):
    first=[1]+list(map(int,bin(x)[3:]))
    flags=first+[mask>>i&1 for i in range(len(rows))]+[1-int(b) for b in reversed(bin(y)[3:])]+[0]
    pos=len(first)
    n=len(flags)
    graph=[0]*n
    for i,j in combinations(range(n),2):
        if not(pos<=i<j<pos+len(rows)) or rows[i-pos]>>(j-pos)&1: graph[i]|=1<<j
        else: graph[j]|=1<<i
    return apex(graph,sum(v<<j for j,v in enumerate(flags)))

def check_host(goal,x,y,rows,mask):
    prof=profile(rows,mask)
    first=[1]+list(map(int,bin(x)[3:]))
    last=[1-int(b) for b in reversed(bin(y)[3:])]+[0]
    p0,p1=profile([0],0),profile([0],1)
    P=(p0,p1)
    count=count_transfer([P[i] for i in first]+[prof]+[P[i] for i in last])
    assert count==goal,(goal,count)
    g=compile_host(x,y,rows,mask)
    if len(g)<=14: assert sum(dp(g)[-1])==goal
    # Every earlier spine vertex reaches final guard directly except internal core arcs,
    # which cannot touch either guard. Initial guard reaches every later spine vertex.
    assert g[-1]&1 and (g[-2]>>(len(g)-1))&1
    assert all(g[0]>>j&1 for j in range(1,len(g)-1))
    assert all(g[j]>>(len(g)-2)&1 for j in range(len(g)-2))
    return len(g)

def main():
    t=1000
    vals=[1]+[t**k+1 for k in range(1,4)]
    assert vals==[1,1001,1000001,1000000001]
    assert vals[3]==vals[1]*(t*t-t+1)
    for k,a in enumerate(vals):
        assert a%4==1
        for b in vals[:k]: assert not sq(a+b)
    for d,ps in [(1,(7,11,13)),(2,(101,9901)),(3,(19,52579))]:
        for p in ps:
            assert factors(p)=={p:1} and factors(t**d+1)[p]==1
            assert next(j for j in range(1,2*d+1) if pow(t,j,p)==1)==2*d
            for n in range(1,31):
                assert (pow(t,n,p)==p-1)==(n%d==0 and n//d%2==1)
    fourth=[1000000001,413350399,107622305,39974047,929673,439227,212853,70955,18037,1079,97,1,0]
    first,second,third=orbit(1001),orbit(1000001),orbit(1000000001)
    assert third==fourth
    assert all(set(a)&set(b)=={0,1} for a,b in combinations((first,second,third),2))
    assert [(a,b) for a,b in zip(third,third[1:]) if b and sq(a+b)]==[(929673,439227)]
    assert 929673+439227==1170**2 and gcd(929673,439227)==9
    tri=[4,1,2]
    assert [check_host(N,x,y,tri,1) for N,x,y in [(1001,4,18),(1000001,1,61481),(1000000001,804,103594)]]==[12,21,31]
    assert check_host(999001,1500,20,[8,1,3,6],8)==21
    assert sum(dp([98,48,107,67,13,88,18])[-1])==101
    assert check_host(9901,178,1,[8,1,3,6],1)==14
    print("Exact prime clocks, fourth aliquot orbit, and literal carrier/ordinal seats passed")

if __name__=='__main__': main()
