from itertools import combinations
from random import Random

def need(ok, *context):
    if not ok:
        raise AssertionError(context)

def table(rows):
    n=len(rows); out=[[0]*n for _ in range(1<<n)]
    for i in range(n): out[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,k in enumerate(out[mask]):
            if k:
                dst=rows[i]&~mask
                while dst:
                    bit=dst&-dst; dst^=bit
                    out[mask|bit][bit.bit_length()-1]+=k
    return out

def prof(rows):
    hp=[1]+[sum(v) for v in table(rows)[1:]]
    M=len(hp)-1
    return hp[M], sum(hp[s]*hp[M^s] for s in range(M+1))

def rows(n, code):
    r=[0]*n
    for b,(i,j) in enumerate(combinations(range(n),2)):
        if code>>b&1: r[i]|=1<<j
        else: r[j]|=1<<i
    return r

def attach(U,R,q):
    # U inherits v's exterior arcs, with v -> q and all others -> v.
    k=len(U); m=len(R); mask=(1<<k)-1
    return [r|1<<(q+k) for r in U]+[(r<<k)|(mask if i!=q else 0) for i,r in enumerate(R)]

def gate(R,q):
    n=len(R); M=(1<<n)-1
    hp=[1]+[sum(v) for v in table(R)[1:]]
    rev=table([M^r^(1<<i) for i,r in enumerate(R)])
    E=O=0
    for s in range(1,M+1):
        if not s&(1<<q) or s.bit_count()==1: continue
        c=rev[s][q]*hp[M^s] # outside q-start nonempty cycle run
        if s.bit_count()%2: E+=c
        else: O+=c
    return hp[M],E,O

def cyclic(U,V,W):
    a,b,c=len(U),len(V),len(W)
    return (
        [r|(((1<<b)-1)<<a) for r in U]
        +[(r<<a)|(((1<<c)-1)<<(a+b)) for r in V]
        +[(r<<(a+b))|((1<<a)-1) for r in W]
    )

def ordered_slots(rows):
    counts=[1]+[sum(v) for v in table(rows)[1:]]
    M=len(counts)-1; n=len(rows)
    out=[[0]*(M+1) for _ in range(n+2)]; out[0][0]=1
    for j in range(1,n+1):
        for s in range(1,M+1):
            a=s
            while a:
                out[j][s]+=counts[a]*out[j-1][s^a]
                a=(a-1)&s
    return [tab[M] for tab in out]

def cycle_transfers(U,V,W):
    data=[ordered_slots(A) for A in (U,V,W)]
    m=min(map(len,data))-2; tally=0
    for r in range(1,m+1):
        A,B,C=(d[r] for d in data)
        X,Y,Z=(d[r+1] for d in data)
        tally+=3*A*B*C+X*B*C+A*Y*C+A*B*Z+A*Y*Z+X*B*Z+X*Y*C
    return tally

def transitive(m): return [(1<<i)-1 for i in range(m)]

def main():
    profiles=set()
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            h,Q=prof(rows(n,code))
            P=Q-h
            need(P>=n*h and P>=2**n-1,n,h,Q)
            if h>1: need(P>=3*2**(n-1)-3,n,h,Q)
            profiles.add((h,Q))
    profiles.add((1,64))
    vals={(Q-h)*(R-g)+2*h*g for h,Q in profiles for g,R in profiles}
    expect=[3,5,9,11,15,17,23,27,33,37,47,51,65,69,91,93,95]
    need(sorted(v for v in vals if v<=98)==expect)
    print('exact cyclic two-core singleton-wrap spectrum through 98:',expect)
    rng=Random(594763)
    for _ in range(300):
        n=rng.randrange(1,6); k=rng.randrange(1,6)
        U=rows(n,rng.randrange(1<<(n*(n-1)//2)))
        R=rows(k,rng.randrange(1<<(k*(k-1)//2)))
        q=rng.randrange(k)
        h,Q=prof(U); A,E,O=gate(R,q)
        need(sum(table(attach(U,R,q))[-1])==(A-2*E)*h+(E+O)*Q,U,R,q)
        qsource=[((1<<k)-1)<<1]+[r<<1 for r in R]
        g,S=prof(R)
        need(sum(table(attach(U,qsource,0))[-1])==(Q-h)*(S-g)+2*h*g,U,R)
    print('300 literal one-outlet and two-core cyclic gluing controls: PASS')
    for sizes,value in [((2,2,2),45),((3,2,2),137)]:
        U,V,W=[transitive(m) for m in sizes]
        need(cycle_transfers(U,V,W)==sum(table(cyclic(U,V,W))[-1])==value)
    for _ in range(100):
        u,v,w=[rng.randrange(1,4) for _ in range(3)]
        U,V,W=[rows(m,rng.randrange(1<<(m*(m-1)//2))) for m in (u,v,w)]
        need(cycle_transfers(U,V,W)==sum(table(cyclic(U,V,W))[-1]))
    print('all three sides finite spectrum through 98:',sorted(expect+[45]))
    cases=[
        ((38,0,2,7,15,30),0,(5,9,10),2,63,(4,5,152,3,11,31,63,123)),
        ((20,1,2,7,14),3,(9,7,11),2,67,(32,33,83,7,11,28,59)),
        ((4,1,2,7),3,(3,3,6),3,69,None),
        ((18,0,3,7,14),2,(5,3,6),3,71,(32,33,35,151,7,24,63,119)),
    ]
    for R,q,aeo,k,N,expected_rows in cases:
        A,E,O=gate(R,q)
        need((A,E,O)==aeo)
        for m in range(1,6):
            outright=attach(transitive(m),R,q)
            direct=sum(table(outright)[-1])
            need(direct==A-2*E+(E+O)*2**m)
            if m==k:
                need(direct==N,N)
                if expected_rows is not None: need(tuple(outright)==expected_rows,N)
    for m in range(2,6):
        triple=[sum(table(attach(transitive(m-(N==67)),R,q))[-1])
                for R,q,_,_,N,_ in cases[1:]]
        need(triple==[9*2**m+i for i in (-5,-3,-1)],m,triple)
    print('four literal outlet arms and shared consecutive-odd triple: PASS')

if __name__=='__main__':
    main()
