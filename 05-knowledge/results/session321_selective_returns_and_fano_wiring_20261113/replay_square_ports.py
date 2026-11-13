"""Literal and independent cut-interface controls at square-sum orders 15..25."""
from math import isqrt

def ck(test, receipt):
    if not test:
        raise AssertionError(receipt)

def graph(n):
    return {u:{v for v in range(1,n+1) if u!=v and isqrt(u+v)**2==u+v}
            for u in range(1,n+1)}

def paths(A):
    out=[]
    def f(p,R):
        if not R:
            if p[0]<p[-1]:
                out.append(p)
            return
        for z in sorted(R&A[p[-1]]):
            f(p+(z,),R-{z})
    for u in A:
        f((u,),set(A)-{u})
    return sorted(out)

def edges(p):
    return {tuple(sorted(e)) for e in zip(p,p[1:])}

def pair(p):
    return p[0],p[-1]

def exposed_covers(A,a,b):
    """Enumerate two-piece path covers by independent DFS on complementary supports."""
    def left(p,R):
        if R and b in R:
            yield p,R
        for z in sorted((R-{b})&A[p[-1]]):
            yield from left(p+(z,),R-{z})
    def right(p,R):
        if not R:
            yield p
        for z in sorted(R&A[p[-1]]):
            yield from right(p+(z,),R-{z})
    return {(P,Q) for P,R in left((a,),set(A)-{a})
                  for Q in right((b,),R-{b})}

def run():
    P={n:paths(graph(n)) for n in range(15,26)}
    counts=[len(P[n]) for n in P]
    ck(counts==[1,1,1,0,0,0,0,0,3,0,10],counts)
    for n in (15,16):
        ck(edges(P[n][0])<edges(P[n+1][0]),n)
    for n in (15,16,17):
        ck(sorted(len(edges(P[n][0])-edges(q)) for q in P[23])==[5,6,6],n)
    A=graph(23)
    channels={}
    for a,b,ear,expected in ((1,12,(24,),0),(11,1,(25,24),4),(11,12,(25,24),6)):
        c=exposed_covers(A,a,b)
        ck(len(c)==expected,(a,b,len(c)))
        def expand(P,Q):
            p=P[::-1]+ear+Q
            return min(p,p[::-1])
        channels[b if a==11 else 24]={expand(P,Q) for P,Q in c}
    ck(channels[24]==set(P[24]),"Q24 channel")
    ck(not (channels[1]&channels[12]),"disjoint ears")
    ck(channels[1]|channels[12]==set(P[25]),"exhaustive Q25 channels")
    print("sizes 15..25:",counts)
    print("channels 1-12, 11-1, 11-12:",[len(channels[t]) for t in (24,1,12)])
    small=sorted(P[23],key=lambda p:(p[0] if p[0]!=18 else 100))
    big=sorted(P[25],key=lambda p:next(v for v in pair(p) if v!=18))
    M=[[len(edges(p)-edges(q)) for q in big] for p in small]
    for p,row in zip(small,M):
        print("old endpoints",pair(p),"lost edges",row)
        ck(min(row)==2,(p,row))
    costs=[next(len(edges(p)-edges(q)) for q in big if pair(q)==pair(p)) for p in small]
    ck(costs==[2,2,3],costs)
    print("same endpoints:",costs)

if __name__=='__main__':
    run()
