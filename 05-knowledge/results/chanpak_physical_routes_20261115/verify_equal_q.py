#!/usr/bin/env python3
"""Independent exact graph checks for the fixed-Q and cyclic routes."""
codes=[1623487625516629146,1695545228144487450,1550304140166016026]
def need(c):
    if not c:raise RuntimeError('fixed-Q route check failed')
def dec(c,n=8):
    return [c>>(n*i)&((1<<n)-1) for i in range(n)]
def trans(n):
    return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def ordinal(*T):
    out=[]; offset=0; n=sum(map(len,T))
    for A in T:
        suffix=((1<<(n-offset-len(A)))-1)<<(offset+len(A))
        out.extend((a<<offset)|suffix for a in A)
        offset+=len(A)
    return out
def cyc(A,B):
    x,y=len(A),len(B)
    return [a|(((1<<y)-1)<<x) for a in A]+[(b<<x)|(1<<(x+y)) for b in B]+[(1<<x)-1]
def counts(A):
    n=len(A); N=1<<n
    d=[[0]*n for _ in range(N)]
    for i in range(n):d[1<<i][i]=1
    h=[1]+[0]*(N-1)
    for m in range(1,N):
        for i,v in enumerate(d[m]):
            if not v:continue
            h[m]+=v
            r=A[i]&~m
            while r:
                j=r&-r;r-=j
                d[m|j][j.bit_length()-1]+=v
    return h[-1],sum(a*b for a,b in zip(h,reversed(h)))

def bouquet(D):
    # T+ on 0,1,2,3,4; replace vertex 2 by D.
    rows=[0]*5
    for i,j in ((0,1),(1,2),(2,0),(0,3),(3,4),(4,0),(1,3),(4,1),(3,2),(4,2)):rows[i]|=1<<j
    B=[(1 if i<2 else len(D) if i==2 else 1) for i in range(5)]
    offsets=[0]
    for s in B:offsets.append(offsets[-1]+s)
    out=[]
    for i,s in enumerate(B):
        direct=sum(((1<<B[j])-1)<<offsets[j] for j in range(5) if rows[i]>>j&1)
        if i==2:out.extend((a<<offsets[i])|direct for a in D)
        else:out.append(direct)
    return out

cores=[dec(c) for c in codes]
need([counts(D) for D in cores]==[(645,8784),(651,8784),(657,8784)])
for A,B,want in zip(codes,codes[1:],[[(0,7),(1,4)],[(1,7),(2,6)]]):
    need([(i,j) for i in range(8) for j in range(i+1,8) if ((A^B)>>(8*i+j))&1]==want)
for t in range(3):
    D=[ordinal(A,trans(t)) for A in cores]
    need([counts(bouquet(A))[0] for A in D]==[52704*(1<<t)+h for h in (645,651,657)])
    for s in range(2,5):
        got=[counts(cyc(A,trans(s)))[0] for A in reversed(D)]
        base=((1<<s)-1)*8784*(1<<t)-((1<<s)-3)*651
        M=6*((1<<s)-3)
        need(got==[base-M,base,base+M])
dual_codes=[1119081703992204326,1119081703992269858,974966515916576802]
dual=[dec(c) for c in dual_codes]
need([counts(D) for D in dual]==[(589,8080),(589,8104),(589,8128)])
for a,b,w in zip(dual_codes,dual_codes[1:],[[(0,2)],[(1,2),(1,7)]]):
    need([(i,j) for i in range(8) for j in range(i+1,8) if ((a^b)>>(8*i+j))&1]==w)
for t in range(3):
    D=[ordinal(A,trans(t)) for A in dual]
    got=[counts(cyc(A,trans(1)))[0] for A in D]
    need(got==[(1<<t)*q+589 for q in (8080,8104,8128)])
print('6 literal fixed-interface cores, 45 whole-tournament route counts, 2-arc edit ledgers')
