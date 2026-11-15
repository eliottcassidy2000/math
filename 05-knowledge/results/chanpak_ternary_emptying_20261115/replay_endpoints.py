from math import gcd
from functools import reduce
from random import Random

def check(p,*msg):
    if not p:raise AssertionError(msg)

def data(rows):
    n=len(rows); full=(1<<n)-1
    H=[1]+[0]*full
    start=[[0]*n for _ in H]; end=[[0]*n for _ in H]
    for A in range(1,full+1):
        for v in range(n):
            if not (A>>v&1):continue
            B=A^(1<<v)
            if B==0:start[A][v]=end[A][v]=1;continue
            start[A][v]=sum(start[B][u] for u in range(n) if rows[v]>>u&1)
            end[A][v]=sum(end[B][u] for u in range(n) if rows[u]>>v&1)
        H[A]=sum(end[A])
    return H,start,end

def slot(H,a):
    Q=[1]+[0]*(len(H)-1)
    for _ in range(abs(a)):
        W=[]
        for A in range(len(H)):
            B=A; x=0
            while True:
                x+=H[B]*Q[A^B]
                if not B:break
                B=(B-1)&A
            W.append(x)
        Q=W
    if a<0:
        Q=[(-1)**A.bit_count()*x for A,x in enumerate(Q)]
    return Q

def obs(paths,Q,v,A=None,b=1):
    if A is None:A=len(paths)-1
    return sum(b**(B.bit_count()-1)*paths[B][v]*Q[A^B] for B in range(1,A+1) if B&A==B and B>>v&1)

def arrows(rows):
    return [(p,q) for p,r in enumerate(rows) for q in range(len(rows)) if (r>>q)&1]

def build(n,arcs):
    rows=[0]*n
    for p,q in arcs:rows[p]|=1<<q
    for i in range(n):
        check(rows[i]>>i&1==0)
        for j in range(i):check((rows[i]>>j&1)+(rows[j]>>i&1)==1)
    return rows

def ordinal(*factors):
    arcs=[]; offset=0; fs=[]
    for g in factors:
        arcs+= [(p+offset,q+offset)for p,q in arrows(g)]
        fs.append((offset,len(g)));offset+=len(g)
    for a,n in fs:
        arcs += [(p,q) for p in range(a,a+n) for q in range(a+n,offset)]
    return build(offset,arcs)

def lifted(core,v):
    n=len(core)
    Arcs=[(u+1,w+1) for u,w in arrows(core)]
    Arcs += [(0,u+1) for u in range(n)]
    Arcs += [(n+1,0),(n+1,v+1)]+[(u+1,n+1) for u in range(n) if u!=v]
    return build(n+2,Arcs)

CORE_SWITCH=((52,17,2,7,12,30),(42,20,1,6,13,30))
for rows,v,want in zip(CORE_SWITCH,(4,3),((50,333),(56,327))):
    H,s,t=data(rows); Q=slot(H,2)
    unit=H[-1],Q[-1],obs(s,H,v),t[-1][v],s[-1][v],obs(t,H,v)
    check(unit==(29,320,43,5,4,want[0]),unit)
    find=data(lifted(rows,v))[0][-1]
    check(find==want[1],find)
    print('seam:',unit,'output', find)

C3=(2,4,1)
FIVE=(6,12,8,1)
FORTY_FIVE_STRONG=(4,1,130,7,143,159,63,75)
FIFTEEN_STRONG=(6,12,24,17,3)
FORTY_FIVE_SMALL=(48,49,3,7,12,28)
SIXTY_THREE=(238,220,216,16,225,206,8,72)
for rows in (FIVE,ordinal(FIVE,C3),ordinal(FIVE,C3,C3),FIFTEEN_STRONG,FORTY_FIVE_SMALL,FORTY_FIVE_STRONG,SIXTY_THREE):
    H,s,t=data(rows); Q=slot(H,2)
    start=[obs(s,H,v) for v in range(len(rows))]; end=[obs(t,H,v) for v in range(len(rows))]
    C=reduce(gcd,start+end)
    val=0
    while C%3==0:val+=1;C//=3
    print('triangle / genuine core:',len(rows),H[-1],Q[-1],val, [w%3 for w in start], [w%3 for w in end])

rng=Random(451563)
checks=0
for n in range(1,8):
    for _ in range(6):
        rows=build(n,[(x,y) if rng.getrandbits(1) else (y,x) for x in range(n) for y in range(x+1,n)])
        H,s,t=data(rows)
        X=slot(H,1);Y=slot(H,2)
        for a in range(-2,n):
            qa=slot(H,a); qu=slot(H,-a-1); qv=slot(H,1-a)
            for A in range(1,1<<n):
                sign=(-1)**(A.bit_count()+1)
                for v in range(n):
                    if not A>>v&1:continue
                    check(obs(t,qa,v,A)==sign*obs(s,qu,v,A))
                    check(obs(t,qa,v,A,-1)==sign*obs(s,qv,v,A,-1))
                    checks+=1
print('exact reflected endpoint identities on induced subsets:',checks)
