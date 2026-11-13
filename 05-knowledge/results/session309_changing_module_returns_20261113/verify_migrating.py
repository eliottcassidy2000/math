from functools import lru_cache
from math import gcd

def need(ok):
    if not ok: raise ValueError('failed clock check')

@lru_cache(None)
def cover(rows,t):
    size=1<<len(rows)
    dp=[[0]*len(rows) for _ in range(size)]
    for i in range(len(rows)):dp[1<<i][i]=1
    for mask in range(1,size):
        for i,c in enumerate(dp[mask]):
            if not c:continue
            rest=rows[i]&~mask
            while rest:
                b=rest&-rest;rest-=b
                dp[mask|b][b.bit_length()-1]+=c
    h=[sum(r) for r in dp];h[0]=1
    if t==1:return h[-1]
    w=h[:]
    for _ in range(1,t):
        u=[0]*size
        for s in range(size):
            r=s
            while True:
                u[s]+=h[r]*w[s^r]
                if not r:break
                r=(r-1)&s
        w=u
    return w[-1]

S=('a','v','c','u','d0','d1','x0','f0','f1')
def graph(m,n,reverse=()):
    items=[('a','a')]+[(x,'v') for x in n]+[('c','c')]+[(x,'u') for x in m]+[(x,x) for x in S[4:]][:5]
    arcs=set()
    for i,(x,l) in enumerate(items):
        for y,r in items[i+1:]:
            if (l,r) in {('a','d0'),('a','d1'),('c','f0'),('c','f1')}:
                arcs.add((y,x))
            else:arcs.add((x,y))
    for x,y in reverse:
        arcs.remove((x,y));arcs.add((y,x))
    return arcs


m=('m0','m1');u=('u0','u1','u2');n0=('n0',);r=('r0','r1')
def transplant(m,n0,u):
    reverse=((u[0],u[-1]),) # the moved triple is a directed triangle
    before=graph(m,n0+u,reverse)
    for x in u:
        before.remove((x,'c'));before.add(('c',x))
    need(before==graph(u+m,n0,reverse))
transplant(m,n0,u)
transplant(m,(),u)
rev=((u[0],u[-1]),)
old=graph(u+m,(),rev)
new=graph(u+m,r,rev)
oldnames={x for e in old for x in e}
need({(x,y) for x,y in new if x in oldnames and y in oldnames}==old)
need(all(('a',x) in new for x in r))
need(all((x,y) in new for x in r for y in oldnames-{'a'}))
U=(40,1,3,6,15,30)
need([cover(U,i) for i in (1,2,3)]==[19,256,1611])
need(gcd(19,15360*1611)==1)
need((1611*pow(19,-1,7))%7==3 and len({pow(3,i,7) for i in range(6)})==6)
need((15360+19)%49==42)
need(pow(1611,42,49)==pow(19,42,49)==1)
need(pow(1611,6,49)==8 and pow(19,6,49)==1)
need((15360*pow(1611,42,13)+pow(19,43,13))%13==1)
print('PASS: literal module crossing/refill and exact (19,256,1611) clock')
