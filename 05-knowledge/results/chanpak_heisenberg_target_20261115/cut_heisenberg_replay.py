from itertools import combinations
from collections import deque
from fractions import Fraction as F
from math import gcd
from random import Random

def all_rows(n):
    edges=list(combinations(range(n),2))
    for mask in range(1<<len(edges)):
        r=[0]*n
        for k,(i,j) in enumerate(edges):
            if mask>>k&1:r[i]|=1<<j
            else:r[j]|=1<<i
        yield tuple(r)

def dp(r):
    n=len(r); full=(1<<n)-1
    start=[[0]*n for _ in range(full+1)]
    end=[[0]*n for _ in range(full+1)]
    h=[1]+[0]*full
    for s in range(1,full+1):
        for i in range(n):
            if s>>i&1:
                t=s^(1<<i)
                start[s][i]=sum(start[t][j] for j in range(n) if r[i]>>j&1) if t else 1
                end[s][i]=sum(end[t][j] for j in range(n) if r[j]>>i&1) if t else 1
        h[s]=sum(start[s])
    return h,start,end

def augment(r,mask):
    n=len(r)
    return tuple(a | (0 if mask>>i&1 else 1<<n) for i,a in enumerate(r))+(mask,)

def ordinal(r,s):
    n=len(r)
    return tuple(a|(((1<<len(s))-1)<<n) for a in r)+tuple(b<<n for b in s)

def state(r,mark):
    h,st,_=dp(r); full=len(h)-1
    starts=[sum(st[s][i] for i in range(len(r)) if mark>>i&1)*h[full^s] for s in range(1,full+1)]
    bias=sum(v*(-1)**(s.bit_count()+1) for s,v in enumerate(starts,1))
    return (h[-1],dp(augment(r,mark))[0][-1],sum(starts),bias)

# (d,u,b,c) acts (x,y) -> (d*x+u,y+b*x+c).
def prod(g,h,m=None):
    d,u,b,c=g; e,v,f,k=h
    out=(d*e,d*v+u,b*e+f,c+k+b*v)
    return tuple(x%m for x in out) if m else out

def inv(g):
    d,u,b,c=g
    return (1/F(d),-u/F(d),-b/F(d),-c+b*u/F(d))

I=(1,0,0,0); A=(2,0,0,0); B=(4,2,-1,0)
def pw(g,k):
    if k<0:return pw(inv(g),-k)
    s=I
    for _ in range(k):s=prod(g,s)
    return s

def main():
    C=prod(pw(A,-2),B)
    C0=C
    C1=prod(prod(A,C),inv(A))
    Cm=prod(prod(inv(A),C),A)
    comm=prod(prod(prod(C1,C0),inv(C1)),inv(C0))
    ctr=prod(prod(pw(Cm,2),pw(C0,-5)),pw(C1,2))
    root=prod(ctr,pw(comm,-4))
    if (comm,ctr,root)!=((1,0,0,F(3,4)),(1,0,0,F(7,2)),(1,0,0,F(1,2))):
        raise ValueError((comm,ctr,root))
    print('Exact central translations:', *(g[-1] for g in (comm,ctr,root)))

    for m in (3,5,7,9,11,15):
        Q=deque([I]); seen={I}
        while Q:
            g=Q.popleft()
            for a in (A,B):
                h=prod(a,g,m)
                if h not in seen:seen.add(h);Q.append(h)
        slopes={g[0] for g in seen}
        if len(seen)!=len(slopes)*m**3//gcd(3,m) or len({(g[1],g[3]) for g in seen})!=m*m:
            raise ValueError((m,len(seen)))
        if m%3==0 and any((u-b)%3 for d,u,b,c in seen if d==1):
            raise ValueError(('3 obstruction',m))
        print('Reduced group / origin orbit:',m,len(seen),m*m)

    count=0
    for n in range(1,6):
        for r in all_rows(n):
            h,st,en=dp(r); full=len(h)-1
            q=sum(h[s]*h[full^s] for s in range(full+1))
            if (q-(-1)**n*h[-1])%3:raise ValueError(('q',r))
            for i in range(n):
                u=sum(st[s][i]*h[full^s] for s in range(1,full+1))
                v=sum(en[s][i]*h[full^s] for s in range(1,full+1))
                if (u-(-1)**(n-1)*v)%3:raise ValueError(('marked',r,i,u,v))
            count+=1
    print('All labeled terminal-block mod-3 sidecars through order 5:',count)

    rr=Random(613)
    pools={n:list(all_rows(n)) for n in range(1,5)}
    for _ in range(400):
        seed=rr.choice(pools[rr.randrange(1,5)])
        mark=rr.randrange(1<<len(seed))
        seed=ordinal(seed,(0,))
        h,n,s,d=state(seed,mark)
        if d:raise ValueError(('bias',seed,mark))
        for word,g in (((0,),A),((1,0),B)):
            z,bits=seed,mark
            for k in word:
                bits|=k<<len(z)
                z=ordinal(z,(0,))
            hh,nn,ss,dd=state(z,bits)
            a,u,b,c=g
            if (hh,dd,F(ss,hh),F(nn-ss,hh))!=(h,0,a*F(s,h)+u,F(n-s,h)+b*F(s,h)+c):
                raise ValueError(('word',word,seed,mark))
    print('Literal append controls for A and B:',800)

if __name__=='__main__':main()
