"""All-length finite gcd certificates, plus independent path DP for the small marked cores."""
from collections import defaultdict, deque
from itertools import product
from math import gcd

bits=tuple(product(range(2),repeat=2))

def require(b,msg=None):
    if not b:
        raise AssertionError(msg)

def core(R,mask):
    f=(1<<len(R))-1
    H=[]
    E=[]
    for rev in (False,True):
        r=tuple(x^(f^(1<<i)) for i,x in enumerate(R)) if rev else R
        A=[[0]*len(r) for _ in range(f+1)]
        B=[1]+[0]*f
        for u in range(1,f+1):
            for i,v in enumerate(r):
                if (u>>i)&1:
                    w=u^(1<<i)
                    A[u][i]=sum(A[w][j] for j in range(len(r)) if (v&w)>>j&1) if w else 1
            B[u]=sum(A[u])
        H.append((B[-1],sum(B[u]*B[f^u] for u in range(f+1))))
        E.append(sum(B[f^u]*sum(A[u][j] for j in range(len(r)) if mask>>j&1)
                 for u in range(1,f+1)))
    require(H[0]==H[1],H)
    return H[0]+tuple(E)

def live(h,Q,O):
    F={}; R=defaultdict(list)
    todo=deque([O]); dist={O:0}
    while todo:
        u=todo.popleft()
        s,r,p,q=u
        F[u]=[]
        for z,w in bits:
            a,b=2*s+h*z-Q*w,r+Q*z-h*w
            if not (-h<=a<=Q and not b%2 and -h<=b//2<=Q):
                continue
            v=(a,b//2,z,w)
            c=(1-z)*s-w*(Q-h-b//2)+h*(p*(1-z)-q*(1-w))
            F[u].append((v,c))
            R[v].append(u)
            if v not in dist:
                dist[v]=dist[u]+c
                todo.append(v)
    todo=deque([O]); C={O}
    while todo:
        for v in R[todo.popleft()]:
            if v not in C:
                C.add(v)
                todo.append(v)
    g=0
    for u in C:
        for v,c in F[u]:
            if v in C:
                g=gcd(g,dist[u]+c-dist[v])
    return F,C,g

def verify(R,mask,H,O,goals):
    h,Q,s,r=core(R,mask)
    require((h,Q,s,r)==H)
    F,C,g=live(h,Q,(s,r,*O))
    output=[]
    for z,w,D in goals:
        u=(s,r,*O); d=0
        for a,b in zip(map(int,z),map(int,w)):
            v,c=next((v,c) for v,c in F[u] if v[2:]==(a,b))
            require(v in C)
            u=v
            d+=c
        require((u,d)==((s,r,*O),D))
        output.append(d)
    print('profile',H,'memory',O,'live vertices',len(C),'gcd',g,'verified shifts',output)
    return g

require(verify((6,12,24,1,11),2,(13,104,18,18),(0,1),[
  ('101001000010','010000001101',10),('111011000010','010010011101',-10)])==10)
require(verify((4,1,130,7,143,159,63,75),35,(45,1440,517,818),(1,1),[
  ('000100011001','010111000111',2),('011100011001','010111110111',-16)])==2)
