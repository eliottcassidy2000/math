#!/usr/bin/env python3
"""Literal controls for a changed-core pair with independently scalable dyadic wings."""
ROWS = (((0,5,9,19,39,15,63),81),((8,5,33,22,7,27),45))

def dp(r):
    n=len(r); f=(1<<n)-1
    a=[[0]*n for _ in range(f+1)]
    for i in range(n):a[1<<i][i]=1
    for s,w in enumerate(a):
        for i,val in enumerate(w):
            if not val:continue
            b=r[i]&~s
            while b:
                c=b&-b;b^=c
                a[s|c][c.bit_length()-1]+=val
    return a

def hc(r):return sum(dp(r)[-1])

def connect(r,u,x='',y=''):
    nu=len(r); n=len(x)+nu+len(y)
    S=((1<<n)-1)
    c=x+''.join(str(u>>i&1) for i in range(nu))+y
    a=[(S^((1<<(i+1))-1)) for i in range(n)]
    for i in range(nu):a[len(x)+i]=(r[i]<<len(x))|(S^((1<<(len(x)+nu))-1))
    return tuple(v|((1-int(c[i]))<<n) for i,v in enumerate(a))+(sum(int(v)<<i for i,v in enumerate(c)),)
def profile(r,U):
    n=len(r);F=(1<<n)-1
    a=dp(r)
    b=dp(tuple(F^x^(1<<i) for i,x in enumerate(r)))
    H=[1]+list(map(sum,a[1:]))
    f=sum(v*H[F^s] for s,v in enumerate(H))
    p=sum(sum(w[i] for i in range(n) if U>>i&1)*H[F^s] for s,w in enumerate(b))
    q=sum(sum(w[i] for i in range(n) if not U>>i&1)*H[F^s] for s,w in enumerate(a))
    k=hc(connect(r,U))
    return H[-1],f,p,q,k
def attach(i,s,t):
    if min(s,t)<1:raise ValueError('exponents must be positive')
    r,U=ROWS[i]
    return connect(r,U,'1'+'0'*s,'1'*t+'0')
def solve(N):
    for i,d in enumerate((101,103)):
        for s in range(1,max(1,N.bit_length())):
            u=1<<s
            v,rem=divmod(N-d-141*u,336*u+228)
            if not rem and v>=2 and v&(v-1)==0:yield i,s,v.bit_length()-1
def primary():
    actual=[profile(r,U) for r,U in ROWS]
    print(actual)
    assert actual==[(9,336,219,132,101),(33,336,195,108,103)]
    for s in range(1,4):
        for t in range(1,4):
            n=336*(1<<(s+t))+141*(1<<s)+228*(1<<t)+101
            for i in range(2):
                N=n+2*i
                assert hc(attach(i,s,t))==N
                assert (i,s,t) in list(solve(N))
    print('PASS: two exact profiles and eighteen explicit strong absolute targets.')
if __name__=='__main__':
    import sys
    if len(sys.argv)==1:primary()
    else:
        for i,s,t in solve(int(sys.argv[1])):print(i,s,t,attach(i,s,t))
