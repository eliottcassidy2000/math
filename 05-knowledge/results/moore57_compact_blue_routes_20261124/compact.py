#!/usr/bin/env python3
"""Compact ordinary Moore encodings, and their checks."""

def arcs(N, f):
    out=[0]*N
    for i in range(N):
        for j in range(i+1,N):
            u,v=(i,j) if f(i,j) else (j,i)
            out[u]|=1<<v
    return out

def forward(T,i,j):
    return bool(T[i]>>j&1)

def field_host(q):
    """q=1: Petersen. q=5: Hoffman–Singleton. Coordinate groups over F5."""
    a=5*q
    def coord(v):return divmod(v%a,5)
    def wins(u,v):
        s,t=u//a,v//a
        i,x=coord(u);j,y=coord(v)
        if s==t:
            if i!=j:return (i<j) if s==0 else (i>j)
            return (y-x)%5 in (1,2)
        assert s==0 and t==1
        return (y-x-i*j)%5!=0
    T=arcs(2*a,wins)
    def alpha(u):
        i,x=coord(u)
        return (1-u//a)*a+5*i+((-x)%5)
    first=list(range(a-5))+list(range(a-5,a))
    i,x=coord(first[-1])
    if (2*x+i*i)%5==0:
        first[-5:]=first[-4:]+first[-5:-4]
    spine=first+[alpha(u) for u in first[::-1]]
    assert all(forward(T,u,v) for u,v in zip(spine,spine[1:]))
    assert all(forward(T,u,v)==forward(T,alpha(v),alpha(u))
               for u in range(2*a) for v in range(2*a) if u!=v)
    return T,spine

def decode_field_host(T,q):
    a=5*q
    s=[u.bit_count() for u in T]
    if q==1:
        U={i for i in range(2*a) if s[i]==6}; W=set(range(2*a))-U
    else:
        w0={i for i in range(2*a) if s[i]==7}
        assert len(w0)==5
        U={i for i in range(2*a) if s[i] in (32,37,42) or
           (s[i] in (22,27) and sum(forward(T,i,j) for j in w0)==4)}
        W=set(range(2*a))-U
    assert len(U)==len(W)==a
    # Read each five-set by score within its now identified side.
    cohorts=[]
    for side in [U,W]:
        cohorts.append([sorted(v for v in side if s[v]==z) for z in sorted({s[v] for v in side})])
    E={i:set() for i in range(2*a)}
    def connect(i,j):E[i].add(j);E[j].add(i)
    for side,C in enumerate(cohorts):
        for S in C:
            assert len(S)==5
            def skew(i,j):return 0 if i==j else 1 if forward(T,i,j) else -1
            for t,i in enumerate(S):
                for j in S[t+1:]:
                    score=sum(skew(i,k)*skew(j,k) for k in S)
                    assert score in (1,-3)
                    if score==(1,-3)[side]:connect(i,j)
    for i in U:
        for j in W:
            if forward(T,j,i):connect(i,j)
    return E

def verify_moore(E,d):
    assert len(E)==d*d+1
    assert all(len(v)==d for v in E.values())
    assert all(len(E[i]&E[j])==(j not in E[i]) for i in E for j in E if i!=j)

def beacon_host(E,order):
    """Unique twin markers + binary positions. A Moore Hamilton spine is optional."""
    n=len(E);r=n.bit_length()
    assert n>=5 and len(order)==n and set(order)==set(E)
    def wins(i,j):
        if j<n:return order[j] not in E[order[i]]
        if i>=n+2:return True
        if j<n+2:return True
        if i>=n:return True
        return bool(i>>(j-n-2)&1)
    return arcs(n+r+2,wins)

def read_beacon(T):
    N=len(T)
    twins=[(i,j) for i in range(N) for j in range(i+1,N)
           if (T[i]^T[j])&~((1<<i)|(1<<j))==0]
    assert len(twins)==1
    a,b=twins[0]
    V={i for i in range(N) if i not in (a,b) and forward(T,i,a)}
    W=[i for i in range(N) if i not in (a,b) and forward(T,a,i)]
    W=sorted(W,key=lambda i:-sum(forward(T,i,j) for j in W))
    address={i:sum(1<<j for j,w in enumerate(W) if forward(T,i,w)) for i in V}
    assert set(address.values())==set(range(len(V)))
    return {address[i]:{address[j] for j in V if (address[i]<address[j] and forward(T,j,i))
                        or (address[i]>address[j] and forward(T,i,j))}
            for i in V}

def step_to_check(T):
    mask=(1<<len(T))-1
    for i in range(len(T)):
        U=1<<i
        while True:
            V=U
            for j in range(len(T)):
                if U>>j&1:V|=T[j]
            if V==U:break
            U=V
        assert U==mask

if __name__=='__main__':
    for q,d in ((1,3),(5,7)):
        T,p=field_host(q)
        E=decode_field_host(T,q)
        verify_moore(E,d)
        step_to_check(T)
        assert read_beacon(beacon_host(E,list(E)))=={i:{list(E).index(j) for j in E[v]} for i,v in enumerate(E)}
        print(f"degree {d}: strong blue plain tournament on {len(T)}, graph recovered")
        if d==3:print("Petersen rows:",tuple(T),"blue spine",p)
