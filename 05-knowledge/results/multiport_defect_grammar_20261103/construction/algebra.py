"""Exact two-birth operator in the square-zero subset ring."""
def need(condition, message):
    if not condition: raise RuntimeError(message)

def birth(rows, cut):
    n=len(rows);full=(1<<n)-1
    need(0<=cut<=full,'birth cut')
    return tuple(row|(((cut>>v)&1)<<n) for v,row in enumerate(rows))+(full^cut,)

def tables(rows, start=None):
    n=len(rows);size=1<<n;ends=[[0]*n for _ in range(size)];h=[0]*size;h[0]=int(start is None)
    for mask in range(1,size):
        for v in range(n):
            if not(mask>>v&1):continue
            rest=mask^(1<<v)
            if rest==0:ends[mask][v]=int(start is None or v==start)
            else:ends[mask][v]=sum(ends[rest][u] for u in range(n) if rows[u]>>v&1)
        h[mask]=sum(ends[mask])
    return h,ends

def product(a,b):
    c=[0]*len(a)
    for mask in range(len(c)):
        s=mask
        while True:
            c[mask]+=a[s]*b[mask^s]
            if s==0:break
            s=(s-1)&mask
    return c

def response(rows,w):
    p,e=tables(rows);W=[v[w] for v in e]
    return p[-1],product(W,p)[-1]

def induced(rows, vertices):
    return tuple(sum((rows[u]>>v&1)<<j for j,v in enumerate(vertices)) for u in vertices)

def strong(rows):
    full=(1<<len(rows))-1
    for rev in (False,True):
        reach=1
        while True:
            new=reach
            for v in range(len(rows)):
                if reach>>v&1:new|=(full^rows[v]^(1<<v)) if rev else rows[v]
            if reach==new:break
            reach=new
        if reach!=full:return False
    return True

def operator(rows,w,a):
    n=len(rows);full=(1<<n)-1;p,e=tables(rows);A,ae=tables(rows,a);S,unused=tables(rows,w)
    W=[v[w] for v in e];E=[v[a] for v in e];Q=[v[w] for v in ae]
    pa=product(p,A);wp=product(W,p)
    L=product(wp,A)[-1];K=product(product(p,p),Q)[-1];J=product(W,E)[-1]
    X=birth(rows,full^(1<<a));C2=birth(X,full);C4=birth(X,1<<w)
    hx=p[-1]-E[-1]+pa[-1];fx=wp[-1]-J+L+K
    sx=S[-1]-tables(rows,w)[1][-1][a]+product(S,A)[-1]
    deletion=induced(X,[v for v in range(n+1) if v!=w]);D=tables(deletion)[0][-1]
    predicted2=(p[-1]+2*pa[-1],wp[-1]+2*L+2*K)
    predicted4=(hx-sx+fx,2*fx-D)
    return dict(rows=list(rows),w=w,a=a,h=p[-1],F=wp[-1],S_a_2=pa[-1],
        e_a=E[-1],L=L,K=K,J=J,H_X=hx,F_X=fx,S_w_X=sx,H_X_minus_w=D,
        X=list(X),C2=list(C2),C4=list(C4),predicted2=predicted2,predicted4=predicted4)
