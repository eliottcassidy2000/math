import sys
from math import gcd, isqrt
from itertools import product
x,y=0,1
def submasks(s):
    full=s
    while 1:
        yield s
        if not s: break
        s=(s-1)&full

def counts(arcs):
    n=len(arcs);N=1<<n
    h=[0]*N;end=[[0]*n for _ in range(N)];xy=[[0]*n for _ in range(N)]
    h[0]=1
    for s in range(1,N):
        for j in range(n):
            if s>>j&1:
                d=s^(1<<j)
                if not d:
                    end[s][j]=1
                    xy[s][j]=int(j==y)
                else:
                    end[s][j]=sum(end[d][i] for i in range(n) if d>>i&1 and arcs[i]>>j&1)
                    xy[s][j]=sum(xy[d][i] for i in range(n) if d>>i&1 and arcs[i]>>j&1)
        h[s]=sum(end[s])
    return h,end,xy

def rev(a):
    return tuple(sum(1<<i for i,b in enumerate(a) if b>>j&1) for j in range(len(a)))

def profile(a):
    h,e,xy=counts(a);_,s,_=counts(rev(a));U=len(h)-1
    r=sum(e[V][x]*h[U^V] for V in range(U+1))
    l=sum(s[V][y]*h[U^V] for V in range(U+1))
    z=0;w=0
    for V in range(U+1):
        D=U^V
        z+=xy[V][x]*sum(h[P]*h[D^P] for P in submasks(D))
        w+=e[V][x]*sum(s[P][y]*h[D^P] for P in submasks(D))
    return h[-1],r,l,z,w

def assemble(A,B,C,i,o):
    # A before B before C, reverse only  x->i and o->y.
    ma,mb,mc=len(A),len(B),len(C);n=ma+mb+mc
    U=(1<<n)-1;V=(1<<(ma+mb))-1; W=(1<<ma)-1
    arcs=[u|(U^W) for u in A]+[(u<<ma)|(U^V) for u in B]+[u<<(ma+mb) for u in C]
    for u,v in ((i,ma+x),(ma+y,ma+mb+o)):
        if not arcs[u]>>v&1: raise ArithmeticError("missing initial forward arc")
        arcs[u]^=1<<v;arcs[v]^=1<<u
    return tuple(arcs)

def F(u):
    z=1
    while u>1:
        if not u&1:z+=u//2+int(bool(u&2))
        u//=2
    return z
def migration(a, k,u):
    # Write the old block first, even though its exterior dominates it.
    n=len(a);out=list(a);out[-1]|=1<<(n+k)
    for j in range(k):
        m=((1<<n)-1)+((1<<k)-(1<<(j+1)))*(1<<n)
        if not u>>(k-1-j)&1: m|=1<<(n+k)
        out.append(m)
    out.append(((1<<(n-1))-1)+(u_to_reverse(k,u)<<n))
    return tuple(out)
def u_to_reverse(k,u):
    return int(format(u,'0%db'%k)[::-1],2) if k else 0
def ports(w):
    a=(2,0);h=1
    for c in w:
        for k,u in dict(A=[(6,20)],B=[(4,2),(8,84)],C=[(6,8),(10,340)])[c]:
            a=migration(a,k,u)
        a_m,b_m=dict(A=(82,-21),B=(5976,-1579),C=(96264,-24407))[c]
        h=a_m*h+b_m
    return a,h

def decode(h):
    inverse=[]
    while h!=1:
        if h<=0:return None
        x=27*h-7
        if x%4:return None
        if h%3==2:r=3
        elif h%3==1:r=2
        else:return None
        a=(x&-x).bit_length()-1-r
        if a<0:return None
        inverse+=['A']*a
        for _ in range(a):
            h,z=divmod(h+21,82)
            if z or h<=0:return None
        if h==1:break
        if h%3==2:digit,m,b='B',5976,1579
        elif h%3==1:digit,m,b='C',96264,24407
        else:return None
        h,z=divmod(h+b,m)
        if z or h<=0:return None
        inverse.append(digit)
    return ''.join(reversed(inverse))

raw = [
    (2,1,((6,132,48,71,43,139,55,93),6,7),((10,48,35,134,77,153,47,87),7,4)),
    (2,2,((8,5,137,82,39,79,151,59),2,5),((10,96,11,98,15,85,149,63),5,4)),
    (1,2,((6,40,138,145,39,77,31,115),2,1),((10,144,67,22,69,31,171,61),1,0)),
    (1,3,((8,13,49,164,75,83,143,55),2,1),((32,13,137,97,15,150,55,91),4,2)),
    (3,1,((66,12,41,97,15,83,150,63),4,5),((66,40,11,193,15,29,54,119),5,2)),
]
def renumber(item):
    A,x,y=item
    order=[x,y]+[i for i in range(len(A)) if i not in (x,y)]
    out=tuple(sum(1<<i for i,v in enumerate(order) if A[u]>>v&1) for u in order)
    edges=set((i,j) for i in range(len(out)) for j in range(len(out)) if out[i]>>j&1)
    if not (all(((i,j) in edges)^((j,i) in edges) for i in range(len(out)) for j in range(i))): raise ArithmeticError("unit-switch certificate check")
    reached={1}
    while len(r:=reached|{j for i,j in edges if i in reached})>len(reached):reached=r
    if not (0 in reached): raise ArithmeticError("unit-switch certificate check")
    return out
def coefficients(q,t,g,r,l,m):
    u,v=2*q+1,2*t+1
    return (g*q*t+r*u*t+l*q*v+m*u*v, g*q+r*u+2*l*q+2*m*u, g*t+2*r*t+l*v+2*m*v,g+2*r+2*l+4*m)
def Hplane(A,B,C,d,u,v):return A*u*v+B*u+C*v+d
def T(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def literal(A,B,C):return counts(assemble(rev(A),B,C,len(A)-1,len(C)-1))[0][-1]

def affine(k,u):
    return ((F(u),u+1),(2**k-u_to_reverse(k,u),2**k))
def mat_product(X,Y):
    return tuple(tuple(sum(X[i][v]*Y[v][j] for v in range(2)) for j in range(2)) for i in range(2))

def decode_52(h):
    answer=[]
    while h!=1:
        if h<1:return None
        if h%16==1: d,k,b='P',292,51
        elif h%16==9: d,k,b='Q',5976,1023
        else:return None
        h,r=divmod(h+b,k)
        if r:return None
        answer.append(d)
    return ''.join(reversed(answer))
def addresses(N,a,bcoef,ccoef,d,q,t):
    z=a*(N-d)+bcoef*ccoef
    if z<=0:return set()
    choices=set()
    for f in range(1,isqrt(z)+1):
        if z%f:continue
        for fx,fy in ((f,z//f),(z//f,f)):
            u,r=divmod(fx-ccoef,a);v,s=divmod(fy-bcoef,a)
            if r or s or u<0 or v<0:continue
            h,k=1+q*u,1+t*v
            left,right=decode_52(h) if q==2 else decode(h),decode_52(k) if t==2 else decode(k)
            if left is not None and right is not None:choices.add((left,right))
    return choices

steps=[
    (),
    ((82,60),(5976,4396),(96264,71856)),
    ((292,120),(5976,2476)),
    ((82,24),),
]
def phases(m,j):
    step=lambda s:{(a*v+b)%m for v in s for a,b in steps[j]}
    s={0}
    while (t:=s|step(s))!=s:s=t
    r=s
    while (u:=step(r))!=r:r=u
    return s,r
def targets(m):
    all_targets=set()
    repeated=set()
    for q,t,a,b,c,d in (
        (2,1,20311,14572,8654,6207),
        (2,2,43209,17882,18144,7507),
        (1,3,29214,8633,21116,6239)
    ):
        U,P=phases(m,q)
        V,Q=phases(m,t)
        for u,v,s in product(U,V,(0,2)):
            n=Hplane(a,b,c,d+s,u,v)%m
            all_targets.add(n)
            if u in P or v in Q:repeated.add(n)
    return all_targets,repeated

def finite_only_addresses(m,j):
    stable=phases(m,j)[1]
    active=[0]
    out=[]
    while active:
        active=[x for x in active if x%m not in stable]
        out.extend(active)
        active=[a*x+b for x in active for a,b in steps[j]]
    return out

certs={
    601669: (2, {2: 2, 3: 5, 619: 1}),
    5745938951: (7, {2: 1, 5: 2, 191: 1, 601669: 1}),
    16183: (3, {2: 1, 3: 2, 29: 1, 31: 1}),
    194197: (5, {2: 2, 3: 1, 16183: 1}),
    109527109: (2, {2: 2, 3: 1, 47: 1, 194197: 1}),
    1095271091: (2, {2: 1, 5: 1, 109527109: 1}),
    8762168729: (3, {2: 3, 1095271091: 1}),
    4100694965173: (5, {2: 2, 3: 2, 13: 1, 8762168729: 1}),
    47124685853114258307047: (5, {2: 1, 5745938951: 1, 4100694965173: 1}),
    93559: (3, {2: 1, 3: 1, 31: 1, 503: 1}),
    16811: (7, {2: 1, 5: 1, 41: 2}),
    33623: (5, {2: 1, 16811: 1}),
    25165874057: (3, {2: 3, 93559: 1, 33623: 1}),
    171077611839487: (3, {2: 1, 3: 1, 11: 1, 103: 1, 25165874057: 1}),
    7869570144616403: (2, {2: 1, 23: 1, 171077611839487: 1}),
    47124685853114258307049: (11, {2: 3, 3: 1, 13: 1, 17: 1, 1129: 1, 7869570144616403: 1}),
}
def lucas(n):
    if n<=10000:return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
    a,fs=certs[n]
    f=1
    for q,v in fs.items():
        if not lucas(q): return False
        f*=q**v
    return f==n-1 and pow(a,n-1,n)==1 and all(gcd(pow(a,(n-1)//q,n)-1,n)==1 for q in fs)

if __name__=='__main__':
    P=affine(8,101);Q=mat_product(affine(3,2),affine(9,169))
    D=affine(6,26)
    if not (D==((19,27),(42,64))): raise ArithmeticError("unit-switch certificate check")
    if not (tuple(7*D[0][i]-3*D[1][i] for i in range(2))==(7,-3)): raise ArithmeticError("unit-switch certificate check")
    o,r=targets(256)
    if not (o==r==set(range(1,256,2))): raise ArithmeticError("unit-switch certificate check")
    print('mod 256: all 128 odd residues infinitely often',flush=True)
    o,r=targets(512)
    if not (set(range(1,512,2))-o=={411,413}): raise ArithmeticError("unit-switch certificate check")
    if not (o-r=={479,481}): raise ArithmeticError("unit-switch certificate check")
    print('mod 512: omitted 411,413; finite-only 479,481',flush=True)
    transient=[[],*(finite_only_addresses(512,j) for j in (1,2,3))]
    if not (list(map(len,transient))==[0,11,6,6]): raise ArithmeticError("unit-switch certificate check")
    found=set()
    for q,t,a,b,c,d in (
        (2,1,20311,14572,8654,6207),
        (2,2,43209,17882,18144,7507),
        (1,3,29214,8633,21116,6239)
    ):
        for u,v,s in product(transient[q],transient[t],(0,2)):
            n=Hplane(a,b,c,d+s,u,v)
            if n%512 in (479,481):found.add(n)
    if not (found=={621848031,621848033}): raise ArithmeticError("unit-switch certificate check")
    print('finite-only residues contain exactly:',*sorted(found),flush=True)
    u,v=120,0
    for c in 'AACAAB':
        a,b=steps[1]['ABC'.index(c)];v=a*v+b
    n=Hplane(20311,14572,8654,6207,u,v)
    if not (n==47124685853114258307047 and lucas(n) and lucas(n+2)): raise ArithmeticError("unit-switch certificate check")
    print('Lucas certificates passed:',n,n+2,flush=True)
    if not (P==((37,102),(90,256)) and Q==((861,2046),(2150,5116))): raise ArithmeticError("unit-switch certificate check")
    if not (all(tuple(5*M[0][i]-2*M[1][i] for i in range(2))==(5,-2) for M in (P,Q))): raise ArithmeticError("unit-switch certificate check")
    for i in range(4):
        for w in product('PQ',repeat=i):
            h=1
            for d in w:h=292*h-51 if d=='P' else 5976*h-1023
            if not (decode_52(h)==''.join(w)): raise ArithmeticError("unit-switch certificate check")
    for j,(q,t,first,second) in enumerate(raw):
        tab=[]
        for rawB in (first,second):
            B=renumber(rawB)
            g,r,l,m1,m2=profile(B)
            a,b,c,d=coefficients(q,t,g,r,l,m1+m2)
            tab.append((d,a,b,c,B,(g,r,l,m1,m2)))
        tab.sort()
        if not (tab[0][1:4]==tab[1][1:4] and tab[1][0]-tab[0][0]==2): raise ArithmeticError("unit-switch certificate check")
        print(q,t,tab,flush=True)
        for d,a,b,c,B,sig in tab:
            g,r,l,m1,m2=sig
            for u,v in ((2,2),(2,3),(3,2),(3,3)) if j<2 else ((2,2),):
                if not (literal(T(u),B,T(v))==g+2**(u-1)*r+2**(v-1)*l+2**(u+v-2)*(m1+m2)): raise ArithmeticError("unit-switch certificate check")
        if j<2:
            d,a,b,c,B,sig=tab[0]
            if not (('','') in addresses(d,a,b,c,d,q,t)): raise ArithmeticError("unit-switch certificate check")
            if not (('P','') in addresses(Hplane(a,b,c,d,120,0),a,b,c,d,q,t)): raise ArithmeticError("unit-switch certificate check")
