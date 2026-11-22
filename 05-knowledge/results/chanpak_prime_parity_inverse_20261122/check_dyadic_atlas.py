"""Standalone exact replay for the parity-refined two-port atlas (stdlib only)."""
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
def u_to_reverse(k,u):
    return int(format(u,'0%db'%k)[::-1],2) if k else 0
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


moves=[
 {},
 {'A':(82,60),'B':(5976,4396),'C':(96264,71856),'E':(55833120,41679296)},
 {'P':(292,120),'Q':(5976,2476)},
 {'D':(82,24)},
]
old=[
 (2,1,20311,14572,8654,6207),
 (2,2,43209,17882,18144,7507),
 (1,3,29214,8633,21116,6239),
]
extra=[
 (2,1,18800,13855,7947,5855),
 (2,1,25558,18428,10801,7787),
]
middle=[
 (5855,(208,25,3,21,164,15,62,110),(425,623,570,384,377)),
 (5857,(90,48,35,134,44,137,62,87),(455,611,560,265,500)),
 (7787,(60,77,72,192,14,30,177,55),(459,650,810,661,441)),
 (7789,(142,28,80,52,97,71,139,62),(489,638,800,418,688)),
]
single_arc=[
 (411,(10,192,35,134,15,27,61,117),5,6),
 (479,(24,37,81,38,138,149,59,79),4,5),
 (551,(72,37,145,134,43,77,30,115),5,7),
]
def phase(move,m):
    image=lambda r:{(a*x+b)%m for x in r for a,b in move.values()}
    s={0}
    while (t:=s|image(s))!=s:s=t
    r=s
    while (t:=image(r))!=r:r=t
    return s,r
def images(m,charts,both=False):
    phases=[None]+[phase(move,m) for move in moves[1:]]
    out=set()
    rep=set()
    mask=m-1
    for q,t,a,b,c,d in charts:
        U,P=phases[q]
        V,Q=phases[t]
        for u in (P if both else U):
            slope=a*u+c
            constant=b*u+d
            for v in (Q if both else V):
                w=(slope*v+constant)&mask
                out.add(w)
                out.add((w+2)&mask)
                if both or u in P or v in Q:
                    rep.add(w)
                    rep.add((w+2)&mask)
    return out,rep
def run(i,w,h=0):
    for letter in w:
        a,b=moves[i][letter]
        h=a*h+b
    return h
def multiplier(i,w):
    result=1
    for letter in w:result*=moves[i][letter][0]
    return result
def v2(z):
    return (z&-z).bit_length()-1

def inverse(h):
    rev=[]
    while h!=1:
        if h<=0 or h%4!=1:return None
        z=27*h-7
        n=v2(z)-(3 if h%3==2 else 2)
        if n<0:return None
        rev+=['A']*n
        for _ in range(n):
            h,s=divmod(h+21,82)
            if s or h<1:return None
        if h==1:break
        l='ECB'[h%3]
        a,b=dict(B=(5976,1579),C=(96264,24407),E=(55833120,14153823))[l]
        h,s=divmod(h+b,a)
        if s or h<1:return None
        rev.append(l)
    return ''.join(rev[::-1])

def check():
    E=((1,0),(0,1))
    for k,i in ((6,8),(9,72),(10,340)):
        assert i<1<<(k-1)
        E=mat_product(affine(k,i),E)
    assert E==((13371651,14153823),(40114950,42461470))
    assert tuple(3*E[0][j]-E[1][j] for j in (0,1))==(3,-1)
    assert v2(55833120)==5 and v2(8678612)==2
    for n in range(6):
        for w in map(''.join,product('ABCE',repeat=n)):
            assert inverse(1+run(1,w))==w
    for j,(d,A,want) in enumerate(middle):
        assert all(((A[i]>>j)&1)^((A[j]>>i)&1) for i in range(8) for j in range(i))
        S=1<<1
        while (Z:=S|sum(1<<k for k in range(8) if any((S>>t)&1 and (A[t]>>k)&1 for t in range(8))))!=S:S=Z
        assert S&1
        sig=profile(A)
        assert sig==want
        g,r,l,m1,m2=sig
        assert coefficients(2,1,g,r,l,m1+m2)==(*extra[j//2][2:5],d)
        for u,v in ((2,2),(2,3),(3,2),(3,3)):
            Z=literal(T(u),A,T(v))
            assert Z==g+2**(u-1)*r+2**(v-1)*l+2**(u+v-2)*(m1+m2)
    for n,A,x0,y0 in single_arc:
        B=list(A);B[x0]^=1<<y0;B[y0]^=1<<x0
        assert counts(A)[0][-1]==n and counts(B)[0][-1]==n+2

    for m,i,l,r,n0 in [
        (512,2,'CEA','D',411),
        (512,0,'','AEC',479),
        (1024,3,'QPQP','BACB',411),
        (1024,4,'PPQQ','EBB',551),
    ]:
        q,t,a,b,c,d=(old+extra)[i]
        x,y=run(q,l),run(t,r)
        N=Hplane(a,b,c,d,x,y)
        assert N%m==n0 and (multiplier(q,l)%m==0 or multiplier(t,r)%m==0)
        if m==1024:assert multiplier(q,l)%m==multiplier(t,r)%m==0
        print(m,n0,n0+2,':',N,N+2,flush=True)

    for m,bank,both,excluded in [
        (512,old,False,set()),
        (1024,old,False,{411,413,551,553}),
        (2**18,old+extra,False,set()),
        (2**19,old+extra,False,{141915,141917}),
        (2**13,old[:2]+extra,True,set()),
        (2**14,old[:2]+extra,True,{13915,13917}),
    ]:
        img,rep=images(m,bank,both)
        assert set(range(1,m,2))-rep==excluded
        if m==2**19:assert img==rep
        print('modulus',m,'charts',len(bank),'two freely growing flanks',both,
              'covered',len(rep),'excluded',sorted(excluded),flush=True)

    Fstar=((1,0),(0,1))
    for k,j in ((9,169),(6,20),(3,2)):
        assert j<1<<(k-1)
        Fstar=mat_product(affine(k,j),Fstar)
    assert Fstar==((70557,167790),(176390,419476))
    assert tuple(5*Fstar[0][i]-2*Fstar[1][i] for i in (0,1))==(5,-2)
    assert {(292*h-51)%16 for h in range(1,64,8)}=={1}
    assert {(5976*h-1023)%64 for h in range(1,64,8)}=={25}
    assert {(490032*h-83895)%64 for h in range(1,64,8)}=={57}
    def inverse_pqf(h):
        digits=[]
        while h!=1:
            if h<1 or h%8!=1:return None
            last='P' if h%16==1 else {25:'Q',57:'F'}.get(h%64)
            if not last:return None
            a,b={'P':(292,51),'Q':(5976,1023),'F':(490032,83895)}[last]
            h,r=divmod(h+b,a)
            if r:return None
            digits.append(last)
        return ''.join(digits[::-1])
    moves[2]['F']=(490032,203068)
    for d in range(6):
        for w in map(''.join,product('PQF',repeat=d)):
            assert inverse_pqf(1+2*run(2,w))==w
    m=2**19
    u,v=run(2,'FFPQ'),run(2,'QQFPFQ')
    n=Hplane(43209,17882,18144,7507,u,v)
    assert (u%m,v%m,n%m)==(76140,167500,141915)
    assert multiplier(2,'QQFPFQ')%m==0
    print('F* repairs final mod-2^19 pair:',n,n+2,flush=True)

if __name__=='__main__':check()
