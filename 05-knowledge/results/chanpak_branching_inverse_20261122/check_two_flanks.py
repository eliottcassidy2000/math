from math import isqrt

B0=(64,93,17,101,41,7,52)
B1=(18,116,97,39,12,81,25)
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

def encode(h,k,s):
    return 7464*h*k-2176*(h+k)+631+6*s

def address_pairs(n):
    s={31:0,5:1}.get(n%32)
    if s is None:return set()
    p=7464*(n-631-6*s)+2176**2
    if p<=0:return set()
    out=set()
    for a in range(1,isqrt(p)+1):
        if p%a:continue
        b=p//a
        h,r=divmod(a+2176,7464);k,t=divmod(b+2176,7464)
        if r or t:continue
        u,v=decode(h),decode(k)
        if u is not None and v is not None:out.update(((u,v,s),(v,u,s)))
    return out

if __name__=='__main__':
    for h,k,u,v in ((1,1,'',''),(1,61,'','A'),(1,4397,'','B'),(1,71857,'','C')):
        for sigma in (0,1):assert (u,v,sigma) in address_pairs(encode(h,k,sigma))
    if not (not address_pairs(3745)): raise ArithmeticError("two-flank inverse check")
    for q, expected in ((B0,(87,283,283,310,321)), (B1,(141,265,265,254,383))):
        g,r,l,z,w=profile(q);m=z+w
        if (g,r,l,z,w) != expected: raise ArithmeticError("middle signature")
        aa=g+3*r+3*l+9*m
        print('signature',g,r,l,z,w,m,'absolute',aa,l+3*m,r+3*m,m,'det',aa*m-(l+3*m)*(r+3*m))
        for u,v in ((2,2),(2,3),(3,2),(3,3)):
            T=lambda n: tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
            G=assemble(T(u),q,T(v),0,v-1)
            n=counts(G)[0][-1]
            if n != g+r*2**(u-1)+l*2**(v-1)+m*2**(u+v-2):
                raise ArithmeticError("transitive flank counts")
        p,h=ports('A')
        Q=assemble(rev(p),q,(2,0),len(p)-1,1)
        if counts(Q)[0][-1] != g*h+r*(3*h-1)+l*2*h+m*2*(3*h-1):
            raise ArithmeticError("migrating flank count")
    print("ten literal graft configurations checked")
