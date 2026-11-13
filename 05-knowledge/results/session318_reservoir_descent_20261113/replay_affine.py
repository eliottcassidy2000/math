from itertools import permutations

def tab(rows):
    n=len(rows); N=1<<n
    d=[[0]*n for _ in range(N)]
    h=[1]+[0]*(N-1)
    for i in range(n): d[1<<i][i]=1
    for m in range(1,N):
        h[m]=sum(d[m])
        for i,c in enumerate(d[m]):
            if c:
                a=rows[i]&~m
                while a:
                    b=a&-a; a-=b; d[m|b][b.bit_length()-1]+=c
    return d,h

def params(rows,u):
    n=len(rows); A=(1<<n)-1
    end,h=tab(rows)
    start,_=tab([A^(1<<i)^r for i,r in enumerate(rows)])
    starts=[sum(c for i,c in enumerate(t) if u>>i&1) for t in start]
    ends=[sum(c for i,c in enumerate(t) if not (u>>i&1)) for t in end]
    q=sum(h[i]*h[A^i] for i in range(A+1))
    S=sum(starts[i]*h[A^i] for i in range(A+1))
    E=sum(ends[i]*h[A^i] for i in range(A+1))
    b=sum(ends[i]*starts[A^i] for i in range(1,A))
    c=starts[A]+ends[A]+b
    return h[A],q,starts[A],ends[A],S,E,c

def require(x):
    if not x: raise RuntimeError('replay mismatch')

def add(rows, cut):
    n=len(rows)
    return tuple(r|((not(cut>>i&1))<<n) for i,r in enumerate(rows))+(cut,)

def join(A,B):
    n=len(A)
    return tuple(a|(((1<<len(B))-1)<<n) for a in A)+tuple(b<<n for b in B)

P=(4,65,50,7,11,27,61,127)
for u in (201,177): print('P',u,params(P,u))
for P,u in [((16,1,3,7,14,31,63),91),((0,33,3,7,15,29,63),101)]: print('short',u,params(P,u))

def trans(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def cyc(m):
    # w=0, a=1, b_j=j+1
    return (2,((1<<m)-1)<<2)+tuple((r<<2)|1 for r in trans(m))

def defect(B):
    if not B: return 0
    x,u=divmod(B,2)
    return defect(x)+u*(x+x%2)

def strong(rows):
    full=(1<<len(rows))-1
    for i in range(len(rows)):
        seen=0; todo=1<<i
        while todo:
            bit=todo&-todo;todo-=bit;seen|=bit
            todo|=rows[bit.bit_length()-1]&~seen
        if seen!=full:return False
    return True

def literal_h(P,p,Z,z,B):
    n=B.bit_length()
    bp=int(f"{B:0{n}b}"[::-1],2) if n else 0
    core=join(join(trans(n),P),Z)
    cut=bp|(p<<n)|(z<<(n+len(P)))
    full=add(core,cut)
    require(strong(full))
    return tab(full)[1][-1]

P0=(16,1,3,7,14,31,63); p0=91
P1=(0,33,3,7,15,29,63); p1=101
require(params(P0,p0)==(9,336,9,2,227,82,73))
require(params(P1,p1)==(9,336,9,0,227,82,73))
for m in range(1,8):
    Z=cyc(m); require(strong(Z))
    mask=(1<<len(Z))-1
    for j in range(1,m+1):
        u=mask^(1<<(j+1))
        h,q,s,e,S,E,c=params(Z,u)
        require(h==2**m+1 and h-s==2**(m-j))
        require(E==2*3**(j-1)*2**(m-j)+2**j-(j==m))
    require(params(Z,mask^2)[2:][0]==2**m)
    h,q,s,e,S,E,c=params(Z,mask^2)
    require(E==h and c==2*h-1)
    if m<=3:
        for B in (1,2,5):
            ell=B%2; D=B-defect(B)-ell
            form=(2**m+1)*(309+418*B+9*D)
            require(literal_h(P0,p0,Z,mask^2,B)==form-2)
            require(literal_h(P1,p1,Z,mask^2,B)==form)
        u=mask^(1<<2)^(1<<(m+1)) if m>1 else mask^(1<<2)
        r=params(Z,u)[0]-params(Z,u)[2]
        require(literal_h(P1,p1,Z,u,1)-literal_h(P0,p0,Z,u,1)==2*r)
print('pure end, unit/chain apertures and dyadic formulas verified')

P=(4,65,50,7,11,27,61,127)
require(params(P,201)==(91,2392,91,36,1679,1058,921))
require(params(P,177)==(91,2392,91,39,1682,1058,919))
for m in range(1,5):
    R=cyc(m); Z=add(R,((1<<len(R))-1)^1); require(strong(Z))
    mask=(1<<len(Z))-1; u=mask^(1<<(len(Z)-1)); t=2**m
    h,q,s,e,S,E,c=params(Z,u)
    require((h,s,E,c)==(3*t-1,2*t-1,3*t,5*t-1))
    if m<=3:
        for p,cc in ((201,921),(177,919)):
            W=join(P,Z)
            v=p|(u<<len(P)); H,Q,ss,ee,SS,EE,C=params(W,v)
            require((H,ss,EE,C)==(273*t-91,273*t-91,10350*t-1058,8037*t-cc))
print('three-coordinate cancellation and literal growing terminal counts verified')
