def require(v, msg):
    if not v: raise AssertionError(msg)
I=((1,0),(0,1)); T=((0,1),(8,9))
def prod(A,B,m=None):
    C=tuple(tuple(sum(A[i][k]*B[k][j] for k in (0,1)) for j in (0,1)) for i in (0,1))
    if m is not None: return tuple(tuple(x%m for x in r) for r in C)
    return C
def mpow(A,n,m=None):
    B=I
    while n:
        if n%2: B=prod(B,A,m)
        A=prod(A,A,m);n//=2
    return B
def adjdiff(A):
    return ((A[0][0]-1,A[0][1]),(A[1][0],A[1][1]-1))
def det(A):return A[0][0]*A[1][1]-A[0][1]*A[1][0]
def val(n,p):
    require(n!=0, 'unexpected zero')
    r=0
    while n%p==0: n//=p;r+=1
    return r
def vp_fact(n,p):
    a=0
    while n: n//=p;a+=n
    return a
def mlog_power(Tm,p,digits):
    modulus=p**digits
    B=adjdiff(Tm); L=[[0,0],[0,0]];P=I
    # k-v_p(k) >= digits for k >= 2*digits+1, p odd.
    for k in range(1,2*digits+1):
        j=val(k,p); pk=p**j
        P=prod(P,B,modulus*p**(digits+2))
        q=pow(k//pk,-1,modulus)
        for i in (0,1):
            for l in (0,1):
                require(P[i][l]%pk==0,'invalid log term')
                L[i][l]=(L[i][l]+(-1)**(k+1)*q*(P[i][l]//pk))%modulus
    return tuple(tuple(x for x in r) for r in L)
for p in (3,5,7,11,13,17,19,23,29,31,37,113):
    ret=T; m=1
    while tuple(tuple(x%p for x in row) for row in ret)!=I:
        ret=prod(ret,T,p);m+=1
    M=mpow(T,m)
    d=val(det(adjdiff(M)),p)+val(M[0][1],p)
    L=mlog_power(M,p,d+2)
    require(val(det((L[0],prod(L,L)[0])),p)==d,'wrong integer jet determinant')
    print('p, period mod p, jet determinant valuation:',p,m,d)
p=113
M=mpow(T,56)
require(mpow(T,56,p)==((58,25),(87,57)), 'wrong ramified shear')
require(all(mpow(T,j,p)!=I for j in (56,56*p-1)), 'wrong ramified order controls')
B=adjdiff(M)
for r in (1,2,3,4,5):
    modulus=p**r
    for k in (0,1,r,p-1):
        j=k*p**(r-1)
        require(mpow(M,j,modulus)==tuple(tuple((I[i][l]+j*B[i][l])%modulus
            for l in (0,1)) for i in (0,1)), 'wrong digit jump')
for v in ((9,89),(5,53)):
    for r in (1,2):
        modulus=p**r
        x,y=v
        seen=[bytearray(modulus) for _ in range(56)]
        for j in range(56*modulus):
            phase=seen[j%56]
            require(not phase[x], 'ramified phase collision')
            phase[x]=1
            x,y=y,(9*y+8*x)%modulus
        require((x,y)==v,'ramified period endpoint changed')
        require(all(sum(a)==modulus for a in seen),'ramified phase missing digit')
print('all-prime jet valuations and ramified replay passed')
for p in (5,7,11,13,17,19,23,29,31,37):
    m=1
    while mpow(T,m,p)!=I: m+=1
    E=adjdiff(mpow(T,m))
    require(all(x%p==0 for row in E for x in row),'invalid return')
    D=tuple(tuple(x//p%p for x in row) for row in E)
    require(det((D[0],prod(D,D,p)[0]))%p!=0,'invalid critical line')
    counts=[]
    for x,y in ((9,89),(5,53)):
        c=0
        for i in range(m):
            c+=(D[0][0]*x+D[0][1]*y)%p==0
            x,y=y%p,(9*y+8*x)%p
        counts.append(c)
    print('p, critical phases AB/BA, total:',p,*counts,m)
