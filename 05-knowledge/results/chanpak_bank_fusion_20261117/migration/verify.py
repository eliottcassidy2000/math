from random import Random

def F(u):
    a,z=0,1
    for w in bin(u)[2:]:
        b=int(w)
        if not b:z+=a+a%2
        a=2*a+b
    return z

def rho_num(k,u):
    return int(format(u,'0'+str(k)+'b')[::-1],2)

def update(h,R,S,k,u):
    if not k:return h+R-S,h+R,h-S
    c=(1<<k)-rho_num(k,u)
    return F(u)*h+(u+1)*R-(1-u%2)*S,(1<<k)*R+c*h,(h if u.bit_length()==k else 0)

def grow(A,root,k,u):
    n=len(A);s=n+k
    out=list(A)+[0]*(k+1)
    out[root]|=1<<s
    out[s]=((1<<n)-1)^(1<<root)
    for j in range(k):
        i=n+j
        out[i]=((1<<n)-1)+sum(1<<v for v in range(i+1,s))
        if u>>(k-1-j)&1:out[s]|=1<<i
        else:out[i]|=1<<s
    return tuple(out),s

def table(A):
    n=len(A)
    D=[[0]*n for _ in range(1<<n)]
    for j in range(n):D[1<<j][j]=1
    for S in range(1,1<<n):
        for i,z in enumerate(D[S]):
            if not z:continue
            p=A[i]&~S
            while p:
                bit=p&-p;p-=bit
                D[S|bit][bit.bit_length()-1]+=z
    return D

def state(A,r):
    D=table(A);M=len(D)-1
    h=[1]+[sum(row) for row in D[1:]]
    AT=tuple(sum(((row>>i)&1)<<j for j,row in enumerate(A)) for i in range(len(A)))
    S=table(AT)[M][r]
    return (h[M],sum(D[V][r]*h[M^V] for V in range(1,M+1)),S),D[M][r],sum(h[V]*h[M^V] for V in range(M+1))

def run(word):
    A,r=(2,4,1),0
    v=(3,3,1)
    for k,u in word:A,r=grow(A,r,k,u);v=update(*v,k,u)
    d,e,p=state(A,r)
    if d!=v or e!=1:raise ArithmeticError(('literal',word,d,e,v))
    return v,p,A,r

for w in (((1,1),),((0,0),(0,0))):
    print(w,run(w),run(w+((0,0),)))
r=Random(961)
for _ in range(60):
    w=[]
    # Raw DP without relying on a precomputed catalogue.
    for _ in range(2):
        k=r.randrange(4)
        w.append((k,r.randrange(1<<k)))
    run(w)
for k in range(1,14):
    for u in range(1<<(k-1)):
        det=(1<<k)*F(u)-(u+1)*((1<<k)-rho_num(k,u))
        if (det==0)!=(u&(u+1)==0) or det<0:raise ArithmeticError('rank')
print('PASS migration, DP, unique ends, and determinant controls.')
def slot(A,q):
    d=table(A); h=[1]+[sum(v) for v in d[1:]]
    s=[1]+[0]*(len(h)-1)
    for _ in range(q):
        t=[]
        for mask in range(len(h)):
            b=mask;z=h[b]*s[mask^b]
            while b:
                b=(b-1)&mask
                z+=h[b]*s[mask^b]
            t.append(z)
        s=t
    return s[-1]
X=[]
for w in (((1,1),),((0,0),(0,0))):
    _,_,a,r=run(w);X.append((a,r))
    if any(slot(a,q)!=q**5+6*q**3+2*q for q in range(1,5)):
        raise ArithmeticError('slot polynomial')
for t in range(3):
    ds=[]
    for a,r in X:
        n=len(a)+t
        s=tuple(((1<<n)-1)-((1<<(j+1))-1) for j in range(t))+tuple(b<<t for b in a)
        ds.append(grow(s,r+t,0,0)[0])
    for q in range(1,5):
        if slot(ds[1],q)-slot(ds[0],q)!=q*q*((q+1)**(t+1)-(q-1)**(t+1)):
            raise ArithmeticError('tagged slot difference')
from math import isqrt
def seam(u,w,j):
    l=84*u+21+2*j
    d=20+2*j+9*(F(u)-1)+(66-2*j)*u+6*(u%2)
    return 819+5*(d+9*(F(w)-1))+142*(l+9*w)+256*l*w
for j in range(1,100):
    u,w=512*j+29,10*j
    if seam(u,w,1)-seam(u,w,0)!=4:raise ArithmeticError('fusion')
p=seam(5149,100,1)
if p!=11136264317 or any(a%d==0 for a in (p,p+2) for d in range(3,isqrt(a)+1,2)):
    raise ArithmeticError('twin prime fixture')
print('PASS neutral slot read, tagged polynomial differences, and actual Farey twin-prime fixture.')
