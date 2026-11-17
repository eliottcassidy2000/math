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

from itertools import product
from math import gcd,isqrt,prod

W={'A':((6,20),),'B':((4,2),(8,84)),'C':((6,8),(10,340))}
maps={'A':(82,21),'B':(5976,1579),'C':(96264,24407)}
A=(174,124,120,240,65,80,1,118)
B=(174,124,120,240,65,16,33,118)

def ck(a,b):
    if a!=b:raise ArithmeticError((a,b))

def h_of(w):
    h=1
    for l in w:
        m,b=maps[l];h=m*h-b
    return h

def v2(x):
    return (x&-x).bit_length()-1

def decode(h):
    if h<1:return None
    track=[]
    while h!=1:
        r=h%3
        if r not in (1,2):return None
        s=v2(27*h-7)-(2 if r==1 else 3)
        if s<0:return None
        for _ in range(s):
            h,j=divmod(h+21,82)
            if j or h<1:return None
            track.append('A')
        if h==1:break
        l='C' if r==1 else 'B'
        m,b=maps[l]
        h,j=divmod(h+b,m)
        if j or h<1:return None
        track.append(l)
    return ''.join(reversed(track))

def endpoint(w):
    return 58284*h_of(w)-15505

def run_reader(I,r,C):
    n,k=len(I),len(C); j=n+k; M=(1<<k)-1
    # Incoming I is reversed. All I -> C -> q; I -> q except q -> r.
    T=tuple(sum(((row>>i)&1)<<z for z,row in enumerate(I)) for i in range(n))
    return tuple(r0|(M<<n)|((i!=r)<<j) for i,r0 in enumerate(T)) + tuple(
       (row<<n)|(1<<j) for row in C)+(1<<r,)

def strong(A):
    for i in range(len(A)):
        r=1<<i
        for _ in A:
            for j,x in enumerate(A):
                if r>>j&1:r|=x
        if r!=(1<<len(A))-1:return False
    return True

def lucas(n,g,fs):
    ck(prod(p**e for p,e in fs),n-1);ck(pow(g,n-1,n),1)
    for p,_ in fs:
        ck(p>=2 and all(p%d for d in range(2,isqrt(p)+1)),True)
        ck(gcd(pow(g,(n-1)//p,n)-1,n),1)

def main():
    ck(((A[5]^B[5]),(A[6]^B[6])),(1<<6,1<<5))
    ck(tuple(i for i in range(8) if A[i]!=B[i]),(5,6))
    ck([strong(T) for T in (A,B)],[True,True])
    ck([(state(T,0)[0][0],state(T,0)[2]) for T in (A,B)],[(187,3736),(197,3728)])
    for l,w in W.items():
        for h in (1,5):
            v=(h,3*h-1,0)
            for k,u in w:v=update(*v,k,u)
            m,b=maps[l]
            ck(v,(m*h-b,3*(m*h-b)-1,0))
    for w in ('','A'):
        T,r=(2,0),1
        for k,u in tuple(e for l in w for e in W[l])+((2,2),):
            T,r=grow(T,r,k,u)
        h=h_of(w)
        v,e,p=state(T,r)
        ck((v[0],v[1],e),(12*h-3,15*h-4,1));ck(strong(T),True)
        if not w:
            ck([sum(table(run_reader(T,r,Z))[-1]) for Z in (A,B)],[42779,42781])
    good=set()
    for d in range(9):
        for letters in product('ABC',repeat=d):
            w=''.join(letters);h=h_of(w)
            ck(decode(h),w);good.add(h)
    ck(len(good),sum(3**d for d in range(9)))
    for h in range(1,50000):
        w=decode(h)
        if w is not None:ck(h_of(w),h)
    n=endpoint('CAB')
    ck(n,2052302173293611)
    lucas(n,2,((2,1),(5,1),(7,1),(18169,1),(4679,1),(344873,1)))
    lucas(n+2,2,((2,2),(3,1),(347,1),(1889,1),(260914747,1)))
    for k in range(2,62,2):
        u=(2**k-4)//3
        ck(F(u),(5*2**k+6*k-14)//18)
        ck((3*F(u)-2**k+rho_num(k,u),3*(u+1)-2**k),(k-3,-1))
    print('Verified two literal rooted states, both literal initial readers, 9841 decoded words, and the twin primes by Lucas certificates.')

if __name__=='__main__':main()
