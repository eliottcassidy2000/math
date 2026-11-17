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

from fractions import Fraction as Q
from itertools import product
M={'A':(9,18,15,32),'B':(5,20,7,32),'C':(9,22,11,32)}
U={'A':17,'B':19,'C':21}
def ck(x,y):
    if x!=y:raise ArithmeticError((x,y))
def height(w):
    h,R=1,2
    for l in w:
        a,b,c,d=M[l];h,R=a*h+b*R,c*h+d*R
    return h,R
def word(h,R):
    w=[]
    while (h,R)!=(1,2):
        if h<1 or R<h:return None
        cand=[l for l,(a,b,c,d) in M.items() if Q(c+d,a+b)<=Q(R,h)<Q(d,b)]
        if len(cand)!=1:return None
        l=cand[0];a,b,c,d=M[l];D=a*d-b*c
        hh,s=divmod(d*h-b*R,D);rr,t=divmod(a*R-c*h,D)
        if s or t or hh<1 or hh>=h:return None
        h,R=hh,rr;w.append(l)
    return ''.join(reversed(w))
def number(w):
    h,R=height(w);return h+(1<<(6*len(w)+1))*R
def decode(d,N):
    if d<0 or N<1:return None
    R,h=divmod(N,1<<(6*d+1))
    w=word(h,R)
    return w if w is not None and len(w)==d else None
def pair(w):
    m,n=2,1
    for l in w:
        m,n={'A':(2*m-n,m),'B':(2*m+n,m),'C':(m+2*n,n)}[l]
    return m,n
def primitive_address(m,n):
    w=[]
    while (m,n)!=(2,1):
        if m<1 or n<1 or m<=n:return None
        if n<m<2*n:w.append('A');m,n=n,2*n-m
        elif 2*n<m<3*n:w.append('B');m,n=n,m-2*n
        elif m>3*n:w.append('C');m,n=m-2*n,n
        else:return None
    return ''.join(reversed(w))
def strong(A):
    for i in range(len(A)):
        r=1<<i
        for _ in A:
            for j,x in enumerate(A):
                if r>>j&1:r|=x
        if r!=(1<<len(A))-1:return False
    return True

def main():
    for l,u in U.items():
        k=5;ck(M[l],(F(u),u+1,32-rho_num(k,u),32))
        a,r=grow((2,0),1,k,u);H,R=height(l)
        d,e,p=state(a,r);ck((d[:2],e,strong(a)),((H,R),1,True))
        b,s=grow(a,r,7,127)
        ck(sum(table(b)[-1]),number(l));ck(strong(b),True)
    for d in range(8):
        H=set()
        for s in product(M,repeat=d):
            w=''.join(s)
            ck(primitive_address(*pair(w)),w)
            ck(word(*height(w)),w)
            N=number(w);ck(decode(d,N),w);H.add(N)
        ck(len(H),3**d)
    q=number('B')
    from math import isqrt
    ck(q,9133);ck(all(q%d for d in range(2,isqrt(q)+1)),True)
    for a in range(3,13):
        mod=1<<a; seen=set()
        for v in ((1,2),height('AC')):
            h,R=v
            for _ in range(a+mod//4):
                h,R=(9*h+18*R)%mod,(15*h+32*R)%mod
                if _>=a:seen.add(h)
        ck(seen,set(range(1,mod,2)))
    print('3280 exact ternary inverses; three literal strong order8 and sealed order16 controls; all dyadic phases through 4096; one prime count.')
if __name__=='__main__':main()
