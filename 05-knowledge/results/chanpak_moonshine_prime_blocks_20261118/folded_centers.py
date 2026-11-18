"""Figurate absolute centers and two overlapping protected windows."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'joint_growth'),str(HERE.parent/'half_reflection'),
             str(HERE.parent/'ans_inverse')]
from joint_decode import join,checksum,LOW,HIGH
from reflection_replay import arm,rooted_last,glue,transitive_arm,rev
from fixtures import state,forward
from math import isqrt
from itertools import product

CAP=((0,5,17,7,43,79,159,63),3)
PAIR65=(((8,33,35,70,15,153,55,95),1),((8,33,19,38,75,149,47,95),1))
CS=(29,17,448,91,-17,91,-15)
AS=(65,17,224,72,-6,135,-11)
BS=(201,65,976,206,-28,570,-52)
G=17*65
def converse(S):
    h,g,b,u,v,A,V=S
    return h,g,b,A,V,u,v
def core(m):
    if not m:return join(CS,CS)
    L=join(AS,BS)
    W=L
    for i in range(1,m):W=join(W,L)
    return join(join(join(CS,W),converse(W)),CS)
def step(m):return 2*17**2*G**(2*m-1)
def counts(m,u):
    assert u>0 and u%2==1
    h,g,b,U,V,A,Z=core(m)
    assert U==A and Z-V==34*G**(2*m)
    return h+2*b*u*u+3*U*u+g*(2*checksum(u)+u+1)-Z,Z-V

def digits(m,r):
    # j counts 65 upgrades; k counts 17 upgrades.
    p=2*m
    j=(23*r)%130  # 17*23 = 391 == 1 mod 130
    k=(r-17*j)//130
    a=max(0,(k-p+16)//17)
    j+=130*a;k-=17*a
    return (j,k) if 0<=j<=p and 0<=k<=p else None

def decode(m,N):
    assert m>=1
    h,g,b,*_=core(m)
    if N<h:return ()
    u=isqrt((N-h)//(2*b))
    if u%2==0:return ()
    C,D=counts(m,u)
    sol=[]
    for e in (0,1):
        r,rem=divmod(N-C-e*D,step(m))
        if rem==0:
            z=digits(m,r)
            if z is not None:sol.append((u,e,*z))
    return tuple(sol)

def two_tape(m,x,y):
    assert x%2==1 and y>0 and y%2==0
    h,g,b,U,V,A,Z=core(m)
    assert U==A
    return h+b*x*y+U*(x+y)+g*(checksum(x)+checksum(y))-Z,Z-V

def rect_decode(m,N):
    assert m>=1
    h,g,b,U,V,_,Z=core(m)
    if N<h:return
    for x in range(1,1+isqrt((N-h)//b),2):
        base=h+U*x+g*checksum(x)-Z
        y=(N-base)//(b*x+U)
        if y%2 or not x<=y<=2*x:continue
        C,D=two_tape(m,x,y)
        for e in (0,1):
            q,r=divmod(N-C-e*D,step(m))
            z=digits(m,q)
            if r==0 and z is not None:yield (x,y,e,*z)

def word(m,u,e,j,k,y=None):
    A=tuple(rooted_last(*s) for s in (LOW,HIGH))
    B=tuple(rooted_last(*s) for s in PAIR65)
    K=rooted_last(*CAP)
    if e:K=rev(K)
    pieces=[]
    for i in range(2*m):
        a,b=A[i<k],B[i<j]
        pieces.extend((a,b) if i<m else (rev(b),rev(a)))
    Q=K
    for T in pieces+[K]:Q=glue(Q,T)
    return glue(glue(rev(transitive_arm(u)),Q),transitive_arm(2*u if y is None else y))

def walls():
    for n in range(2,10):
        M=CS
        for _ in range(1,n):M=join(M,CS)
        Q=set()
        for bits in product((0,1),repeat=n-2):
            w=(0,*bits,0);Z=CS
            for b in w[1:]:Z=join(Z,converse(CS) if b else CS)
            k=sum(x!=y for x,y in zip(w,w[1:]))
            assert Z==(M[0]+2*17**(n-2)*k,*M[1:])
            Q.add(k)
        assert Q==set(range(0,n,2))
    print('orientation-wall palette through nine arms verified')
def main():
    walls()
    assert arm(*CAP)==CS and arm(*PAIR65[0])==BS
    assert arm(*PAIR65[1])==(203,*BS[1:])
    assert arm(*LOW)==AS and arm(*HIGH)==(69,*AS[1:])
    C,D=counts(0,1)
    A,B=(word(0,1,e,0,0) for e in (0,1))
    r=len(A)-1
    a,b=arm(A,r),arm(B,r)
    assert a[0]==C and b==(a[0]+D,*a[1:])
    A,B=(forward(T,r,(0,)) for T in (A,B))
    a,b=state(A,r),state(B,r)
    assert a[0]==C and b==(a[0]+D,*a[1:])
    print('literal two centers, same external signature',C,C+D,'shift',D)
    for m in (1,2,64,65,128):
        for u in (1,3,65,1001,(1<<49)+1):
            C,D=counts(m,u)
            assert D==65*step(m)
            for e in (0,1):
                for j,k in ((0,0),(2*m,2*m),(m,m)):
                    N=C+e*D+step(m)*(17*j+130*k)
                    assert any(p==u and s==e and (17*a+130*b)==17*j+130*k
                               for p,s,a,b in decode(m,N))
            if m>=65:
                R=17*129+65
                N=C+R*step(m)
                assert len(decode(m,N))==2
        x,y=9,14
        C,D=two_tape(m,x,y)
        R=17*m+130*m
        N=C+D+step(m)*R
        assert any((a,b,e,17*j+130*k)==(x,y,1,R)
                   for a,b,e,j,k in rect_decode(m,N))
        if m>=65:
            print('overlap',m,17*129+65,294*m-17*129,'two inverse choices at every odd address')
    m=65;u=1;C,D=counts(m,u)
    for s in decode(m,C+(17*129+65)*step(m)):
        assert s[0]==u
        W=word(m,*s)
        assert all(not ((r>>i)&1) and r.bit_length()<=len(W) for i,r in enumerate(W))
    print('compiled both physical large graphs; large counts certified by transfer')
if __name__=='__main__':main()
