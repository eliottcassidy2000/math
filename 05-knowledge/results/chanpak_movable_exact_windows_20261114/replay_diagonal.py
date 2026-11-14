from collections import defaultdict

def ending(A):
    n=len(A); F=(1<<n)-1
    out=[[0]*n for _ in range(F+1)]
    for i in range(n):out[1<<i][i]=1
    for S in range(1,F+1):
        for i,x in enumerate(out[S]):
            if x:
                B=A[i]&~S
                while B:
                    b=B&-B;B^=b
                    out[S|b][b.bit_length()-1]+=x
    return out
def count(A):
    return sum(ending(A)[-1]) if A else 1
def rev(x,n):
    return int(f'{x:0{n}b}'[::-1],2)
def g(x):
    h=1
    while x:
        x,b=divmod(x,2)
        if not b:h+=x+(x&1)
    return h

def profile(A,U):
    n=len(A);F=(1<<n)-1
    ee=ending(A)
    ss=ending([r^F^(1<<i) for i,r in enumerate(A)])
    h=[1]+[sum(row) for row in ee[1:]]
    f=sum(h[s]*h[F^s] for s in range(F+1))
    P=sum(ss[s][i]*h[F^s] for s in range(F+1) for i in range(n) if U>>i&1)
    Q=sum(ee[s][i]*h[F^s] for s in range(F+1) for i in range(n) if not U>>i&1)
    S=sum(ee[F][i] for i in range(n) if U>>i&1)
    T=sum(ss[F][i] for i in range(n) if not U>>i&1)
    K=count(tuple(r|((((U>>i)&1)^1)<<n) for i,r in enumerate(A))+(U,))
    E=sum(ee[s][i]*ss[F^s][j] for s in range(F+1) for i in range(n) for j in range(n) if U>>i&1 and U>>j&1)
    assert K-P+S==h[F]-E
    return h[F],f,P,Q,S,T,K

def physical(A,U,x,y):
    left=f'{x:b}'
    right=''.join('1' if b=='0' else '0' for b in f'{y:b}'[::-1])
    l,n,r=len(left),len(A),len(right); F=(1<<(l+n+r))-1
    C=[]
    for i,b in enumerate(left):
        C.append((F^((2<<i)-1))|(((b=='0')*1)<<(l+n+r)))
    for i,R in enumerate(A):
        C.append((R<<l)|(((1<<r)-1)<<(l+n))|((((U>>i)&1)^1)<<(l+n+r)))
    for j,b in enumerate(right,l+n):
        C.append((F^((2<<j)-1))|(((b=='0')*1)<<(l+n+r)))
    return C+[sum((b=='1')<<i for i,b in enumerate(left))|(U<<l)|sum((b=='1')<<(l+n+i) for i,b in enumerate(right))]

A=(2,4,9,3,15);B=(0,5,9,19,7)
banks=[(35,[(A,11),(A,7),(B,3),(B,5)]),(83,[(B,26),(B,28),(A,24),(A,20)])]
for W,choices in banks:
    for j,(C,U) in enumerate(choices):
        h,f,P,Q,S,T,K=profile(C,U)
        assert (h,f,P+Q,K-P+S)==(5,64,W,-7+2*j)
        for M in [2,4,8,16,32,64,256]:
            for u in [1,2,3,19]:
                A0=M*u
                expect=f*A0*A0+A0*(W-f)+h*((M-1)*u+g(u)+g(u-1)+(u&1)-2)-7+2*j
                assert expect==f*A0*(A0-1)+Q*A0+P*(A0-1)+h*(g(A0)+g(A0-1)-2)+K+S
                if (M,u) in [(2,1),(2,2),(4,1),(8,1),(8,2),(16,1)]:
                    assert count(physical(C,U,A0,A0-1))==expect
    print('Verified full literal four-odd windows at W',W)
for u in range(1,513):
    A0=8*u; top=A0*A0-u+g(u)+g(u-1)+(u&1)-1
    for a in range(3):
        x,y=A0+a,A0-1-a
        assert x*y+g(x)+g(y)-1==top-2*a
        if u<=2:assert count(physical((),0,x,y))==top-2*a
print('Verified near-square three-odd windows and direct literal controls')

tilt=[((0,133,9,147,135,159,191,5),46),
      ((0,5,9,131,15,159,191,23),172),
      ((0,1,3,23,39,143,63,95),146),
      ((0,5,9,131,15,159,191,23),60),
      ((0,5,9,131,15,159,191,23),58)]
for j,(C,U) in enumerate(tilt):
    h,f,P,Q,S,T,K=profile(C,U)
    assert (h,f,2*P+Q,K-P+S)==(9,672,981,-51-2*j)
    for M in [2,4,8,16,32,64,256]:
        for u in [1,2,3,19]:
            x,y=M*u,2*M*u-1
            expect=1344*x*x+309*x+9*((M-1)*u+g(u)+g(2*u-1)+(u&1)-2)-51-2*j
            assert expect==f*x*y+Q*x+P*y+h*(g(x)+g(y)-2)+K+S
            if (M,u) in [(2,1),(2,2),(4,1),(8,1)]:
                assert count(physical(C,U,x,y))==expect
print('Verified full literal five-odd tilted windows at 2P+Q=981')

tilt7=[((2,8,3,21,7,31,63),46),
       ((2,4,9,19,7,31,63),46),
       ((0,1,11,19,39,79,31),41),
       ((2,4,9,19,7,31,63),45),
       ((2,4,9,19,7,31,63),53)]
for j,(C,U) in enumerate(tilt7):
    h,f,P,Q,S,T,K=profile(C,U)
    assert (h,f,4*P+Q,K-P+S)==(9,336,645,-51-2*j)
    for M,u in [(2,1),(2,2),(4,1),(8,1),(32,19)]:
        x,y=M*u,4*M*u-1
        expect=1344*x*x+309*x+9*((M-1)*u+g(u)+g(u-1)+(u&1)-2)-51-2*j
        assert expect==f*x*y+Q*x+P*y+h*(g(x)+g(y)-2)+K+S
        if (M,u) in [(2,1),(2,2),(4,1),(8,1)]:
            assert count(physical(C,U,x,y))==expect
print('Verified the order-seven tilted bank directly, including literal tournaments')
