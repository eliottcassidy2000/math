"""Literal central modules, movable wings, and exact inverse. Python stdlib."""
from dataclasses import dataclass
from math import isqrt
from random import Random

def binary_h(x):
    a,h=0,1
    for b in bin(x)[2:]:
        d=int(b)
        if not d: h+=a+(a&1)
        a=2*a+d
    return h

def rev(x,m): return int(format(x, f'0{m}b')[::-1], 2) if m else 0

@dataclass(frozen=True)
class Core:
    masks: tuple
    out: int
    h: int
    f: int
    P: int
    Q: int
    S: int
    T: int
    K: int

k41=Core((18,20,9,3,12),3,15,120,42,63,6,9,41)
k43=Core((18,20,9,3,12),5,15,120,42,63,6,9,43)
k65=Core((56,57,51,4,40,8),17,17,224,72,135,6,11,65)
k69=Core((54,32,26,51,2,20),10,17,224,72,135,6,11,69)
k373=Core((100,13,104,1,47,74,26),29,95,1276,694,445,62,39,373)
k375=Core((100,13,104,1,47,74,26),57,95,1276,694,445,62,39,375)

def value(C,x,y):
    return (C.K+C.h*(binary_h(x)+binary_h(y)-2)+C.Q*x+C.P*y
            +C.f*x*y+C.T*(x&1)+C.S*(y&1))

def rooted(C,x,y):
    m=y.bit_length()
    return ((C.f*x+C.P+C.h)<<m)-C.h*rev(y,m)

def graph(C,x,y):
    left, right=x.bit_length(),y.bit_length()
    n=left+len(C.masks)+right
    out=[0]*(n+1)
    apex_out=(x and rev(x,left))|C.out<<left
    if right:
        apex_out|=(((1<<right)-1)^y)<<(left+len(C.masks))
    out[n]=apex_out
    for i in range(n):
        if not apex_out>>i&1:out[i]|=1<<n
        for j in range(i+1,n):
            if left<=i<j<left+len(C.masks) and not C.masks[i-left]>>(j-left)&1:
                out[j]|=1<<i
            else: out[i]|=1<<j
    return out

def count(rows):
    n=len(rows); m=1<<n; H=[0]*m; H[0]=1
    S=[[0]*n for _ in range(m)]
    for z in range(1,m):
        rem=z
        while rem:
            v=(rem&-rem).bit_length()-1
            rem&=rem-1
            t=z^(1<<v)
            if not t: S[z][v]=1
            else:
                d=t&rows[v]
                while d:
                    w=(d&-d).bit_length()-1
                    d&=d-1
                    S[z][v]+=S[t][w]
        H[z]=sum(S[z])
    return H,S

def check_profile(C):
    H,A=count(C.masks)
    n=len(C.masks); m=1<<n
    conv=tuple((m-1)^(1<<i)^v for i,v in enumerate(C.masks))
    _,B=count(conv)
    U=[i for i in range(n) if C.out>>i&1]
    V=[i for i in range(n) if not C.out>>i&1]
    P=sum(A[z][i]*H[(m-1)^z] for i in U for z in range(1,m))
    Q=sum(B[z][i]*H[(m-1)^z] for i in V for z in range(1,m))
    S=sum(B[-1][i] for i in U)
    T=sum(A[-1][i] for i in V)
    assert (H[-1],sum(H[z]*H[(m-1)^z] for z in range(m)),P,Q,S,T,count(graph(C,0,0))[0][-1])==(C.h,C.f,C.P,C.Q,C.S,C.T,C.K)

def inverse_pair(C,N,R):
    # H≥P*y, so all relevant normalized right lengths are bounded.
    if N<C.K or R<=0:return
    for m in range((N//C.P).bit_length()+1):
        D=1<<m
        v=(-R*pow(C.h,-1,D))%D
        y=rev(v,m)
        if y.bit_length()!=m:continue
        q,t=divmod(R+C.h*v,D)
        if t:continue
        x,r=divmod(q-C.P-C.h,C.f)
        if x>=0 and not r and value(C,x,y)==N:yield x,y

def inverse_compiler(C,N):
    """Report all normalized witnesses, including both possible wing imbalances."""
    if N<C.K:return
    seen=set()
    for flip in (False,True):
        P,Q,S,T=(C.P,C.Q,C.S,C.T) if not flip else (C.Q,C.P,C.T,C.S)
        for u in range(1+isqrt(N//C.f)):
            base=C.K+C.h*(binary_h(u)-1)+Q*u+T*(u&1)
            A=P+C.f*u
            Y=(N-base)//A
            if Y<u:continue
            r,D=0,1
            while A*D*D<=2*C.h*Y or A*D<=2*C.h:
                r,D=r+1,D*2
            for j in range(min(D,Y+1)):
                a,b,k=A,S,N-base+C.h
                for i in range(r):
                    if j>>i&1:
                        k-=a+b
                        a,b=2*a,0
                    else:
                        a,b=2*a+C.h,C.h
                lo=max(0,(u-j+D-1)//D)
                hi=(Y-j)//D
                while lo<=hi:
                    z=(lo+hi)//2
                    trial=a*z+C.h*binary_h(z)+b*(z&1)
                    if trial==k:
                        ans=(u,D*z+j) if not flip else (D*z+j,u)
                        assert value(C,*ans)==N
                        if ans not in seen:
                            yield ans
                            seen.add(ans)
                        break
                    if trial<k:lo=z+1
                    else:hi=z-1

def literal(C,x,y):
    rows=graph(C,x,y)
    A,S=count(rows); m=len(A)-1
    assert A[-1]==value(C,x,y)
    assert sum(row[-1]*A[m^i] for i,row in enumerate(S))==rooted(C,x,y)
    if y:
        # Count paths ending at the apex by deleting their last vertex.
        old=m>>1
        assert not rows[-1]>>(len(rows)-2)&1
        assert A[old]==C.h

def main():
    for C in (k41,k43,k65,k69,k373,k375):
        check_profile(C)
    for C in (k41,k43,k65,k69):
        for x in range(16):
            for y in range(16):
                literal(C,x,y)
    for C in (k41,k65,k373):
        sums={}
        for x in range(100):
            for y in range(100):
                N=value(C,x,y)
                sums.setdefault(N,set()).add((x,y))
        for N in range(1, min(C.P,C.Q)*100):
            assert set(inverse_compiler(C,N))==sums.get(N,set())
        rng=Random(6145)
        for _ in range(30):
            x,y=rng.choice(((1,18),(24,1),(9,9),(0,20),(18,0)))
            x,y=rng.randrange(1<<x),rng.randrange(1<<y)
            N,R=value(C,x,y),rooted(C,x,y)
            assert (x,y) in set(inverse_compiler(C,N))
            assert (x,y) in set(inverse_pair(C,N,R))
    for r in range(50):
        for C,D,x,y in ((k65,k69,1+17*r,36+448*r),(k373,k375,7+95*r,204+2552*r)):
            N=value(C,x,y)
            assert N==value(D,x,y+1)
            assert rooted(C,x,y)-rooted(D,x,y+1)==C.h*(1<<(y.bit_length()-1))
            if r==0:
                literal(C,x,y)
                literal(D,x,y+1)
                print('exact neutral return:',N,'root responses',rooted(C,x,y),rooted(D,x,y+1))
    for C in (k41,k43):
        for X in range(1,20):
            for Y in range(1,20):
                A,B=(16*X+7,Y),(X,16*Y+8)
                assert rooted(C,*A)==rooted(C,*B)
                assert value(C,*A)-value(C,*B)==6*((X&1)+(Y&1)-1)
        for d in range(8):
            seen=set()
            for i in range(d+1):
                x=(16**(d-i))*3+7*(16**(d-i)-1)//15
                y=(16**i)*10+8*(16**i-1)//15
                seen.add((value(C,x,y),rooted(C,x,y),-C.h))
            assert len(seen)==1
        for x,y in ((23,2),(1,40),(23,20),(1,328)):
            literal(C,x,y)
        assert value(C,23,2)==C.K+7122 and rooted(C,23,2)==11253
        assert (23,2) in set(inverse_pair(C,C.K+7122,11253))
        assert (1,40) in set(inverse_pair(C,C.K+7122,11253))
        assert sorted(r.bit_count() for r in graph(C,23,2))!=sorted(r.bit_count() for r in graph(C,1,40))
        assert value(C,23,20)==value(C,1,328)==C.K+57798
    assert 57841==7*8263
    assert all(n%k for n in (57839,8263) for k in range(2,isqrt(n)+1))
    print("all checks passed")

if __name__=="__main__": main()
