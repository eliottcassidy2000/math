"""Target-first all-width compiler for the ordered 5/11 two-wing bank."""
from math import isqrt
from dataclasses import replace
from replay_f1 import core5,core11,inverse,make

def binary_h(x):
    h=1
    while x:
        x,p=divmod(x,2)
        if not p:h+=x+x%2
    return h

def joined(C,D):
    # C,D are nonempty ordinal subblocks, with one common apex.
    return type(C)((),0,C.h*D.h,C.f*D.f,
               D.f*C.P+C.h*D.P,C.f*D.Q+D.h*C.Q,
               C.h*D.S,D.h*C.T,
               D.h*C.K+C.h*(D.K-D.h)+C.P*D.Q+C.S*D.T)

def value(C,x,y):
    return C.K+C.h*(binary_h(x)+binary_h(y)-2)+C.Q*x+C.P*y+C.f*x*y+C.T*(x%2)+C.S*(y%2)

def wings(C,N):
    """All x,y>0 yielding exact baseline N, including asymmetric widths."""
    if N<C.f:return
    for converse in (False,True):
        D=(replace(C,P=C.Q,Q=C.P,S=C.T,T=C.S) if converse else C)
        for u in range(1,1+isqrt(N//C.f)):
            if converse and u*u>=N//C.f+1:break
            H0=D.K+D.h*(binary_h(u)-1)+D.Q*u+D.T*(u%2)
            A=D.P+D.f*u
            upper=(N-H0)//A
            if upper<u:continue
            bits=0;scale=1
            while A*scale*scale<=2*D.h*upper or A*scale<=2*D.h:
                bits+=1;scale*=2
            for tail in range(min(scale,upper+1)):
                p,q,target=A,D.S,N-H0+D.h
                for i in range(bits):
                    if tail>>i&1:
                        target-=p+q
                        p,q=2*p,0
                    else:
                        p,q=2*p+D.h,D.h
                lo=max(0,(u+int(converse)-tail+scale-1)//scale)
                hi=(upper-tail)//scale
                while lo<=hi:
                    z=(lo+hi)//2
                    test=p*z+D.h*binary_h(z)+q*(z%2)
                    if test==target:
                        y=scale*z+tail
                        ans=(y,u) if converse else (u,y)
                        assert value(C,*ans)==N
                        yield ans
                        break
                    if test<target:lo=z+1
                    else:hi=z-1

def decode(N,saturated=False):
    """Yield (a,b,x,y,i,j). Empty output is an exact bank rejection."""
    A=core5[0];a=1
    while A.f*core11[0].f<=N:
        C=joined(A,core11[0]);b=1
        while C.f<=N:
            if not saturated or (a>=9 and b>=10):
                L=(5**(a-1))*(11**(b-1))
                for d in range(11*a+10*b+1):
                    if saturated and not 90<=d<=11*a+10*b-90:continue
                    bits=inverse(d,a,b,11,10)
                    if bits is None:continue
                    for x,y in wings(C,N-2*L*d):
                        i,j=bits
                        yield a,b,x,y,i,j
            C=joined(C,core11[0]);b+=1
        A=joined(A,core5[0]);a+=1

def literal_spec(spec):
    a,b,x,y,i,j=spec
    return make([core5[1]]*i+[core5[0]]*(a-i)+[core11[1]]*j+[core11[0]]*(b-j),x,y)

if __name__=='__main__':
    from sys import argv
    if len(argv)>1:
        N=int(argv[1]);w=next(decode(N,'--saturated' in argv),None)
        print("reject" if w is None else (w,literal_spec(w)))
    else:
        C=joined(core5[0],core11[0])
        for x,y in [(1,1),(2,3),(8,2),(1,1000),(331,7),(20117,2),(2,30009)]:
            N=value(C,x,y)
            assert (x,y) in set(wings(C,N))
            d=11
            assert (1,1,x,y,1,0) in set(decode(N+2*d))
        N=value(C,1,1)
        assert list(decode(23))==[]
        v=next(decode(N+22))
        assert v==(1,1,1,1,1,0)
        print("PASS: exact absolute 5/11 bank including both asymmetric wing orientations")
