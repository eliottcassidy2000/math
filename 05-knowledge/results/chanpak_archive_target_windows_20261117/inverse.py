from math import gcd,isqrt
from functools import lru_cache
@lru_cache(None)
def F(v):
    if not v:return 1
    t=v//2
    return F(t)+(0 if v%2 else t+t%2)
def H(u,w,side):
    return 137+2*side+432*u*w+258*u+129*w+45*(F(u)+F(w)-2)
def ceil(a,b):return -((-a)//b)
def inverse(N,side):
    c=137+2*side
    if N<c or N%3!=c%3 or N%2==0:return []
    out=set()
    for orient in (0,1):
        # Try fixing the smaller actual flank. Every solution has such an orientation.
        for v in range(0,1+isqrt((N-c)//432),2):
            D=432*v+(129 if orient==0 else 258)
            C=c+v*(258 if orient==0 else 129)+45*(F(v)-2)
            if C+45>N:continue
            L=1;ell=0
            while D*D*L*L<64*N:L*=2;ell+=1
            for r in range(0,L,2 if L>1 else 1):
                s=f'{r:0{ell}b}' if ell else ''
                X=int(''.join('1' if z=='0' else '0' for z in reversed(s)),2) if s else 0
                B=D*L+45*X
                Y=N-C-D*r-45*(F(r)-1)
                g=gcd(B,45)
                if Y%g:continue
                m=45//g
                low=max(0,ceil(Y-90,B+45));high=(Y-45)//B
                q0=(Y//g)*pow(B//g,-1,m)%m
                q=low+(q0-low)%m
                if q>high:continue
                x=L*q+r
                if x%2 or x<v:continue
                u,w=(v,x) if not orient else (x,v)
                if H(u,w,side)==N:out.add((u,w))
    return sorted(out)
if __name__=='__main__':
    for side in (0,1):
        for u,w in ((0,0),(2,0),(0,8),(16,12),(102,0),(0,50000),(246,1296),(2,8320000)):
            if (u,w) not in inverse(H(u,w,side),side):raise ValueError((u,w))
    small={}
    for u in range(0,42,2):
        for w in range(0,84,2):
            N=H(u,w,0)
            if N<=10000:small.setdefault(N,[]).append((u,w))
    for N in range(137,10001,6):
        if inverse(N,0)!=sorted(small.get(N,[])):raise ValueError(N)
    print("exact target inverses agree with direct scan and selected larger witnesses")

    u=3533710
    if inverse(1000000007,0)!=[(u,0)] or inverse(1000000009,1)!=[(u,0)]:
        raise ValueError("prescribed twin targets")
    print("exact 1000000007 / 1000000009 on a single contact switch; u=3533710, w=0")
