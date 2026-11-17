from functools import lru_cache
from math import isqrt,log2
@lru_cache(None)
def F(x):
    if not x:return 1
    u=x//2
    return F(u) if x%2 else F(u)+u+u%2
def read(m,u,w):
    q=1<<m;h=q+1;P=4*q+2*(3**m-1)
    if w%2==0:
        e=1;L=2*q-1;d=2*q-2;R=P-2*h;rs=-q
    else:
        e=q;L=P-h-(2*q-1);d=R0_comp(m);R=h;rs=-1
    ell=P*u+L
    dx=d+h*(F(u)-1)+u*R-rs*(u%2)
    first=91*h+5*(dx+h*(F(w)-1))+142*(ell+h*w)+256*ell*w+2*h*(w%2)
    return first+min(0,2*e-2*h*(w%2))
@lru_cache(None)
def R0_comp(m):
    # Even chosen-slot marked paths from nonzero starts to the only possible end 0.
    # For all nonempty subsets ending 0, the even count excluding singleton is
    # (weighted total + signed total)/2, with signed total=-the original full start0=-1.
    return (1<<m)
def select(N):
    ans=[]
    for m in range(1,N.bit_length()):
        P=4*(1<<m)+2*(3**m-1)
        for w in range(1,N//(256*((1<<(m+1))-1))+1):
            for u in range(0,N//(256*P*w)+1):
                H=read(m,u,w)
                if N in (H,H+2):ans.append((m,u,w,H))
    return ans
if __name__=='__main__':
    for N in (8191,65537,458759,57337,3599,53549,53551):
        print(N,select(N))
