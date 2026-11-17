from functools import lru_cache
from math import isqrt
@lru_cache(None)
def F(x):
    if not x:return 1
    q=x//2
    return F(q) if x%2 else F(q)+q+(q%2)

@lru_cache(None)
def rc(r,t):
    j=0
    for i in range(t):
        j=(j<<1)+(1-((r>>i)&1))
    return j

def ceildiv(x,y):return -((-x)//y)

def select(N,h,P,L,R,Ls,Rs,k):
    if N<k:return None
    rem=N-k
    small=isqrt(rem//P)
    for left,right,ls,rs,rev in [(R,L,Rs,Ls,False),(L,R,Ls,Rs,True)]:
        for u in range(small+1):
            T=rem-h*(F(u)-2)-left*u+ls*(u%2)
            A=right+P*u
            if T<h-max(0,rs):continue
            t=0
            while A*(1<<t)<isqrt(N): t+=1
            y=1<<t
            for r in range(y):
                C=A*y+h*rc(r,t)
                D=T-A*r-h*(F(r)-1)+(rs*(r%2) if t else 0)
                g=t>0 and (r>>(t-1))==0
                tail=rs if not t else 0
                for q in range(max(0,ceildiv(D-2*h-max(0,-tail),C+h)),(D-h+max(0,tail))//C+1):
                    w=y*q+r
                    if w<u:continue
                    if C*q+h*F(q)+(h*g-tail)*(q%2)==D:
                        return (w,u) if rev else (u,w)
    return None

if __name__=='__main__':
    for n in (257,631,8191,65537,458759,1000000007,7000000049):
        for k in (137,139): print(n,k,select(n,45,432,129,258,-15,-30,k))
