from math import isqrt
from random import Random

def ck(x,y):
    if x!=y:raise AssertionError((x,y))
def rev(r,j):
    return int(format(r,f'0{j}b')[::-1],2) if j else 0
def U(u):
    b=[i for i in range(u.bit_length()) if u>>i&1]
    return sum((1<<(i-k-1))+(i==k+1) for i in b for k in b if i>k)
def phi(u,h,p):
    return u*p+(1-U(u))*h
def pair(h,p,j,rho):
    q=1<<j
    M,d=divmod(q*p-2,h)
    if d or M%2 or not (0<=rho<q and rho%2==0):return None
    A=M-rho
    if A<=0:return None
    r=rev(rho,j)
    u=2*q*A+r
    N=2*h*A*A+4*A+r*p+h*(1-U(A)-U(r))
    return N,u
def decode(h,p,N):
    if N<1:return []
    answers=[]
    for j in range(N.bit_length()+1):
        M,d=divmod((1<<j)*p-2,h)
        if d or M<1 or M&1:continue
        R=isqrt(N//(2*h))
        for A in range(max(1,R-9),R+10):
            rho=M-A
            a=pair(h,p,j,rho)
            if a and a[0]==N:answers.append((j,rho))
    return answers
def TT(n):return tuple((1<<n)-(1<<(j+1)) for j in range(n))
def literal(u,B):
    n=u.bit_length();m=len(B);g=n+m
    a=[r|(((1<<m)-1)<<n) for r in TT(n)]+[(r<<n)|(1<<g) for r in B]+[0]
    for i,b in enumerate(bin(u)[2:]):
        if b=='1':a[g]|=1<<i
        else:a[i]|=1<<g
    return tuple(a)
def table(a):
    d=[[0]*len(a) for _ in range(1<<len(a))];h=[1]+[0]*(len(d)-1)
    for i in range(len(a)):d[1<<i][i]=1
    for m in range(1,len(d)):
        h[m]=sum(d[m])
        for i,c in enumerate(d[m]):
            if not c:continue
            z=a[i]&~m
            while z:
                b=z&-z;z-=b;d[m|b][b.bit_length()-1]+=c
    return h
def profile(B):
    h=table(B);f=len(h)-1
    return h[-1],sum(h[s]*h[f^s] for s in range(f+1))
def strong(a):
    for i in range(len(a)):
        r=1<<i
        for _ in a:
            for j,x in enumerate(a):
                if r>>j&1:r|=x
        if r!=(1<<len(a))-1:return False
    return True
def main():
    B=(2,4,9,3);ck(profile(B),(5,32))
    for T in (TT(1),TT(2),TT(3),B):
        h,p=profile(T)
        for j in range(3):
            for rho in range(0,1<<j,2):
                a=pair(h,p,j,rho)
                if a is None:continue
                N,u=a
                ck(phi(u,h,p),N);ck(phi(u+(1<<j),h,p),N+2)
                old,new=literal(u,T),literal(u+(1<<j),T)
                ck(strong(old) and strong(new),True)
                ck(sum(x!=y for x,y in zip(old,new)),2)
                ck((table(old)[-1],table(new)[-1]),(N,N+2))
    rng=Random(2182);cnt=0
    for a,H,P,L in ((0,1,4,1),(1,5,32,4),(2,25,1024,20),(3,125,32768,100)):
        for j in range(60):
            q=1<<j
            if (q*P-2)%H:continue
            for rho in set([0]+[2*rng.randrange(max(1,q//2)) for _ in range(12)]):
                N,u=pair(H,P,j,rho)
                ck(phi(u,H,P),N);ck(phi(u+q,H,P),N+2)
                ck((j,rho) in decode(H,P,N),True);cnt+=1
    # Test a full source-height residue output and monotonicity at one nontrivial h.
    h,p,j=5,32,12
    V=[pair(h,p,j,r)[0] for r in range(0,1<<j,2)]
    ck(all(x>y for x,y in zip(V,V[1:])),True)
    ck({x%h for x in V},set(range(h)))
    ck(pair(1,4,0,0)[0],17)
    ck(pair(5,32,0,0)[0],379)
    for p in (17,19,379):
        ck(all(p%d for d in range(2,isqrt(p)+1)),True)
    print(cnt,'formula/inverse checks, literal final graph counts, source 5-adic phase and 2048 strict/modular controls pass.')
if __name__=='__main__':main()
