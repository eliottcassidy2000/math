"""Construct a literal reader word near an odd target, with two exact residues."""
from fractions import Fraction
from itertools import combinations

def K(d):
    return (1<<(d-1))+(d==1)

def rev(b,n):
    if n==0:return 0
    return int(f'{b:0{n}b}'[::-1],2)

def H(b,n):
    sites=[i for i in range(n) if b>>i&1]
    J=sum(K(j-i) for i,j in combinations(sites,2))
    return -11+82*((1<<(n+1))-rev(b,n))+9*(2*b-J)

def decode(N,k,q=1,a=0):
    if N<89 or not N%2 or q<1 or not q%2 or k<1:
        raise ValueError('require large odd target, positive dyadic depth and odd address modulus')
    n=(N//88).bit_length()-1
    D=1<<n
    L0=5*k+4*(q-1).bit_length()+9
    if n<L0:raise ValueError('not enough free bits at this target scale for the proven guarantee')
    e=(n-L0)//2
    L=n-2*e
    z=Fraction(N,D)
    x,y=((0,2-z/82) if z<=164 else ((z/2-82)/9,0))
    B=1<<e
    u,v=(min(B-1,int(B*w)) for w in (x,y))
    fixed=(u<<(n-e))|rev(v,e)
    exterior=[j for j in range(n) if fixed>>j&1]
    mask=(1<<(k+1))-1; mod=1<<k
    layers=[{(0,(H(fixed,n)-1)//2 % mod, fixed%q):None}]
    for i in range(e,e+L):
        weight=9*(1<<i)-41*(1<<(n-1-i))
        weight-=9*sum(K(abs(i-j))//2 for j in exterior if abs(i-j)<=k+1)
        residues=[]
        for m in range(mask+1):
            residues.append((weight-9*sum(K(d)//2 for d in range(1,k+2) if m>>(d-1)&1))%mod)
        nxt={}
        qi=pow(2,i,q)
        for state in layers[-1]:
            m,r,u0=state
            nm=(m<<1)&mask
            nxt.setdefault((nm,r,u0),(state,0))
            nxt.setdefault((nm|1,(r+residues[m])%mod,(u0+qi)%q),(state,1))
        layers.append(nxt)
    w=next(s for s in layers[-1] if s[1:]==((N-1)//2%mod,a%q))
    b=fixed
    for j in range(L,0,-1):
        w,digit=layers[j][w]
        b|=digit<<(e+j-1)
    hw=H(b,n)
    assert hw%(2*mod)==N%(2*mod) and b%q==a%q
    return n,b,hw

if __name__=='__main__':
    N=(109<<45)+369
    for a in range(3):
        n,b,h=decode(N,4,3,a)
        print('target',N,'count',h,'signed residual',h-N,'address mod 3',b%3,'bits',n)
