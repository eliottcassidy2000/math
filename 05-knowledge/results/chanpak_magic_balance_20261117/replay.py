from math import gcd,isqrt
from itertools import product
from slot_profiles import profile,append,external,ordinal,strong,end_profiles
A=(158,164,56,242,98,1,167,52)
B0=(254,240,162,54,4,144,60,88)
B1=(142,136,18,20,131,31,191,44)
R=(52,365,176,69,266,280,437,315,13)
S=[profile(R,m) for m in (113,43,39)]
if S != [(1521,26412,10422,11997,750,631,k) for k in (-2387,-2385,-2383)]:raise ValueError("rotor profile")
for j in range(3):
    t=external(R,(113,43,39)[j],1)
    if not strong(t) or end_profiles(t)[0][-1]!=18457+2*j:raise ValueError("real rotor")
def stack(tiles):
    state=(1,0,1,0);Z=0;P=1
    for r in tiles:
        state=append(state,r)
        Z=Z*r[0]+P*r[3];P*=r[1]
    J,W,K,_=state
    return J,P,W,Z,K
def magic(c,a,b):
    return [c+a,c-a-b,c+b,c-a+b,c,c+a-b,c-b,c+a+b,c-a]
L=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def square_check(vals,l):
    if len(set(vals))!=9 or any(sum(vals[i] for i in s)!=l for s in L):raise ValueError(vals)
for a,b in ((1,3),(9,2)):
    J,P,W,Z,K=stack([S[1]]*(a+b));T=2*(J//S[0][0]);c=2*W+K
    X=magic(c,a*T,b*T)
    for i,j in product(range(3),repeat=2):
        r=stack([S[i]]*a+[S[j]]*b)
        if r[:-1]!=(J,P,W,Z) or 2*r[2]+r[-1]!=c+T*(a*(i-1)+b*(j-1)):
            raise ValueError("rotor coefficients")
    square_check(X,3*c)
    if (a,b)==(9,2) and sorted((x-c)//T for x in X)!=sorted((x**3-15*x)//2 for x in range(-4,5)):
        raise ValueError("cubic")
    print("rotor",a,b,"center",c,"step",T)
T=[profile(A,144),profile(A,6),profile(B0,61),profile(B1,79)]
J,P,W,Z,K=stack(T)
if (J,P,W,Z,K)!=(1423779289,205520896000000,46594073974899,142977375089376,-32411071075897):
    raise ValueError("center fixture")
alpha=2*P-J-Z
C0=2*W+K-J
def F(v):
    if v==0:return 1
    t=v//2
    return F(t)+(0 if v%2 else t+t%2)
def center(v):
    return C0+alpha*v+J*F(v)
def ceil(a,b):return -((-a)//b)
def candidates(N):
    if N<center(0) or not N%2:return []
    L=1;ell=0
    while alpha*alpha*L*L<4*N:
        L*=2;ell+=1
    hits=[]
    for r in range(0,L,2 if L>1 else 1):
        X=int(''.join('1' if b=='0' else '0' for b in reversed(f'{r:0{ell}b}')),2) if ell else 0
        B=alpha*L+J*X
        Y=N-C0-alpha*r-J*(F(r)-1)
        low=max(0,ceil(Y-2*J,B+J));high=(Y-J)//B
        vq=(Y*pow(B,-1,J))%J
        q=low+(vq-low)%J
        if q<=high:
            v=L*q+r
            if not v%2 and center(v)==N:hits.append(v)
    return hits
if gcd(alpha,J)!=1:raise ValueError("not primitive")
for v in (0,2,14,10000,128,10**18):
    x=center(v)
    if v not in candidates(x) or candidates(x+2):
        raise ValueError((v,x))
print("exact center inverse",center(0),center(128))
