from math import gcd, isqrt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import canon, extend, strong, hp_all, rooted, ends, glue

def require(x):
    if not x:raise RuntimeError("all-port replay failed")

def state(T,v):
    a=list(rooted(T))[v]
    H,L=a[:2]
    F=hp_all(T)
    m=len(F)-1
    return H,F[m^(1<<v)],L,ends(T)[m][v]

def step(X,j,k):
    H,h,L,E=X; d=1<<j; F=H+d*L
    if k==0:return F,h,2*d*L-h,h
    if k==j+1:return F,H,F,H
    t=j-k
    return F,H+(d//2)*L if j>=2 else H+L-E,(2<<t)*H+(1<<k)*3**t*L-(E if k==1 else 0),(1<<(k-1))*L

def decorated_word(T,v,w):
    n=len(T)
    for j,p in w:
        if p=="O":k=0
        if p=="P":k=j
        if p=="T":k=j+1
        T=glue(T,v,j)
        if k:v=n+k-1
        n=len(T)
    return T,v

def word(X,w):
    for j,p in w:
        X=step(X,j,{"O":0,"P":j,"T":j+1}[p])
    return X

checks=0;bank={(0,)}
for n in range(1,6):
    if n>1:bank={canon(extend(A,c)) for A in bank for c in range(1<<(n-1))}
    for T in bank:
        if not strong(T):continue
        for v in range(n):
            X=state(T,v)
            for j in (1,2,4):
                A=glue(T,v,j)
                require(strong(A))
                for k in range(j+2):
                    Z=state(A,n+k-1 if k else v)
                    if Z!=step(X,j,k): raise RuntimeError((T,v,X,j,k,Z,step(X,j,k)))
                    checks+=1
print(f"all old/fresh formula cases passed: {checks}")

starters=[((0,),0,[(1,"T")],109,31),((2,4,9,3),0,[(1,"T")],623,177)]
for T,v,init,N,e in starters:
    w=init+[(1,"P"),(1,"O"),(1,"T")]
    A,v=decorated_word(T,v,w)
    require(state(A,v)==(N,e,N,e) and strong(A))
print("literal strong third-port outputs 109 and 623")

sources=[(2,8,3,21,39,15),(2,16,3,7,45,15),(2,4,41,19,7,27)]
for A,e in zip(sources,[15,17,19]):
    B,v=decorated_word(A,2,[(1,"T")])
    require(state(B,v)==(73,e,73,e))
    C,u=decorated_word(B,v,[(4,"O"),(1,"T")])
    H=365*(1<<4)+73-2*e
    d=73*(1+(1<<4))
    require(state(C,u)==(H,d,H,d) and gcd(H,d)==1 and strong(C))
    # Further (unbuilt) full graph states use proved letter recurrence.
    for j in (4,8,16,32,64):
        N=(1+5*(1<<j))*73-2*e
        ul=(1+(1<<j))*73
        require(word((73,e,73,e),[(j,"O"),(1,"T")])==(N,ul,N,ul) and gcd(N,ul)==1)
print("all three strong ready seeds and exact gap-four rails passed")

for H in range(3,70,2):
 for e in range(1,H+1,2):
    if gcd(H,e)!=1: continue
    X=(H,e,H,e)
    for j in (1,2,4,8):
        d=1<<j
        for r in range(4 if j>=2 else 1):
            b=word(X,[(j,"O")]+[(j,"P")]*r+[(j,"T")])
            expected=gcd(1+d,2*H+e)
            require(gcd(*b)==expected)
            if expected==1:
                require(gcd(1+d,2*b[0]+b[1])==1)
        for l in (1,2):
            b=word(X,[(1,"P"),(j,"O"),(l,"T")])
            require(gcd(*b)==gcd(4*d-3,8*H-e))
print("exact and renewable odd gcd gates certified on independent recurrence controls")

# Fixed O2 cheap binary; the cheap branch with I_{2,1} is independently
# tested against literal DPs above. These controls use both O2 suffixes.
def mm(A,B,p):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))%p for j in range(2)) for i in range(2))
def mv(A,x,p):
    return tuple(sum(a*b for a,b in zip(row,x))%p for row in A)
def power(A,n,p):
    R=((1,0),(0,1))
    while n:
        if n%2:R=mm(R,A,p)
        A=mm(A,A,p); n//=2
    return R
A=((2201,-212),(709,-68)); B=((2321,-272),(685,-80))
W=[[(2,"O"),(2,"P"),(4,"P"),(1,"T")],[(2,"O"),(4,"P"),(2,"P"),(1,"T")]]
for M,w in zip([A,B],W):
    for H,e in [(5,1),(73,15),(73,17)]:
        N,D=mv(M,(H,e),10**20)
        require(word((H,e,H,e),w)==(N,D,N,D))
    require(M[0][0]*M[1][1]-M[0][1]*M[1][0]==640)
seed=(5,1)
require(len({mv(power(A,n,7),seed,7) for n in range(48)})==48)
require(power(A,48,49)==((43,14),(35,1)))
C=((6,2),(5,0));D=((4,1),(4,2))
require(mm(B,C,7)==mm(D,B,7))
require(mv(C,seed,7)==(4,4) and mv(D,seed,7)==(0,1))
require(mv(power(A,336,49),seed,49)==seed)
require(all(mv(power(A,n,49),seed,49)!=seed for n in range(1,336)))
require(mv(B,seed,10**10)[0]==11333 and 11333//7==1619)
require(all(1619%i for i in range(2,isqrt(1619)+1)))
# Independent path DP for exact order-17 realization.
G,w=decorated_word((0,),0,[(2,"T")]+W[1])
require(len(G)==17 and strong(G))
require(hp_all(G)[-1]==11333)
print("fixed O2 binary, local 7-adic controls, and 7*1619 order-17 witness passed")
