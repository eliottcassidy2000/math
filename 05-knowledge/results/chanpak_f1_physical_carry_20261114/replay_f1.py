from dataclasses import dataclass
from itertools import product
from math import gcd, prod

@dataclass(frozen=True)
class Core:
    masks:tuple
    U:int
    h:int
    f:int
    P:int
    Q:int
    S:int
    T:int
    K:int

core15=(
    Core((18,20,9,3,12),3,15,120,42,63,6,9,41),
    Core((18,20,9,3,12),5,15,120,42,63,6,9,43),
)
core29=(
    Core((10,44,57,48,35,1),24,29,320,93,198,9,20,91),
    Core((12,45,56,48,35,1),12,29,320,93,198,9,20,95),
)
core77=(
    Core((58,32,115,6,10,88,27,127),250,77,2096,1622,380,62,0,419),
    Core((254,12,40,128,174,74,30,102),111,77,2096,1622,380,62,0,421),
)

core5=(
    Core(tuple(x&511 for x in (510,1008,994,486,1004,448,384,768,512)),105,5,1024,716,810,0,0,579),
    Core(tuple(x&511 for x in (510,1020,504,496,896,272,816,864,512)),45,5,1024,716,810,0,0,581),
)
core11=(
    Core((1022,1020,1016,704,520,536,560,624,760,0),812,11,3200,671,929,11,11,209),
    Core((1022,996,896,902,910,924,956,768,512,0),678,11,3200,671,929,11,11,213),
)
core17=(
    Core((56,57,51,4,40,8),17,17,224,72,135,6,11,65),
    Core((54,32,26,51,2,20),10,17,224,72,135,6,11,69),
)

def path_arrays(rows):
    m=1<<len(rows)
    K=[1]+[0]*(m-1)
    D=[[0]*len(rows) for i in range(m)]
    for s in range(1,m):
        rem=s
        while rem:
            b=rem&-rem; rem^=b
            v=b.bit_length()-1
            t=s^b
            if not t: D[s][v]=1
            else:
                u=t&rows[v]
                while u:
                    a=u&-u;u^=a
                    D[s][v]+=D[t][a.bit_length()-1]
        K[s]=sum(D[s])
    return K,D

def make(sources,x=0,y=0):
    A=bin(x)[2:] if x else ''
    B=''.join('10'[int(ch)] for ch in reversed(bin(y)[2:])) if y else ''
    labels=[((0,),int(ch)) for ch in A]+[(c.masks,c.U) for c in sources]+[((0,),int(ch)) for ch in B]
    n=sum(len(R) for R,U in labels)
    rows=[]; apex=0
    for R,U in labels:
        j=len(rows); m=len(R)
        apex|=U<<j
        across=(1<<n)-(1<<(j+m))
        for k in range(m):
            rows.append((R[k]<<j)|across|((not (U>>k&1))<<n))
    return rows+[apex]

def audit(C):
    K,D=path_arrays(C.masks)
    m=len(K)-1
    R=tuple((m^(1<<i)^a) for i,a in enumerate(C.masks))
    _,E=path_arrays(R)
    U=[i for i in range(len(R)) if C.U>>i&1]
    V=[i for i in range(len(R)) if not(C.U>>i&1)]
    Y=lambda X,W: sum(X[s][v]*K[m^s] for v in W for s in range(1,m+1))
    out=(K[-1],sum(K[s]*K[m^s] for s in range(m+1)),Y(D,U),Y(E,V),sum(E[-1][v] for v in U),sum(D[-1][v] for v in V))
    assert out == (C.h,C.f,C.P,C.Q,C.S,C.T)
    assert path_arrays(make([C]))[0][-1]==C.K

def inverse(k,a,b,A,B):
    # The permissible j interval, with a single residue condition.
    lo=max(0,(k-A*a+B-1)//B)
    hi=min(b,k//B)
    r=k*pow(B,-1,A)%A
    j=lo+(r-lo)%A
    if j>hi:return None
    i=(k-B*j)//A
    assert 0<=i<=a and 0<=j<=b and A*i+B*j==k
    return i,j

if __name__=='__main__':
    for C in core5+core11+core15+core17+core29+core77:audit(C)
    Ns={}
    for i,j in product(range(2),repeat=2):
        Ns[i,j]=path_arrays(make([core15[i],core29[j]],2,2))[0][-1]
    assert all(N==Ns[0,0]+2*(29*i+30*j) for (i,j),N in Ns.items())
    for a,b,A,B in [(29,28,29,30),(50,70,29,30),(153,28,29,154),(200,50,29,154),(9,16,17,10),(9,10,11,10)]:
        C=(A-1)*(B-1)
        S={A*i+B*j for i in range(a+1) for j in range(b+1)}
        for k in range(A*a+B*b+1):
            assert (inverse(k,a,b,A,B) is not None)==(k in S)
        assert all(k in S for k in range(C,A*a+B*b-C+1))
        assert C-1 not in S and A*a+B*b-C+1 not in S
    assert gcd(2096*320,77*29)==1
    for a in range(1,15):
        for b in range(1,15):
            A=[core5[0]]*a+[core11[0]]*b
            Q=sum(C.Q*prod(X.f for X in A[:i])*prod(X.h for X in A[i+1:]) for i,C in enumerate(A))
            assert Q%5==pow(4,a-1,5)
    print("matched literal profiles, two-wing neutral unit, and conductor inverses passed", Ns)
