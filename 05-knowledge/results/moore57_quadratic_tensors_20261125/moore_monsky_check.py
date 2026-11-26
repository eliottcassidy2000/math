"""Exact binary Moore-to-Monsky lift, with small graph and arithmetic controls."""
from itertools import combinations
from fractions import Fraction
from math import isqrt, prod
from random import Random


def rank(A):
    U={}
    for x in A:
        while x:
            h=x.bit_length()
            if h not in U:
                U[h]=x
                break
            x^=U[h]
    return len(U)


def multiply(A,B):
    C=[]
    for x in A:
        r=0
        while x:
            b=x&-x
            r^=B[b.bit_length()-1]
            x^=b
        C.append(r)
    return C


def graph(n,edges):
    A=[0]*n
    for i,j in edges:
        assert i!=j
        A[i]|=1<<j
        A[j]|=1<<i
    return A


def MooreDefect(A):
    return [y^x for y,x in zip(multiply(A,A),A)]


def fiber_lift(A):
    n=len(A)
    L=[0]*(3*n)
    for i in range(n):
        step=1 if (n-1-i)%2 else -1
        for x in range(3):
            a=3*i+x
            L[a]|=1<<(3*i+(x+step)%3)
            for j in range(i+1,n):
                for y in range(3):
                    b=3*j+y
                    if not (A[i]>>j&1) or x==y:L[a]|=1<<b
                    else:L[b]|=1<<a
    assert all(((L[i]>>j)^(L[j]>>i))&1 for i,j in combinations(range(3*n),2))
    assert all(not (x>>i&1) for i,x in enumerate(L))
    return L


def selmer(T):
    n=len(T)
    L=[r|((r.bit_count()%2)<<i) for i,r in enumerate(T)]
    F=[x^y^(1<<i) for i,(x,y) in enumerate(zip(multiply(L,L),L))]
    M=[(x^(1<<i))|(1<<(i+n)) for i,x in enumerate(L)]+[(1<<i)|(x<<n) for i,x in enumerate(L)]
    assert 2*n-rank(M)==n-rank(F)
    return n-rank(F)


def hoffman_singleton():
    E=[]
    for i in range(5):
        for x in range(5):
            E.extend(((5*i+x,5*i+(x+1)%5),(25+5*i+x,25+5*i+(x+2)%5)))
            E.extend((5*i+x,25+5*j+(x+i*j)%5) for j in range(5))
    return graph(50,E)


def check(A,claim=None,moore=False):
    n=len(A)
    B=MooreDefect(A)
    r=selmer(fiber_lift(A))
    assert r==2*(n-rank(B))
    if claim is not None:assert r==claim
    if moore:
        assert all((((A[i]>>j)&1)+(A[i]&A[j]).bit_count()==1)
                   for i,j in combinations(range(n),2))
    return r


def paths(T):
    n=len(T)
    allbits=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            if not dp[s][i]:continue
            q=T[i]&(allbits^s)
            while q:
                b=q&-q
                q-=b
                dp[s|b][b.bit_length()-1]+=dp[s][i]
    return sum(dp[-1])


def apex(word):
    m=len(word)
    T=[sum(1<<j for j in range(i+1,m))+(1<<m)*int(b) for i,b in enumerate(word)]
    T.append(sum(1<<i for i,b in enumerate(word) if b=='0'))
    return T


def catalan():
    P=(Fraction(2),Fraction(3))
    def add(P,Q):
        x,y=P
        a,b=Q
        if x==a and b==-y:return None
        m=(3*x*x/(2*y)) if P==Q else ((b-y)/(a-x))
        u=m*m-x-a
        return u,m*(x-u)-y
    Q=add(P,P)
    assert Q==(0,1) and add(Q,Q)==(0,-1)
    R=add(P,Q)
    assert R==(-1,0) and add(R,R) is None
    assert 3**2-2*2**2==1


def c5():
    return graph(5,((i,(i+1)%5) for i in range(5)))


def check_primes(T,P):
    assert len(T)==len(P)==len(set(P))
    for i,p in enumerate(P):
        assert p%8==3 and all(p%d for d in range(2,isqrt(p)+1))
        for j,q in enumerate(P):
            if j!=i:assert (pow(q,(p-1)//2,p)==p-1)==bool(T[i]>>j&1)
    return prod(P)


if __name__=='__main__':
    petersen=list(combinations(range(5),2))
    P=graph(10,((i,j) for i,j in combinations(range(10),2)
      if not (set(petersen[i])&set(petersen[j]))))
    R=graph(9,((i,j) for i,j in combinations(range(9),2) if i//3==j//3 or i%3==j%3))
    print("C5",check(c5(),2,True),"Petersen",check(P,18,True),
          "Hoffman–Singleton",check(hoffman_singleton(),98,True),"Paley9",check(R,18))
    rng=Random(2429)
    for n in range(1,21):
        for _ in range(5):check(graph(n,((i,j) for i,j in combinations(range(n),2) if rng.getrandbits(1))))
    bad=graph(10, [*combinations(range(4),2),*((i,j) for i in range(4,7) for j in range(7,10))])
    assert check(bad,16)==16
    assert rank(MooreDefect(bad))==2
    primes=(3,19,83,131,523,1051,1307,12539,13451,128987,155723,205427,280883,612331,1708411)
    T=fiber_lift(c5())
    m=check_primes(T,primes)
    assert m==91050429687411693345720737007317859408027933788445680739
    print("Explicit pentagon twist:",m,"pure 2-Selmer dimension",selmer(T),"Hamiltonian paths",paths(T))
    experiments=[(check(A), paths(fiber_lift(A))) for A in (graph(2,[]),graph(2,[(0,1)]))]
    assert experiments==[(4,9),(2,45)]
    print("Two-fibre controls:",experiments)
    witnesses=[(selmer(T),paths(T)) for T in
               ([*(r|(1<<5) for r in apex('0101')),0],
                [*(r|(1<<5) for r in apex('0001')),0],
                fiber_lift(graph(2,[])))]
    assert witnesses==[(0,9),(2,9),(4,9)]
    print("Six vertices, H=9:",witnesses)
    catalan()
