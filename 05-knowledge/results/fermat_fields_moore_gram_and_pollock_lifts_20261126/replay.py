from itertools import combinations, product
from collections import Counter
from fractions import Fraction as F
from math import prod

def comp(a,b):return tuple(b[i] for i in a) # route a then b
def inv(a):return tuple(a.index(i) for i in range(len(a)))
def matrix(A,root):
    d=len(A[root]);B=[]
    for n in sorted(A[root]):B.append(sorted(A[n]-{root}))
    assert len(set().union(*map(set,B)))==d*(d-1)
    X={(i,j):tuple(B[j].index(next(iter(A[u]&set(B[j])))) for u in B[i])
       for i in range(d) for j in range(d) if i!=j}
    return X,d
def graph(n,edges):
    A=list(map(set,[[]]*n))
    for i,j in edges:A[i].add(j);A[j].add(i)
    return A
def petersen():
    lab=list(combinations(range(5),2))
    return graph(10,((i,j) for i,j in combinations(range(10),2) if not set(lab[i])&set(lab[j])))
def hs():
    e=[]
    for i,x in product(range(5),repeat=2):
        e.extend([(5*i+x,5*i+(x+1)%5),(25+5*i+x,25+5*i+(x+2)%5)])
        e.extend((5*i+x,25+5*j+(x+i*j)%5) for j in range(5))
    return graph(50,e)
def residual(X,i,j,O,b):
    p=[X[i,j]]+[comp(X[i,k],X[k,j]) for k in O if k not in (i,j)]
    return tuple(tuple(sorted(set(range(b))-{v[u] for v in p})) for u in range(b))
def singleton(R):
    assert all(len(r)==1 for r in R)
    x=tuple(r[0] for r in R);assert len(set(x))==len(x);return x

def test(A,nroots):
    tests=Counter()
    for root in range(nroots):
        X,d=matrix(A,root);b=d-1
        for k in (1,2):
            for N in combinations(range(d),k):
                O=[i for i in range(d) if i not in N]
                Rs={(i,j):residual(X,i,j,O,b) for i,j in permutations2(O)}
                for i,j in permutations2(O):
                    for u in range(b):
                        assert set(Rs[i,j][u])=={comp(X[i,a],X[a,j])[u] for a in N}
                        assert len(Rs[i,j][u])==k
                if k==1:
                    R={ij:singleton(x) for ij,x in Rs.items()}
                    for i,j,l in combinations(O,3):assert comp(R[i,j],R[j,l])==R[i,l]
                    anchor=O[0];a=N[0];ren=X[anchor,a]
                    for i in O:
                        if i!=anchor:assert comp(R[i,anchor],ren)==X[i,a]
                else:
                    a,c=N
                    for i in O:
                        p=[X[i,a]]+[comp(X[i,j],X[j,a]) for j in O if i!=j]
                        D=singleton(tuple(tuple(set(range(b))-{v[u] for v in p}) for u in range(b)))
                        assert comp(inv(X[i,c]),D)==X[c,a]
                tests[k]+=1
    return tests
def permutations2(O):
    return ((i,j) for i in O for j in O if i!=j)

from math import isqrt
def order(x,p):
    if x%p==0:return None
    f=p-1
    while f%2==0 and pow(x,f//2,p)==1:f//=2
    assert pow(x,f,p)==1
    return f
u=(-17,-5,-4,12,14);v=(-16,-10,3,7,16)
def moments(y):return [sum(z**i for z in y) for i in range(1,5)]
assert moments(u)==moments(v)
for m in range(5):
    p=2**(2**m)+1
    assert all(p%d for d in range(2,isqrt(p)+1))
    l=3*pow(2,-1,p)%p
    a=order(2,p),order(l,p),order(l*l,p)
    assert a[0]==2**(m+1)
    if m>=2:assert a[1:]==(p-1,(p-1)//2)
    q=(2**(2**m+2)+2*(3**(2**m)-1))%p
    if m>=1:assert q
    ranks=tuple(min(4,len({t%p for t in w})) for w in (u,v))
    print('Fermat',p,'orders (2,O,OO)',a,'two-color witness',q,'Pollock ranks',ranks,'Gram determinant',178158400%p)

def mul(A,B):
    return [[sum(x*z for x,z in zip(row,col)) for col in zip(*B)] for row in A]
def Kold(X,O,b,k):
    n=len(O)*b;K=[[0]*n for _ in range(n)]
    for i,x in enumerate(O):
        for t in range(b):K[i*b+t][i*b+t]=k
        for j,y in enumerate(O):
            if i!=j:
                for row,cols in enumerate(residual(X,x,y,O,b)):
                    assert len(cols)==k
                    for c in cols:K[i*b+row][j*b+c]=1
    return K

def certificate(K,q,k):
    n=len(K);A=mul(K,K)
    if k==1:return all(A[i][j]==q*K[i][j] for i in range(n) for j in range(n))
    B=mul(A,K)
    return all(B[i][j]-2*q*A[i][j]+(q*q-1)*K[i][j]==2*(q-1) for i in range(n) for j in range(n))
for name,A in [('Petersen',petersen()),('Hoffman-Singleton',hs())]:
    X,d=matrix(A,0);b=d-1;count=0
    for k in (1,2):
        q=d-k
        for N in combinations(range(d),k):
            O=[i for i in range(d) if i not in N]
            K=Kold(X,O,b,k)
            assert certificate(K,q,k)
            if k==2:
                a,c=N
                for i in range(b):
                    used={comp(X[a,j],X[j,c])[i] for j in O}
                    assert used==set(range(b))-{X[a,c][i]}
            count+=1
    print(name,'exact spectral controls',count)

I=(0,1,2);P=(1,2,0);Q=inv(P)
X={(0,1):I,(1,0):I,(0,2):I,(2,0):I,(1,2):P,(2,1):Q}
assert not certificate(Kold(X,list(range(3)),3,1),3,1)
print('width-three locally valid false star rejected')

def solve(A,b,p):
    M=[[(x%p) for x in r]+[y%p] for r,y in zip(A,b)]
    for i in range(len(A)):
        k=next(k for k in range(i,len(A)) if M[k][i])
        M[i],M[k]=M[k],M[i]
        q=pow(M[i][i],-1,p);M[i]=[x*q%p for x in M[i]]
        for j in range(len(A)):
            if i!=j:
                q=M[j][i];M[j]=[(a-q*c)%p for a,c in zip(M[j],M[i])]
    return [r[-1] for r in M]

for p in [5,17,257,65537]:
    assert sorted(z%(p*p) for z in u)!=sorted(z%(p*p) for z in v)
    goal=[a+p*p*b for a,b in zip(moments(u),[1,-2,3,-4])]
    for w in [u,v]:
        J=next(J for J in combinations(range(5),4) if len({w[j]%p for j in J})==4)
        A=[[d*w[j]**(d-1) for j in J] for d in range(1,5)]
        r=solve(A,[(a-b)//(p*p) for a,b in zip(goal,moments(w))],p)
        x=list(w)
        for j,t in zip(J,r):x[j]+=p*p*t
        assert all((a-b)%(p**3)==0 for a,b in zip(moments(x),goal))
print('four-moment collision: two distinct branches lifted at every Fermat prime >=5 in the list')
