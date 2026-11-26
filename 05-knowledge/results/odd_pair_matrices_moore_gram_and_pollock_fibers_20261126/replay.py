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
print('Petersen controls',dict(test(petersen(),10)))
print('Hoffman-Singleton controls',dict(test(hs(),50)))
I=(0,1,2);P=(1,2,0);Q=inv(P)
X={(0,1):I,(1,0):I,(0,2):I,(2,0):I,(1,2):P,(2,1):Q}
R={ij:singleton(residual(X,*ij,range(3),3)) for ij in permutations2(range(3))}
assert comp(R[0,1],R[1,2])!=R[0,2]
print('false local completion: old pairs disjoint; last-star cocycle fails',R[0,1],R[1,2],R[0,2])

def T(x):return (3*x+1)//2 if x&1 else x//2
count=0
for x in range(1,20000,2):
    z=x;r=s=0
    while z&1:r+=1;z=T(z)
    while not z&1:s+=1;z=T(z)
    D=2**(r+s);A=3**r
    assert D*z==A*(x+1)-2**r
    assert z+1==(A*(x+1)+D-1)//D
    assert (z>x)==(A>D)
    count+=1
print('Collatz run blocks checked',count)

def det3(X):
    a,b,c=X
    return sum(a[i]*(b[(i+1)%3]*c[(i+2)%3]-b[(i+2)%3]*c[(i+1)%3]) for i in range(3))
def gram(y):
    q,b,v=[sum(t**k for t in y) for k in (2,3,4)]
    G=((5,0,q),(0,q,b),(q,b,v))
    assert det3(G)==sum(((a-b)*(a-c)*(b-c))**2 for a,b,c in combinations(y,3))
    return (q//2,b//6,v),det3(G)
from itertools import combinations_with_replacement
least=10**9;hits=[]
for z in combinations_with_replacement(range(-8,9),4):
    f=-sum(z)
    if not z[-1]<=f<=8:continue
    w=(*z,f);p,D=gram(w)
    if D and D<least:least,hits=D,[]
    if D==least:hits.append(w)
assert least==12 and hits==[(-1,0,0,0,1)]
u=(-17,-5,-4,12,14);v=(-16,-10,3,7,16)
assert gram(u)==gram(v) and prod(u)!=prod(v)
print('Pollock collision:',gram(u),'constant terms',-prod(u),-prod(v))
