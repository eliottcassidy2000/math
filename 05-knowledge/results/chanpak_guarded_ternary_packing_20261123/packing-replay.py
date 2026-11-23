from itertools import permutations,combinations

p=5
H=(19,61,63,65,23);depth=(1,2,2,2,1)
seeds=[h*pow(3,-r,p)%p for h,r in zip(H,depth)]
assert seeds==[3,4,2,0,1]
assert [ (4-x)%5 for x in seeds]==seeds[::-1]

def legal(W):
    return all((i+1)%5 not in W for i in W)
def hop(i,j):
    assert (i-j)%5 in (2,3)
    return [frozenset((i,j)),frozenset((j,))]
def compile(i,j):
    d=(j-i)%5
    if d==0:return hop(i,(i+2)%5)+hop((i+2)%5,i)
    if d in (2,3):return [
        frozenset((i,j)),frozenset((j,)),frozenset((j,(j+2)%5)),frozenset((j,))]
    k=(i+(-2 if d==1 else 2))%5
    return hop(i,k)+hop(k,j)
for r in range(4):
    for h in range(5):
        i=h*pow(3,-r,5)%5
        for d in (-1,0,1):
            j=(3*h+2*d)*pow(3,-r-1,5)%5
            R=[frozenset((i,))]+compile(i,j)
            assert len(R)==5 and R[-1]==frozenset((j,))
            assert all(W and legal(W) for W in R)
            assert all(len(U^V)==1 for U,V in zip(R,R[1:]))
            assert (pow(4,i,11)*pow(4,2*pow(3,-r-1,5)*d,11))%11==pow(4,j,11)
assert next(j for j in range(1,10) if pow(3,j,5)==1)==4

QR={pow(4,i,11) for i in range(5)}
assert QR=={1,3,4,5,9}
assert {sum(pow(4,j,11) for j in W)%11 for i in range(5)
         if legal(W:={i,(i+2)%5})}==set(range(1,11))-QR
def inv(x):return 3*pow(x,-1,11)%11
assert [x for x in sorted(QR) if inv(x)==x]==[5]
assert all(((y-x)%11 in QR)==((inv(x)-inv(y))%11 in QR)
           for x,y in permutations(QR,2))
pal=lambda N:sum(all((b-a)%11 in QR for a,b in zip(P,P[1:])) for P in permutations(N))
assert pal(QR)==pal(QR|{0})==15

# Golden W5 exact normalization.
q=[[4,1,7],[3,1,6],[6,1,14]]
c=[[1,0,4],[1,0,3],[3,1,6]]
s=[[0,1,1],[0,0,1],[1,0,1]]
M=[[int(i==j) for j in range(3)] for i in range(3)]
for Z in (q,c,c,q,s,s,s,s,q):
    M=[[sum(M[i][k]*Z[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
assert all(z%5==0 for row in M for z in row)
A,B,C=[tuple(x//5 for x in row) for row in M]
assert A[0]*(B[1]*C[2]-B[2]*C[1])-A[1]*(B[0]*C[2]-B[2]*C[0])+A[2]*(B[0]*C[1]-B[1]*C[0])==1

def dp(A):
    n=len(A);N=1<<n;E=[[0]*n for _ in range(N)]
    for i in range(n):E[1<<i][i]=1
    for T in range(1,N):
        for u,h in enumerate(E[T]):
            if not h:continue
            e=A[u]&~T
            while e:
                k=e&-e;e-=k;E[T|k][k.bit_length()-1]+=h
    return E
def rev(A):
    return tuple(((1<<len(A))-1)^r^(1<<i) for i,r in enumerate(A))
def profile(A):
    U=dp(A);V=dp(rev(A));m=(1<<len(A))-1
    h=[1]+list(map(sum,U[1:]))
    return h[-1],U[-1],V[-1],[sum(V[S][i]*h[S^m] for S in range(m+1)) for i in range(len(A))],[
        sum(U[S][i]*h[S^m] for S in range(m+1)) for i in range(len(A))]
def bridge(A,B,a,b):
    k=len(A);n=k+len(B)
    mask=((1<<n)-1)^((1<<k)-1)
    R=[r|mask for r in A]+[r<<k for r in B]
    R[a]^=1<<(k+b);R[k+b]^=1<<a
    return tuple(R)
def tans(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
for A in ((2,8,3,5),(2,4,17,7,11)):
    h,E,S,C,D=profile(A)
    for i in range(len(A)):
        got=set()
        for B in tuple(tans(j) for j in range(1,5))+((2,4,1),(14,4,8,2),(2,8,3,5),(0,)):
            hb,Eb,Sb,Cb,Db=profile(B)
            for j in range(len(B)):
                n=profile(bridge(A,B,i,j))[0]
                assert n==h*hb-E[i]*Sb[j]+C[i]*Db[j]
                if n in (3*h-2,3*h+2):got.add(n)
        assert len(got)<2
from collections import Counter
assert Counter((a+2*b)%5 for a in (-1,0,1) for b in (-1,0,1))=={0:1,1:2,2:2,3:2,4:2}
print("five normalized states, four-phase legal toggles, eleven-state reader, and bridge controls pass")

from collections import Counter
assert Counter((x+y)%5 for x in (-1,0,1) for y in (-2,0,2))=={0:1,1:2,2:2,3:2,4:2}
print('golden two-phase mixing counts pass')
