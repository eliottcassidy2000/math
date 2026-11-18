"""Independent subset counts of physical coding moves, no repo imports."""
from itertools import combinations, product
from math import isqrt
from ans_replay import run
def rev(A):
    return tuple(sum(((A[j]>>i)&1)<<j for j in range(len(A))) for i in range(len(A)))
def masks(E):
    A=[0]*(1+max(max(u,v) for u,v in E))
    for u,v in E:A[u]|=1<<v
    return tuple(A)
def rename(A,P):
    return masks(tuple((P[u],P[v]) for u in range(len(A)) for v in range(len(A)) if A[u]>>v&1))
def dp(A):
    n=len(A);M=(1<<n)-1
    assert all((A[u]>>v&1)+(A[v]>>u&1)==1 for u,v in combinations(range(n),2))
    E=[[0]*n for _ in range(M+1)]
    for v in range(n):E[1<<v][v]=1
    for X in range(1,M+1):
        for u,c in enumerate(E[X]):
            if c:
                Y=A[u]&~X
                while Y:
                    b=Y&-Y;Y^=b
                    E[X|b][b.bit_length()-1]+=c
    return E
def H(A):return sum(dp(A)[-1])
def state(A,r):
    E=dp(A);F=dp(rev(A));M=len(E)-1;h=[1]+[sum(x) for x in E[1:]]
    return h[M],E[M][r],sum(F[X][r]*h[M^X] for X in range(1,M+1)),h[M^(1<<r)]
def forward(A,r,w):
    for b in w:
        n=len(A)
        A=tuple((a|(1<<n))^(int(i==r)*b<<n) for i,a in enumerate(A))+(b<<r,)
    return A
def flip(A,u,v):
    B=list(A);B[u]^=1<<v;B[v]^=1<<u;return tuple(B)
def frameflip(A,P):
    assert len(P)==len(A) and len(set(P))==len(A)
    B=A
    for u,v in zip(P,P[1:]):
        assert A[u]>>v&1
        B=flip(B,u,v)
    return B
def wrap(D,I,w):
    n=len(D);v=n+w
    return tuple(a|(((I>>i)&1)<<v) for i,a in enumerate(D))+tuple(
        ((1<<n)-1)|sum(1<<j for j in range(i+1,v)) for i in range(n,v)
    )+(((1<<v)-1)^I,)
def coeff(D,I):
    E=dp(D);F=dp(rev(D));M=len(E)-1
    h=[1]+[sum(x) for x in E[1:]]
    e=[sum(z for i,z in enumerate(row) if I>>i&1) for row in E]
    s=[sum(z for i,z in enumerate(row) if I>>i&1) for row in F]
    return sum(e[X]*h[M^X] for X in range(1,M+1)),h[M]-sum(e[X]*s[M^X] for X in range(1,M))
def lanes():
    rows=(
       ((56,17,3,6,12,30),25, (0,5), (0,5,2,1,4,3), (168,-17),(126,-13)),
       ((36,1,34,7,15,26),29, (2,5), (1,0,2,5,4,3), (172,-79),(129,-59)))
    for i,(D,I,(u,v),P,X,Y) in enumerate(rows):
        E=flip(D,u,v)
        assert coeff(D,I)==X and coeff(E,I)==Y
        for w in range(1,6):
            A,B=wrap(D,I,w),wrap(E,I,w)
            assert H(A)==X[0]*(1<<w)+X[1] and H(B)==Y[0]*(1<<w)+Y[1]
            assert 4*H(B)==3*H(A)+(1 if i else -1)
            if w>=2:
                R=(6,*range(8,w+6),*P,w+6,7)
                C=frameflip(A,R)
                Z=R[::-1];j=Z.index(u)
                Z=Z[j:]+Z[:j]
                assert Z[-1]==v
                perm=list(range(w+7));perm[6],perm[7]=7,6
                assert rename(frameflip(C,Z),perm)==B
                if w==2: print('rim representatives',H(A),H(C),H(B))
            elif i:
                assert rename(frameflip(A,(7,6,5,4,3,2,1,0)),(3,1,0,5,4,7,6,2))==B
            else:
                C=frameflip(A,(0,3,2,1,4,7,6,5))
                D1=frameflip(B,(0,7,6,5,3,1,4,2))
                assert H(C)==255 and rename(C,(3,1,0,4,2,5,6,7))==D1
    print('width 1 boundaries, affine coefficients and two-line lifts verified')

def main():
    for A,scale in (((4,9,18,21,35,15,63,127),1),((4,9,66,53,39,135,59,95),7)):
        S=tuple(scale*x for x in (29,11,31,13))
        assert state(A,0)==S
        for w in ((1,0,1,0,0,1,1),(1,1,1,1,1,1,1)):
            for tail in ((),(1,)):
                V=state(forward(A,0,w+tail),0)
                assert V==run(S,w+tail)
                print('physical continuation',V)
    A=(126,108,56,112,98,64,4)
    B=flip(A,2,5);R=(23,2,50,9)
    assert state(A,3)==R and state(B,3)==(25,*R[1:])
    for k in range(5):
        for w in product((0,1),repeat=k):
            a,b=map(lambda G:state(forward(G,3,w),3),(A,B))
            assert a==run(R,w) and b==(a[0]+2,*a[1:])
    w=(1,0,0,0,1,0,1,1,0,1)
    for T,p in ((A,32561),(B,32563)):
        assert H(forward(T,3,w))==p
        assert all(p%j for j in range(2,isqrt(p)+1))
    print('same arc and same continuation give twin primes 32561,32563')
    lanes()
if __name__=='__main__':main()
