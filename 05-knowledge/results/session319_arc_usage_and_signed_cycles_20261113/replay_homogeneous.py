from itertools import permutations, combinations
from collections import Counter
from functools import cache
A=(54,64,122,115,98,66,1)
B=(54,32,122,115,98,64,3)
R=(0,1,2,4,5,6)

def require(x):
    if not x: raise RuntimeError('replay mismatch')

def arc(rows,i,j): return rows[i]>>j&1

def good(rows,p): return all(arc(rows,i,j) for i,j in zip(p,p[1:]))

def cyc(rows,p): return good(rows,p) and arc(rows,p[-1],p[0])

p=(0,5,2,3,4,1,6)
C=tuple(sum(1<<j for j in range(7) if arc(A,p[i],p[j])) for i in range(7))
require(C==(54,96,122,115,98,64,1))
require(tuple((u^64 if i==1 else u^2 if i==6 else u) for i,u in enumerate(C))==B)
for rows in (A,B):
    for j,want in ((1,(27,9,1,18)),(2,(9,9,1,9))):
        ps=[p for p in set(permutations(R+(3,)*j)) if good(rows,p)]
        obs=(len(ps),sum(p[0]==3 for p in ps),sum(p[-1]==3 for p in ps),
             sum(sum(i==2 and k==3 for i,k in zip(p,p[1:])) for p in ps))
        require(obs==want)

def polys(rows,z):
    @cache
    def pack(X):
        if not X:return 1
        v,*T=X
        return z*pack(tuple(T))+sum(2*z*pack(tuple(i for i in T if i not in p))
              for l in range(2,len(T)+1,2) for p in permutations(T,l) if cyc(rows,(v,)+p))
    G=[0,0]; fixed=Counter(); mixed=Counter(); pair=[]
    for l in range(2,len(R)+1):
        for p in permutations(R,l):
            if cyc(rows,(3,)+p):
                G[l%2]+=pack(tuple(i for i in R if i not in p))
                mixed[l+1]+=1
                for q in permutations(set(R)-set(p),3):
                    if q[0]==min(q) and cyc(rows,q):pair.append(((3,)+p,q))
        if l>=3:
            for p in permutations(R,l):
                if p[0]==min(p) and cyc(rows,p):fixed[l]+=1
    return (pack(R),*G),fixed,mixed,pair

for z in range(7):
    U,fA,mA,pA=polys(A,z)
    V,fB,mB,pB=polys(B,z)
    require(U==(z**6+8*z**4+8*z**2,z**4+3*z**2+1,z**3+3*z))
    require(V==(z**6+8*z**4+6*z**2,z**4+4*z**2+1,z**3+2*z))
    require(fA-fB==Counter({4:2,5:1}) and not fB-fA)
    require(mA-mB==Counter({5:1,6:1}) and not mB-mA)
    require(mA[7]==mB[7]==1 and pA==[] and pB==[((3,0,2),(1,5,6))])

def expand(rows,U):
    n=len(U); slots=R+(3,)*n
    def edge(i,j):
        a,b=slots[i],slots[j]
        return arc(U,i-6,j-6) if a==b==3 else arc(rows,a,b)
    return tuple(sum(1<<j for j in range(6+n) if edge(i,j)) for i in range(6+n))

def tab(rows):
    n=len(rows);N=1<<n
    d=[[0]*n for _ in range(N)]
    h=[1]+[0]*(N-1)
    for i in range(n):d[1<<i][i]=1
    for m in range(1,N):
        h[m]=sum(d[m])
        for i,v in enumerate(d[m]):
            if v:
                b=rows[i]&~m
                while b:
                    q=b&-b;b-=q;d[m|q][q.bit_length()-1]+=v
    return h,d

def all_q(h,limit):
    N=len(h); q=[0]*N;q[0]=1; vals=[0]
    for _ in range(limit):
        nextq=[0]*N
        for M in range(N):
            S=M
            while True:
                nextq[M]+=q[M^S]*h[S]
                if not S:break
                S=(S-1)&M
        q=nextq; vals.append(q[-1])
    return vals

def arc_use(hs,he,i,j):
    N=len(hs); M=N-1
    return sum(he[S][i]*hs[M^S][j] for S in range(N) if S>>i&1 and not(S>>j&1))

def converse(rows):
    mask=(1<<len(rows))-1
    return tuple(mask^(1<<i)^v for i,v in enumerate(rows))

for U in ((0,),(2,0),(6,4,0),(2,4,1),(14,28,24,16,1)):
    n=len(U); A1=expand(A,U); B1=expand(B,U)
    Hu,_=tab(U); Hu=all_q(Hu,5)
    a,ae=tab(A1);b,be=tab(B1)
    va,vb=all_q(a,4),all_q(b,4)
    require(va[1]==vb[1]==9*(Hu[1]+Hu[2]))
    for z in (2,3,4):
        require(vb[z]-va[z]==z*z*((z-1)*Hu[z+1]-(z+1)*Hu[z-1]))
    _,aS=tab(converse(A1));_,bS=tab(converse(B1))
    for i in range(6,6+n):
        require(ae[-1][i]==be[-1][i] and aS[-1][i]==bS[-1][i])
        # x is outside index 2 in this expansion.
        require(arc_use(aS,ae,2,i)==arc_use(bS,be,2,i))
        for j in range(6,6+n):
            if arc(A1,i,j):require(arc_use(aS,ae,i,j)==arc_use(bS,be,i,j))
    if U in ((0,),(2,0),(6,4,0)):
        for i in range(6,6+n):
            j=i-5
            require(ae[-1][i]==2**(j-1) and aS[-1][i]==9*2**(n-j))
            require(arc_use(aS,ae,2,i)==9*2**(n-j)+9*(j==1))
            for k in range(i+1,6+n):
                require(arc_use(aS,ae,i,k)==9*2**(n-k+i)+9*(k==i+1))
print('quotient and one-edge relabelling verified; five expanded controls; full difference and named marginals pass')
