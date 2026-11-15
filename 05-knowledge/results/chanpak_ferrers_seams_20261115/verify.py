#!/usr/bin/env python3
"""Exhaustive small seam controls, independent word enumeration and tournament DP."""
from itertools import product,combinations
from math import comb

def walk_counts(D):
    p,q=len(D),len(D[0])
    X=[[0]*(q+1) for _ in range(p+1)]
    Y=[[0]*(q+1) for _ in range(p+1)]
    for i in range(1,p+1): X[i][0]=1
    for j in range(1,q+1): Y[0][j]=1
    for i,j in product(range(1,p+1),range(1,q+1)):
        X[i][j]=X[i-1][j]+(1-D[i-1][j-1])*Y[i-1][j]
        Y[i][j]=Y[i][j-1]+D[i-1][j-1]*X[i][j-1]
    F=[[x+y for x,y in zip(xx,yy)] for xx,yy in zip(X,Y)]
    F[0][0]=1
    return F

def counter(D,aa,bb):
    out=0
    for choose in combinations(range(aa+bb),aa):
        choose=set(choose)
        v=[int(t not in choose) for t in range(aa+bb)]
        i=j=-1
        for k,side in enumerate(v):
            if side: j+=1
            else: i+=1
            if k and side!=v[k-1] and D[i][j]!=side: break
        else: out+=1
    return out

def check_seam(D):
    p,q=len(D),len(D[0])
    R=[[D[p-i-1][q-j-1] for i in range(p)] for j in range(q)]
    F,G=walk_counts(D),walk_counts(R)
    for i,j in product(range(p+1),range(q+1)):
        assert F[i][j]==counter(D,i,j) and G[j][i]==counter(R,j,i)
    L=sum(F[i][j]*G[q-j][p-i]-1 for i,j in product(range(p+1),range(q+1)))
    P={(i,j) for i,j in product(range(p),range(q)) if F[i+1][j+1]>1}
    S={(i,j) for i,j in product(range(p),range(q)) if G[q-j][p-i]>1}
    assert len(P)+len(S)<=L
    U=P|S
    M=U.copy()
    for i,j in product(range(p),range(q)):
        if not D[i][j] or (i,j) in U: continue
        for u,w in product(range(i+1),range(j,q)):
            if D[u][w] or (u,w) in U: continue
            assert (i,w) in P and (u,j) in S
            M.add((i,j))
    assert len(M)<=L+L*L//4
    E=[[int(any(D[u][w] and (u,w) not in M for u,w in product(range(i,p),range(j+1)))) for j in range(q)] for i in range(p)]
    assert sum(x!=y for u,v in zip(D,E) for x,y in zip(u,v))<=len(M)
    assert all(E[i][j]<=E[i-1][j] for i,j in product(range(1,p),range(q)))
    assert all(E[i][j]>=E[i][j-1] for i,j in product(range(p),range(1,q)))
    badones=[(i,j) for i,j in product(range(p),range(q)) if D[i][j]==0]
    if len(badones)==1:
        i,j=badones[0]; i+=1; j+=1
        bound=(p-i+(j>1))*(q-j+1)+(j-1+(i<p))*i
        assert L>=bound and ((i,j)==(p,1) or bound>=p+q-1)
    return L

def hp(A):
    p=len(A); W=1<<p
    d=[[0]*p for _ in range(W)]
    for i in range(p):d[1<<i][i]=1
    for S in range(1,W):
        for i in range(p):
            if d[S][i]:
                for j in range(p):
                    if not(S>>j&1) and A[i][j]: d[S|1<<j][j]+=d[S][i]
    return sum(d[-1])

def chain(parts):
    n=1+sum(s-1 for s in parts); A=[[int(i<j) for j in range(n)] for i in range(n)]
    i=0
    for s in parts:
        j=i+s-1; A[i][j]=0; A[j][i]=1; i=j
    return A

def run():
    num=0
    for p,q in [(2,2),(2,3),(3,2),(3,3),(3,4),(4,3)]:
        for bits in product(range(2),repeat=p*q):
            check_seam([bits[i*q:(i+1)*q] for i in range(p)])
            num+=1
    print(f'{num} interleaving matrices: literal-word enumeration, edit certificate and single reversal bound')
    corner=[]
    for a in (3,4,5):
        line=[]
        for b in (3,4,5):
            A=chain([a,b]); i,j=a-2,a
            A[i][j]=0; A[j][i]=1
            line.append(hp(A)-hp(chain([a,b])))
        corner.append(line)
    assert corner==[[2,4,6],[4,6,8],[6,8,10]],corner
    for a,b in product((3,4,5),repeat=2):
        for left,right in product((False,True),repeat=2):
            parts=([4] if left else [])+[a,b]+([4] if right else [])
            A=chain(parts); i,j=(3*left+a-2,3*left+a)
            base=hp(A); A[i][j]=0; A[j][i]=1
            assert hp(A)>=base+corner[a-3][b-3]*5**(left+right)
    print(f'corner surplus {corner}; 36 whole-tournament DP controls')
if __name__=='__main__': run()
