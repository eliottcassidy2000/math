#!/usr/bin/env python3
"""Exact, bounded parity/interface checks. Standard library only."""
from itertools import product
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import random, json, hashlib, subprocess

@lru_cache(None)
def pf(M):
    n=len(M)
    if n%2: return 0
    if not n: return 1
    return sum((-1)**(j+1)*M[0][j]*pf(tuple(tuple(M[a][b] for b in range(1,n) if b!=j) for a in range(1,n) if a!=j)) for j in range(1,n))
def mat(M): return tuple(map(tuple,M))
def cof(M): return [(-1)**i*pf(tuple(tuple(M[a][b] for b in range(len(M)) if b!=i) for a in range(len(M)) if a!=i)) for i in range(len(M))]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def tournament(n,bits):
    M=[[0]*n for _ in range(n)]; k=0
    for i in range(n):
        for j in range(i+1,n): M[i][j]=1 if bits>>k&1 else -1; M[j][i]=-M[i][j]; k+=1
    return mat(M)
def join(A,B,u,v):
    m,n=len(A),len(B)
    return mat([list(A[i])+[u[i]*v[j] for j in range(n)] for i in range(m)]+[[-v[i]*u[j] for j in range(m)]+list(B[i]) for i in range(n)])
def border(A,x): return mat([list(r)+[x[i]] for i,r in enumerate(A)]+[[-t for t in x]+[0]])
def inv(A):
    n=len(A); R=[[F(v) for v in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for c in range(n):
        p=next(i for i in range(c,n) if R[i][c]); R[c],R[p]=R[p],R[c]
        d=R[c][c]; R[c]=[x/d for x in R[c]]
        for i in range(n):
            if i!=c:
                d=R[i][c]; R[i]=[x-d*y for x,y in zip(R[i],R[c])]
    return [r[n:] for r in R]
def bil(a,I,b): return sum(a[i]*I[i][j]*b[j] for i in range(len(a)) for j in range(len(b)))
def det(A):
    R=[[F(x) for x in row] for row in A]; out=F(1); n=len(R)
    for c in range(n):
        nz=[i for i in range(c,n) if R[i][c]]
        if not nz: return F(0)
        p=nz[0]
        if p!=c: R[c],R[p]=R[p],R[c]; out=-out
        d=R[c][c]; out*=d
        for i in range(c+1,n):
            q=R[i][c]/d
            for j in range(c+1,n): R[i][j]-=q*R[c][j]
    return out

def main():
    counts={'rank_one_even_total':0,'mixed_rooted':0,'multiport_even':0,'multiport_odd':0,'determinant_crosschecks':0,'weighted_singular_controls':0,'cofactor_covariance':0}
    # All ordered positive block sizes of even total <=6; all tournament blocks;
    # all rank-one sign matrices (u[0]=1 removes their double representation).
    for total in (2,4,6):
        for m in range(1,total):
            n=total-m
            for ai in range(1<<(m*(m-1)//2)):
                A=tournament(m,ai)
                for bi in range(1<<(n*(n-1)//2)):
                    B=tournament(n,bi)
                    for uv in product((-1,1),repeat=total-1):
                        u=(1,)+uv[:m-1]; v=uv[m-1:]; S=join(A,B,u,v)
                        expected=pf(A)*pf(B) if m%2==0 else dot(cof(A),u)*dot(cof(B),v)
                        if pf(S)!=expected: raise AssertionError(('rankone',A,B,u,v))
                        d=[1 if i%3==0 else -1 for i in range(total)]
                        SD=mat([[d[i]*S[i][j]*d[j] for j in range(total)] for i in range(total)])
                        if pf(SD)!=__import__('math').prod(d)*expected: raise AssertionError('switch')
                        if pf(mat([[-x for x in row] for row in S]))!=(-1)**(total//2)*expected: raise AssertionError('reversal')
                        if m%2:
                            DA=d[:m]; Ap=mat([[DA[i]*A[i][j]*DA[j] for j in range(m)] for i in range(m)])
                            if cof(Ap)!=[__import__('math').prod(DA)*DA[i]*cof(A)[i] for i in range(m)]: raise AssertionError('cofactor covariance')
                            counts['cofactor_covariance']+=1
                        counts['rank_one_even_total']+=1
    # Complete 2+1, 2+3, 4+1 mixed parity attachment universes, all roots.
    for m,n in ((2,1),(2,3),(4,1)):
        for ai in range(1<<(m*(m-1)//2)):
            A=tournament(m,ai); I=inv(A)
            for bi in range(1<<(n*(n-1)//2)):
                B=tournament(n,bi); p=cof(B)
                for uv in product((-1,1),repeat=m+n-1):
                    u=(1,)+uv[:m-1]; v=uv[m-1:]; S=join(A,B,u,v)
                    for w in product((-1,1),repeat=m+n):
                        value=pf(A)*(dot(p,w[m:])+dot(p,v)*bil(u,I,w[:m]))
                        if pf(border(S,w))!=value: raise AssertionError('mixed')
                        counts['mixed_rooted']+=1
    rng=random.Random(102)
    for trial in range(300):
        m=rng.randrange(1,6); n=rng.choice([j for j in range(1,6) if (m+j)%2==0])
        def weighted(q):
            M=[[0]*q for _ in range(q)]
            for i in range(q):
                for j in range(i+1,q): M[i][j]=rng.randrange(-2,3); M[j][i]=-M[i][j]
            return mat(M)
        A=weighted(m); B=weighted(n); u=[rng.randrange(-2,3) for _ in range(m)]; v=[rng.randrange(-2,3) for _ in range(n)]
        S=join(A,B,u,v)
        expected=pf(A)*pf(B) if m%2==0 else dot(cof(A),u)*dot(cof(B),v)
        if pf(S)!=expected: raise AssertionError('weighted rank-one')
        counts['weighted_singular_controls']+=1
    for n in range(1,8):
        for k in range(1,5):
            if (n+k)%2: continue
            for trial in range(25):
                A=tournament(n,rng.randrange(1<<(n*(n-1)//2)))
                U=[[rng.choice((-1,1)) for _ in range(k)] for _ in range(n)]
                E=tournament(k,rng.randrange(1<<(k*(k-1)//2)))
                whole=mat([list(A[i])+U[i] for i in range(n)]+[[-U[j][i] for j in range(n)]+list(E[i]) for i in range(k)])
                if n%2==0:
                    I=inv(A); cols=list(zip(*U)); K=mat([[E[i][j]+bil(cols[i],I,cols[j]) for j in range(k)] for i in range(k)])
                    expected=pf(A)*pf(K); counts['multiport_even']+=1
                else:
                    # Every possible pivot, with permutation sign retained.
                    for a in range(n):
                        inds=[i for i in range(n) if i!=a]; C=mat([[A[i][j] for j in inds] for i in inds]); I=inv(C)
                        b=[A[i][a] for i in inds]; V=[[U[i][j] for i in inds] for j in range(k)]
                        gamma=[U[a][j]+bil(b,I,V[j]) for j in range(k)]
                        K=[[E[i][j]+bil(V[i],I,V[j]) for j in range(k)] for i in range(k)]
                        residual=mat([[0]+gamma]+[[-gamma[i]]+K[i] for i in range(k)])
                        expected=(-1)**a*pf(C)*pf(residual)
                        if pf(whole)!=expected: raise AssertionError('odd multiport')
                        counts['multiport_odd']+=1
                if pf(whole)!=expected: raise AssertionError('multiport')
                if det(whole)!=pf(whole)**2: raise AssertionError('independent determinant')
                counts['determinant_crosschecks']+=1
    A=mat([[0,1],[-1,0]]); u=(1,1); B=((0,),); S=join(A,B,u,(1,)); Iplus=inv([[1,1],[-1,1]])
    hostiles=[]
    for x in ((1,1),(1,-1)):
        W=border(S,x+(1,))
        hostiles.append({'x':x,'matrix':W,'pfaffian':pf(W),'determinant':int(det(W)), 'individual_normalized_root_energy':str(bil(x,Iplus,x)), 'mixed_response_u_Ainv_x':str(bil(u,inv(A),x))})
    out={'status':'PASS','universe':counts,'hostile':hostiles,'claims':'bounded exhaustive rank-one and mixed rooted checks; seeded multiport checks; all-size proofs separate','seed':102}
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
