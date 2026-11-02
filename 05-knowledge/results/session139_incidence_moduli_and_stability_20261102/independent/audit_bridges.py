#!/usr/bin/env python3
"""Independent exact audit: no SymPy; bivariate rational polynomial arithmetic.
Run normally and with python3 -O. Claims are checked with explicit exceptions.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from collections import Counter
import argparse

class P:
    def __init__(self,x=0):
        self.d={k:F(v) for k,v in x.items() if v} if isinstance(x,dict) else ({(0,0):F(x)} if x else {})
    def __add__(self,o):
        o=o if isinstance(o,P) else P(o); d=self.d.copy()
        for k,v in o.d.items(): d[k]=d.get(k,0)+v
        return P(d)
    __radd__=__add__
    def __neg__(self): return P({k:-v for k,v in self.d.items()})
    def __sub__(self,o): return self+-as_p(o)
    def __rsub__(self,o): return as_p(o)+-self
    def __mul__(self,o):
        o=as_p(o); d={}
        for (i,j),v in self.d.items():
            for (k,l),w in o.d.items(): d[i+k,j+l]=d.get((i+k,j+l),0)+v*w
        return P(d)
    __rmul__=__mul__
    def __eq__(self,o): return self.d==as_p(o).d
    def __repr__(self): return repr({str(k):str(v) for k,v in sorted(self.d.items())})
    def at(self,a,b): return sum(v*a**i*b**j for (i,j),v in self.d.items())
def as_p(x): return x if isinstance(x,P) else P(x)
def check(ok,label):
    if not ok: raise RuntimeError(label)
def matmul(A,B): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def transpose(A): return list(map(list,zip(*A)))
def eye(n): return [[P(int(i==j)) for j in range(n)] for i in range(n)]
def det(A):
    n=len(A); out=P()
    for p in permutations(range(n)):
        term=P((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))
        for i in range(n): term=term*A[i][p[i]]
        out+=term
    return out
def characteristic(A):
    n=len(A); return [sum(((-1)**k)*det([[A[i][j] for j in s] for i in s]) for s in combinations(range(n),k)) for k in range(n+1)]
def triangle(A,s):
    i,j,k=s; return A[i][j]*A[j][k]*A[k][i]+A[i][k]*A[k][j]*A[j][i]
def inv(A):
    n=len(A); H=P()
    for p in permutations(range(n)):
        t=P(1)
        for i,j in zip(p,p[1:]): t*=A[i][j]
        H+=t
    V=sum(A[i][j]*A[j][i] for i,j in combinations(range(n),2))
    tau=sum(triangle(A,s) for s in combinations(range(n),3))
    M23=sum(A[i][j]*A[j][i]*triangle(A,[k for k in range(n) if k not in (i,j)]) for i,j in combinations(range(n),2)) if n==5 else P()
    u=det([[A[i][j]+int(i==j) for j in range(n)] for i in range(n)])
    v=-det(A)
    return H,V,tau,M23,u,v,characteristic(A)

def cospectral_audit(mutate):
    a=P({(1,0):1});b=P({(0,1):1}); n=5
    K=[[P() for j in range(n)] for i in range(n)]
    K[0][3]=a;K[3][0]=-a;K[1][4]=b;K[4][1]=-b
    A=[[(int(i!=j)+K[i][j])*F(1,2) for j in range(n)] for i in range(n)]
    w=[1,1,-2,0,0]
    Q=[[P(int(i==j)-F(w[i]*w[j],3)) for j in range(n)] for i in range(n)]
    check(matmul(Q,transpose(Q))==eye(n),'orthogonality')
    check(all(sum(row)==1 for row in Q),'Q fixes all-one root')
    B=matmul(matmul(Q,A),transpose(Q))
    check(all(B[i][j]+B[j][i]==int(i!=j) for i in range(n) for j in range(n)),'weighted complement relation')
    for row in B:
        for p in row:
            constant=p.d.get((0,0),0); radius=sum(abs(v) for k,v in p.d.items() if k!=(0,0))
            check(0<=constant-radius<=constant+radius<=1,'weighted domain all |a|,|b|<=1')
    HA,VA,tA,MA,uA,vA,pA=inv(A); HB,VB,tB,MB,uB,vB,pB=inv(B)
    check(pA==pB,'equal determinant-zeta polynomial')
    check(VA==VB and tA==tB and uA==uB and vA==vB,'scalar invariants')
    check(HA==F(15,2),'constant source Hamiltonian')
    coeff=F(1,55) if mutate=='hamiltonian' else F(1,54)
    check(HB-HA==P({(2,2):coeff}),'H difference a^2 b^2/54')
    check(MB-MA==P({(2,2):F(1,108)}),'M23 difference a^2 b^2/108')
    check(MA==(10-a*a-b*b)*F(1,16),'M23 source formula')
    point=(F(1,4),F(1,2))
    exact=[p.at(*point) for p in [HA,HB,VA,tA,MA,MB,uA,vA]]
    check(exact[2:4]==[F(155,64),F(5,2)],'V,tau numerical witness')
    check(exact[-2:]==[F(233,1024),F(-37,256)],'u,v numerical witness')
    check(exact[1]-exact[0]==F(1,3456),'numerical H/Delta difference')
    print('Symbolic family PASS; polynomial H difference',HB-HA,'M23 difference',MB-MA)
    print('Point values H_A,H_B,V,tau,M23_A,M23_B,u,v:',','.join(map(str,exact)))
    print('Common det(I-zA) coefficients:',[str(p.at(*point)) for p in pA])
    print('Endpoints a=0 or b=0 retain H; nonzero ab loses H despite full spectrum and all-one-root conjugacy.')

# Binary checks independently compute cyclic triples via outdegrees, not
# Pfaffian signs; exhaustive n<=6 and sharp ordinal C3 witnesses through n=20.
def binary_stats(A):
    n=len(A); e={}
    for s in combinations(range(n),4):
        cyc=sum(all(sum(A[i][j] for j in t if j!=i)==1 for i in t) for t in combinations(s,3))
        i,j,k,l=s
        pf=(2*A[i][j]-1)*(2*A[k][l]-1)-(2*A[i][k]-1)*(2*A[j][l]-1)+(2*A[i][l]-1)*(2*A[j][k]-1)
        check(cyc%2==(pf*pf-1)//8,'cyclic/Pf parity')
        e[s]=cyc%2
    return e

def cocycle_audit(mutate):
    distributions={}
    for n in range(4,7):
        pairs=list(combinations(range(n),2)); dist=Counter()
        for mask in range(1<<len(pairs)):
            A=[[0]*n for _ in range(n)]
            for b,(i,j) in enumerate(pairs): A[i][j]=(mask>>b)&1;A[j][i]=1-A[i][j]
            e=binary_stats(A);k=sum(e.values());dist[k]+=1
            for s in combinations(range(n),5): check(sum(e[t] for t in combinations(s,4))%2==0,'5-set cocycle identity')
            check(k==0 or k>=n-3,'nonzero support lower bound')
        distributions[n]=dict(sorted(dist.items()))
    check(distributions[5]=={0:384,2:640},'old five-set census recovered')
    print('Full labelled binary universes n=4,5,6 PASS:',distributions)
    for n in range(4,21):
        A=[[int(i<j) for j in range(n)] for i in range(n)]
        A[0][2]=0;A[2][0]=1
        check(sum(binary_stats(A).values())==n-3,'sharp C3+singletons witness')
    n=7; residues={1,2,4}; A=[[int(i!=j and (j-i)%n in residues) for j in range(n)] for i in range(n)]
    e=binary_stats(A); k=sum(e.values())
    check(k==(35 if mutate=='paley' else 14),'Paley7 has 14, not all 35, odd foursets')
    print('Ordinal C3 witness n=4..20 PASS; Paley7 odd-fourset count',k,'of',len(e))
    # A single nonzero 4-face is NOT a cocycle on 5 vertices: explicit hostile.
    fake={s:int(s==(0,1,2,3)) for s in combinations(range(5),4)}
    check(sum(fake.values())%2==1,'arbitrary-support negative control is not cocycle')
    print('Hostile arbitrary singleton 4-support fails cocycle, so domain condition is essential.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mutate',choices=['hamiltonian','paley']);args=ap.parse_args()
    cospectral_audit(args.mutate);cocycle_audit(args.mutate)
    print('INDEPENDENT BRIDGE AUDIT PASS')
