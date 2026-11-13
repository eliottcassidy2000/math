#!/usr/bin/env python3
"""Small independent exact controls for O1-O3 and C1-C3; no producer imports."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PINS = {
    'ROOT_OBSERVERS.md': 'bbe31d87ba4a7d328add13194af9dde80f68a2d6d03bcadb5de3e47bce4063c7',
    'CLIFFORD_BOUNDARY.md': '6615b63e6e7c1bb692facbe99ade8a925b96cebb963ecdd8305d56f13ca7fa10',
}
COUNT = 0


def check(ok, label):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(label)


def matrix(a): return [[Q(x) for x in row] for row in a]
def eye(n): return matrix([[int(i == j) for j in range(n)] for i in range(n)])
def transpose(a): return list(map(list, zip(*a)))
def mul(a,b): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def add(a,b,scale=1): return [[x+scale*y for x,y in zip(r,s)] for r,s in zip(a,b)]
def zero(a): return all(not x for r in a for x in r)


def path_controls():
    out=[]
    for m in range(2,15):
        n=m+1
        degree=[1]+[2]*(m-1)+[1]
        P=matrix([[Q(int(abs(i-j)==1),degree[i]) for j in range(n)] for i in range(n)])
        D=matrix([[degree[i]*int(i==j) for j in range(n)] for i in range(n)])
        L=matrix([[int(i==0),int(i>0)] for i in range(n)])
        M=matrix([[1,0],[0,2*m-1]])
        E=matrix([[int(i==0) for i in range(n)],[Q(degree[i],2*m-1) if i else 0 for i in range(n)]])
        R=add(eye(n),mul(L,E),-1)
        B=mul(mul(R,P),L)
        gram=mul(mul(transpose(B),D),B)
        EP2L=mul(mul(mul(E,P),P),L)
        EPL=mul(mul(E,P),L)
        defect=mul(M,add(EP2L,mul(EPL,EPL),-1))
        check(defect==gram and not zero(gram),'reversible Gram identity and positive hostile')
        G=add(P,eye(n),-1)
        first=mul(mul(E,G),L)
        second=mul(mul(mul(E,G),G),L)
        check(mul(M,add(second,mul(first,first),-1))==gram,'derivative normalization')
        half_second=[[x/2 for x in r] for r in second]
        check(mul(M,add(half_second,mul(first,first),-1))!=gram,'Taylor coefficient is not second derivative')
        check(EP2L[0][0]==Q(1,2) and mul(EPL,EPL)[0][0]==Q(1,2*m-1),'return probability discrepancy')
        colors=[0]+[1]*m
        counts=[2]
        while len(set(colors))<n:
            blocks=sorted(set(colors))
            signatures=[(colors[i],tuple(sum(P[i][j] for j in range(n) if colors[j]==c) for c in blocks)) for i in range(n)]
            keys=sorted(set(signatures))
            colors=[keys.index(s) for s in signatures]
            counts.append(len(keys))
        check(counts==list(range(2,m+2)),'exact m-minus-one strict rounds')
        # Every equality class after k rounds is the stated singleton prefix plus tail.
        out.append({'m':m,'strict_rounds':len(counts)-1,'actual_return':'1/2','quotient_return':str(Q(1,2*m-1))})
    return out


def partitions(n):
    if n==0:
        yield []
        return
    for p in partitions(n-1):
        for k in range(max(p,default=-1)+2):
            yield p+[k]


def algebra_controls():
    tested=0
    for p in partitions(3):
        labels=sorted(set(p))
        L=matrix([[int(p[i]==c) for c in labels] for i in range(3)])
        # Lagrange products reconstruct all class indicators from one class-separating function.
        for c in labels:
            values=[]
            for x in p:
                z=Q(1)
                for d in labels:
                    if d!=c: z*=Q(x-d,c-d)
                values.append(z)
            check(values==[int(x==c) for x in p],'finite algebra indicator interpolation')
        for bits in product((0,1),repeat=9):
            A=matrix([bits[3*i:3*i+3] for i in range(3)])
            AL=mul(A,L)
            B=[AL[p.index(c)] for c in labels]
            stable=all(AL[i]==AL[j] for i in range(3) for j in range(3) if p[i]==p[j])
            check((AL==mul(L,B))==stable,'row sums iff invariant block algebra')
            if stable:
                check(mul(A,AL)==mul(L,mul(B,B)),'word intertwining')
            tested+=1
    # Original insufficient hypothesis: span(1,f) contains a coarse partition's indicators.
    f=[0,1,2]
    check([x*x for x in f] != [0,1,2], 'square not in affine span: endpoints force coefficients 0,1')
    check([f[0]**2,f[1]**2,f[2]**2]==[0,1,4], 'three distinct levels produce discrete algebra')
    return {'binary_3_by_3_matrix_partition_pairs':tested,
            'repaired_hypothesis':'V0 is exactly the block-constant algebra of Pi0'}


def clifford_controls():
    balances=[]
    for p,q in product(range(5),repeat=2):
        n=p+q
        if not n: continue
        R=matrix([[int(i==j)*(1 if i<p else -1) for j in range(n)] for i in range(n)])
        J=matrix([[0]*n for _ in range(n)])
        for i in range(min(p,q)):
            J[p+i][i]=1
            J[i][p+i]=-1
        check(zero(add(mul(J,R),mul(R,J))),'odd block map anticommutes')
        check((mul(J,J)==[[-x for x in row] for row in eye(n)])==(p==q),'negative-square extension iff balanced')
        balances.append([p,q,p==q])
    minima=[]
    for s in range(9):
        feasible=[(a+b,a,b) for a,b in product(range(18),repeat=2) if b-a==s]
        check(min(feasible)==(s,0,s),'least stabilization is s negative dimensions')
        minima.append(s)
    # Balanced, incidence-compatible positive control with an explicit orthogonal intertwiner.
    R=matrix([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]])
    B=matrix([[0,1],[-1,0]])
    A=matrix([[1,1],[1,2]])
    D=matrix([[2,-1],[-1,1]])
    check(mul(B,A)==mul(D,B),'explicit sector intertwiner')
    J=matrix([[0,0,0,1],[0,0,-1,0],[0,1,0,0],[-1,0,0,0]])
    C=matrix([[1,1,0,0],[1,2,0,0],[0,0,2,-1],[0,0,-1,1]])
    check(zero(add(mul(J,R),mul(R,J))) and mul(J,J)==[[-x for x in row] for row in eye(4)],'Cl11 relation positive control')
    check(mul(J,C)==mul(C,J),'incidence-compatible extension')
    for a,b in product(range(-3,4),repeat=2):
        J=matrix([[0,a],[b,0]])
        C=matrix([[1,0],[0,2]])
        check((mul(J,C)==mul(C,J))==(a==b==0),'balanced spectral hostile')
    # Original class basis: enumerate all nonnegative binary 3x3 matrices.
    R=matrix([[0,1,0],[1,0,0],[0,0,1]])
    survivors=0
    for bits in product((0,1),repeat=9):
        J=matrix([bits[3*i:3*i+3] for i in range(3)])
        if zero(add(mul(J,R),mul(R,J))):
            check(zero(J),'nonnegative anticommuter vanishes')
            survivors+=1
    check(survivors==1,'only zero nonnegative binary matrix')
    return {'sector_dimensions':balances,'minimal_additions':minima,
            'nonnegative_matrices_tested':512,'nonnegative_survivors':survivors,
            'scope':'finite linear controls, not extra tournament classes or a Bott recurrence'}


def main():
    actual={p:hashlib.sha256((HERE.parent/p).read_bytes()).hexdigest() for p in PINS}
    check(actual==PINS,'root reviewed bytes unchanged')
    result={'status':'FINITE-EXACT independent controls; root proof status unchanged',
            'reviewed_sha256':actual,'algebra':algebra_controls(),'paths':path_controls(),'clifford':clifford_controls()}
    result['checks']=COUNT
    result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(result,sort_keys=True,indent=2))


if __name__=='__main__': main()
