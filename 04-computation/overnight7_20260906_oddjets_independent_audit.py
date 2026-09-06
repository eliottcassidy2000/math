#!/usr/bin/env python3
"""Independent odd-prime three-jet Smith audit, standard-library only.

Every residual minor is reconstructed by integer determinant interpolation,
not permutation expansion. Full Smith exponents use p-local pivot reduction.
No primary mathematical producer, SymPy or integer Smith routine is imported.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb,factorial
from pathlib import Path
import argparse,hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)
def vp(n,p):
    n=F(n)
    if not n:return float('inf')
    a,b=abs(n.numerator),n.denominator; v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v
def det_int(a):
    a=[row[:] for row in a];n=len(a);previous=1;sign=1
    for k in range(n-1):
        nonzero=next((i for i in range(k,n) if a[i][k]),None)
        if nonzero is None:return 0
        if nonzero!=k:a[k],a[nonzero]=a[nonzero],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=a[i][j]*pivot-a[i][k]*a[k][j]
                if value%previous:raise RuntimeError('Bareiss nonintegral division')
                a[i][j]=value//previous
        for i in range(k+1,n):a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]
def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
def interpolate(values):
    n=len(values);delta=[F(x) for x in values]
    answer=[F(0)]*n;falling=[F(1)]
    for k in range(n):
        scale=delta[0]/factorial(k)
        for j,c in enumerate(falling):answer[j]+=scale*c
        delta=[delta[j+1]-delta[j] for j in range(len(delta)-1)]
        falling=mul(falling,[-k,1])
    need(all(c.denominator==1 for c in answer),'integer polynomial from finite differences')
    return {j:int(c) for j,c in enumerate(answer) if c}
def poly_minor(I,J):
    near=[i for i in I if i>=3]
    degree=sum(sorted(J,reverse=True)[:len(near)])-sum(i%3 for i in near)
    def value(a):
        return det_int([[comb(j,i%3)*(a**(j-i%3) if i>=3 else 1) for j in J] for i in I])
    polynomial=interpolate([value(a) for a in range(degree+1)])
    need(sum(c*(-2)**j for j,c in polynomial.items())==value(-2),'independent polynomial identity control')
    return sum(J)-sum(i%3 for i in I),polynomial
def sub23(poly):
    out={}
    for j,c in poly.items():
        for k in range(j+1):out[k]=out.get(k,0)+c*comb(j,k)*2**(j-k)*3**k
    return {j:c for j,c in out.items() if c}
FRONTS={
0:{1:{(1,0,0)},2:{(3,1,1),(4,0,0)},3:{(7,1,1),(9,0,0)},4:{(12,4,4),(13,1,1)}},
3:{1:{(1,0,1),(3,0,0)},2:{(3,1,3),(4,0,1),(6,0,0)},3:{(7,1,2),(9,0,0)},4:{(12,4,6),(13,1,2)}},
5:{1:{(1,0,0)},2:{(3,1,1),(4,0,0)},3:{(7,1,2),(8,1,1),(9,0,0)},4:{(12,4,5),(13,1,1)}}}
def smith_local(nodes,p):
    a=[[F(comb(j,r)*x**(j-r) if j>=r else 0) for j in range(9)] for x in nodes for r in range(3)]
    answer=[]
    for k in range(9):
        v,i,j=min((vp(a[i][j],p),i,j) for i in range(k,9) for j in range(k,9))
        need(v!=float('inf') and v>=0,'nonnegative finite p-local pivot')
        a[k],a[i]=a[i],a[k]
        for row in a:row[k],row[j]=row[j],row[k]
        pivot=a[k][k];answer.append(v)
        for i in range(k+1,9):
            ratio=a[i][k]/pivot
            need(vp(ratio,p)>=0,'p-integral elimination multiplier')
            for j in range(k+1,9):a[i][j]-=ratio*a[k][j]
            a[i][k]=F(0)
        # Corresponding column eliminations clear the pivot row and do not
        # change the remaining submatrix because its pivot column is zero.
        for j in range(k+1,9):a[k][j]=F(0)
    need(answer==sorted(answer),'ordered p-local Smith exponents')
    return tuple(answer)
def predicted(p,e,d):
    delta=[0]*4
    if not d:
        if not e:return (0,)*9
        extra=(1,2,2,2) if p==3 else (0,0,1,1) if p==5 else (0,0,0,0)
        delta += [s*e+b for s,b in zip((1,3,7,12),extra)]
    elif p==3:
        delta += [min(e+1,3*e),min(3*e+d+2,4*e+1,6*e),
                  min(7*e+d+1,9*e),min(12*e+4*d+2,13*e+d+1)]
    else:
        delta += [e,min(3*e+d,4*e),min(7*e+d+(p==5),9*e),
                  min(12*e+4*d+(p==5),13*e+d)]
    largest=max(6*e,8*e-(p==3)) if not d else 8*e+5*d-(p==3)
    delta += [27*e+9*d-largest,27*e+9*d]
    return tuple(delta[i+1]-delta[i] for i in range(9))
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--primary-certificate',type=Path)
    args=parser.parse_args()
    records=[];bykey={}
    for r in range(1,5):
        for I in combinations(range(6),r):
            for J in combinations(range(3,9),r):
                W,poly=poly_minor(I,J)
                need(bool(poly),'nonzero symbolic residual minor')
                row={'rank':r,'rows':list(I),'degrees':list(J),'weight':W,
                     'polynomial':[[j,c] for j,c in poly.items()]}
                records.append(row);bykey[I,J]=(W,poly)
    need(len(records)==886 and sum(len(row['polynomial']) for row in records)==2623,'complete polynomial universe')
    for p in (0,3,5):
        for r in range(1,5):
            points={(row['weight'],j,j+(vp(c,p) if p else 0)) for row in records if row['rank']==r for j,c in row['polynomial']}
            front={a for a in points if not any(a!=b and all(y<=x for x,y in zip(a,b)) for b in points)}
            need(front==FRONTS[p][r],'all close-case coefficient frontiers')
        shifts=(1,2,2,2) if p==3 else (0,0,1,1) if p==5 else (0,0,0,0)
        termcount=0
        for row in records:
            r=row['rank'];w=row['weight'];s=(1,3,7,12)[r-1]
            poly=sub23(dict(row['polynomial'])) if p==3 else dict(row['polynomial'])
            for c in poly.values():
                need(w>=s and w+(vp(c,p) if p else 0)>=s+shifts[r-1],'all equilateral coefficient floors')
                termcount+=1
        print('COEFFICIENT_FLOOR',p or 'generic',termcount)
    # Identified low-weight factorizations, reconstructed by interpolation.
    witnesses=[((2,),(3,),[3]),((0,),(3,),[1]),((1,2),(3,4),[6]),
        ((0,1),(3,4),[1]),((2,5),(3,4),[0,-18,18]),
        ((0,1,2),(3,4,5),[1]),((1,2,5),(3,4,5),[0,30,-90,60]),
        ((2,4,5),(3,4,5),[0,0,0,0,60,-90,30]),
        ((0,2,5),(3,4,5),[0,12,-42,30]),
        ((0,1,2,5),(3,4,5,6),[0,-3,18,-30,15]),
        ((1,2,4,5),(3,4,5,6),[0,0,0,0,90,-360,540,-360,90])]
    for I,J,coeff in witnesses:need(bykey[I,J][1]=={j:c for j,c in enumerate(coeff) if c},'exact attaining witness')
    # Convex domination uses d>=1, not integrality of e,d.
    for e,d in ((F(1,2),F(1)),(F(5,4),F(3,2)),(F(7,3),F(1))):
        need(8*e+d >= ((9*e)+(7*e+d+1))/2,'real-cone quinary redundancy')
    e,d=F(1,2),F(1,4)
    need(8*e+d<min(9*e,7*e+d+1),'d>=1 scope hostile')
    literal=[]
    for p in (3,5,7,11,13):
        for e,d in ((0,1),(1,1),(2,1),(4,1),(1,3),(3,2)):
            for u in (1,-1,p+2):
                nodes=(0,p**e,p**(e+d)*u)
                actual=smith_local(nodes,p)
                need(actual==predicted(p,e,d),'fresh close full p-local Smith control')
                literal.append([p,e,d,u,list(actual)])
        for a in range(2,p):
            for e in (0,1,3):
                actual=smith_local((0,p**e,p**e*a),p)
                need(actual==predicted(p,e,0),'complete residue equilateral full Smith control')
                literal.append([p,e,0,a,list(actual)])
    need(smith_local((0,4,8),2)[-1]==18 and smith_local((0,4,12),2)[-1]==19,'dyadic outside-scope hostile')
    if args.primary_certificate:
        primary=json.loads(args.primary_certificate.read_text())
        # The primary wrapper may include the certificate beside other metadata.
        rows=primary if isinstance(primary,list) else primary['minors']
        need(rows==records,'all primary minor polynomials independently reconstructed')
    print('FULL_LITERAL_SMITH_CONTROLS',len(literal),'plus two dyadic hostiles')
    print('POLYNOMIAL_SHA256',hashlib.sha256(json.dumps(records,sort_keys=True,separators=(',',':')).encode()).hexdigest())
    print('LITERAL_SHA256',hashlib.sha256(json.dumps(literal,separators=(',',':')).encode()).hexdigest())
    print('PASS',CHECKS,'explicit gates; no primary imports, SymPy or integer Smith routines')
if __name__=='__main__':main()
