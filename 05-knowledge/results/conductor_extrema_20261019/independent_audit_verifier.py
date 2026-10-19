#!/usr/bin/env python3
"""Independent modular-rank certificates and exact small controls.
No producer code is imported; no assert gates or floating-point arithmetic.
A maximal modular rank is an exact rational-rank lower bound.
The matching upper bounds are the proved dimension/parametrization bounds.
"""
from math import comb, isqrt
from itertools import combinations
import json
import sympy as S


def require(x,label):
    if not x: raise RuntimeError(label)


PRIMES=(1000003,1000033)
for p in PRIMES:
    require(all(p%d for d in range(2,isqrt(p)+1)), 'rank modulus must be prime')


def rankmod(A,p):
    A=[[x%p for x in row] for row in A];r=0
    for j in range(len(A[0])):
        k=next((k for k in range(r,len(A)) if A[k][j]),None)
        if k is None: continue
        A[k],A[r]=A[r],A[k]
        inv=pow(A[r][j],-1,p);A[r]=[x*inv%p for x in A[r]]
        for k in range(r+1,len(A)):
            c=A[k][j]
            if c:A[k]=[(x-c*y)%p for x,y in zip(A[k],A[r])]
        r+=1
        if r==len(A):break
    return r


def evaluation(nodes,power,q):
    exponents=sorted({a+power*b for a in range(q+1) for b in range(q-a+1)})
    return [[t**e for e in exponents] for t in nodes]


def certificate(A,target):
    for p in PRIMES:
        r=rankmod(A,p)
        require(r<=target,'rank exceeds claimed upper bound')
        if r==target:return p
    raise RuntimeError('No full-rank modular certificate in the stated two-prime budget')


def thresholds(r):
    s=0
    while comb(s+2,2)<r:s+=1
    return s,r//2


def bounds(r):
    s,c=thresholds(r)
    return (s-1)*r-comb(s+2,3)+1,(r-2)**2//4


def main():
    for r in range(3,10001):
        s,c=thresholds(r);lo,hi=bounds(r)
        require(lo==sum(r-comb(q+2,2) for q in range(1,s)),'minimum sum')
        require(hi==sum(r-2*q-1 for q in range(1,c)),'maximum sum')
        require(lo<=hi,'bound order')
    matrix_count=deleted_count=0;prime_use={p:0 for p in PRIMES}
    for r in range(3,31):
        s,c=thresholds(r);nodes=list(range(1,r+1))
        for power,d in ((r,s),(2,c)):
            for q in range(c+2):
                target=min(r,comb(q+2,2)) if power==r else min(r,2*q+1)
                A=evaluation(nodes,power,q)
                prime_use[certificate(A,target)]+=1;matrix_count+=1
            before=evaluation(nodes,power,d-1)
            target=min(r,comb(d+1,2)) if power==r else min(r,2*d-1)
            for i in range(r):
                prime_use[certificate(before[:i]+before[i+1:],target)]+=1;deleted_count+=1
    # Full exact-rational evaluation rather than modular certificates for the
    # equal-Hilbert/different-conductor hostile and small mixed controls.
    def general_eval(points,q):
        return S.Matrix([[x**a*y**b*z**(q-a-b) for a in range(q+1) for b in range(q-a+1)]
                         for x,y,z in points])
    def ranks_and_separators(points):
        r=len(points);hs=[];ds=[None]*r
        for q in range(r//2+1):
            A=general_eval(points,q);h=A.rank();hs.append(h)
            for i in range(r):
                if ds[i] is None:
                    rows=[j for j in range(r) if j!=i]
                    if h-A.extract(rows,range(A.cols)).rank()==1:ds[i]=q
        return hs,ds
    D7=[(1,t,t*t) for t in range(6)]+[(1,6,37)]
    G7=[(1,t,t**7) for t in range(1,8)]
    controls=[]
    for name,points,expected in [('D7',D7,[3]*6+[2]),('G7',G7,[3]*7)]:
        require(all(S.Matrix(triple).det()!=0 for triple in combinations(points,3)),'hostile matroid')
        h,d=ranks_and_separators(points)
        require(h==[1,3,6,7] and d==expected,'same Hilbert/different conductors')
        controls.append({'name':name,'ranks':h,'conductors':d})
    # A third, non-endpoint family prevents testing only the two attainers.
    intermediate=[]
    for r in (6,7,10,13):
        points=[(1,t,t**3) for t in range(1,r+1)]
        require(all(S.Matrix(triple).det()!=0 for triple in combinations(points,3)),'cubic control matroid')
        hs=[general_eval(points,q).rank() for q in range(r//2+1)]
        require(all(h==min(r,3*q) for q,h in enumerate(hs) if q),'cubic Hilbert control')
        delta=sum(r-h for h in hs[1:]);lo,hi=bounds(r)
        require(lo<=delta<=hi,'intermediate bounds')
        require((delta==hi)==(hs[2]==5),'upper conic equality control')
        intermediate.append({'r':r,'ranks':hs,'defect':delta,'bounds':[lo,hi]})
    print(json.dumps({'status':'PASS','formula_r_range':[3,10000],
                      'witness_r_range':[3,30],'full_modular_rank_certificates':matrix_count,
                      'row_deletion_modular_certificates':deleted_count,'prime_use':prime_use,
                      'same_Hilbert_hostile':controls,'third_family_controls':intermediate,
                      'scope':'Exact finite algebra only; unrestricted proof reviewed separately.'},indent=2,sort_keys=True))


if __name__=='__main__':main()
