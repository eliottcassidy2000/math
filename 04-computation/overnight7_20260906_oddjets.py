#!/usr/bin/env python3
"""Complete odd-prime three-node three-Hasse-jet Smith certificate.

Proof universe: every residual rank1..4 minor (886 total), expanded by
determinant permutations. Close-pair lower bounds keep both depths.
Equilateral lower bounds retain the p=3 residue substitution explicitly.
The d>=1 convex simplification at p=5 is separated from monomial dominance.
No repository imports; normal/-O gates remain active. A JSON polynomial
manifest is written beside this source, or at --certificate PATH.
"""
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations,permutations,product
from math import comb
from pathlib import Path
import argparse
import json
import sys
from sympy import Matrix,ZZ
from sympy.matrices.normalforms import smith_normal_form

sys.stdout.reconfigure(newline="\n")
GATES=0
ROWS=tuple((node,r) for node in (0,1) for r in range(3))


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def vp(value,p):
    value=Q(value)
    if not value:
        return float('inf')
    a,b=abs(value.numerator),value.denominator
    out=0
    while a%p==0:
        a//=p;out+=1
    while b%p==0:
        b//=p;out-=1
    return out


def symbolic_minor(I,J):
    coefficients=defaultdict(int)
    r=len(I)
    for sigma in permutations(range(r)):
        coefficient=(-1)**sum(sigma[i]>sigma[j] for i in range(r) for j in range(i+1,r))
        power=0
        for i,j in enumerate(sigma):
            node,order=ROWS[I[i]]
            degree=J[j]
            coefficient*=comb(degree,order)
            power+=node*(degree-order)
        coefficients[power]+=coefficient
    W=sum(J)-sum(ROWS[i][1] for i in I)
    return W,{b:c for b,c in sorted(coefficients.items()) if c}


def pmul(a,b):
    out=defaultdict(int)
    for i,x in a.items():
        for j,y in b.items():
            out[i+j]+=x*y
    return {i:c for i,c in sorted(out.items()) if c}


def pprod(*factors):
    out={0:1}
    for factor in factors:
        out=pmul(out,factor)
    return out


def substitute(poly,center,scale):
    out=defaultdict(int)
    for b,c in poly.items():
        for r in range(b+1):
            out[r]+=c*comb(b,r)*center**(b-r)*scale**r
    return {b:c for b,c in sorted(out.items()) if c}


def symbolic_universe():
    records=[]
    counts=[]
    for r in range(1,5):
        count=0
        for I in combinations(range(6),r):
            for J in combinations(range(3,9),r):
                W,poly=symbolic_minor(I,J)
                need(bool(poly),('nonzero residual minor',I,J))
                # A separate literal integer determinant catches sign/index errors.
                at=-2
                literal=Matrix([[comb(j,ROWS[i][1])*at**(ROWS[i][0]*(j-ROWS[i][1])) for j in J] for i in I]).det()
                need(sum(c*at**b for b,c in poly.items())==literal,('literal minor control',I,J))
                records.append({'rank':r,'rows':list(I),'degrees':list(J),'weight':W,
                                'polynomial':[[b,c] for b,c in poly.items()]})
                count+=1
        counts.append(count)
    need(counts==[36,225,400,225] and len(records)==886,'complete residual minor universe')
    print('SYMBOLIC_UNIVERSE ranks1..4',counts,'total886 monomials',sum(len(r['polynomial']) for r in records))
    return records


FRONTS={
    'generic':{1:{(1,0,0)},2:{(3,1,1),(4,0,0)},3:{(7,1,1),(9,0,0)},4:{(12,4,4),(13,1,1)}},
    'p3':{1:{(1,0,1),(3,0,0)},2:{(3,1,3),(4,0,1),(6,0,0)},3:{(7,1,2),(9,0,0)},4:{(12,4,6),(13,1,2)}},
    'p5':{1:{(1,0,0)},2:{(3,1,1),(4,0,0)},3:{(7,1,2),(8,1,1),(9,0,0)},4:{(12,4,5),(13,1,1)}}}


def lower_certificates(records):
    # generic means coefficient valuation bounded below by zero, not one sampled prime.
    for name,p in (('generic',None),('p3',3),('p5',5)):
        for r in range(1,5):
            costs={(row['weight'],b,b+(vp(c,p) if p else 0))
                   for row in records if row['rank']==r for b,c in row['polynomial']}
            front={cost for cost in costs if not any(other!=cost and all(x<=y for x,y in zip(other,cost)) for other in costs)}
            need(front==FRONTS[name][r],('complete close-pair Pareto set',name,r,front))
            for cost in costs:
                need(any(all(x<=y for x,y in zip(witness,cost)) for witness in front),
                     ('every close monomial dominated on e,d-1>=0',name,r,cost))
            print('CLOSE_FRONT',name,'rank',r,sorted(front))
    # e>=1 and a,1-a units; at p3 the entire allowed residue class is a=2+3z.
    minW={1:1,2:3,3:7,4:12}
    beta={'generic':(0,0,0,0),'p3':(1,2,2,2),'p5':(0,0,1,1)}
    for name,p in (('generic',None),('p3',3),('p5',5)):
        terms=0
        for row in records:
            r=row['rank'];W=row['weight'];s=minW[r];b0=beta[name][r-1]
            poly=dict(row['polynomial'])
            if p==3:
                poly=substitute(poly,2,3)
            for power,c in poly.items():
                cost=vp(c,p) if p else 0
                need(W>=s and W+cost>=s+b0,
                     ('equilateral unbounded e=1+E coefficient inequality',name,row['rows'],row['degrees'],power,c))
                terms+=1
        print('EQUILATERAL_LOWER',name,'e>=1 slopes',(1,3,7,12),'intercepts',beta[name],'all_terms',terms)


def attaining_identities():
    A={1:1};Am1={1:1,0:-1}
    witnesses=[
      ('unit-entry',(0,),(3,),{0:1}),
      ('three-entry',(2,),(3,),{0:3}),
      ('unit-rank2',(0,1),(3,4),{0:1}),
      ('six-rank2',(1,2),(3,4),{0:6}),
      ('matching-rank2',(2,5),(3,4),pprod({0:18},A,Am1)),
      ('unit-rank3',(0,1,2),(3,4,5),{0:1}),
      ('first-rank3',(1,2,5),(3,4,5),pprod({0:30},A,Am1,{1:2,0:-1})),
      ('second-rank3',(2,4,5),(3,4,5),pprod({0:30},{4:1},Am1,{1:1,0:-2})),
      ('p5-redundant-rank3',(0,2,5),(3,4,5),pprod({0:6},A,Am1,{1:5,0:-2})),
      ('linear-rank4',(0,1,2,5),(3,4,5,6),pprod({0:3},A,Am1,{2:5,1:-5,0:1})),
      ('quartic-rank4',(1,2,4,5),(3,4,5,6),{4:90,5:-360,6:540,7:-360,8:90})]
    out=[]
    for name,I,J,factored in witnesses:
        W,actual=symbolic_minor(I,J)
        need(actual==factored,('exact attaining factorization',name))
        out.append({'name':name,'rows':list(I),'degrees':list(J),'weight':W,'polynomial':[[b,c] for b,c in actual.items()]})
        print('ATTAIN',name,'rows',I,'degrees',J,'weight',W,'polynomial',actual)
    # The two equilateral rank3 linear factors cannot simultaneously gain extra valuation.
    need({1:2,0:-1}=={1:2,0:2*(-2)+3},'(2a-1)-2(a-2)=3 polynomial identity')
    for a in range(2,29,3):
        need(min(vp(2*a-1,3),vp(a-2,3))==1,'ternary complete residue rank3 attainment control')
    for p in (5,7,11,13):
        for a in range(2,p):
            need(min(vp(2*a-1,p),vp(a-2,p))==0,'odd nonternary equilateral rank3 control')
    return out


def largest(p,e,d=0):
    if p==2:
        raise ValueError('odd-prime theorem; dyadic unit sidecar is different')
    if d:
        return 8*e+5*d-(p==3)
    return max(6*e,8*e-(p==3))


def predicted(p,e,d):
    if p==2:
        raise ValueError('odd prime required')
    if d==0:
        if e==0:
            return (0,)*9
        if p==3:
            return (0,0,0,e+1,2*e+1,4*e,5*e,7*e-1,8*e-1)
        if p==5:
            return (0,0,0,e,2*e,4*e+1,5*e,7*e-1,8*e)
        return (0,0,0,e,2*e,4*e,5*e,7*e,8*e)
    if p==3:
        D4=min(3*e,e+1)
        D5=min(6*e,4*e+1,3*e+d+2)
        D6=min(9*e,7*e+d+1)
        D7=min(13*e+d+1,12*e+4*d+2)
    else:
        D4=e
        D5=min(4*e,3*e+d)
        D6=min(9*e,7*e+d+(p==5))
        D7=min(13*e+d,12*e+4*d+(p==5))
    D=[0,0,0,0,D4,D5,D6,D7,19*e+4*d+(p==3),27*e+9*d]
    return tuple(D[j+1]-D[j] for j in range(9))


def literal_smith(nodes,p):
    H=Matrix([[comb(j,r)*x**(j-r) if j>=r else 0 for j in range(9)] for x in nodes for r in range(3)])
    S=smith_normal_form(H,domain=ZZ)
    return tuple(vp(S[i,i],p) for i in range(9))


def precision_identities():
    N0={0:2,1:3,2:2}
    N1={0:7,1:-7,2:2}
    Nt={0:2,1:-7,2:7}
    diff={b:N0.get(b,0)-N1.get(b,0) for b in set(N0)|set(N1)}
    need({b:c for b,c in diff.items() if c}=={0:-5,1:10},'quadratic common-zero elimination')
    need(sum(Q(c)*Q(1,2)**b for b,c in N0.items())==4,'common root at half excluded at odd nonfive primes')
    need(3 not in {r*r%5 for r in range(5)},'quinary quadratic discriminant nonsquare')
    need(all(sum(c*r**b for b,c in N0.items())%5 for r in range(5)),'complete p5 unit numerator check')
    need(all(sum(c*2**b for b,c in P.items())%3 for P in (N0,N1,Nt)),'ternary equilateral numerator unit control')
    # Coordinatewise Pareto does not remove this cost, but convex dominance does.
    for e,d in product(range(10),range(1,10)):
        need(min(9*e,8*e+d,7*e+d+1)==min(9*e,7*e+d+1),'integer-depth p5 rank3 redundancy control')
    e=Q(1,2);d=Q(1,4)
    need(min(9*e,8*e+d,7*e+d+1)!=min(9*e,7*e+d+1),
         'dropping d>=1 invalidates the frontier simplification')
    need(all((Q(a)+Q(b))/2<=c for a,b,c in zip((9,0,0),(7,1,2),(8,1,1))),
         'symbolic convex domination on e,d-1>=0')


def literal_controls():
    count=0
    for p,e,d in product((3,5,7,11),(0,1,2,4,7),(1,2,4)):
        for unit in (1,1+p,-1):
            nodes=(0,p**e,p**(e+d)*unit)
            actual=literal_smith(nodes,p)
            expected=predicted(p,e,d)
            need(actual==expected,('full close-pair Smith',p,e,d,unit,actual,expected))
            need(actual[-1]==largest(p,e,d) and tuple(sorted(actual))==actual,'close precision and monotone factors')
            count+=1
    close_count=count
    for p in (3,5,7,11):
        stop=p*p if p<=7 else p
        for a in range(2,stop):
            if a%p in (0,1):
                continue
            for e in (0,1,2,5):
                nodes=(0,p**e,p**e*a)
                actual=literal_smith(nodes,p)
                expected=predicted(p,e,0)
                need(actual==expected,('full equilateral Smith',p,e,a,actual,expected))
                need(actual[-1]==largest(p,e,0) and tuple(sorted(actual))==actual,'equilateral precision and monotone factors')
                count+=1
    for p,e,d in ((3,2,1),(5,1,3),(7,2,0),(11,1,1)):
        a=2 if d==0 else p**d
        moved=tuple(-5-2*x for x in (0,p**e,p**e*a))
        need(literal_smith(moved,p)==predicted(p,e,d),('signed unit-affine normalization',p,e,d))
        count+=1
    need(literal_smith((0,4,8),2)[-1]==18 and literal_smith((0,4,12),2)[-1]==19,
         'dyadic hostile preserved outside odd-prime scope')
    try:
        predicted(2,2,1)
    except ValueError:
        need(True,'odd-prime API rejects dyadic specialization')
    else:
        need(False,'dyadic theorem accidentally extended')
    print('LITERAL_CONTROLS close',close_count,'equilateral',count-close_count-4,'affine4 total',count,'plus dyadic18/19 hostile')
    print('ODD_PRIME_LARGEST close=8e+5d-[p3]; equilateral=max(6e,8e-[p3])')
    print('EQUILATERAL_E1 p3',predicted(3,1,0),'p5',predicted(5,1,0),'p>=7',predicted(7,1,0))
    print('SHALLOW_SLOPE_HOSTILE p7 e4 d1 D7',sum(predicted(7,4,1)[:7]),'versus13e+d=53')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args()
    records=symbolic_universe()
    lower_certificates(records)
    witnesses=attaining_identities()
    precision_identities()
    literal_controls()
    certificate={'scope':'all886 residual rank1..4 integer polynomial identities; odd-prime close/equilateral laws',
                 'rows':ROWS,'minors':records,'witnesses':witnesses,
                 'close_frontiers':{name:{str(r):sorted(front) for r,front in ranks.items()} for name,ranks in FRONTS.items()}}
    raw=(json.dumps(certificate,sort_keys=True,indent=2)+'\n').encode('utf-8')
    args.certificate.write_bytes(raw)
    print('POLYNOMIAL_CERTIFICATE_SHA256',sha256(raw).hexdigest())
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)


if __name__=='__main__':
    main()
