"""Exact enumeration formulas and independent cyclic-composition controls.

Burnside's general tournament formula is classical, not a new algorithm.
The sampling-family generating functions follow from SAMPLING.md.
"""
from pathlib import Path
from fractions import Fraction
from itertools import combinations
from math import gcd,factorial,comb
from collections import Counter
import json,sympy as s
HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def odd_partitions(n,lo=1):
    if n==0:yield ()
    for a in range(lo,n+1,2):
        for tail in odd_partitions(n-a,a):yield (a,)+tail
def tournaments(n):
    value=Fraction(0)
    for part in odd_partitions(n):
        e=sum((a-1)//2 for a in part)+sum(gcd(part[i],part[j]) for i in range(len(part)) for j in range(i))
        z=1
        for a,m in Counter(part).items():z*=a**m*factorial(m)
        value+=Fraction(2**e,z)
    need(value.denominator==1,'Burnside integrality')
    return value.numerator
def L(n,k):
    if n<k:return 0
    return (comb(n-1,k-1)+(k-1)*(n%k==0))//k
def literal(n,k):
    reps=set()
    for cuts in combinations(range(1,n),k-1):
        ends=(0,)+cuts+(n,);v=tuple(ends[i+1]-ends[i] for i in range(k))
        reps.add(min(v[j:]+v[:j] for j in range(k)))
    return len(reps)
def run():
    z=s.Symbol('z')
    F3=z/(1-z)+(z**3/(1-z)**3+2*z**3/(1-z**3))/3
    F5=F3+(z**5/(1-z)**5+4*z**5/(1-z**5))/5
    series3=s.series(F3,z,0,65).removeO().expand()
    series5=s.series(F5,z,0,65).removeO().expand()
    rows=[]
    for n in range(1,65):
        a=1+L(n,3);b=a+L(n,5)
        need(a==series3.coeff(z,n) and b==series5.coeff(z,n),'rational GF coefficients')
        if n<=18:
            need(L(n,3)==literal(n,3) and L(n,5)==literal(n,5),'literal cyclic compositions')
        rows.append({'n':n,'parabola_classes':a,'power_curve_m_ge3_classes':b})
    old=json.loads((HERE/'tournament_receipt.json').read_text())
    for lev in old['levels']:
        r=rows[lev['n']-1]
        need(r['parabola_classes']==lev['parabola_types'],'independent sign-word parabola census')
        need(r['power_curve_m_ge3_classes']==lev['cubic_types'],'independent sign-word cubic census')
    burn=[tournaments(n) for n in range(31)]
    need(burn[:9]==[1,1,1,2,4,12,56,456,6880],'classical small tournament controls')
    out={'scope':'Classical Burnside formula n0..30; proved sampling-family formulas n1..64; independent cyclic-composition enumeration through18 and original sign-word census through14. No new OEIS terms or general complexity improvement claimed.','A000568_existing_formula':burn,'parabola_GF':str(s.factor(F3)),'power_curve_GF':str(s.factor(F5)),'sampling_counts':rows}
    tag='.optimized' if not __debug__ else ''
    (HERE/('aspect_sequences'+tag+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'burnside_n30':burn[-1],'sampling_n18':rows[17],'sampling_n64':rows[-1]}))
if __name__=='__main__':run()
