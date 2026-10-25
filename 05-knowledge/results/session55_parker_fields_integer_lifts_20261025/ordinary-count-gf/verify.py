#!/usr/bin/env python3
"""Exact integer-polynomial replay; Python standard library only."""
from fractions import Fraction
from pathlib import Path
import math,json,hashlib,argparse
P=Path(__file__).resolve().parent
PHI={1:[1,-1],2:[1,1],3:[1,1,1],4:[1,0,1],5:[1,1,1,1,1],6:[1,-1,1]}
ORDERS={1:4,2:3,3:2,4:2,5:2,6:2}
def check(t,msg):
    if not t:raise RuntimeError(msg)
def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def add(a,b):
    c=[0]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return trim(c)
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)
def power(a,k):
    r=[1]
    for _ in range(k):r=mul(r,a)
    return r
def divrem(a,b):
    r=list(map(Fraction,a)); q=[Fraction(0)]*max(1,len(a)-len(b)+1)
    while len(r)>=len(b) and any(r):
        k=len(r)-len(b); x=r[-1]/b[-1];q[k]=x
        for i,y in enumerate(b):r[i+k]-=x*y
        trim(r)
    return trim(q),trim(r)
def C(n):
    g=lambda d:math.gcd(n,d)
    return 4*g(4)+4*g(6)+2*g(2)**2+2*g(5)+g(2)*g(4)-2-4*g(2)
def indicatorC(n):return 7+sum(a*(n%d==0) for d,a in ((2,13),(3,8),(4,12),(5,8),(6,8)))
def count(n):return n**3-(4*math.gcd(n,2)+4)*n*n+C(n)*n
# Each tuple is scalar, power j, substitution z^d.
TERMS=((1,3,1),(-8,2,1),(7,1,1),(-16,2,2),(26,1,2),(24,1,3),(48,1,4),(40,1,5),(48,1,6))
def run():
    D=[1]
    for d,k in ORDERS.items():D=mul(D,power(PHI[d],k))
    N=[0]
    for a,j,d in TERMS:
        base=[1]+[0]*(d-1)+[-1];den=power(base,j+1)
        num0={1:[0,1],2:[0,1,1],3:[0,1,4,1]}[j]
        num=[0]*(d*(len(num0)-1)+1)
        for k,c in enumerate(num0):num[k*d]=a*c
        q,r=divrem(D,den);check(not any(r),'denominator divides D')
        N=add(N,mul(num,q))
    N=[int(x) for x in N]
    Q=[27,37,119,201,346,438,639,649,729,635,576,400,300,160,94,34,14,2]
    check(N==[0]*9+[8*x for x in Q],'exact numerator polynomial')
    residues={}
    for d,phi in PHI.items():
        _,r=divrem(N,phi);check(any(r),'no cyclotomic cancellation')
        residues[d]=[str(x) for x in r]
    coeff=[0]+[count(n) for n in range(1,361)]
    check(all(C(n)==indicatorC(n) for n in range(1,61)),'indicator decomposition')
    product=mul(coeff,D)
    check(product[:361]==N+[0]*(361-len(N)),'coefficient convolution through n360')
    reconstructed=[]
    for n in range(361):
        value=sum(a*(n//d)**j for a,j,d in TERMS if n>0 and n%d==0)
        check(value==coeff[n],'Lambert summand coefficient')
        reconstructed.append(value)
    return {'status':'PASS','denominator_ascending':D,'numerator_ascending':N,'numerator_quotient_8z9_ascending':Q,'cyclotomic_pole_orders':ORDERS,'numerator_cyclotomic_remainders':residues,'coefficients_n0_to_60':coeff[:61],'series_verified_through':360,'minimal_period':60,'maximal_proper_period_counterexamples':{'C60':C(60),'C72':C(72),'C80':C(80),'C90':C(90)}}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();out=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if a.write:(P/'certificate.json').write_text(out)
    else:check((P/'certificate.json').read_text()==out,'persisted certificate')
    print('PASS: exact rational identity, noncancellation, indicator decomposition, coefficients through 360')
