"""Exact phase symmetry, drift sum, and an H25/H3 eleven-clock family."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import subprocess
import types

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def g(n):
    result=1
    while n>=2:
        q=n//2
        if n%2==0:result+=q+(q&1)
        n=q
    return result


def order2(p):
    k=1
    while pow(2,k,p)!=1:k+=1
    return k


def phase_data(p):
    k=order2(p);Q=2**k-1;h=Q//p;data={}
    for i in range(1,p):
        r=h*i;rev=int(f'{r:0{k}b}'[::-1],2)
        a=Fraction(p*(Q-rev),Q);d=int(2*i<p)
        drift=p*(g(r)-1+d*(r&1))-a*r
        if drift.denominator!=1:raise RuntimeError(('nonintegral drift',p,i))
        data[i]=dict(r=r,a=a,d=d,drift=int(drift))
    for i,z in data.items():
        doubled=data[2*i%p];opposite=data[p-i]
        if doubled['a']!=(z['a']+p*z['d'])/2 or opposite['a']!=p-z['a']:
            raise RuntimeError(('phase symmetry',p,i))
        if doubled['drift']!=z['drift'] or opposite['drift']!=z['drift']:
            raise RuntimeError(('drift symmetry',p,i))
    chi=1 if p%4==1 else -1
    if sum(z['drift'] for z in data.values())!=Fraction(k*(p-chi)**2,8):
        raise RuntimeError(('drift sum',p))
    return k,data


def M(B):return 25*g(B)+122*B+16*(B&1)+48
def N(C):return 3*g(C)+12*C+10
def E(n):return 25*g(n)-3*g(9*n)+3*n+25*(n&1)+3206


def main():
    catalog=[]
    for p in range(3,44,2):
        if any(p%d==0 for d in range(2,int(p**.5)+1)):continue
        k,data=phase_data(p)
        denominator=1
        from math import lcm
        for z in data.values():denominator=lcm(denominator,z['a'].denominator)
        catalog.append(dict(prime=p,period=k,integral_coefficients=denominator==1,
                            coefficient_denominator=denominator,a1=str(data[1]['a']),
                            drift_values=sorted({z['drift'] for z in data.values()})))
    for m in range(1,9):
        n=2**m+1;k,data=phase_data(n)
        if k!=2*m:raise RuntimeError(('Fermat period',m))
        for i,z in data.items():
            expected=1+int(f'{i-1:0{m}b}'[::-1],2)
            if z['a']!=expected:raise RuntimeError(('Fermat coefficient',m,i))
        if data[1]['drift']!=2:raise RuntimeError(('Fermat first drift',m))
        if m>=3 and data[3]['drift']!=6*2**(m-2)+6:raise RuntimeError(('Fermat split',m))
    for k in range(2,9):
        n=2**k-1;period,data=phase_data(n)
        if period!=k or data[1]['drift']!=2**(k-1):raise RuntimeError(('Mersenne period',k))
        if k>=3 and data[3]['drift']!=2**(k-2)+2:raise RuntimeError(('Mersenne split',k))
    period,eleven=phase_data(11)
    for n in range(128):
        if E(n)!=M(64*n+31)-N(576*n+46):raise RuntimeError(('physical eleven defect',n))
        if E(1024*n+93)!=E(n)+36:raise RuntimeError(('eleven loop',n))
    families=[]
    for k in(7,8,9):
        q=50*2**k-7;m=125*2**(k-7)-88
        if E(q)!=-36*m:raise RuntimeError(('eleven seed',k))
        n=1024**m*q+93*(1024**m-1)//1023
        B,C=64*n+31,576*n+46
        if E(n)!=0 or M(B)!=N(C):raise RuntimeError(('eleven family',k))
        families.append(dict(k=k,q=q,loops=m,n_bits=n.bit_length(),H_bits=M(B).bit_length()))
    path='05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
    source=subprocess.check_output(['git','show','HEAD:'+path],cwd=ROOT)
    module=types.ModuleType('inherited_compiler');exec(compile(source,path,'exec'),module.__dict__)
    literal=[]
    for n in(0,1):
        B,C=64*n+31,576*n+46
        charts=[([50,8,3,5,46,14],38,B,M(B)),([50,40,35,5,46,8],35,B,M(B)+2),
                ([2,4,1],3,2*C+1,N(C)),([2,4,1],1,4*C+3,N(C)+2)]
        saved=[]
        for core,cut,address,expected in charts:
            rows=module.construct(core,cut,0,address)
            if module.hp_subsets(rows)[-1]!=expected:raise RuntimeError(('literal',n,cut))
            saved.append(dict(core=core,cut=cut,address=address,vertices=len(rows),H=expected))
        literal.append(dict(n=n,charts=saved))
    report=dict(status='PROVED phase criteria and named families; finite catalog is not exhaustive prime classification',
                catalog_universe='all odd primes3..43',catalog=catalog,
                family_controls='Fermat-shaped2^m+1,m1..8; Mersenne-shaped2^k−1,k2..8; all eleven phases; actuator n0..127; actual infinite-family members k7,8,9; four literal charts at n0,1',
                eleven_phases=[dict(i=i,word=f"{z['r']:010b}",a=int(z['a']),d=z['d'],drift=z['drift']) for i,z in eleven.items()],
                H25_H3_family=families,literal_controls=literal,
                inherited_compiler_sha256=hashlib.sha256(source).hexdigest())
    (HERE/'phase_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(catalog=catalog,H25_H3_family=families),indent=2))


if __name__=='__main__':main()
