#!/usr/bin/env python3
"""Exact finite controls for growing-alphabet information bounds. Standard library only."""
import argparse,json,math
from collections import Counter
from decimal import Decimal,localcontext
from fractions import Fraction
from itertools import product
from pathlib import Path
HERE=Path(__file__).resolve().parent
CHECKS=Counter()
def need(ok,label):
    CHECKS[label]+=1
    if not ok:raise RuntimeError(label)
def entropy_compare(a,target):
    """Sign of f(a)-target via exact integer powers for rational arguments."""
    p,q=a.numerator,a.denominator;u,v=target.numerator,target.denominator
    left=(p+q)**(v*(p+q));right=2**(q*u)*p**(v*p)*q**(v*q)
    return (left>right)-(left<right)
def numerical_root(target):
    with localcontext() as c:
        c.prec=70;lo=Decimal(0);hi=Decimal(1);ln2=Decimal(2).ln()
        for _ in range(210):
            a=(lo+hi)/2;f=((1+a)*(1+a).ln()-a*a.ln())/ln2
            if f<target:lo=a
            else:hi=a
        return (lo+hi)/2

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    need(17**17<2**70,'row_entropy_exact_power_bound')
    need(65**128<2**15*64**128,'strong_lower_bound_induction_step')
    need(64**128==2**(15*64-192),'strong_lower_bound_base_case')
    # Independent enumeration of weak count vectors, including total-mass arrays.
    for n in range(1,7):
        for L in range(1,5):
            rows=sum(sum(x)<=n for x in product(range(n+1),repeat=L))
            need(rows==math.comb(n+L,L),'row_stars_and_bars')
    for N in range(7):
        for K in range(1,6):
            arrays=sum(sum(x)==N for x in product(range(N+1),repeat=K))
            need(arrays==math.comb(N+K-1,K-1),'fixed_mass_stars_and_bars')
    tested=[]
    for n in range(64,513):
        L=n//16;N=n*(n-1)//2;R=math.comb(n+L,L)
        need(R**16<2**(6*n),'explicit_row_upper_bound')
        need(n**128<=2**(15*n-192),'explicit_strong_lower_bound')
        need(R**n*math.factorial(n)<2**(n*(n-3)//2),'exact_pigeonhole_failure')
        need(math.comb(N+n*L-1,n*L-1)*math.factorial(n)<2**(n*(n-3)//2),'exact_fixed_mass_failure')
        tested.append([n,L])
    # Rigorous rational brackets. Approximate decimals below are explanatory only.
    rowlo,rowhi=Fraction(1048,10000),Fraction(1049,10000)
    masslo,masshi=Fraction(1469,10000),Fraction(147,1000)
    need(entropy_compare(rowlo,Fraction(1,2))<0,'row_root_lower_rational_bracket')
    need(entropy_compare(rowhi,Fraction(1,2))>0,'row_root_upper_rational_bracket')
    need(entropy_compare(2*masslo,Fraction(1))<0,'mass_root_lower_rational_bracket')
    need(entropy_compare(2*masshi,Fraction(1))>0,'mass_root_upper_rational_bracket')
    ar=numerical_root(Decimal('0.5'));am=numerical_root(Decimal(1))/2
    out={'status':'PASS','scope':'Class separation by one outgoing-row multiset signature; not authenticated global wiring reconstruction.',
         'explicit_obstruction':{'all_order_theorem':'Every integer n>=64 and 1<=L<=n/16 fails','exact_finite_controls':[64,512],'maximal_integer_L_tested':'floor(n/16)','orders_tested':len(tested)},
         'necessary_asymptotic_constants':{'row_bound':{'equation':'f(a)=1/2','approximate':format(ar,'.24f'),'exact_rational_bracket':['0.1048','0.1049']},'fixed_mass_bound':{'equation':'f(2a)=1','approximate':format(am,'.24f'),'exact_rational_bracket':['0.1469','0.1470']}},
         'integer_certificate':{'17^17':str(17**17),'2^70':str(2**70),'65^128_lt_2^15_times_64^128':True},'checks':dict(CHECKS)}
    dest=HERE/'growing_alphabet_receipt.json'
    if args.write:dest.write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads(dest.read_text())==out,'frozen_receipt_matches')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
