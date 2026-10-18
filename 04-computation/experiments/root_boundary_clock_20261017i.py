#!/usr/bin/env python3
"""Evaluate dyadic boundary forcing without enumerating a huge source interval.

Exact integer arithmetic; uses R8's proved band period, the independently
checked K recursion, and the sparse support theorem of this session.
Periodicity of the resulting error is tested, not assumed.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path

import numpy as np
from root_full_support_20261017i import K0, ideal_step, boundary_seeds, require


def build_tables(depth):
    tabs=[K0]
    coefficient_bound=3
    for d in range(1,depth+1):
        bands=2*3**(d+1)
        coefficient_bound*=bands
        require(coefficient_bound<2**63,'int64 table safety from R4')
        prev=tabs[-1][np.arange(bands)%len(tabs[-1])]
        tabs.append(ideal_step(prev))
    return tabs


def kval(k,d,s,tabs):
    b=((4*3**d*(1<<k))//((1<<d)*s)).bit_length()-1
    if b<0:return 0
    tab=tabs[d]
    return int(tab[b%len(tab),s%tab.shape[1]])


def push_one(s):
    a=2 if s%8==1 else 1
    t=(3*s+1)//(1<<a)
    while t%8==5:t=(t-1)//4
    return t,(-1)**a


def forcing(k,d,t,tabs):
    """Coefficient of L k_(d-1)-k_d at literal primitive t."""
    upper=4*3**(d-1)*(1<<k)//(1<<(d-1))
    total=0
    for A in range(1,(3*upper+3).bit_length()+2):
        a=1 if A%2 else 2
        ca=(3+(1<<a)-(1<<A))//3
        num=(1<<A)*t-ca
        if num%3:continue
        s=num//3
        if s<=upper:
            require(s>0 and s%2 and s%8!=5,'native inverse')
            require(3*(3*s+1)==(1<<a)*((1<<(A-a))*(3*t+1)-1),'inverse braid identity')
            total+=(-1)**a*kval(k,d-1,s,tabs)
    return total-kval(k,d,t,tabs)


def log2_floor(x):
    k=x.numerator.bit_length()-x.denominator.bit_length()
    if x < Q(2)**k:k-=1
    require(Q(2)**k<=x<Q(2)**(k+1),'rational binary logarithm')
    return k


def forcing_threshold(d,t):
    """A sufficient exact k threshold in the eventual-periodicity proof."""
    C=Q(36*3**(d-1),2**(d-1))
    D=3*t+1
    kappa=log2_floor(C/D)
    gap=D-C/(Q(2)**(kappa+1))
    require(gap>0,'non-dyadic boundary')
    A0=1
    while (1<<A0)*gap<=7:A0+=1
    threshold=max(0,A0-kappa-1,-log2_floor(Q(4*3**d,(2**d)*t)))
    for A in range(1,A0):
        a=1 if A%2 else 2
        num=(1<<A)*D-(3+(1<<a))
        if num%9==0:
            s=num//9
            threshold=max(threshold,-log2_floor(Q(4*3**(d-1),(2**(d-1))*s)))
    return threshold


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--depth',type=int,default=5)
    ap.add_argument('--max-k',type=int,default=1600)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    tabs=build_tables(args.depth)
    vectors=[[{} for _ in range(args.max_k+1)] for _ in range(args.depth+1)]
    sums=[]
    for k in range(args.max_k+1):
        e={1:(-2,-1,1,-1,1,2)[k%6]}
        vectors[0][k]=e
        for d in range(1,args.depth+1):
            nxt={}
            for s,c in e.items():
                t,sgn=push_one(s);nxt[t]=nxt.get(t,0)+sgn*c
            for t in boundary_seeds(1,d):
                nxt[t]=nxt.get(t,0)+forcing(k,d,t,tabs)
            e={t:c for t,c in nxt.items() if c}
            vectors[d][k]=e
    threshold=0
    for d in range(args.depth+1):
        if d:
            threshold=max(threshold,max(forcing_threshold(d,t) for t in boundary_seeds(1,d)))
        period=2*3**(d+1)
        failures=[k for k in range(args.max_k-period+1) if vectors[d][k]!=vectors[d][k+period]]
        energies=[sum(c*c for c in e.values()) for e in vectors[d]]
        coefficient_bound=2**d*3**(1+d*(d+3)//2)
        error_l1_bound=2+6*d*(d+1)*coefficient_bound
        observed_l1=max(sum(abs(c) for c in e.values()) for e in vectors[d])
        observed_coefficient=int(np.abs(tabs[d]).max())
        require(observed_l1<=error_l1_bound,'R4 error bound')
        require(observed_coefficient<=coefficient_bound,'R4 coefficient bound')
        require(threshold+period<=args.max_k+1,'complete finite bank for all-k consequence')
        require(not any(k>=threshold for k in failures),'proved eventual period')
        sums.append(dict(d=d,tested_period=period,proved_period_threshold=threshold,period_failures=failures,
                         all_k_max_energy=max(energies[:threshold+period]),
                         max_energy=max(energies),maximizing_k=energies.index(max(energies)),
                         coefficient_bound=coefficient_bound,observed_max_coefficient=observed_coefficient,
                         error_l1_bound=error_l1_bound,observed_max_error_l1=observed_l1))
    # Existing complete-support observations are independent finite controls.
    whole=Path(__file__).with_name('energy_extraction_20261017i')/'root_full_support.json'
    checked=sample_checks=0
    ideal_energies=[]
    for d,tab in enumerate(tabs):
        period,q=tab.shape
        # Python integers: square and accumulate without int64 overflow.
        numerator=sum(sum(int(x)**2 for x in row)*2**(period-1-b) for b,row in enumerate(tab))
        ideal_energies.append(Q(numerator,2**(d+1)*(2**period-1)))
    if whole.exists():
        for row in json.loads(whole.read_text())['rows']:
            n,d=row['N'],row['d']
            if n&(n-1) or d>args.depth:continue
            k=n.bit_length()-1
            e=vectors[d][k]
            require(sum(c*c for c in e.values())==row['full_error'],'independent full-support energy')
            require(sorted(e)==row['error_seeds'],'independent full-support set')
            require([[t,c] for t,c in sorted(e.items())]==row['error_vector'],'independent signed error vector')
            L=log2_floor(Q(4*3**d*n,2**d))
            B=2**d*3**(1+d*(d+3)//2)
            require(abs(Q(row['ideal_lattice_energy'])-n*ideal_energies[d])<=B**2*(3*3**(d+1)*(L+1)+1),
                    'R4 exact sampling bound')
            checked+=1
            sample_checks+=1
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(dict(status='FINITE-EXACT',max_k=args.max_k,
        full_support_controls=checked,sampling_bound_controls=sample_checks,
        ideal_energies=[str(x) for x in ideal_energies],summary=sums,
        vectors=[[[[t,c] for t,c in sorted(e.items())] for e in layer] for layer in vectors]),separators=(',',':'))+'\n')
    print(json.dumps(sums,indent=2))


if __name__=='__main__':main()
