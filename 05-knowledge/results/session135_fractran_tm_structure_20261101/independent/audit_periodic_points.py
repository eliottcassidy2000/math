#!/usr/bin/env python3
"""Finite-horizon caps and exact periodic-point densities, independently checked."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from math import prod
import json,sys
sys.dont_write_bytecode=True
from audit_early import F as PRIMEGAME, P as PRIMES
ROOT=Path(__file__).resolve().parents[1]/'acceleration'
def check(ok,msg):
    if not ok:raise RuntimeError(msg)
def raw(n,rows):
    for i,(a,b) in enumerate(rows):
        if n*a%b==0:return n*a//b,i
    return n,None
def observe(n,rows,T):
    trace=[];start=n;period=None
    for t in range(1,T+1):
        n,i=raw(n,rows)
        if i is None:trace.append('halt');break
        trace.append(i)
        if n==start and period is None:period=t
    return trace,period

def density(rows,primes,caps,L):
    denom=prod(p**c for p,c in zip(primes,caps));nums={};counts={}
    weights=[[(p-1)*p**(c-e-1) if e<c else 1 for e in range(c+1)] for p,c in zip(primes,caps)]
    cases=0
    for x in product(*(range(c+1) for c in caps)):
        n=prod(p**e for p,e in zip(primes,x));_,period=observe(n,rows,L)
        if period is not None:
            nums[period]=nums.get(period,0)+prod(w[e] for w,e in zip(weights,x));counts[period]=counts.get(period,0)+1
        cases+=1
    return {k:str(F(v,denom)) for k,v in nums.items()},counts,cases

def main():
    pool=[(2**max(a,0)*3**max(b,0),2**max(-a,0)*3**max(-b,0)) for a,b in product(range(-1,2),repeat=2)]
    cases=0
    for rows in product(pool,repeat=2):
        M=[]
        for p in (2,3):
            es=[]
            for a,b in rows:
                e=0
                while b%p==0:b//=p;e+=1
                es.append(e)
            M.append(max(es))
        for T in range(1,5):
            caps=[T*m for m in M]
            for x in product(range(11),repeat=2):
                n=2**x[0]*3**x[1];y=[min(e,c) for e,c in zip(x,caps)];capped=2**y[0]*3**y[1]
                check(observe(n,rows,T)==observe(capped,rows,T),('finite cap',rows,T,x,caps));cases+=1
    # No post-horizon halt claim: both trace prefixes have two decrements,
    # but only the capped representative has halted immediately after them.
    check(observe(4,[(1,2)],2)[0]==observe(8,[(1,2)],2)[0],'T observations')
    check(raw(1,[(1,2)])[1] is None and raw(2,[(1,2)])[1]==0,'post-T hostile')
    toy,tcounts,_=density([(3,2),(2,3)],(2,3),(2,2),2)
    reverse,_,_=density([(2,3),(3,2)],(2,3),(2,2),2)
    prime,counts,total=density(PRIMEGAME,PRIMES,(2,)*10,2)
    check(toy=={2:'5/12'} and reverse=={2:'5/9'},'toy densities')
    analytic=F(4,7)*(F(12,143)+F(12,169))*prod(F(p-1,p) for p in (17,19,23,29))
    check(prime=={2:str(analytic)} and analytic==F(2654208,36409529),'PRIMEGAME point density')
    producer=json.loads((ROOT/'periodic.normal.json').read_text())
    for key,expected in [('toy',toy),('reversed_priority_toy',reverse),('primegame_period2',prime)]:
        check(producer[key]['density_by_least_period']=={str(k):v for k,v in expected.items()},('producer density',key))
    # Independently execute both literal ordered lists under actual prime
    # permutations, including spectator factors and a non-involutive 3-cycle.
    def rename_integer(n,mapping):
        remaining=n;renamed=1
        for p,q in mapping.items():
            while remaining%p==0:remaining//=p;renamed*=q
        return remaining*renamed
    renaming_cases=renaming_steps=readout_cases=0
    experiments=[([(3,2),(2,3)],{2:3,3:2}),
                 ([(3,2),(2,3)],{2:3,3:5,5:2}),
                 ([(1,10),(5,2)],{2:3,3:5,5:2}),
                 (PRIMEGAME,{2:3,3:5,5:2})]
    for rows,mapping in experiments:
        check(set(mapping)==set(mapping.values()),'finite permutation')
        inverse={q:p for p,q in mapping.items()}
        renamed_rows=[(rename_integer(a,mapping),rename_integer(b,mapping)) for a,b in rows]
        for n in range(1,513):
            m=rename_integer(n,mapping)
            check(rename_integer(m,inverse)==n,'inverse relabeling')
            for p,q in mapping.items():
                a=n;b=m;u=v=0
                while a%p==0:a//=p;u+=1
                while b%q==0:b//=q;v+=1
                check(u==v,'renamed output valuation');readout_cases+=1
            left=n;right=m
            for _ in range(8):
                left,i=raw(left,rows);right,j=raw(right,renamed_rows)
                check(i==j and rename_integer(left,mapping)==right,'literal conjugacy')
                if i is None:break
                renaming_steps+=1
            renaming_cases+=1
    polynomial=0
    for K in range(1,8):
        for x,y in product(range(12),repeat=2):
            for a,b,c,d in product(range(3),repeat=4):
                f=lambda u,v:a+b*u+c*v+d*u*v+u*u*v
                check(min(f(min(x,K),min(y,K)),K)==min(f(x,y),K),'positive polynomial cap');polynomial+=1
    # Affine phase costs over [5/14,7/5,11/7], input7*2^17*3.
    value=7*2**17*3
    for i in [0,1]*17:
        value,row=raw(value,[(5,14),(7,5),(11,7)]);check(row==i,'cost example word')
    check(value==21,'cost example endpoint')
    R=17;x=[17,1];D=[-1,0];prefixes=[[0,0],[-1,0]];coeffs=[([2,3],5),([4,1],7)]
    direct=sum(sum(a*z for a,z in zip(coef,[u+v+t*d for u,v,d in zip(x,offset,D)]))+b for t in range(R) for offset,(coef,b) in zip(prefixes,coeffs))
    C0=sum(sum(a*(u+v) for a,u,v in zip(coef,x,offset))+b for offset,(coef,b) in zip(prefixes,coeffs));C1=sum(sum(a*d for a,d in zip(coef,D)) for coef,b in coeffs)
    check(direct==R*C0+R*(R-1)*C1//2,'affine cost sum')
    print(json.dumps({'status':'PASS','finite_horizon_rule_and_return_cases':cases,'horizons':[1,4],'post_horizon_halt_counterexample':{'inputs':[4,8],'T':2},'primegame_period_cutoff':2,'primegame_capped_states':total,'primegame_representative_period_counts':counts,'primegame_periodic_point_density':str(analytic),'toy_periodic_point_density':toy[2],'reversed_toy_periodic_point_density':reverse[2],'toy_eventual_cycle_basin_density':'2/3','positive_polynomial_cap_cases':polynomial,'affine_cost_example':direct,'prime_permutation_input_cases':renaming_cases,'prime_permutation_literal_steps':renaming_steps,'renamed_readout_checks':readout_cases},indent=2))
if __name__=='__main__':main()
