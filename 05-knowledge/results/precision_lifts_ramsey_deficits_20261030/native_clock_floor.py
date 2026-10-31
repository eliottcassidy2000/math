#!/usr/bin/env python3
"""Exact necessary (clock,length) filters; not a positive-word enumerator."""
from math import isqrt
from pathlib import Path
import json


def need(x):
    if not x:raise ValueError('clock/lattice control failed')


def floor_sum(n,m,a,b):
    ans=0
    while True:
        qa,a=divmod(a,m);qb,b=divmod(b,m)
        ans+=qa*n*(n-1)//2+qb*n
        y=a*n+b
        if y<m:return ans
        n,b=divmod(y,m);m,a=a,m


def count(N,M,a,b):
    return floor_sum(N,M,a,a+b-1)-floor_sum(N,M,a-1,a-1+b)


tests=0
for M in range(2,38):
    for a in (0,1,M//2,M-1,M+3):
        for b in (-M-1,-1,0,1,M//2,M+1):
            for N in (0,1,M//2,M-1):
                literal=sum(1 for A in range(1,N+1) if 1<=(a*A+b)%M<A)
                need(count(N,M,a,b)==literal);tests+=1

q=5**14;M=4*5**13;step=isqrt(M)+1
need(pow(3,M,q)==1 and pow(3,M//2,q)!=1 and pow(3,M//5,q)!=1)
baby={};x=1
for j in range(step):
    baby.setdefault(x,j);x=3*x%q
factor=pow(pow(3,step,q),-1,q)
def log3(target):
    x=target
    for i in range(step+1):
        if x in baby:
            value=(i*step+baby[x])%M
            need(pow(3,value,q)==target)
            return value
        x=x*factor%q
    raise ValueError('target not found in unit group')

a=log3(2)
source=json.loads(Path('native14_controls.json').read_text())
records=[]
for pattern in source['patterns']:
    target_m,target_t=pattern['required_core_map_mod_5_14']
    b=log3(target_m);cap=pattern['max_integer_core_clock']
    need(cap<M)
    total=count(cap,M,a,b)
    def select(rank):
        lo,hi=1,cap
        while lo<hi:
            mid=(lo+hi)//2
            if count(mid,M,a,b)>=rank:hi=mid
            else:lo=mid+1
        length=(a*lo+b)%M
        need(1<=length<lo)
        need(pow(3,length,q)==target_m*pow(2,lo,q)%q)
        need(count(lo,M,a,b)==rank and count(lo-1,M,a,b)==rank-1)
        return [lo,length]
    row=dict(pattern,b=b,necessary_clock_length_pairs=total,
        first_pair=select(1) if total else None)
    if total<=100:row['all_necessary_clock_length_pairs']=[select(i) for i in range(1,total+1)]
    records.append(row)

print(json.dumps(dict(status='FINITE-EXACT necessary slope constraints; translation, word order and native filters still required',
    q=q,unit_group_order=M,log_base3_of2=a,floor_sum_literal_checks=tests,patterns=records,
    pass_checks=True),indent=2))
