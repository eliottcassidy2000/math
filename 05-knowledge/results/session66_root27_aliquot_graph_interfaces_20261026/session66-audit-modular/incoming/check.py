#!/usr/bin/env python3
"""Independent table-shift iff and prime-power compatibility controls."""
from collections import Counter
from math import lcm
import json
import hashlib
from pathlib import Path

checks=Counter()
def need(ok,label,witness=None):
    checks[label]+=1
    if not ok:
        raise RuntimeError((label,witness))
def order2(m):
    if m<=1 or m%2==0:
        raise ValueError('odd modulus at least three')
    x=1
    for j in range(1,m+1):
        x=2*x%m
        if x==1:
            return j
    raise RuntimeError('order bound failed')
def factor_powers(m):
    result=[]
    p=3
    while p*p<=m:
        if m%p==0:
            power=1
            while m%p==0:
                m//=p
                power*=p
            result.append((p,power))
        p+=2
    if m>1:
        result.append((m,m))
    return result
def v2(x):
    return (x & -x).bit_length()-1
def table(m,k,b):
    inv=pow(pow(2,k,m),-1,m)
    # Odd representatives, rather than implicitly assuming source coverage.
    return {(x%m,((3*x+b)*inv)%m) for x in range(1,2*m,2)}

base=Path(__file__).resolve().parent
source_manifest=json.loads((base/'SOURCE_MANIFEST.json').read_text())
for row in source_manifest['sources']:
    data=(base/row['local']).read_bytes()
    blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    need(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256']
         and blob==row['git_blob'],'source_integrity')

records=[]
for m in range(3,5002,2):
    d=order2(m)
    orbit=[]
    x=1
    for t in range(d):
        orbit.append(x)
        x=2*x%m
    shifts=[t for t,x in enumerate(orbit) if x==m-1]
    pp=factor_powers(m)
    base_orders=[order2(p) for p,power in pp]
    lifted_orders=[order2(power) for p,power in pp]
    ss=[v2(o) for o in base_orders]
    criterion=len(set(ss))==1 and ss[0]>0
    need(bool(shifts)==criterion,'prime_criterion',m)
    for (p,power),base,lifted in zip(pp,base_orders,lifted_orders):
        need(v2(base)==v2(lifted),'prime_power_preserves_two_part',(p,power))
        quotient=lifted//base
        while quotient%p==0:
            quotient//=p
        need(lifted%base==0 and quotient==1,'prime_power_order_ratio',(p,power))
    need(d==lcm(*lifted_orders),'CRT_order',m)
    if criterion:
        s=ss[0]
        U=lcm(*(o//2**s for o in lifted_orders))
        t=2**(s-1)*U
        need(shifts==[t] and t==d//2 and pow(2,t,m)==m-1,'least_shift_constructor',m)
        need(all((t//(o//2))%2==1 and t%(o//2)==0 for o in lifted_orders),'odd_multiple_local_halves',m)
    if m<=101:
        records.append({'m':m,'order':d,'shift':shifts[0] if shifts else None,'prime_order_two_parts':ss})

for m in range(3,102,2):
    d=order2(m)
    for k in range(1,7):
        plus=table(m,k,1)
        minus=table(m,k,-1)
        need(minus=={((-x)%m,(-y)%m) for x,y in plus},'negation_tables',(m,k))
        need({x for x,y in plus}==set(range(m)),'odd_sources_cover_residues',(m,k))
        for t in range(d):
            shifted=table(m,k+t,-1)
            source_negated={((-x)%m,y) for x,y in plus}
            equal=shifted==source_negated
            need(equal==(pow(2,t,m)==m-1),'full_table_shift_iff',(m,k,t))
            zero_plus=next(y for x,y in plus if x==0)
            zero_minus=next(y for x,y in shifted if x==0)
            need((zero_plus==zero_minus)==(pow(2,t,m)==m-1),'constant_term_necessity',(m,k,t))

# Slopes alone lose the exact obstruction when 3 is not invertible.
m,k,t=3,1,0
inv=pow(pow(2,k,m),-1,m)
inv_shift=pow(pow(2,k+t,m),-1,m)
slope_equal=(3*inv_shift)%m==(-3*inv)%m
need(slope_equal and table(m,k+t,-1)!={((-x)%m,y) for x,y in table(m,k,1)},'slope_cancellation_hostile')
# Unit translation is not even a self-map of the odd source domain.
need(all((x-1)%2==0 for x in range(-101,102,2)),'unit_translation_domain_hostile')
need(sum(1 for x in range(1,2001,2) if (x-1)%2==1)==0,'producer_T1b_empty_filter')
# The excluded modulus one is trivial, without an even-order requirement.
need(all(table(1,k,b)=={(0,0)} for k in range(1,7) for b in (-1,1)),'modulus_one_boundary')
print(json.dumps({'status':'PASS','checks':dict(sorted(checks.items())),
                  'total_checks':sum(checks.values()),'moduli_classified':2500,
                  'full_table_moduli':50,'full_table_k_range':[1,6],
                  'records_through_101':records,
                  'proof_repairs':['use source residue zero, not slope cancellation',
                                   'state odd modulus m>=3; m=1 is trivial separately',
                                   'T1b passing counterexample tests translation by 2, not 1']},
                 indent=2,sort_keys=True))
