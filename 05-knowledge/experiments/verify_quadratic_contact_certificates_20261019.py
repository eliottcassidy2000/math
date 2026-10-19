#!/usr/bin/env python3
"""Replay all three quadratic-contact certificates using only Python stdlib."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from collections import Counter
from math import isqrt
import json
from rational_contact_certificate_core import Poly,variables,decode,system
checks=Counter()
def check(b,name):
    checks[name]+=1
    if not b:raise RuntimeError(name)
a,b,d,v,m=variables(5)
patterns=['adjacent_plus','adjacent_minus','opposite'];records=[]
for pattern in patterns:
    path=Path(__file__).resolve().parent.parent/'results'/('quadratic_membership_'+pattern+'.json')
    if not path.exists(): path=Path(__file__).with_name('quadratic_membership_'+pattern+'.json')
    payload=path.read_bytes();data=json.loads(payload);inputs=system(pattern)
    check(data['variables']==['a','b','d','v','m'],'variable_order')
    check([decode(p) for p in data['inputs']]==inputs,'original_coefficient_system')
    vals=[]
    for i,node in enumerate(data['nodes']):
        if 'input' in node:
            check(0<=node['input']<len(inputs),'input_reference');z=inputs[node['input']]
        else:
            z=Poly(5)
            for j,p in node['terms']:
                check(0<=j<i,'acyclic_reference');z+=decode(p)*vals[j]
        vals.append(z)
    target=(m*m+4 if pattern=='adjacent_plus' else
            (m-2)**3*(m+2)*(m*m+20*m+52) if pattern=='adjacent_minus' else m*m+16)
    check(decode(data['target'])==target,'intended_target')
    check(vals[data['target_node']]==target,'exact_polynomial_combination')
    records.append({'pattern':pattern,'input_equations':len(inputs),'certificate_nodes':len(vals),
                    'certificate_sha256':sha256(payload).hexdigest()})
# Standard-library checks of the elementary adjacent-minus exceptional cases.
check(20*20-4*52==192 and isqrt(192)**2!=192,'quadratic_has_no_rational_roots')
q=a*a+b*b;c=m*a/2;rho=1+m*(q-1)/2
check((c-a).substitute({4:2})==0 and (rho-q).substitute({4:2})==0,'m_two_coincident_poles')
ps=system('adjacent_minus')
for sigma in [1,-1]:
    # First substitute m,v. Substitution of the polynomial d is implemented
    # directly as an exact sum of monomials, independent of any CAS.
    def replace_d(P):
        out=Poly(5)
        for mon,c0 in P.terms.items():
            term=Poly.scalar(5,c0)
            for i,power in enumerate(mon):term*=((-2*sigma*a-b) if i==2 else variables(5)[i])**power
            out+=term
        return out
    e4=replace_d(ps[4].substitute({4:-2,3:2*sigma}))
    e6=replace_d(ps[6].substitute({4:-2,3:2*sigma}))
    check(e4==-16*(3*a*a+2*a-b*b+1),'specialized_coefficient')
    check(e6==32*(3*a*a+2*sigma*a*b+b*b-1),'specialized_norm')
    check((3*a*a+2*a-b*b+1)+(3*a*a+2*sigma*a*b+b*b-1)==2*a*(3*a+1+sigma*b),'elementary_case_split')
    bb=-sigma*(3*a+1)
    check(3*a*a+2*sigma*a*bb+bb*bb-1==2*a*(3*a+2),'nonzero_a_value')
    check(ps[1].substitute({0:F(-2,3),1:sigma,2:F(sigma,3),3:2*sigma,4:-2})==F(128,27),'exceptional_branch_contradiction')
print(json.dumps({'status':'PROVED: exact certificate replay, standard library only',
 'checks':dict(checks),'certificates':records},indent=2))
