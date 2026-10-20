#!/usr/bin/env python3
"""Tiny exact state/algebra controls; standard library; no repository writes."""
from fractions import Fraction as F
from itertools import product
import json

checks={}
def check(k,b):
    if not b: raise RuntimeError(k)
    checks[k]=checks.get(k,0)+1

def basis_add(basis,v):
    v=list(map(F,v))
    for pivot,row in basis:
        scale=v[pivot]
        if scale: v=[x-scale*y for x,y in zip(v,row)]
    pivot=next((i for i,x in enumerate(v) if x),None)
    if pivot is None: return False
    scale=v[pivot];v=[x/scale for x in v]
    basis.append((pivot,v));basis.sort(key=lambda z:z[0])
    return True

def closure(size,maps,obs,algebra):
    basis=[]
    for v in [[1]*size]+obs: basis_add(basis,v)
    while True:
        old=[row for _,row in basis];before=len(basis)
        for f in old:
            for t in maps: basis_add(basis,[f[t[x]] for x in range(size)])
        if algebra:
            for f,g in product(old,repeat=2): basis_add(basis,[a*b for a,b in zip(f,g)])
        if len(basis)==before:return basis

def canonical_labels(values):
    seen={};return [seen.setdefault(v,len(seen)) for v in values]

def partition(size,maps,obs):
    labels=canonical_labels([tuple(f[x] for f in obs) for x in range(size)])
    depth=0
    while True:
        new=canonical_labels([(labels[x],)+tuple(labels[t[x]] for t in maps) for x in range(size)])
        if new==labels:return labels,depth
        labels=new;depth+=1

def chi(a,x):return (-1)**((a&x).bit_count()%2)
def apply_matrix(rows,x):return sum(((r&x).bit_count()%2)<<i for i,r in enumerate(rows))
def pull_label(rows,a):
    v=0
    for i,r in enumerate(rows):
        if (a>>i)&1:v^=r
    return v

def span2(values):
    out={0}
    for v in values:out|={x^v for x in list(out)}
    return out

def label_closure(rows_list,observed):
    out=span2(observed)
    while True:
        new=span2(out|{pull_label(rows,a) for rows in rows_list for a in out})
        if new==out:return out
        out=new

def orbit_labels(rows_list,observed):
    out={0}|set(observed)
    while True:
        new=out|{pull_label(rows,a) for rows in rows_list for a in out}
        if new==out:return out
        out=new

cases=[]
def run_case(name,size,maps,obs,expected_linear,expected_algebra):
    lin=closure(size,maps,obs,False);alg=closure(size,maps,obs,True)
    labels,depth=partition(size,maps,obs);q=len(set(labels))
    check('linear_dimension',len(lin)==expected_linear)
    check('algebra_dimension',len(alg)==expected_algebra==q)
    check('word_horizon',depth<=size-1)
    for c in range(q):
        trial=[(p,row[:]) for p,row in alg]
        check('quotient_indicator',not basis_add(trial,[int(labels[x]==c) for x in range(size)]))
    for t in maps:
        for x,y in product(range(size),repeat=2):
            if labels[x]==labels[y]:check('stable_partition',labels[t[x]]==labels[t[y]])
    case={'name':name,'states':size,'linear_dimension':len(lin),'algebra_dimension':len(alg),'classes':q,'refinement_depth':depth}
    cases.append(case)
    return labels

run_case('nonlinear_five_state',5,[[0,2,3,3,4]],[[0,0,0,1,1]],4,4)
run_case('distinct_scalar_identity',3,[list(range(3))],[[0,1,2]],2,3)
run_case('summed_character_hostile',4,[list(range(4))],[[chi(1,x)+chi(2,x) for x in range(4)]],2,3)

binary=[
 ('identity_two_bits',2,[(1,2)],[0],[1,2],3,4),
 ('cyclic_three_register',3,[(2,4,1)],[0],[1],4,8),
 ('affine_two_map_quotient',3,[(2,1,4),(3,2,4)],[1,4],[1],4,4),
 ('zero_map_forgets_future',3,[(0,0,0)],[0],[1],2,2),
 ('no_observed_character',3,[(2,4,1)],[0],[],1,1),
]
for name,n,rows_list,shifts,obs_masks,ld,ad in binary:
    maps=[[apply_matrix(rows,x)^b for x in range(1<<n)] for rows,b in zip(rows_list,shifts)]
    obs=[[chi(a,x) for x in range(1<<n)] for a in obs_masks]
    labels=run_case(name,1<<n,maps,obs,ld,ad)
    L=label_closure(rows_list,obs_masks);K=orbit_labels(rows_list,obs_masks)
    check('binary_dimension_formula',len(L)==ad and len(K)==ld)
    for x,y in product(range(1<<n),repeat=2):
        check('annihilator_equivalence',(labels[x]==labels[y])==all(chi(a,x^y)==1 for a in L))
    for rows,b in zip(rows_list,shifts):
        for a,x in product(L,range(1<<n)):
            check('affine_pullback',chi(a,apply_matrix(rows,x)^b)==chi(a,b)*chi(pull_label(rows,a),x))
    # Arbitrary quotient law recovered from full L-moments and used for dynamics.
    q=len(L);representatives=[labels.index(c) for c in range(q)]
    raw=list(range(1,q+1));lam=[F(v,sum(raw)) for v in raw]
    mu={a:sum(lam[c]*chi(a,x) for c,x in enumerate(representatives)) for a in L}
    for c,x in enumerate(representatives):
        inverse=sum(mu[a]*chi(a,x) for a in L)/q
        check('quotient_fourier_inverse',inverse==lam[c])
    for t,rows,b in zip(maps,rows_list,shifts):
        pushed=[F(0)]*q
        for c,x in enumerate(representatives):pushed[labels[t[x]]]+=lam[c]
        for a in L:
            actual=sum(pushed[c]*chi(a,x) for c,x in enumerate(representatives))
            check('moment_transition',actual==chi(a,b)*mu[pull_label(rows,a)])

# Fully observed trajectories, indistinguishable separate time expectations.
n=3;plus=[F(1+chi(7,x),8) for x in range(8)];minus=[F(1-chi(7,x),8) for x in range(8)]
for a in [0,1,2,4]:check('cyclic_register_mean_hostile',sum((p-q)*chi(a,x) for x,(p,q) in enumerate(zip(plus,minus)))==0)
check('cyclic_register_joint_repair',sum((p-q)*chi(7,x) for x,(p,q) in enumerate(zip(plus,minus)))==2)

print(json.dumps({'status':'PASS','exact_checks':sum(checks.values()),'checks':checks,'cases':cases,
 'scope':'Eight declared tiny systems; exact rational closure and independent partition refinement. No broad census or conjecture closure.'},indent=2,sort_keys=True))
