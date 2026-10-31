#!/usr/bin/env python3
"""Extra own-stdlib finite state audit and exact integer lattice determinant."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import json
from verify_primegame import FRACS, HAND, PRIMES, step_integer, step_exponent, step_control, value

ROOT=Path(__file__).resolve().parent
checks=0

def require(x,msg):
    global checks
    checks+=1
    if not x: raise RuntimeError(msg)

def determinant(matrix):
    a=[[Fraction(x) for x in row] for row in matrix];det=Fraction(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None: return Fraction(0)
        if pivot!=i: a[i],a[pivot]=a[pivot],a[i];det=-det
        v=a[i][i];det*=v;a[i]=[x/v for x in a[i]]
        for j in range(i+1,len(a)):
            v=a[j][i];a[j]=[x-v*y for x,y in zip(a[j],a[i])]
    return det

small_onehot_states=0
for counters in product(range(4),repeat=4):
    for t in (1,11,13,17,19,23,29):
        e=counters+tuple(int(t==p) for p in PRIMES[4:])
        n,j=step_integer(value(e));ee,k=step_exponent(e);cc,l=step_control(e)
        require(j==k==l and n==value(ee) and ee==cc,'onehot control mismatch')
        require(sum(ee[4:])<=1,'onehot invariant mismatch')
        small_onehot_states+=1

small_general_states=0
for e in product(range(2),repeat=10):
    n,j=step_integer(value(e));ee,k=step_exponent(e)
    require(j==k and n==value(ee),'general vector mismatch')
    small_general_states+=1

cols=(8,12,9,0,13,4,5,11,2,3)
vectors=[]
for j in cols:
    take,give=HAND[j]
    vectors.append([give.count(i)-take.count(i) for i in range(10)])
matrix=[[vectors[j][i] for j in range(10)] for i in range(10)]
det=determinant(matrix)
require(det==-1,'lattice minor not unimodular')

# Our own elementary factorization for F,U,S, independent of the fixed-basis factorizer.
def fus(n):
    require(n>=2,'FUS scope requires n>=2')
    remaining=n;exponents=[];p=2
    while p*p<=remaining:
        a=0
        while remaining%p==0:remaining//=p;a+=1
        if a:exponents.append(a)
        p+=1
    if remaining>1:exponents.append(1)
    tau=1
    for a in exponents:tau*=a+1
    prime=(len(exponents)==1 and exponents[0]==1)
    squarefree=all(a==1 for a in exponents)
    F=tau-2;U=len(exponents)-int(prime);S=2**len(exponents)-1-int(squarefree)
    return {'N':n,'factor_exponents':exponents,'F':F,'U':U,'S':S,'defect':F-U-S}

samples=[fus(n) for n in (2,7,17,55,15,825,343,60)]+[fus(2**p) for p in PRIMES]
for p,r in zip(PRIMES,samples[-10:]):
    require((r['F'],r['U'],r['S'],r['defect'])==(p-1,1,1,p-3),'output FUS mismatch')
require(step_integer(91)==(17,0),'priority witness mismatch')
report={'onehot_states_tested':small_onehot_states,'onehot_counter_range':[0,3],'general_binary_vectors_tested':small_general_states,'unimodular_minor_instruction_indices_1based':[x+1 for x in cols],'minor_determinant':int(det),'FUS_samples':samples,'priority_witness':{'integer':91,'chosen_fraction':'17/91','chosen_next':17,'later_enabled_fraction':'11/13','wrong_later_next':77},'explicit_checks':checks}
(ROOT/'structural_audit_results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
