#!/usr/bin/env python3
"""Exact verification for the follow-on forbidden-pattern classification."""
from itertools import combinations
from math import gcd
from fractions import Fraction
from pathlib import Path
import json

def invariants(q,B):
    J=[sum((b-t)%q in B for b in B) for t in range(q)]
    E=sum(j*j for j in J)
    T=sum((a+b)%q in B for a in B for b in B)
    s=len(B);D=s**3-q*T
    return {'size':s,'energy':E,'triples':T,'surplus_numerator':D,'conditional_mean':str(1+Fraction(D,(q-s)**3)),'J':J}

def primes(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,int(p**.5)+1))]

def symsets(q):
    orbits=[{a,(-a)%q} for a in range(1,q) if a<=(-a)%q]
    for bits in range(1<<len(orbits)):
        B=set()
        for i,O in enumerate(orbits):
            if bits>>i&1:B|=O
        yield B

out={};n=0
for p in primes(499):
    if p<11:continue
    for a in range(2,p-1):
        B={1,p-1,a,p-a};I=invariants(p,B)
        triangle={2,p-2,pow(2,-1,p),(-pow(2,-1,p))%p}
        energy={3,p-3,pow(3,-1,p),(-pow(3,-1,p))%p}
        assert I['triples']==(6 if a in triangle else 0)
        assert I['energy']==(44 if a in energy else 36)
        assert not triangle&energy
        n+=1
out['prime_ratio_classification_cases']=n
out['prime_ratio_range']='every prime 11 <= p <= 499, every a not 0,+1,-1'
# All symmetric loopless connected cyclic patterns of every size, through order12.
census=[];first=None
for q in range(2,13):
    seen={};count=0
    pairs=set()
    for B in symsets(q):
        if not B or gcd(q,*B)!=1:continue
        I=invariants(q,B);count+=1
        if len(B)==4:pairs.add((I['energy'],I['triples'],I['surplus_numerator']))
        key=(len(B),I['energy'])
        for prevD,prevB in seen.get(key,[]):
            if prevD*I['surplus_numerator']<0 and first is None:
                first={'q':q,'B1':sorted(prevB),'B2':sorted(B),'I1':invariants(q,prevB),'I2':I}
        seen.setdefault(key,[]).append((I['surplus_numerator'],B))
    census.append({'q':q,'connected_patterns_all_sizes':count,'size4_energy_triples_surplus':sorted(pairs)})
assert first['q']==12
out['connected_small_census']=census
out['first_connected_cyclic_example']=first
out['general_cyclic_order6']=[{'B':sorted(B),**invariants(6,B)} for B in ({1,5},{2,4})]
assert [v['energy'] for v in out['general_cyclic_order6']]==[6,6]
assert [v['triples'] for v in out['general_cyclic_order6']]==[0,2]
assert out['general_cyclic_order6'][0]['J']==out['general_cyclic_order6'][1]['J']
# All symmetric sets, including those containing0, on odd cyclic groups <=27.
recon=0
for q in range(3,28,2):
    for base_B in symsets(q):
        for B in (base_B,base_B|{0}):
            J=[sum((b-t)%q in B for b in B) for t in range(q)]
            assert all(J[(2*x)%q]%2==int(x in B) for x in range(q))
            recon+=1
out['odd_order_parity_reconstruction_cases']=recon
out['all_checks_passed']=True
base=Path(__file__).resolve().parent
(base/'forbidden_pattern_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('ALL CHECKS PASSED')
