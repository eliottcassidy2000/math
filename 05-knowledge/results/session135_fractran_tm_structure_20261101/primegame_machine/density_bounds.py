#!/usr/bin/env python3
"""Exact lower density bounds from complete, unbounded onehot and T=0 domains."""
from fractions import Fraction as F
from math import prod
import json

CONTROL=(11,13,17,19,23,29)
C0=prod(F(p-1,p) for p in CONTROL)
s=sum((F(1,p) for p in CONTROL),F(0))
onehot_mass=C0*(1+s)
prime_onehot=C0*(F(2,3)+(F(1,11)+F(1,13))*F(9,35)
                 +(F(1,17)+F(1,19)+F(1,29))*F(3,5)+F(1,23))
cycle2_onehot=onehot_mass-prime_onehot
# Additional disjoint closed domain: k>=2,T=a+c+i=0.
cycle2_T0_multi=F(2,5)*(F(22,23)-C0*(1+sum((F(1,p) for p in CONTROL if p!=23),F(0))))
cycle2_lower=cycle2_onehot+cycle2_T0_multi
# Inherited periodic-point density gives a separate cycle6-basin lower bound.
cycle6_lower=F(37396469504,151504041213525)
lower={'prime':prime_onehot,'cycle2':cycle2_lower,'cycle6':cycle6_lower}
if not all(x>0 for x in lower.values()) or sum(lower.values())>1:
    raise RuntimeError('density bounds inconsistent')

def exact(x):
    return {'numerator':x.numerator,'denominator':x.denominator,'decimal_for_display_only':float(x)}

report={'status':'EXACT bounds; displayed decimals are approximations, not limits from a census',
        'onehot_mass':exact(onehot_mass),'onehot_prime_density':exact(prime_onehot),
        'onehot_cycle2_density':exact(cycle2_onehot),'T0_multicontrol_cycle2_density':exact(cycle2_T0_multi),
        'basin_intervals':{key:{'lower':exact(value),'upper':exact(1-sum(v for k,v in lower.items() if k!=key))}
                           for key,value in lower.items()},
        'unassigned_mass':exact(1-sum(lower.values()))}
print(json.dumps(report,sort_keys=True,indent=2))
