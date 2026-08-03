#!/usr/bin/env python3
"""Exact certificate and hostile controls for THM-3211.

The analytical proof is in the theorem file.  This companion checks the
Bernoulli formula on 27,342 exact lanes, exhausts the four small primitive
channels, and proves the three sharp-ray numerator identities by a pinned
independent THM-3171 affine-branch engine plus every finite head.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import gcd
import subprocess
from lrc_uniform_channel_limit_engine_thm3211 import (
    H, EDGES, L, mass, cleared_num, ray_limit, bernoulli_limit,
    exact_period_scout,
)

FLOOR=F(1,105)
records=[]
primitive=0
formula_checks=0
for Q in range(2,101):
    for P in range((Q+1)//2,Q):
        if gcd(P,Q)>1 or P+Q<8:
            continue
        primitive += 1
        for i,j in EDGES:
            for e,f in ((H[i],H[j]),(H[j],H[i])):
                got=ray_limit(e,P,f,Q)
                closed=bernoulli_limit(e,P,f,Q)
                if got!=closed:
                    raise RuntimeError(('limit formula',P,Q,e,f,got,closed))
                records.append((P,Q,e,f,closed))
                formula_checks += 1

small=[]
for P,Q in ((3,5),(4,5),(4,7),(5,6)):
    values=[]
    for i,j in EDGES:
        for e,f in ((H[i],H[j]),(H[j],H[i])):
            values.append((bernoulli_limit(e,P,f,Q),e,f))
    if any(x<FLOOR for x,e,f in values):
        raise RuntimeError(('small floor',P,Q,values))
    equal=tuple((e,f) for x,e,f in values if x==FLOOR)
    small.append((P,Q,min(x for x,e,f in values),equal))

expected_equal=((1,2),(3,1),(2,3))
if small[0][3]!=expected_equal or any(row[3] for row in small[1:]):
    raise RuntimeError(('equality classification',small))

engine_bytes=subprocess.check_output([
    'git','show',
    '75d0c078d2c2:04-computation/lrc_small_channel_symbolic_tail_thm3171.py'
])
engine_hash=sha256(engine_bytes).hexdigest()
if engine_hash!='d73273a4cf4b88bea2890e001166d96cb07dd9b61f3a248ff1538ec44579796a':
    raise RuntimeError(('engine hash',engine_hash))
engine={'__name__':'lrc3171_engine_for_independent_certificate'}
exec(compile(engine_bytes,'canonical_thm3171_engine','exec'),engine)

ray_specs=((0,1,2,4032,96),(1,3,1,4032,744),(5,2,3,4032,48))
ray_polys=[]
ray_stable_max=0
for edge,e,f,a,b in ray_specs:
    M=abs(5*e-3*f) or 1
    for residue in range(1,M+1):
        if not engine['ray_is_stable'](3,5,e,f,residue,M,1):
            raise RuntimeError(('ray not stable',edge,e,f,residue,M))
        ray_stable_max=max(ray_stable_max,residue+M)
        for j in range(4):
            g=residue+M*(1+j)
            if cleared_num(e,3*g,f,5*g)!=a*g*g+b*g:
                raise RuntimeError(('sharp tail polynomial',e,f,g))
    for g in range(1,ray_stable_max+1):
        if cleared_num(e,3*g,f,5*g)!=a*g*g+b*g:
            raise RuntimeError(('sharp head polynomial',e,f,g))
    ray_polys.append((e,f,a,b))

period_rows=[]
for Q in range(5,31):
    P=Q-1
    M,d,_=exact_period_scout(12,P,1,Q)
    period_rows.append((Q,M,d))

payload='\n'.join(map(str,records)).encode()
print('LRC REFLECTED-CHANNEL RAY-LIMIT BERNOULLI CERTIFICATE')
print(f'primitive_channels_Q_le_100={primitive}')
print(f'exact_limit_formula_checks={formula_checks}')
print(f'formula_record_sha256={sha256(payload).hexdigest()}')
print('coarse_tail=limit >= 1/49-1/(3PQ); PQ>=31 gives >=1/105+1/7595')
print(f'small_channel_minima={tuple(small)}')
print(f'equality_ray_polynomials={tuple(ray_polys)}')
print(f'equality_ray_branch_stability_max_g={ray_stable_max};engine_sha256={engine_hash}')
print('exception_hostile=(P,Q;e,f;g)=(3,5;6,1;2), mass=2030/280393, limit=17/1680')
print(f'period_scout_P_eq_Q_minus_1_e12_f1={tuple(period_rows)}')

print('scope=uniform_ray_limits_and_three_exact_sharp_rays_not_LRC14_or_finite_head_replacement')
print('FAILED_CHECKS=NONE')
