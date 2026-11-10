#!/usr/bin/env python3
"""Independent subset-cover / whole-tournament DP for sharp cyclic inflation gaps."""
import json
from collections import Counter
from itertools import permutations
from verify import fixed_masks, rows, path_count, require, scc_count


def restrict(out, subset):
    vs=[i for i in range(len(out)) if (subset>>i)&1]
    return [sum(((out[a]>>b)&1)<<j for j,b in enumerate(vs)) for a in vs]


def inflate(out):
    n=len(out); result=[r|(1<<n) for r in out]
    result += [1<<(n+1), (1<<n)-1]
    return result


def two_covers(out):
    full=(1<<len(out))-1
    return sum(path_count(restrict(out,s))*path_count(restrict(out,full^s)) for s in range(1,full))


def permutation_paths(out):
    return sum(all((out[a]>>b)&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(out))))

result={}
for p in (3,5,7):
    cases=[]
    for mask in fixed_masks((p,)):
        out=rows(p,mask); h=path_count(out); f2=two_covers(out); big=inflate(out)
        direct=path_count(big)
        require(direct==3*h+f2,('ordered-cover formula',p,mask))
        require(scc_count(big)==1,'sharp witness strong')
        require(permutation_paths(out)==h,'independent permutation path count')
        cases.append({'mask':mask,'H':h,'F2':f2,'inflation_H':direct})
    mu=min(c['H'] for c in cases); gamma=min(c['inflation_H'] for c in cases)
    result[p]={'mu':mu,'gamma':gamma,'multi_orbit_lower':mu*mu+2*p,'nontransitive_action_lower':min(gamma,mu*mu+2*p),'cases':cases}
require(result[7]['gamma']==2163 and result[7]['multi_orbit_lower']>2163,'sharp7 gap')
# Second independent whole-witness check (9! permutations; no subset DP).
witness=min(result[7]['cases'],key=lambda c:c['inflation_H'])
require(permutation_paths(inflate(rows(7,witness['mask'])))==2163,'independent9vertex brute path check')
print(json.dumps({'status':'PASS','prime_action_gap':result},indent=2,sort_keys=True))
