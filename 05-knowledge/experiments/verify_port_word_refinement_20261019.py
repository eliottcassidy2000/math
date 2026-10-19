#!/usr/bin/env python3
"""Independent orbit-symmetry refinement; leaves the base verifier unchanged."""
from functools import reduce
from math import gcd
from fractions import Fraction as F
import json
from verify_port_word_seams_20261019 import (
    ALPHA,BETA,GRAPHS,CORES,shuffle_kernel,profile_tables,scalar,
    assemble,explicit_extensions,profiles,
)


def primitive(v):
    d=reduce(gcd,v)
    assert d>0
    return tuple(x//d for x in v)


def refined_generators(kernels):
    left={(1,)*5};right={(1,)*4}
    for K in kernels:
        left.update(primitive(K[i]) for i in (1,2))
        left.add(primitive(tuple(a+b for a,b in zip(K[3],K[4]))))
        right.update(primitive(tuple(K[i][j] for i in range(1,5))) for j in (0,1,2))
        right.add(primitive(tuple(K[i][3]+K[i][4] for i in range(1,5))))
    return sorted(left),sorted(right)


def run():
    kernels=[shuffle_kernel(g) for g in GRAPHS]
    left,right=refined_generators(kernels)
    assert len(left)==8 and len(right)==6
    expected_left=[(1,1,1,1,1),(1,1,2,1,2),(1,2,1,2,1),(1,2,2,3,3),
                   (2,3,3,3,3),(2,3,3,4,4),(2,3,6,4,8),(2,6,3,8,4)]
    expected_right=[(1,1,1,1),(1,2,2,1),(1,3,3,1),
                    (2,1,1,2),(2,2,3,3),(3,6,8,4)]
    assert left==expected_left and right==expected_right
    reports={}
    for name,core in CORES.items():
        N,C=profile_tables(core)
        # Endpoint permutation involutions preserve this particular marked event.
        assert N[3]==N[4] and C[3]==C[4]
        assert all(row[2]==row[3] for row in N+C)
        checks=[]
        for i,l in enumerate(left):
            for j,r in enumerate(right):
                den=scalar(l,N,r);num=scalar(l,C,r)
                low=525*num-236*den;high=289*den-525*num
                assert den>0 and low>=0 and high>=0
                checks.append({'left_generator':i,'right_generator':j,
                    'numerator':num,'denominator':den,'probability':str(F(num,den)),
                    'lower_deficit_525C_minus_236N':low,
                    'upper_deficit_289N_minus_525C':high})
        probs=[F(t['numerator'],t['denominator']) for t in checks]
        expected=(F(1095,2372),F(600,1153)) if name=='B9' else (F(236,525),F(18,35))
        assert (min(probs),max(probs))==expected
        reports[name]={'probability_range':list(map(str,expected)),
             'minimum_lower_deficit':min(t['lower_deficit_525C_minus_236N'] for t in checks),
             'minimum_upper_deficit':min(t['upper_deficit_289N_minus_525C'] for t in checks),
             'all_48_generator_checks':checks}
    # Fixed direct involution controls on a full mixed exterior: actual extensions
    # are deliberately not enumerated here because the analytic swap is global.
    # Verify instead that both designated outgoing points remain maximal and both
    # incoming points remain minimal for two nontrivial seam choices.
    for names,graphs in [(('B8','B9'),(3,)),(('B9','B8'),(4,))]:
        p,offsets=assemble(names,graphs)
        last=CORES[names[-1]];first=CORES[names[0]]
        for u in (6,8):
            x=offsets[-1]+last['labels'].index(u)
            assert not any(ds>>x&1 for ds in p)
        for v in (0,3):
            x=first['labels'].index(v)
            assert p[x]==0
    assert F(236,525)==F(9,20)-F(1,2100)
    return {'status':'all 48 refined generator tests per core passed; actual-context sharpness not claimed',
            'alpha_profiles_nonempty':ALPHA[1:],'beta_profiles':BETA,
            'primitive_left_generators':left,'primitive_right_generators':right,
            'common_balance_bound':'236/525','improvement_over_base_bound':str(F(236,525)-F(24,55)),
            'distance_below_matching_bound_9_over_20':'1/2100','cores':reports,
            'scope':'Permutation-orbit outer cone for uniform exterior-poset laws, including whole prefixes/suffixes. No claim of exact realization of its individual rays.'}


if __name__=='__main__':
    print(json.dumps(run(),indent=2))
