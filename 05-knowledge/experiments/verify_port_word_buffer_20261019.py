#!/usr/bin/env python3
"""One-core buffer certificate only; no iterations or word census."""
from fractions import Fraction as F
from functools import reduce
from math import gcd
import json
from verify_port_word_seams_20261019 import (
    CORES,GRAPHS,shuffle_kernel,profile_tables,matmul,scalar,
)


def primitive(v):
    d=reduce(gcd,v)
    assert d
    return tuple(x//d for x in v)


def run():
    kernels=[shuffle_kernel(g)[1:] for g in GRAPHS]  # four suffixes by five prefixes
    left={(1,)*5};right={(1,)*4};lsource={};rsource={}
    for name,core in CORES.items():
        core['N'],core['C']=profile_tables(core)
        for g,K in enumerate(kernels):
            NK=matmul(core['N'],K);KN=matmul(K,core['N'])
            for i,row in enumerate(NK):
                if any(row):
                    ray=primitive(row);left.add(ray)
                    lsource.setdefault(ray,[]).append({'core':name,'seam_index':g,'row':i})
            for j in range(4):
                col=tuple(row[j] for row in KN)
                if any(col):
                    ray=primitive(col);right.add(ray)
                    rsource.setdefault(ray,[]).append({'core':name,'seam_index':g,'column':j})
    left=sorted(left);right=sorted(right)
    assert (len(left),len(right))==(36,15)
    result={}
    for name,core in CORES.items():
        checks=[]
        for i,l in enumerate(left):
            for j,r in enumerate(right):
                den=scalar(l,core['N'],r);num=scalar(l,core['C'],r)
                lower=12875*num-5872*den;upper=7003*den-12875*num
                assert den>0 and lower>=0 and upper>=0
                checks.append({'left_ray':i,'right_ray':j,'numerator':num,'denominator':den,
                    'probability':str(F(num,den)),
                    'lower_deficit_12875C_minus_5872N':lower,
                    'upper_deficit_7003N_minus_12875C':upper})
        lo=min(checks,key=lambda x:F(x['numerator'],x['denominator']))
        hi=max(checks,key=lambda x:F(x['numerator'],x['denominator']))
        expected=(F(759,1622),F(33199,64887)) if name=='B9' else (F(5872,12875),F(58366,115755))
        assert (F(lo['probability']),F(hi['probability']))==expected
        result[name]={'probability_range':list(map(str,expected)),
            'minimizing_generator_pair':{'left_index':lo['left_ray'],'right_index':lo['right_ray']},
            'maximizing_generator_pair':{'left_index':hi['left_ray'],'right_index':hi['right_ray']},
            'minimum_lower_deficit':min(x['lower_deficit_12875C_minus_5872N'] for x in checks),
            'minimum_upper_deficit':min(x['upper_deficit_7003N_minus_12875C'] for x in checks),
            'all_540_generator_checks':checks}
    beta=F(5872,12875)
    assert beta>F(9,20)>F(236,525)>F(24,55)
    return {'status':'all 540 one-buffer generator tests per core passed; no additional buffer iteration',
            'left_ray_count':len(left),'right_ray_count':len(right),
            'primitive_left_rays':left,'primitive_right_rays':right,
            'left_ray_sources':[lsource.get(x,[{'boundary':'empty'}]) for x in left],
            'right_ray_sources':[rsource.get(x,[{'boundary':'empty'}]) for x in right],
            'common_balance_bound':str(beta),'improvement_over_9_over_20':str(beta-F(9,20)),
            'cores':result,
            'scope':'Exterior is empty or contains a full adjacent B8/B9 core. Arbitrary earlier/later nonnegative profile weights are allowed. Cone sharpness is not actual-context sharpness.'}


if __name__=='__main__':
    print(json.dumps(run(),indent=2))
