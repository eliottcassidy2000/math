#!/usr/bin/env python3
"""Provable rational basin-density intervals. Polynomial, but NOT fast in practice.

--bits m uses cap=m+4, (m+5)^10 states. Use --cap 0/1 for small smoke tests.
The decimal display is absent: integer numerators/denominator are the result.
"""
from itertools import product
from math import prod
import argparse,json
from primegame_basin import PRIMES,classify_vector,require


def approximate(cap):
    if not isinstance(cap,int) or cap<0: raise ValueError('nonnegative cap required')
    denom=prod(p**(cap+1) for p in PRIMES)
    weights=tuple(tuple((p-1)*p**(cap-e) for e in range(cap+1)) for p in PRIMES)
    sums={'prime':0,'cycle2':0,'cycle6':0}; counts=dict.fromkeys(sums,0)
    for x in product(range(cap+1),repeat=10):
        label=classify_vector(x)['basin']
        sums[label]+=prod(weights[i][e] for i,e in enumerate(x)); counts[label]+=1
    retained=prod(p**(cap+1)-1 for p in PRIMES)
    require(sum(sums.values())==retained,'valuation-cylinder mass identity')
    error=denom-retained
    return dict(cap=cap,states=(cap+1)**10,common_denominator=denom,
                error_numerator=error,census=counts,
                intervals={k:{'lower_numerator':v,'upper_numerator':v+error} for k,v in sums.items()})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--cap',type=int)
    group.add_argument('--bits',type=int)
    a=parser.parse_args()
    if a.bits is not None and a.bits<0: parser.error('--bits must be nonnegative')
    result=approximate(a.cap if a.cap is not None else a.bits+4)
    if a.bits is not None:
        require(result['error_numerator']*2**a.bits<result['common_denominator'],'requested precision')
        result['requested_bits']=a.bits
    print(json.dumps(result,sort_keys=True,indent=2))
