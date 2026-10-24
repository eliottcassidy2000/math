#!/usr/bin/env python3
"""Hostile controls; no assert statements, hence also effective under -O."""
from copy import deepcopy
import json
from pathlib import Path
from check import verify_arithmetic, verify_cover, graph, require, lengths, PATTERN, mr_failure
ROOT=Path(__file__).resolve().parent
original=json.loads((ROOT/'inputs/prime_offsets_certificates.json').read_text())

def reject(name,operation):
    try:operation()
    except (ValueError,KeyError,TypeError,IndexError):print('REJECTED:',name);return
    raise RuntimeError('accepted hostile: '+name)

a=deepcopy(original);a['rows'][3]['plus']['offset']=11
reject('wrong Fortunate endpoint/minimality length',lambda:verify_arithmetic(a))
a=deepcopy(original);a['primality_certificates']['5']['base']=1
reject('invalid Lucas base',lambda:verify_arithmetic(a))
a=deepcopy(original);a['primality_certificates']['5']['factors']=[[2,1]]
reject('incomplete Lucas factorization',lambda:verify_arithmetic(a))
a=deepcopy(original);a['rows'][3]['plus']['witnesses'][9]={'prime':True}
reject('composite 221 passed as prime in the minimality interval',lambda:verify_arithmetic(a))
a=deepcopy(original);a['rows'][0]['P']=3
reject('shifted primorial indexing',lambda:verify_arithmetic(a))
reject('duplicated vertex in a spanning cover',lambda:verify_cover(graph(range(1,4)),[[1,3],[3]]))
reject('nonsquare adjacency in a cover',lambda:verify_cover(graph(range(1,4)),[[1,2,3]]))
reject('raw numerical gaps misread as the run-length word',lambda:require([0,0,0,2,10,0,0]==PATTERN,'raw-gap mismatch'))
reject('ordinary Hamiltonicity misread as the complete graph word',lambda:require(lengths([False]*5+[True]*3+[False]*5+[True,False,True,True])==PATTERN,'ordinary-path mismatch'))
reject('Miller-Rabin base equal to the tested prime',lambda:require(mr_failure(5,5),'invalid-base failure is not a compositeness witness'))
print('PASS: 10 hostile controls rejected')
