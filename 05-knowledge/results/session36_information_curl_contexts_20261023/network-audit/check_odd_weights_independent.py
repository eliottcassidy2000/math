#!/usr/bin/env python3
"""Full dyadic profile for every odd global color-weight ratio, finite n<=7."""
import json
from pathlib import Path
from collections import Counter
from check_metagraph_independent import check,laplacian
from symmetric_dyadic import valuations
HERE=Path(__file__).resolve().parent

def run():
    source=json.loads((HERE.parent/'metagraph/dyadic_network_receipt.json').read_text())
    expected=json.loads((HERE.parent/'metagraph/dyadic_lift_receipt.json').read_text());out=[]
    for s,e in zip(source,expected):
        n=s['nodes'];B=Counter({(u,v):w for u,v,w in s['blue_edges']});K=Counter({(u,v):w for u,v,w in s['paired_black_edges']});L=laplacian(n,B);M=laplacian(n,K)
        tests=[]
        for t in (1,3,5,7):
            A=[[t*L[i][j]+2*M[i][j] for j in range(n-1)] for i in range(n-1)]
            vals=valuations(A,precision=3)
            hist={str(k):vals.count(k) for k in sorted(set(vals))}
            check(hist==e['valuation_histogram'],'odd ratio dyadic profile changed')
            tests.append(dict(ratio_mod8=t,histogram=hist,precision_bits=3,precision_exhausted=False))
        out.append(dict(n=s['n'],residue_tests=tests))
    ans=dict(status='PASS',scope='Actual full networks n=3..7; odd positive global blue/literal-black conductances only',rows=out)
    (HERE/'odd_weights_independent_receipt.json').write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2));print('ALL INDEPENDENT ODD-WEIGHT CHECKS PASSED')
if __name__=='__main__':run()
