"""Independent exhaustive check of the modular kernel engine used in exclusions."""
from itertools import product
from pathlib import Path
import json

from primary_native_clocks import modular_kernel


def run():
    records=[]
    for p,k in ((2,2),(3,2)):
        q=p**k
        cases=0
        for entries in product(range(q),repeat=4):
            A=[list(entries[:2]),list(entries[2:])]
            basis,vals=modular_kernel(A,p,k)
            actual={v for v in product(range(q),repeat=2)
                    if all(sum(x*y for x,y in zip(row,v))%q==0 for row in A)}
            generated={(0,0)}
            for b in basis:
                generated={tuple((v[i]+t*b[i])%q for i in range(2))
                           for v in generated for t in range(q)}
            predicted=p**sum(vals)*q**(2-len(vals))
            if actual!=generated or len(actual)!=predicted:
                raise ValueError(('kernel completeness failure',p,A))
            cases+=1
        records.append({'prime':p,'precision':k,'complete_2_by_2_matrices':cases,
                        'all_source_vectors_per_matrix':q*q})
    out={'status':'FINITE-EXACT independent full kernel enumeration', 'records':records}
    (Path(__file__).resolve().parent/'clock_kernel_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':run()
