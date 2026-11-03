#!/usr/bin/env python3
"""Keep the exact named next-operation boundary, not just its scalar H."""
from pathlib import Path
from math import comb,factorial
import json
from scout import endpoint,clone,require

HERE=Path(__file__).resolve().parent


def defect(rows):
    d=endpoint(tuple(rows));h=[sum(x)for x in d];h[0]=1;full=len(d)-1
    return [sum(d[s][v]*h[full^s]for s in range(1,full+1))-h[full]for v in range(len(rows))]


if __name__=='__main__':
    sources=[r['target']for r in json.loads((HERE/'returns.json').read_text())]
    sources.append(json.loads((HERE/'next_0.json').read_text())['programmes'][0]['target'])
    result=[]
    for t in sources:
        v=3;observations=[]
        for m in range(1,6):
            d=defect(clone(tuple(t),v,m));observations.append({q:d[q+(m-1 if q>v else 0)]for q in range(len(t))if q!=v})
        cs={}
        for q in observations[0]:
            coefficients=[]
            for m in range(1,5):
                lower=sum(c*sum((-1)**(r-j)*comb(r,j)*j**m for j in range(r+1))for r,c in enumerate(coefficients,1))
                num=observations[m-1][q]-lower
                require(num%factorial(m)==0,'integral ordered-cover coefficient')
                coefficients.append(num//factorial(m))
            want=sum(c*sum((-1)**(r-j)*comb(r,j)*j**5 for j in range(r+1))for r,c in enumerate(coefficients,1))
            require(want==observations[4][q],'extra boundary sample')
            cs[q]=coefficients
        result.append({'source':t,'H':sum(endpoint(tuple(t))[-1]),'clone':v,
                       'exterior_defect_run_coefficients':cs,
                       'constant_exterior_defects':{q:c[0]for q,c in cs.items()if c[1:]==[0,0,0]},
                       'values':observations})
    (HERE/'boundary.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in result:print(r['H'],r['constant_exterior_defects'])
