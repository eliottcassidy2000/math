#!/usr/bin/env python3
"""Declared finite neighbourhood of E_n; not all sources or a readiness proof."""
import argparse,itertools,json
from pathlib import Path
from verify import en,profile
P=Path(__file__).resolve().parent
def strong(rows):
    full=(1<<len(rows))-1
    for q in range(len(rows)):
        seen=1<<q
        while True:
            nxt=seen
            for i in range(len(rows)):
                if seen>>i&1:nxt|=rows[i]
            if nxt==seen:break
            seen=nxt
        if seen!=full:return False
    return True
def main():
    a=argparse.ArgumentParser();a.add_argument('--distance',type=int,default=2);a.add_argument('--max-order',type=int,default=10);opts=a.parse_args()
    report=[]
    for n in range(7,opts.max_order+1):
        edges=list(itertools.combinations(range(n),2));base=en(n);found=[];tests=0
        for es in itertools.combinations(edges,opts.distance):
            r=base.copy()
            for i,j in es:r[i]^=1<<j;r[j]^=1<<i
            h,k=profile(r);tests+=1
            if h==67:found.append({'J':[z-h for z in k],'flips':es,'rows':r,'strong':strong(r)})
        report.append({'n':n,'flips_per_test':opts.distance,'tested':tests,'count67':found})
        print(n,tests,'hits',len(found),'least strong J',min((min(f['J']) for f in found if f['strong']),default=None),flush=True)
    (P/f'nearby67_d{opts.distance}.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
