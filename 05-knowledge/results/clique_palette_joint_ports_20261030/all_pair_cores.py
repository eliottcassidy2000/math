"""Apply the incoming two-tail conductor to every eleven-subset of 1..13."""
from fractions import Fraction as F
from itertools import combinations
from math import ceil, floor
from pathlib import Path
import argparse
import json
import lrc_banks as native

HERE=Path(__file__).resolve().parent

def one_tail(intervals, radius, speed):
    for left,right in intervals:
        lo=ceil(speed*left+radius-1)
        hi=floor(speed*right-radius)
        if lo<=hi:
            return max(left,(lo+radius)/speed)
    return None

def run():
    radius=F(1,14); records=[]; queries=0; references=0
    for core in combinations(range(1,14),11):
        intervals=native.safe_set(core,radius)
        width=max((b-a for a,b in intervals),default=F(0))
        if width<=0: raise RuntimeError('No positive core width')
        cutoff=ceil(1/width); branches=[]
        for u in range(1,cutoff):
            joint=native.intersection(intervals,native.safe_set((u,),radius))
            w2=max((b-a for a,b in joint),default=F(0))
            if w2<=0: raise RuntimeError(('Point/empty branch',core,u))
            tail=ceil(2*radius/w2); rejected=[]
            for v in range(u,tail):
                t=one_tail(joint,radius,v)
                if t is None: rejected.append(v)
                elif not all(min((s*t)%1,1-(s*t)%1)>=radius for s in (*core,u,v)):
                    raise RuntimeError('Invalid native common-time witness')
                # Independent full intersection at every finite-head query.
                truth=native.intersection(joint,native.safe_set((v,),radius))
                if bool(truth)!=(t is not None): raise RuntimeError('Literal-comb disagreement')
                queries+=1; references+=1
            branches.append({'smaller_speed':u,'larger_tail':tail,'rejected_larger_speeds':rejected})
        records.append({'core':list(core),'both_large_cutoff':cutoff,'branches':branches,
                        'all_pairs_accepted':all(not b['rejected_larger_speeds'] for b in branches)})
    return {'families':records,'core_count':len(records),'finite_queries':queries,
            'independent_literal_comparisons':references,
            'all_pairs_accepted':all(x['all_pairs_accepted'] for x in records),
            'maximum_both_large_cutoff':max(x['both_large_cutoff'] for x in records),
            'scope':'Incoming two-tail theorem plus independently computed finite certificates for 78 cores.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path); args=ap.parse_args()
    result=run(); text=json.dumps(result,separators=(',',':'),sort_keys=True)+'\n'
    if args.output: args.output.write_text(text,encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k!='families'},sort_keys=True))
