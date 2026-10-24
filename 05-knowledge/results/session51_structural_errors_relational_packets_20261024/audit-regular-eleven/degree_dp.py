#!/usr/bin/env python3
"""Independent degree coefficient recursion using literal indexed subsets.
No equal-degree binomial compression or Landau oracle; delete a maximum-degree
vertex and choose its in-neighbors, decrementing exactly those residuals.
"""
import argparse, functools, itertools, json, math
@functools.lru_cache(None)
def coefficient(d):
    n=len(d)
    if not n:return 1
    if d[0]<0 or d[-1]>=n or sum(d)!=n*(n-1)//2:return 0
    remaining=d[:-1]; losses=n-1-d[-1]
    result=0
    for chosen in itertools.combinations(range(n-1),losses):
        next_d=list(remaining)
        for j in chosen:next_d[j]-=1
        result+=coefficient(tuple(sorted(next_d)))
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    rows=[]
    for n in (3,5,7,9,11):
        h=(n-1)//2
        fixed=coefficient(tuple([h-1]*h+[h]*h))
        almost=math.comb(n-1,h)*fixed
        regular=coefficient(tuple([h]*n))
        if almost!=regular:raise RuntimeError('independent labeled bijection fails')
        rows.append(dict(regular_order=n,fixed_assignment_almost=fixed,all_labelled_almost=almost,all_labelled_regular=regular))
    # Entire small labeled universes, independent of the recurrence.
    brute=[]
    for n in (3,5):
        edges=list(itertools.combinations(range(n),2));hist={}
        for mask in range(1<<len(edges)):
            d=[0]*n
            for b,(i,j) in enumerate(edges):d[i if (mask>>b)&1 else j]+=1
            key=tuple(d);hist[key]=hist.get(key,0)+1
        for key,value in hist.items():
            if coefficient(tuple(sorted(key)))!=value:raise RuntimeError('literal labeled degree histogram differs')
        brute.append(dict(order=n,labelled_tournaments=1<<len(edges),distinct_labelled_degree_vectors=len(hist)))
    out=dict(status='PASS',method='maximum-degree deletion; literal indexed incoming-neighbor subsets',rows=rows,brute_controls=brute,cache=coefficient.cache_info()._asdict())
    with open(a.output,'w') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
