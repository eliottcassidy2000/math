#!/usr/bin/env python3
"""Independent basis audit: fixed-Hamilton-cycle covers, recursive counts,
plus literal permutation verification of every published obstruction/card.
No import from either producer.
"""
import json,math
from pathlib import Path
from itertools import combinations,permutations,product
from functools import lru_cache
from collections import defaultdict,Counter

def require(p,msg):
    if not p:raise RuntimeError(msg)
def induced(a,vertices):
    return [[bool(a[u]>>v&1) for v in vertices] for u in vertices]
def direct(matrix):
    n=len(matrix)
    return sum(all(matrix[p[k]][p[k+1]] for k in range(n-1)) for p in permutations(range(n)))
def count(a):
    n=len(a);incoming=[sum(1<<u for u in range(n) if a[u]>>v&1) for v in range(n)]
    @lru_cache(None)
    def end(mask,v):
        rest=mask^(1<<v)
        if not rest:return 1
        options=rest&incoming[v]; answer=0
        while options:
            b=options&-options; options-=b;answer+=end(rest,b.bit_length()-1)
        return answer
    return sum(end((1<<n)-1,v) for v in range(n))
def remove(a,v):
    keep=[u for u in range(len(a)) if u!=v]
    return [sum(((a[u]>>w)&1)<<j for j,w in enumerate(keep)) for u in keep]
def normal(a):
    n=len(a);groups=defaultdict(list)
    for v in range(n):groups[a[v].bit_count()].append(v)
    best=None;edges=list(combinations(range(n),2))
    for parts in product(*(list(permutations(groups[d])) for d in sorted(groups))):
        p=sum(parts,()); code=0
        for k,(u,v) in enumerate(edges):code|=((a[p[u]]>>p[v])&1)<<k
        if best is None or code<best:best=code
    return best

def main():
    root=Path(__file__).resolve().parent
    data=json.loads((root/'basis21.json').read_text());records=data['strong_obstructions']+data['nonstrong_obstructions']
    direct_checks=0
    for rec in records:
        a=rec['rows'];n=len(a)
        got=direct(induced(a,list(range(n))))
        require(got==rec['H'] and got>21,'literal obstruction count fails');direct_checks+=1
        for v in range(n):
            h=direct(induced(a,[u for u in range(n) if u!=v]))
            require(h==rec['deletion_H'][v] and h<=21,'literal card count fails');direct_checks+=1
    cover=[]
    for n in (6,7):
        base=[0]*n;fixed=set()
        for u in range(n):
            v=(u+1)%n;base[u]|=1<<v;fixed.add(tuple(sorted((u,v))))
        free=[e for e in combinations(range(n),2) if e not in fixed]
        found={}; total=1<<len(free)
        for code in range(total):
            a=base[:]
            for k,(u,v) in enumerate(free):
                if code>>k&1:a[u]|=1<<v
                else:a[v]|=1<<u
            h=count(a)
            if h<=21:continue
            if any(count(remove(a,v))>21 for v in range(n)):continue
            found[normal(a)]=h
        expected={rec['degree_code']:rec['H'] for rec in data['strong_obstructions'] if rec['n']==n}
        require(found==expected,'fixed-Hamilton-cycle cover mismatch')
        cover.append({'n':n,'labelled_cycle_completions':total,'minimal_bad_classes':len(found)})
    # Two boundary mutations must fail, independently of producer assertions.
    row=records[0];bad=dict(row,H=row['H']+2)
    require(direct(induced(bad['rows'],list(range(bad['n']))))!=bad['H'],'mutated H accepted')
    padded=[x|(1<<len(row['rows'])) for x in row['rows']]+[0]
    require(count(padded)>21 and count(remove(padded,len(padded)-1))>21,'padded nonminimal hostile failed')
    print(json.dumps({'status':'PASS','threshold':21,'literal_permutation_checks':direct_checks,
        'independent_fixed_cycle_covers':cover,'mutations_rejected':['wrong H','bad tournament with singleton padding'],
        'basis_counts':data['counts']},indent=2,sort_keys=True))
if __name__=='__main__':main()
