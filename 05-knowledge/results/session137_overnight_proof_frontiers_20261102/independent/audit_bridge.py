"""Independent all-labelled order6 two-edge H25/H27 bridge census."""
from tournament_kernel import *
import json

def run():
    n=6;pairs=list(combinations(range(n),2));hs={};sources=[]
    for c in range(1<<15):
        a=matrix(n,c)
        if strong(a):
            h=hp(a);hs[c]=h
            if h==25:sources.append(c)
    counts={'bridges':0,'disjoint_edges':0,'strong_pivot_parent':0,'opposite_direction':0}
    witnesses=[]
    for c in sources:
        a=matrix(n,c)
        for i,j in combinations(range(15),2):
            d=c^(1<<i)^(1<<j)
            if hs.get(d)!=27:continue
            counts['bridges']+=1
            pivot=set(pairs[i])&set(pairs[j])
            if not pivot:counts['disjoint_edges']+=1;continue
            r=pivot.pop();v=next(x for x in pairs[i] if x!=r);w=next(x for x in pairs[j] if x!=r)
            counts['strong_pivot_parent']+=strong(induced(a,[x for x in range(n) if x!=r]))
            counts['opposite_direction']+=a[r][v]!=a[r][w]
            if len(witnesses)<3:witnesses.append({'H25_code':c,'H27_code':d,'pivot':r,'flipped_vertices':[v,w]})
    require(counts=={'bridges':1440,'disjoint_edges':0,'strong_pivot_parent':0,'opposite_direction':0},'bridge census')
    masks=[4,1,2,7,15];a=tuple(tuple((masks[i]>>j)&1 for j in range(5)) for i in range(5))
    require(hp(append(a,25))==25 and hp(append(a,16))==27,'explicit bypass')
    print(json.dumps({'scope':'FINITE-EXACT all labelled strong6 H25 to H27 pairs differing in exactly two arcs',
                      'counts':counts,'first_witnesses':witnesses,'explicit_base_masks':masks,
                      'explicit_cuts':[25,16],'explicit_H':[25,27]},sort_keys=True,indent=2))
if __name__=='__main__':run()
