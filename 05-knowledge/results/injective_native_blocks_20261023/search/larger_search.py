"""Bounded deterministic search for larger primitive source-charge circuits."""
import json,random,time,sys
from collections import defaultdict
from itertools import combinations
import explore as e
SOURCES=[176,1394,1425,1457,1477,1480,1532,1549,1638,1691,1841,1921,1923,1974,2071,2087,2176,2185,2553,2581,2835,2837,2905,2945,3000,3015,3360,3369]
TARGETS=(1965,2234)
BIT={s:sum((e.Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in SOURCES}
def key(ids,mod=False):
    x=sum(e.X[s] for s in ids);y=0
    for s in ids:y^=BIT[s]
    return (x%128 if mod else x,y)
def primitive(a,b):
    for k in range(1,len(a)):
        left={key(ids,True) for ids in combinations(a,k)}
        if any(key(ids,True) in left for ids in combinations(b,k)):return False
    return True
def main():
    rng=random.Random(20261023);start=time.time();found=[]
    sizes=(8,) if '16' in sys.argv else (4,5,6)
    for k in sizes:
        buckets=defaultdict(list);seen=set();tested=0;answer=None
        for sample in range(150000):
            if time.time()-start>55:break
            a=tuple(sorted(rng.sample(SOURCES,k)))
            if a in seen:continue
            seen.add(a);ak=key(a)
            for b in buckets[ak]:
                if set(a)&set(b):continue
                tested+=1
                if primitive(a,b):answer=(a,b);break
            if answer:break
            buckets[ak].append(a)
        if answer:
            aa,bb=answer;slots=[]
            for s,sgn in [(s,1) for s in aa]+[(s,-1) for s in bb]:
                rows=e.paths(s);cost,pins,i,j=e.menu(rows)[TARGETS]
                if sgn<0:i,j=j,i
                slots.append(dict(source=s,sign=sgn,canon=int(e.NODES[s]['canon_a']),cost=cost,pins=pins,
                                  before=rows[i],after=rows[j],all_unchanged_clocks=[[h,u] for h,(u,v) in enumerate(zip(rows[i]['path'],rows[j]['path'])) if u==v]))
            found.append(dict(size=2*k,positive=aa,negative=bb,positive_X=[e.X[s] for s in aa],negative_X=[e.X[s] for s in bb],
                              exact_X_sum=sum(e.X[s] for s in aa),sampled_subsets=len(seen),balanced_disjoint_candidates_tested=tested,slots=slots))
            print('FOUND',2*k,aa,bb,'sum',sum(e.X[s] for s in aa),'sampled',len(seen),'tested',tested,flush=True)
        else:print('NOT FOUND',2*k,'sampled',len(seen),'tested',tested,flush=True)
    out={'seed':20261023,'sample_cap_per_size':150000,'source_universe':SOURCES,'targets':TARGETS,'circuits':found}
    (e.HERE/('larger16_selected.json' if '16' in sys.argv else 'larger_selected.json')).write_text(json.dumps(out,indent=2));print('elapsed',time.time()-start)
if __name__=='__main__':main()
