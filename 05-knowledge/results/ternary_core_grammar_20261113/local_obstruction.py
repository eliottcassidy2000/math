"""Literal local certificate for 21 and the merge/orphan decomposition."""
from itertools import permutations
from collections import Counter
import numpy as np
from triangle_probe import attach,ends,need,save

def paths(T):
    return [p for p in permutations(range(len(T))) if all(T[a]>>b&1 for a,b in zip(p,p[1:]))]

def valid(T,p):return all(T[a]>>b&1 for a,b in zip(p,p[1:]))

def cycle_count(T):
    # Independent odd-cycle collection formula, at orders <=6.
    n=len(T);cycles=[]
    for size in (3,5):
        if size>n:continue
        for p in permutations(range(n),size):
            if p[0]!=min(p):continue
            if valid(T,p) and T[p[-1]]>>p[0]&1:cycles.append(sum(1<<v for v in p))
    disjoint=sum(a&b==0 for i,a in enumerate(cycles) for b in cycles[i+1:])
    return 1+2*len(cycles)+4*disjoint

def run():
    T=np.array([6,12,8,1],dtype=np.int64);tables=[];breakdowns=[];total=0
    for r in range(4):
        table=[]
        for first in range(8):
            row=[]
            for second in range(8):
                mask=first+(second<<3);U=tuple(map(int,attach(T,r,mask)))
                ps=paths(U);h=len(ps)
                need(h==int(ends(np.array(U)).sum())==cycle_count(U),'three independent counts')
                merge=Counter();orphan=Counter()
                for p in ps:
                    a=tuple(v for v in p if v<4);b=tuple(v for v in p if v in (r,4,5))
                    va=valid(U,a);vb=valid(U,b)
                    if va and vb:merge[a,b]+=1
                    else:orphan[('A' if not va else '')+('B' if not vb else '')]+=1
                need(len(merge)==15,'all 5x3 pairs have a merge')
                extra=sum(c-1 for c in merge.values());lost=sum(orphan.values())
                need(h==15+extra+lost,'exact deficit')
                breakdowns.append(dict(root=r,mask=mask,H=h,surplus=extra,orphans=dict(orphan)))
                row.append(h);total+=1
            table.append(row)
        tables.append(dict(root=r,table=table,values=sorted({v for row in table for v in row})))
    small=[];C=np.array([2,4,1],dtype=np.int64)
    for mask in range(16):
        U=tuple(map(int,attach(C,0,mask)));p=paths(U)
        need(len(p)==cycle_count(U)==int(ends(np.array(U)).sum()),'triangle/triangle counts')
        small.append(len(p))
    need(all(21 not in t['values'] for t in tables),'local omission21')
    save('local_obstruction',dict(scope='All 256 rooted C4/triangle cross matrices and all16 C3/triangle matrices; literal paths, endpoint DP and OCF agree',tables=tables,triangle_triangle=small,breakdowns=breakdowns))
    print('C3/C3',sorted(set(small)))
    for t in tables:print('ROOT',t['root'],t['values']);print(*t['table'],sep='\n')

if __name__=='__main__':run()
