#!/usr/bin/env python3
"""Independent cap BFS and physical subset-DP audit. Uses inherited profile
matrix generator as theorem input, not the construction author's BFS."""
import json,collections,time
from functools import lru_cache
import ribbon_inherited as ri

def independent_h(adj):
    n=len(adj)
    incoming=[sum(1<<v for v in range(n) if adj[v]>>j&1) for j in range(n)]
    @lru_cache(None)
    def f(mask,last):
        rest=mask^(1<<last)
        if not rest:return 1
        prev=rest&incoming[last];total=0
        while prev:
            bit=prev&-prev;prev-=bit;total+=f(rest,bit.bit_length()-1)
        return total
    return sum(f((1<<n)-1,v) for v in range(n))

def bfs(w, cap, budget=200000):
    order,idx,mats,final=ri.transfer_system(w)
    alphabet=ri.letters(w)
    updates=[]
    for a in alphabet:
        row=[[] for _ in order]
        for (dest,src),mul in mats[a].items(): row[dest].append((src,mul))
        updates.append(row)
    start=(1,)+(0,)*(len(order)-1)
    seen={start:()}; queue=collections.deque([start]); hits={};edges=0
    while queue:
        v=queue.popleft()
        for ai,row in enumerate(updates):
            vv=tuple(min(cap+1,sum(m*v[s] for s,m in es)) for es in row)
            edges+=1
            h=min(cap+1,sum(vv[i] for i in final))
            if h>cap:continue
            word=seen[v]+(ai,)
            if h not in hits:hits[h]=word
            if vv not in seen:
                seen[vv]=word;queue.append(vv)
                if len(seen)>budget:raise RuntimeError('budget')
    for h,word in hits.items():
        adj=ri.tournament_from_word(w,[alphabet[i] for i in word])
        if independent_h(adj)!=h:raise RuntimeError(('physical',w,cap,h,word))
    return {'span':w,'cap':cap,'profile_dimension':len(order),'reachable_pruned_vectors':len(seen),
        'edges_checked':edges,'max_representative_length':max(map(len,seen.values())),
        'values':sorted(hits),'witness_lengths':{h:len(x) for h,x in hits.items()},
        'witness_47':list(hits[47]) if 47 in hits else None,'alphabet':[sorted(a) for a in alphabet]}

if __name__=='__main__':
    out=[]
    for w,c in [(1,99),(2,99),(3,47)]:
        t=time.monotonic();v=bfs(w,c);v['seconds']=round(time.monotonic()-t,3)
        print(json.dumps(v),flush=True);out.append(v)
    open('cap_receipt.json','w').write(json.dumps(out,indent=2)+'\n')
