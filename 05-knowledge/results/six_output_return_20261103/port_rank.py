"""Independent triangle matchings, collapsed walks and clone-run degree."""
import json,itertools
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok:raise ValueError(msg)
def tournament(n,mask):
    t=[0]*n;k=0
    for i in range(n):
      for j in range(i+1,n):
        if mask>>k&1:t[i]|=1<<j
        else:t[j]|=1<<i
        k+=1
    return tuple(t)
def triangle_matching(T,q):
    left=[v for v in range(len(T)) if T[q]>>v&1]
    right=[v for v in range(len(T)) if T[v]>>q&1];owner={}
    def augment(a,seen):
        for b in right:
            if T[a]>>b&1 and b not in seen:
                seen.add(b)
                if b not in owner or augment(owner[b],seen):owner[b]=a;return True
        return False
    for a in left:augment(a,set())
    return len(owner),sorted((a,b) for b,a in owner.items())
def walk_coefficients(T,q):
    outside=[v for v in range(len(T)) if v!=q];out=[0]*(len(T)+1)
    for perm in itertools.permutations(outside):
      for mask in range(1,1<<len(T)):
        path=[]
        for i in range(len(T)):
            if mask>>i&1:path.append(q)
            if i<len(perm):path.append(perm[i])
        if all(T[a]>>b&1 for a,b in zip(path,path[1:])):out[mask.bit_count()]+=1
    return out
def one_separator(T,q):
    for w in range(len(T)):
        if w==q:continue
        keep=((1<<len(T))-1)^(1<<w);reach=1<<q;old=-1
        while old!=reach:
            old=reach
            for v in range(len(T)):
                if reach>>v&1:reach|=T[v]&keep
        if not any((reach>>v&1) and (T[v]>>q&1) for v in range(len(T)) if v!=q):return w
    return None

def minimum_cycle_guard(T,q):
    """Independent deletion/reachability algorithm, with no triangle graph."""
    outside=[v for v in range(len(T)) if v!=q]
    for size in range(len(outside)+1):
      for deleted in itertools.combinations(outside,size):
        keep=((1<<len(T))-1)^sum(1<<v for v in deleted)
        reach=1<<q;old=-1
        while old!=reach:
            old=reach
            for v in range(len(T)):
                if reach>>v&1:reach|=T[v]&keep
        if not any(v!=q and (reach>>v&1) and (T[v]>>q&1) for v in range(len(T))):
            return deleted
    raise ValueError('deleting the full exterior always removes every q-cycle')
def strong(T):
    for q in range(len(T)):
        seen=1<<q;old=-1
        while old!=seen:
            old=seen
            for v in range(len(T)):
                if seen>>v&1:seen|=T[v]
        if seen!=(1<<len(T))-1:return False
    return True
def run():
    hist=Counter();checks=0;expanded=None
    for n in range(1,6):
      for mask in range(1<<(n*(n-1)//2)):
        T=tournament(n,mask)
        for q in range(n):
            rank,matching=triangle_matching(T,q);c=walk_coefficients(T,q)
            actual=max(i for i,x in enumerate(c) if x)
            need(actual==rank+1,'matching exactly controls largest run coefficient')
            guard=minimum_cycle_guard(T,q)
            need(len(guard)==rank,'minimum exterior cycle deletion equals matching rank')
            if rank==1:
                w=one_separator(T,q);need(w is not None,'one cycle guardian exists')
                if expanded is None and strong(T) and T[q].bit_count()>1 and len(T)-1-T[q].bit_count()>1:
                    expanded=dict(T=T,q=q,guardian=w,matching=matching,coefficients=c)
            checks+=1;hist[n,rank]+=1
    inputs=[]
    for d in json.loads((HERE/'sources.json').read_text()):
        T=tuple(d['target']);q=d['port'];rank,matching=triangle_matching(T,q);w=one_separator(T,q)
        need(strong(T) and rank==1 and w is not None,'all six actual outputs have the common input type')
        inputs.append(dict(head=3*d['source_H']+d['offset'],port=q,guardian=w,degree=T[q].bit_count(),matching=matching))
    out=dict(checks=checks,minimum_cycle_guard_checks=checks,histogram=[dict(n=n,rank=r,count=v) for (n,r),v in sorted(hist.items())],expanded_port=expanded,six_inputs=inputs)
    (HERE/'port_rank.json').write_text(json.dumps(out,indent=2)+'\n');print('port rank',checks,expanded,flush=True);return out
if __name__=='__main__':run()
