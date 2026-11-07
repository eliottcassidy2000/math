#!/usr/bin/env python3
"""Independent square-sum path census and endpoint-sensitive onset certificates.
No external modules, no inherited solver, and no unsound degree-two pruning.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations
from math import isqrt
import json
from pathlib import Path

def sq_graph(n):
    return tuple(tuple(j for j in range(1,n+1) if j!=i and isqrt(i+j)**2==i+j)
                 for i in range(1,n+1))

def paths(n):
    a=sq_graph(n)
    leaves=[i+1 for i,x in enumerate(a) if len(x)==1]
    # A specified leaf is an endpoint of every spanning path; fixing it
    # as the first vertex counts each reversal orbit exactly once.
    starts=leaves[:1] if leaves else range(1,n+1)
    found=[]
    def dfs(p,used):
        if len(p)==n:
            if leaves or n==1 or p[0]<p[-1]: found.append(tuple(p))
            return
        for v in a[p[-1]-1]:
            if not (used>>(v-1)&1): dfs(p+[v],used|(1<<(v-1)))
    for s in starts: dfs([s],1<<(s-1))
    return found

def count_dp(n):
    a=sq_graph(n)
    leaf=next((i+1 for i,x in enumerate(a) if len(x)==1),None)
    @lru_cache(None)
    def count(mask,last):
        if mask==0:return 1
        return sum(count(mask^(1<<(v-1)),v) for v in a[last-1] if mask>>(v-1)&1)
    allbits=(1<<n)-1
    if leaf:return count(allbits^(1<<(leaf-1)),leaf)
    directed=sum(count(allbits^(1<<(v-1)),v) for v in range(1,n+1))
    return directed if n==1 else directed//2

def two_port_forests(n,u,v):
    """Count spanning two-path forests with u,v in distinct components.
    Each component has its designated port at an endpoint. Singletons count once.
    Independent of a Hamilton-path-on-Q_(n+1) search.
    """
    a=sq_graph(n)
    def prefix_counts(root,excluded):
        counts=Counter()
        def extend(last,mask):
            counts[mask]+=1
            for w in a[last-1]:
                if w!=excluded and not(mask>>(w-1)&1):
                    extend(w,mask|(1<<(w-1)))
        extend(root,1<<(root-1))
        return counts
    left=prefix_counts(u,v);right=prefix_counts(v,u);full=(1<<n)-1
    return sum(c*right.get(full^mask,0) for mask,c in left.items())

def need(cond,msg):
    if not cond: raise RuntimeError(msg)

def main():
    expected={1:1,15:1,16:1,17:1,23:3,25:10}
    rows=[]
    for n in range(1,26):
        pp=paths(n)
        need(len(pp)==expected.get(n,0),f'census n={n}')
        need(count_dp(n)==len(pp),f'independent DP n={n}')
        for p in pp:
            need(sorted(p)==list(range(1,n+1)),f'permutation n={n}')
            need(all(isqrt(x+y)**2==x+y for x,y in zip(p,p[1:])),f'edge n={n}')
        rows.append({'n':n,'unoriented_paths':len(pp),'endpoint_pairs':
                     [[list(k),v] for k,v in sorted(Counter(tuple(sorted((p[0],p[-1]))) for p in pp).items())]})
    onsets=[]
    for n in (15,23,25):
        pp=paths(n); neighbors=sq_graph(n)[n-1]
        need(len(neighbors)==2,'onset degree two')
        need(two_port_forests(n-1,*neighbors)==len(pp),'two-port forest check')
        forests=[]
        for p in pp:
            k=p.index(n)
            need(0<k<n-1,'onset must be interior')
            l,r=p[:k],p[k+1:]
            need({l[-1],r[0]}==set(neighbors),'designated ports')
            forests.append([list(l),list(r)])
        onsets.append({'n':n,'ports':neighbors,'count':len(pp),'paths':[list(p) for p in pp],
                       'two_path_forests_after_deletion':forests})
    # The local count is not determined by the old scalar p=0 and birth degree2:
    # all these vertex additions have those same inputs, but differing outputs.
    comparison=[]
    for n in (13,14,15,19,20,21,22,23,25):
        need(not paths(n-1),f'old count n={n}')
        need(len(sq_graph(n)[n-1])==2,f'birth degree n={n}')
        need(two_port_forests(n-1,*sq_graph(n)[n-1])==len(paths(n)),'hostile forest control')
        comparison.append([n,len(paths(n))])
    out={'status':'PASS','universe':'every Q_n for 1<=n<=25; all paths up to reversal',
         'algorithms':['simple-path DFS with justified leaf-root symmetry reduction','memoized subset endpoint recurrence','independent designated-port two-path forest convolution'],
         'census':rows,'onsets':onsets,'zero_count_degree_two_birth_controls':comparison,
         'no_affine_transport_witnesses':[{'n':23,'p':3,'2p+1':7},{'n':25,'p':10,'2p+1':21}]}
    target=Path(__file__).with_name('RESULT.json');target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','counts':[[r['n'],r['unoriented_paths']] for r in rows],
                     'onset_ports':[[r['n'],r['ports'],r['count']] for r in onsets],
                     'same_scalar_birth_controls':comparison},indent=2))
if __name__=='__main__':main()
