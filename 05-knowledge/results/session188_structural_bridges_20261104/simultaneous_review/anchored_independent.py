#!/usr/bin/env python3
from itertools import combinations,permutations,product
from functools import lru_cache
from pathlib import Path
import json
from independent import polynomial_bank,rows,modules,flip,category,need


def truth_bank(n,B,bank):
    output=[]
    for vector in bank:
        a=[0]*(1<<(2*n))
        for s,c in zip(B,vector):a[s]=c
        for bit in range(2*n):
            half=1<<bit
            for start in range(0,len(a),2*half):
                for i in range(half):a[start+half+i]+=a[start+i]
        output.append(a)
    return output


def submasks(s):
    t=s
    while True:
        yield t
        if not t:break
        t=(t-1)&s


def main():
    root=Path(__file__).resolve().parent.parent
    published=json.loads((root/'simultaneous_response'/'anchored.json').read_text())
    records={}
    for r in published['records']:
        key=(r['n'],r['source_code'],tuple(r['module_masks']),r['chosen_outer_module'],r['anchor'])
        need(key not in records,'duplicate anchor record')
        records[key]=r
    families=anchors=covered=average=boundary=0;seen=set()
    for n in range(1,6):
        B,bank=polynomial_bank(n);V=truth_bank(n,B,bank);full=(1<<n)-1
        def value(code,a,c):return V[code][a|(c<<n)]
        for code in range(len(bank)):
            adj=rows(n,code);mods=modules(adj);flips={M:flip(n,M) for M in mods}
            @lru_cache(None)
            def endpoint(M,v):
                names=[u for u in range(n) if M>>u&1]
                count=0
                for word in permutations(names):
                    if (word[0]==v or word[-1]==v) and all(adj[u]>>z&1 for u,z in zip(word,word[1:])):
                        count+=int(word[0]==v)+int(word[-1]==v)
                return count%2
            for mask in range(1,1<<len(mods)):
                chosen=[M for j,M in enumerate(mods) if mask>>j&1]
                if category(chosen)=='overlapping':continue
                maxima=[M for M in chosen if not any(M!=N and M&N==M for N in chosen)]
                exposed=[]
                for M in maxima:
                    hidden=0
                    for N in chosen:
                        if N!=M and N&M==N:hidden|=N
                    for v in range(n):
                        if M>>v&1 and not hidden>>v&1:exposed.append((M,v))
                if not exposed:covered+=1;continue
                families+=1;target=code
                for M in chosen:target^=flips[M]
                for M,v in exposed:
                    names=[u for u in range(n) if M>>u&1];m=len(names)
                    ext=sum(1<<u for u in range(n) if not M>>u&1 and adj[u]>>names[0]&1)
                    d=sum(adj[v]>>u&1 for u in names if u!=v);b=endpoint(M,v)
                    queries=[(ext|(1<<u)|((1<<v) if adj[v]>>u&1 else 0),ext|(1<<v)) for u in names if u!=v]
                    if d%2:queries.append((ext|(1<<v),ext|(1<<v)))
                    if (m-1+d)%2:queries.append((ext,ext|(1<<v)))
                    original=sum(value(code,a,c) for a,c in queries)
                    only=sum(value(code^flips[M],a,c) for a,c in queries)
                    every=sum(value(target,a,c) for a,c in queries)
                    need(only==every,'aggregate invariance')
                    if m>=3:need((original-every)%4==2*(1+b)%4,'anchor congruence')
                    else:
                        need(b==1 and (original-every)%4==2,'n=2 E1 boundary control')
                        boundary+=1
                    q=(ext|(1<<v),ext)
                    counts=[value(t,*q) for t in (code,code^flips[M],target)]
                    need(counts[1]==counts[2],'singleton invariance')
                    if b:need((counts[0]-counts[2])%4==2,'endpoint separator')
                    key=(n,code,tuple(chosen),M,v);need(key in records,'missing anchor record');rec=records[key]
                    need(rec['target_code']==target and rec['query_masks']==[list(q) for q in queries],'record query semantics')
                    need([rec['source_query_sum'],rec['single_outer_query_sum'],rec['all_reversals_query_sum']]==[original,only,every],'record exact sums')
                    need(rec['endpoint_counts_source_single_all']==counts and rec['endpoint_parity']==b,'record endpoint data')
                    seen.add(key);anchors+=1
            if n<=4:
                for M in mods:
                    ext=full^M;names=[u for u in range(n) if M>>u&1]
                    for a,c,cy in product(submasks(ext),submasks(ext),(0,M)):
                        left=sum(value(code,a|(1<<u),c|cy) for u in names)
                        right=sum(value(code^flips[M],a|(1<<u),c|cy) for u in names)
                        need(left==right,'exact port averaging');average+=1
    need(seen==set(records),'record coverage')
    need((families,anchors,covered,average)==(15818,33808,1064,3828),'independent scope counts')
    return dict(status='PASS',exposed_laminar_families=families,anchor_certificates=anchors,
                fully_covered_out_of_scope_families=covered,exact_averaging_checks=average,
                two_vertex_boundary_certificates=boundary,
                source='Literal-word full banks independently reconstructed; no producer response functions imported.')


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
