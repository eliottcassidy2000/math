"""Literal-permutation audit of eleven H71 cores and exact ray reconstruction.

This imports no birth-budget or repository counting implementation and does
not regenerate or certify the inherited 6880-class isomorphism census.
"""
from itertools import permutations
from pathlib import Path
import hashlib,json
from neutral_scc import marked, sccs


def need(value,message):
    if not value:raise RuntimeError(message)


def literal(rows):
    n=len(rows);full=(1<<n)-1
    h=[0]*(full+1);h[0]=1
    end=[[0]*n for _ in range(full+1)]
    starts=[0]*n;attempts=0
    for mask in range(1,full+1):
        vertices=[v for v in range(n) if mask>>v&1]
        for path in permutations(vertices):
            attempts+=1
            if all(rows[a]>>b&1 for a,b in zip(path,path[1:])):
                h[mask]+=1;end[mask][path[-1]]+=1
                if mask==full:starts[path[0]]+=1
    P2=sum(h[mask]*h[full^mask] for mask in range(full+1))
    fields=[]
    for root in range(n):
        F=sum(end[mask][root]*h[full^mask] for mask in range(1,full+1) if mask>>root&1)
        q=h[full^(1<<root)];s=starts[root]
        d=sum(bool(row>>root&1) for row in rows)
        fields.append(dict(root=root,start=s,deletion=q,F=F,indegree=d,delta=F-h[full]-q+s))
    return h[full],P2,fields,attempts


def oddpart(n):
    while n%2==0:n//=2
    return n


def main():
    directory=Path(__file__).resolve().parent
    input_path=directory.parent/'returns/birth_budget_receipt.json'
    receipt=json.loads(input_path.read_text(encoding='utf-8'))
    census=directory.parents[3]/receipt['source']['path']
    census_sha=hashlib.sha256(census.read_bytes()).hexdigest()
    need(census_sha==receipt['source']['sha256'],'inherited census source hash')
    need(sum(bool(line.strip()) for line in census.read_text().splitlines())==6880,
         'inherited census line count only')
    results=[];attempts=0;support={}
    for core in receipt['strong_core_presentations']:
        rows=core['rows'];n=len(rows);full=(1<<n)-1
        need(all(not (rows[v]>>v&1) and rows[v]<=full for v in range(n)),
             'loopless in-range rows')
        need(all(bool(rows[a]>>b&1)!=bool(rows[b]>>a&1)
                 for a in range(n) for b in range(a+1,n)), 'tournament pairs')
        reach=[row|(1<<v) for v,row in enumerate(rows)]
        for k in range(n):
            for v in range(n):
                if reach[v]>>k&1:reach[v]|=reach[k]
        need(all(row==full for row in reach),'strong core')
        h,p2,fields,count=literal(rows);attempts+=count
        need(h==core['H']==71,'literal H')
        need(p2==core['P2'],'literal P2')
        need(fields==core['roots'],'all literal marked profiles')
        need(sum(r['F'] for r in fields)==p2-h,'complete cover partition')
        need(oddpart(p2)!=405 and all(oddpart(r['F'])!=405 for r in fields),'odd-part exclusion')
        for amplitude,extra in [(p2,0)]+[(r['F'],1) for r in fields]:
            odd=oddpart(amplitude);exponent=(amplitude//odd).bit_length()-1+extra
            support[odd]=min(support.get(odd,exponent),exponent)
        results.append(dict(rows=rows,order=len(rows),P2=p2,F=[r['F'] for r in fields]))
    need(len(results)==11 and sum(r['order'] for r in results)==81,'declared finite universe')
    rays=[dict(oddpart=o,minimum_exponent=e) for o,e in sorted(support.items())]
    need(rays==receipt['balanced_t71_ray_support'] and len(rays)==62,'complete ray compression')
    need(405 not in support and support[101]==3,'405 gap and 404 family threshold')
    neighbor=receipt['neighboring_ray'];previous=tuple(neighbor['parent_rows']);oldroot=neighbor['core_root']
    seed=marked(previous,oldroot)
    need((seed[0],seed[1],seed[3])==(71,0,808),'literal weak seed channels')
    levels=[]
    for level in neighbor['literal_levels']:
        rows=tuple(level['rows']);root=level['root'];n=len(previous)
        expected=list(previous);expected[oldroot]|=1<<n;expected.append(((1<<n)-1)^(1<<oldroot))
        need(tuple(expected)==rows,'literal N/R birth adjacency')
        z=marked(rows,root)
        need(z==tuple(level['state'])==(404*2**level['k']+71,71,71,404*2**level['k']+71),
             'independent full neighboring ray state')
        need(len(sccs(rows))==1,'neighboring ray strong')
        levels.append(dict(k=level['k'],state=z,order=len(rows)))
        previous=rows;oldroot=root
    result={'status':'FINITE-EXACT','optimized':not __debug__,
            'scope':'Eleven supplied literal strong H71 cores, 62-ray compression, and four actual neighboring-ray graphs; inherited isomorphism census completeness not regenerated.',
            'literal_permutations_tested':attempts,'root_profiles':81,'cores':results,
            'exact_ray_support':rays,'neighboring_levels':levels,
            'inherited_census_sha256':census_sha,
            'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest()}
    name='birth_budget_audit.optimized.json' if not __debug__ else 'birth_budget_audit.json'
    (directory/name).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','core_presentations':11,'root_profiles':81,'literal_permutations_tested':attempts}))


if __name__=='__main__':main()
