"""Independent named-endpoint subset convolution for vertex arrival.

For an unordered neighbor pair a,b of the newcomer, split the old vertex set
as A containing a and not b. Fixed-start simple paths on A and its complement
are exactly the two marked components. Their OTHER endpoints remain recorded.
"""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import combinations
from math import isqrt
import argparse,json
HERE=Path(__file__).resolve().parent;checks=0
def need(x,why):
    global checks
    checks+=1
    if not x:raise RuntimeError(why)
def prefix_paths(m,start):
    adj={v:[w for w in range(1,m+1) if v!=w and isqrt(v+w)**2==v+w] for v in range(1,m+1)};result=defaultdict(Counter);nodes=0
    def dfs(mask,last):
        nonlocal nodes
        nodes+=1;result[mask][last]+=1
        for w in adj[last]:
            if not mask>>(w-1)&1:dfs(mask|(1<<(w-1)),w)
    dfs(1<<(start-1),start);return result,nodes
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args();paths={r['n']:[tuple(p) for p in r['paths']] for r in json.loads((HERE/'paths.json').read_text())['orders']};out={'status':'FINITE-EXACT named joining-port and other-endpoint convolution for every arrival3..32','orders':[]}
    for n in range(3,33):
        neighbors=[a for a in range(1,n) if isqrt(a+n)**2==a+n];full=(1<<(n-1))-1;bank={};stateinfo=[]
        for a in neighbors:
            bank[a],nodes=prefix_paths(n-1,a);stateinfo.append({'start':a,'prefixes':nodes,'vertex_subsets':len(bank[a])})
        prediction=Counter();extension=0;paircounts=[];edgeuses=Counter()
        for a in neighbors:
            for p in paths[n-1]:
                if a in (p[0],p[-1]):
                    other=p[-1] if p[0]==a else p[0];prediction[tuple(sorted((n,other)))]+=1;extension+=1;edgeuses[a]+=1
        for a,b in combinations(neighbors,2):
            ends=Counter();component_splits=0;witness=None
            for mask,byend in bank[a].items():
                if mask>>(b-1)&1:continue
                complement=full^mask
                if complement not in bank[b]:continue
                component_splits+=1
                for u,cu in byend.items():
                    for v,cv in bank[b][complement].items():
                        need(u!=v,'endpoints of disjoint components overlap');ends[tuple(sorted((u,v)))]+=cu*cv
                        if witness is None:witness={'component_a':[i+1 for i in range(n-1) if mask>>i&1],'component_b':[i+1 for i in range(n-1) if complement>>i&1],'other_endpoints':[u,v]}
            amount=sum(ends.values());prediction+=ends;edgeuses[a]+=amount;edgeuses[b]+=amount;paircounts.append({'joining_ports':[a,b],'two_cover_count':amount,'component_vertex_splits':component_splits,'other_endpoint_counts':[[list(e),k] for e,k in sorted(ends.items())],'first_split':witness})
        actual=Counter(tuple(sorted((p[0],p[-1]))) for p in paths[n]);need(prediction==actual,('full endpoint recurrence',n))
        for a in neighbors:need(edgeuses[a]==sum(any({u,v}=={a,n} for u,v in zip(p,p[1:])) for p in paths[n]),('new edge incidence recurrence',n,a))
        row={'n':n,'new_neighbors':neighbors,'endpoint_extension_count':extension,'interior_count':sum(r['two_cover_count'] for r in paircounts),'joining_pairs':paircounts,'new_edge_uses':[[a,k] for a,k in sorted(edgeuses.items())],'prefix_state_sizes':stateinfo};out['orders'].append(row);print(n,'endpoint',extension,'interior',row['interior_count'],'prefixes',sum(x['prefixes'] for x in stateinfo),flush=True)
    out['checks']=checks
    if args.write:(HERE/'recurrence.json').write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads((HERE/'recurrence.json').read_text())==out,'frozen recurrence receipt')
    print('checks',checks,flush=True)
if __name__=='__main__':main()
