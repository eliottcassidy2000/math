"""Independent endpoint-marked edge-subset census and incidence inventories."""
from pathlib import Path
from itertools import combinations
from collections import Counter,defaultdict
from math import isqrt,comb
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];checks=0
def need(x,why):
    global checks
    checks+=1
    if not x:raise RuntimeError(why)
def graph(n):return [(a,b) for a in range(1,n+1) for b in range(a+1,n+1) if isqrt(a+b)**2==a+b]
def edge_paths(n,edges):
    if n==1:return [(1,)],1
    incidence=[0]*(n+1)
    for k,(a,b) in enumerate(edges):incidence[a]|=1<<k;incidence[b]|=1<<k
    if any(not incidence[v] for v in range(1,n+1)):return [],0
    leaves={v for v in range(1,n+1) if incidence[v].bit_count()==1};paths=[];nodes=0
    for ends in combinations(range(1,n+1),2):
        if not leaves<=set(ends):continue
        degree=[0]+[1 if v in ends else 2 for v in range(1,n+1)]
        def search(selected,excluded):
            nonlocal nodes
            nodes+=1
            while True:
                old=(selected,excluded)
                for v in range(1,n+1):
                    present=(incidence[v]&selected).bit_count();unknown=incidence[v]&~(selected|excluded);available=unknown.bit_count();d=degree[v]-present
                    if d<0 or d>available:return
                    if d==0:excluded|=unknown
                    elif d==available:selected|=unknown
                if selected&excluded:return
                if old==(selected,excluded):break
            parent=list(range(n+1))
            def find(v):
                while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
                return v
            for k,(a,b) in enumerate(edges):
                if selected>>k&1:
                    x,y=find(a),find(b)
                    if x==y:return
                    parent[x]=y
            if selected.bit_count()==n-1:
                adjacency=defaultdict(list)
                for k,(a,b) in enumerate(edges):
                    if selected>>k&1:adjacency[a].append(b);adjacency[b].append(a)
                path=[ends[0]]
                while len(path)<n:
                    nxt=[v for v in adjacency[path[-1]] if len(path)==1 or v!=path[-2]];need(len(nxt)==1,'accepted edge set not a path');path.append(nxt[0])
                need(path[-1]==ends[1] and len(set(path))==n,'endpoint completeness');paths.append(tuple(path));return
            # Excluded edges may disconnect the remaining allowable graph.
            for k,(a,b) in enumerate(edges):
                if not excluded>>k&1:parent[find(a)]=find(b)
            if len({find(v) for v in range(1,n+1)})>1:return
            choices=[]
            for v in range(1,n+1):
                u=incidence[v]&~(selected|excluded);d=degree[v]-(incidence[v]&selected).bit_count()
                if u:choices.append((comb(u.bit_count(),d),u.bit_count(),v,u))
            if not choices:return
            u=min(choices)[3];bit=u&-u;search(selected|bit,excluded);search(selected,excluded|bit)
        search(0,0)
    need(len(paths)==len(set(paths)),'edge-subset duplicate');return sorted(paths),nodes
USER=[[(23,13),(23,2),(25,11),(25,24),(15,1),(15,10),(17,8),(16,9)],[(13,3),(13,12),(14,2),(14,11),(19,6),(19,17),(20,5),(20,16),(21,4),(21,15),(22,3),(22,14),(24,1),(24,12)],[(18,7),(1,3),(4,12),(7,2)]]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args();literal=json.loads((HERE/'paths.json').read_text());out={'status':'FINITE-EXACT complete square-sum paths up to reversal, 1<=n<=32; independent endpoint edge-subset enumeration','orders':[],'user_partition':[],'user_correction':{'original_pair':[3,12],'confirmed_pair':[4,12],'status':'User confirmed the original3:12 was a typo; use4:12.'}};by={}
    for item in literal['orders']:
        n=item['n'];edges=graph(n);expected=sorted(tuple(p) for p in item['paths']);actual,nodes=edge_paths(n,edges);need(actual==expected,('independent complete path set',n));count=len(actual);use=Counter(tuple(sorted((a,b))) for p in actual for a,b in zip(p,p[1:]));end=Counter(tuple(sorted((p[0],p[-1]))) for p in actual);vends=Counter(v for p in actual for v in (p[0],p[-1]));E=[{'edge':list(e),'uses':use[e],'status':'no_paths' if not count else 'never' if not use[e] else 'forced' if use[e]==count else 'sometimes'} for e in edges]
        need(sum(use.values())==count*(n-1),'edge incidence sum')
        for v in range(1,n+1):need(sum(use[e] for e in edges if v in e)==2*count-vends[v],('endpoint incidence identity',n,v))
        arrival=[e for e in edges if n in e];extension=sum(n in (p[0],p[-1]) for p in actual);interior=count-extension
        if n>=2:need(sum(use[e] for e in arrival)==extension+2*interior,('arrival handshake',n))
        row={'n':n,'path_count':count,'edge_subset_nodes':nodes,'endpoint_pairs':[[list(e),c] for e,c in sorted(end.items())],'edges':E,'arrival_edges':[list(e) for e in arrival],'new_vertex_endpoint_paths':extension,'new_vertex_interior_paths':interior}
        by[n]=row;out['orders'].append(row);print(n,'paths',count,'edge-subset nodes',nodes,flush=True)
    need([by[n]['path_count'] for n in range(23,33)]==[3,0,10,12,35,52,19,20,349,392],'known counts')
    for group in USER:
        records=[]
        for stated in group:
            e=tuple(sorted(stated));birth=max(e);valid=e in graph(25);timeline=[]
            for n in range(birth,33):
                if e in graph(n):timeline.append({'n':n,**next(x for x in by[n]['edges'] if tuple(x['edge'])==e)})
            records.append({'stated_pair':list(stated),'edge':list(e),'valid_square_edge':valid,'sum':sum(e),'birth':birth,'arrival_degree':isqrt(2*birth-1)-isqrt(birth),'birth_path_count':by[birth]['path_count'],'birth_predecessor_path_count':by[birth-1]['path_count'] if birth>1 else None,'timeline':timeline})
        out['user_partition'].append(records)
    selected={tuple(sorted(e)) for g in USER for e in g};out['unlisted_G25_edges']=[list(e) for e in graph(25) if e not in selected]
    out['input_paths_hash_protocol']='SHA-256 after replacing CRLF with LF; no other normalization.'
    out['input_paths_sha256']=hashlib.sha256((HERE/'paths.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest();out['checks']=checks
    if args.write:(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads((HERE/'verification.json').read_text())==out,'frozen verification')
    print('checks',checks,flush=True)
if __name__=='__main__':main()
