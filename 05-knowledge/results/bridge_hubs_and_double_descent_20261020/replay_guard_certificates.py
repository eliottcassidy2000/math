#!/usr/bin/env python3
"""Read-only literal replay of the 65 saved rigid-guard edge certificates.
Packaging convenience, not a new enumeration of eligible source classes.
Imports only this audit's independently implemented tournament helpers.
"""
from pathlib import Path
from collections import defaultdict
import json
from independent_controls import check, decode, encode, reflect, edge_key, paths, triangle, connected
ROOT=Path(__file__).resolve().parent

def K(rows,path):
    full=(1<<len(rows))-1
    out=[full^r^(1<<i) for i,r in enumerate(rows)]
    for a,b in zip(path,path[1:]):out[a]^=1<<b;out[b]^=1<<a
    return out

def main():
    certs=json.loads((ROOT/'controls'/'boundary_hub_certificates.json').read_text())
    covered=defaultdict(set);sources={};counts=defaultdict(int)
    for c in certs:
        n=c['n'];rows=decode(c['source_representative'],n)
        if c['source_conversed']:rows=[((1<<n)-1)^r^(1<<i) for i,r in enumerate(rows)]
        check(not connected(rows),('source must be nonstrong',c))
        p=tuple(c['first_path']);pp=tuple(c['other_path'])
        check((p[0],p[1],p[-1])==(pp[0],pp[1],pp[-1]),('shared triple',c))
        a,b,z=p[0],p[1],p[-1];d=[r.bit_count() for r in rows];delta=2*d[a]-d[b]-d[z]-3
        check(d[a]+d[b]!=n and delta!=0,('hub tests',c))
        check(delta==c['hub_c3_change'],('stored potential',c))
        expected=rows[:]
        for u,v in ((a,b),(a,z)):expected[u]^=1<<v;expected[v]^=1<<u
        check(triangle(expected)-triangle(rows)==delta,('literal potential',c))
        edge_ids=[];rotation_ids=[]
        for path in (p,pp):
            check(set(path)==set(range(n)) and len(path)==n,('permutation',c))
            check(all((rows[u]>>v)&1 for u,v in zip(path,path[1:])),('source path arcs',c))
            t=encode(rows,path);check(t!=reflect(t,n),('source blackness',c));edge_ids.append(edge_key(t,n))
            partner=K(rows,path);check(connected(partner),('strong partner',c));q=path[1:]+path[:1]
            check(all((partner[u]>>v)&1 for u,v in zip(q,q[1:])),('rotation arcs',c))
            qt=encode(partner,q);check(qt!=reflect(qt,n),('rotation blackness',c));rotation_ids.append(edge_key(qt,n))
            check(K(partner,q)==expected,('literal common hub',c))
        check(edge_ids==[c['target_edge'],c['other_edge']],('source edge IDs',c))
        check(rotation_ids==c['rotation_edges'],('rotation edge IDs',c))
        check(edge_ids[0] not in (edge_ids[1],*rotation_ids),('target-edge avoidance',c))
        key=(n,c['source_node']);covered[key].add(edge_ids[0]);sources[key]=rows;counts[n]+=1
    for (n,u),rows in sources.items():
        all_edges={edge_key(encode(rows,p),n) for p in paths(rows) if encode(rows,p)!=reflect(encode(rows,p),n)}
        check(covered[n,u]==all_edges,('all Q edges at listed source',n,u))
    result={'status':'FINITE-EXACT; all 65 saved literal certificates passed','scope':'All incident Q edges at the nine listed eligible source nodes; not an independent enumeration of eligible n8 classes.','n7_routes':counts[7],'n8_routes':counts[8],'listed_source_nodes':len(sources),'total_routes':len(certs)}
    check(result==json.loads((ROOT/'results'/'guard_certificate_replay.json').read_text()),'Bundled guard replay result mismatch')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
