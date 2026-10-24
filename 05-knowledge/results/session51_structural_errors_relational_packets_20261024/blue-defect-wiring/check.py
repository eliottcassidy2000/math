#!/usr/bin/env python3
"""Actual complement-line pairing on the 4200 blue ports of 24 defect nodes.
No target tournament canonicalization; boundary targets are not class-merged.
"""
from pathlib import Path
from itertools import combinations,permutations,product
from collections import Counter,defaultdict
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(HERE/'normal.json'));args=ap.parse_args()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for pin in pins['files']:need(hashlib.sha256((HERE/pin['path']).read_bytes()).hexdigest()==pin['sha256'],'frozen source SHA256')
    need(hashlib.sha256((HERE/'sources/converse-MANIFEST.json').read_bytes()).hexdigest()=='30b62a8fdf93065b19f071050d9fcfe5ad08c5f79cfcfbf30afcdae53babde5a','exact audited converse packet')
    base_manifest=json.loads((HERE/'sources/converse-MANIFEST.json').read_text())
    base_hashes={r['path']:r['sha256'] for r in base_manifest['files']}
    for original,archived in [('normal.json','converse-receipt.json'),('port-normal.json','blue-port-receipt.json'),('THEOREM.md','converse-THEOREM.md')]:
        need(hashlib.sha256((HERE/'sources'/archived).read_bytes()).hexdigest()==base_hashes[original],
             'archived source is exactly the audited base artifact')
    records=json.loads((HERE/'sources/converse-receipt.json').read_text())['representatives']
    ports=json.loads((HERE/'sources/blue-port-receipt.json').read_text())['representatives']
    expected={r['representative_index']:r['blue_ports'] for r in ports}
    pairs=[(i,j) for i in range(10) for j in range(i+2,10)];nbits=len(pairs);full=(1<<nbits)-1
    position={pair:k for k,pair in enumerate(pairs)}
    reflection=tuple(position[(9-j,9-i)] for i,j in pairs)
    word_owner={};word_path={};count=Counter();all_rows={};cross_degrees={}
    for rec in records:
        if not rec['self_converse']:continue
        idx=rec['index'];T=tuple(rec['rows']);q=tuple(rec['anti_certificates'][0]['anti_involution'])
        need(rec['automorphism_order']==1 and rec['error']==40,'complete positive self-converse slice rigid at error40')
        cross_degrees[idx]=min(r.bit_count() for r in T)-2
        orbits=[(v,q[v]) for v in range(10) if v<q[v]]
        for selections in product(*orbits):
            for half in permutations(selections):
                path=half+tuple(q[v] for v in half[::-1])
                if not all(T[a]>>b&1 for a,b in zip(path,path[1:])):continue
                word=sum(((T[path[i]]>>path[j])&1)<<k for k,(i,j) in enumerate(pairs))
                need(all(((word>>k)&1)==((word>>reflection[k])&1) for k in range(nbits)),'every represented word literally blue')
                need(word not in word_owner,'rigidity and ordinary class distinction give unique words')
                word_owner[word]=idx;word_path[word]=path;count[idx]+=1
                rows=[0]*10
                for i in range(9):rows[i]|=1<<(i+1)
                for k,(i,j) in enumerate(pairs):
                    a,b=(i,j) if word>>k&1 else (j,i);rows[a]|=1<<b
                need(all(((rows[i]>>j)&1)==((T[path[i]]>>path[j])&1) for i in range(10) for j in range(10)),'literal normalized rows match source path')
                all_rows[word]=tuple(rows)
    need(dict(count)==expected and len(word_owner)==4200,'all and only 4200 previously counted blue ports')
    adjacency=Counter();internal_ports=Counter();boundary_ports=Counter();transport=Counter();edges=[];boundary=[]
    for word,idx in sorted(word_owner.items()):
        target=full^word
        need(target!=word and (full^target)==word,'complement is fixed-point-free involution')
        need(all(((target>>k)&1)==((target>>reflection[k])&1) for k in range(nbits)),'complement partner is blue')
        # Independently construct the partner by flipping actual non-path arcs
        # in the normalized row tuple and verify the encoded complement.
        flipped=list(all_rows[word])
        for i,j in pairs:flipped[i]^=1<<j;flipped[j]^=1<<i
        encoded=sum(((flipped[i]>>j)&1)<<k for k,(i,j) in enumerate(pairs))
        need(encoded==target,'literal non-path arc flip equals complement encoding')
        source_degrees=[r.bit_count() for r in all_rows[word]]
        target_degrees=[r.bit_count() for r in flipped]
        need(all(target_degrees[i]==9-source_degrees[i]+int(i==0)-int(i==9) for i in range(10)),
             'inherited endpoint-score transport on every blue port')
        ds,dt=source_degrees[0],source_degrees[9]
        moment_delta=sum(d*d for d in target_degrees)-sum(d*d for d in source_degrees)
        need(moment_delta==2*(dt-ds+1),'squared-score endpoint transport')
        def cyclic_triangles(rows):
            return sum(((rows[a]>>b)&1)==((rows[b]>>c)&1)==((rows[c]>>a)&1)
                       for a,b,c in combinations(range(10),3))
        triangle_delta=cyclic_triangles(flipped)-cyclic_triangles(all_rows[word])
        need(triangle_delta==ds-dt-1,'literal cyclic-triangle endpoint transport')
        transport[cross_degrees[idx],ds,dt,triangle_delta,target in word_owner]+=1
        if target in word_owner:
            other=word_owner[target];internal_ports[idx]+=1
            need(cross_degrees[idx]==cross_degrees[other]==2 and ds==5 and dt==4,'all internal edges obey structural degree obstruction')
            need(all_rows[target]==tuple(flipped),'internal partner actual tournament rows exact')
            if word<target:
                a,b=sorted((idx,other));adjacency[a,b]+=1
                edges.append({'source_index':idx,'target_index':other,'word':word,'complement_word':target,'source_path':word_path[word],'target_path':word_path[target]})
        else:
            boundary_ports[idx]+=1;boundary.append({'source_index':idx,'word':word,'complement_word':target})
    need(sum(internal_ports.values())==2*len(edges),'internal port count even and twice edge count')
    need(2*len(edges)+len(boundary)==4200,'port-edge-boundary conservation')
    need(all(internal_ports[i]+boundary_ports[i]==count[i] for i in count),'per-node port conservation')
    degree=Counter()
    for (a,b),m in adjacency.items():degree[a]+=m;degree[b]+=m
    need(degree==internal_ports,'loops contribute twice to induced degree')
    # Connected components use the actual positive-slice adjacency, ignoring
    # multiplicities but retaining isolated vertices; loops do not connect.
    neighbors={i:set() for i in count}
    for a,b in adjacency:neighbors[a].add(b);neighbors[b].add(a)
    todo=set(count);components=[]
    while todo:
        seed=min(todo);seen={seed};queue=[seed]
        while queue:
            v=queue.pop()
            for w in neighbors[v]-seen:seen.add(w);queue.append(w)
        todo-=seen;components.append(sorted(seen))
    def component_count(omit_edge=None,omit_vertex=None):
        vertices=set(count)-({omit_vertex} if omit_vertex is not None else set())
        graph={v:set() for v in vertices}
        for eid,edge in enumerate(edges):
            if eid==omit_edge:continue
            a,b=edge['source_index'],edge['target_index']
            if a not in vertices or b not in vertices:continue
            graph[a].add(b);graph[b].add(a)
        unseen=set(vertices);cc=0
        while unseen:
            cc+=1;todo=[unseen.pop()]
            while todo:
                for w in graph[todo.pop()]&unseen:unseen.remove(w);todo.append(w)
        return cc
    baseline=component_count();need(baseline==len(components),'independent component count')
    bridges=[eid for eid in range(len(edges)) if component_count(omit_edge=eid)>baseline]
    articulation=[v for v in sorted(count) if component_count(omit_vertex=v)>baseline]
    cycle_rank=len(edges)-len(count)+baseline
    result={'status':'PASS','scope':'actual induced blue multigraph on the 24 positive-error self-converse order-ten nodes; boundary target classes are not identified','base_converse_manifest_sha256':pins['base_converse_manifest_sha256'],'normalized_nonpath_pair_order':pairs,'positive_nodes':len(count),'blue_ports':len(word_owner),'induced_blue_edges':len(edges),'induced_blue_loop_edges':sum(m for (a,b),m in adjacency.items() if a==b),'induced_nonloop_blue_edges':sum(m for (a,b),m in adjacency.items() if a!=b),'internal_blue_ports':sum(internal_ports.values()),'boundary_blue_edges':len(boundary),'adjacency_multiplicities':[[a,b,m] for (a,b),m in sorted(adjacency.items())],'induced_cycle_rank':cycle_rank,'bridge_edge_indices':bridges,'articulation_node_indices':articulation,'induced_degree_histogram':sorted(Counter(internal_ports[i] for i in count).items()),'endpoint_transport_histogram':[[list(k),v] for k,v in sorted(transport.items())],'adjacency_multiplicity_histogram':sorted(Counter(adjacency.values()).items()),'components':components,'nodes':[{'source_representative_index':i,'cross_degree':cross_degrees[i],'blue_ports':count[i],'internal_ports':internal_ports[i],'boundary_edges':boundary_ports[i]} for i in sorted(count)],'internal_edges':edges,'boundary_edges':boundary}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('internal_edges','boundary_edges')},indent=2,sort_keys=True))
if __name__=='__main__':main()
