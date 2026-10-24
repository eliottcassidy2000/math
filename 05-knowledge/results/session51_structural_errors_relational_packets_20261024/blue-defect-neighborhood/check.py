#!/usr/bin/env python3
"""Complete blue neighborhood of the 24 positive-error order-ten nodes.
Full cross-matrix canonicalization by actual stable colors and intrinsic cycles.
No permutation search beyond at most25 cell rotations on any target.
"""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import product
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
COUNTS=Counter()
def need(ok,msg):
    COUNTS[msg]+=1
    if not ok:raise RuntimeError(msg)
def stable(T):
    colors=(0,)*len(T)
    while True:
        signatures=tuple((colors[v],)+tuple(sum((T[v]>>w)&1 for w in range(len(T)) if colors[w]==c) for c in range(max(colors)+1)) for v in range(len(T)))
        names={s:i for i,s in enumerate(sorted(set(signatures)))}
        nxt=tuple(names[s] for s in signatures)
        if len(set(nxt))==len(set(colors)):return tuple(tuple(v for v in range(len(T)) if nxt[v]==c) for c in range(max(nxt)+1))
        colors=nxt

def intrinsic_cycles(T):
    C=stable(T);out=[]
    need(sorted(map(len,C)) in ([1]*10,[1]*4+[3]*2,[5]*2),'self-converse even ten stable-size menu')
    for cell in C:
        if len(cell)==1:out.append(cell);continue
        succ={}
        for v in cell:
            targets=[w for w in cell if T[v]>>w&1]
            need(len(targets)==(len(cell)-1)//2,'regular internal cell')
            if len(cell)==3:succ[v]=targets[0]
            else:
                x,y=targets;succ[v]=x if T[x]>>y&1 else y
        v=min(cell);cycle=[]
        for _ in cell:cycle.append(v);v=succ[v]
        need(v==cycle[0] and set(cycle)==set(cell),'intrinsic successor is one cycle')
        need(all(((T[cycle[i]]>>cycle[j])&1)==int(0<(j-i)%len(cell)<=(len(cell)-1)//2) for i in range(len(cell)) for j in range(len(cell))),'intrinsic coordinates are oriented cyclic tournament')
        out.append(tuple(cycle))
    return tuple(out)
def normalize(T,order):return tuple(sum(((T[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)
def canonical(T):
    C=intrinsic_cycles(T);best=None;bestorder=None;rotations=1
    for cell in C:rotations*=len(cell)
    need(rotations<=25,'at most twenty-five complete cross-matrix normalizations')
    for shifts in product(*(range(len(cell)) for cell in C)):
        order=tuple(v for cell,k in zip(C,shifts) for v in cell[k:]+cell[:k]);candidate=normalize(T,order)
        if best is None or candidate<best:best,bestorder=candidate,order
    return best,bestorder

def reverse(T):return tuple(((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T))
def error(T):
    C=intrinsic_cycles(T);group=[]
    for shifts in product(*(range(len(cell)) for cell in C)):
        p=list(range(len(T)))
        for cell,k in zip(C,shifts):
            for i,v in enumerate(cell):p[v]=cell[(i+k)%len(cell)]
        if all(((T[i]>>j)&1)==((T[p[i]]>>p[j])&1) for i in range(len(T)) for j in range(len(T))):group.append(tuple(p))
    todo=set(range(len(T)));sizes={v:len(c) for c in C for v in c};orbits=[]
    while todo:
        v=min(todo);orb={g[v] for g in group};need(orb<=todo,'full automorphism orbits partition');todo-=orb;orbits.append(sorted(orb))
    return sum(sizes[o[0]]-len(o) for o in orbits),group,orbits

def decode(word,pairs):
    T=[0]*10
    for i in range(9):T[i]|=1<<(i+1)
    for k,(i,j) in enumerate(pairs):
        a,b=(i,j) if word>>k&1 else (j,i);T[a]|=1<<b
    return tuple(T)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(HERE/'normal.json'));args=ap.parse_args()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for pin in pins['files']:need(hashlib.sha256((HERE/pin['path']).read_bytes()).hexdigest()==pin['sha256'],'frozen input SHA256')
    need(hashlib.sha256((HERE/'sources/wiring-MANIFEST.json').read_bytes()).hexdigest()=='cd89cad144786c67f076c235409c90e3b4de21e0a5f81926139184068f1d676d','exact wiring freeze')
    base=json.loads((HERE/'sources/wiring-MANIFEST.json').read_text());hashes={r['path']:r['sha256'] for r in base['files']}
    for original,archived in [('normal.json','wiring-receipt.json'),('THEOREM.md','wiring-THEOREM.md')]:need(hashlib.sha256((HERE/'sources'/archived).read_bytes()).hexdigest()==hashes[original],'exact archived wiring artifact')
    wiring=json.loads((HERE/'sources/wiring-receipt.json').read_text());pairs=[tuple(p) for p in wiring['normalized_nonpath_pair_order']]
    need(len({e['word'] for e in wiring['boundary_edges']})==4020 and
         len({e['complement_word'] for e in wiring['boundary_edges']})==4020,
         'boundary source and target presentation words are distinct')
    target_records=[];key_counts=Counter();patterns=Counter()
    for edge in wiring['boundary_edges']:
        T=decode(edge['complement_word'],pairs)
        need(all(((T[i]>>j)&1)==((T[9-j]>>(9-i))&1) for i in range(10) for j in range(10)),'boundary target self-converse by literal reversal')
        key,order=canonical(T);key_counts[key]+=1;patterns[tuple(sorted(map(len,stable(T))))]+=1
        need(normalize(T,order)==key,'every target canonical witness retains full matrix')
        target_records.append({'source_index':edge['source_index'],'source_word':edge['word'],'target_word':edge['complement_word'],'canonical_key':key,'canonical_vertex_order':order})
    keys=sorted(key_counts);ids={key:i for i,key in enumerate(keys)};targets=[];class_patterns=Counter();aut_hist=Counter()
    for i,T in enumerate(keys):
        need(canonical(T)[0]==T,'canonical target idempotence')
        need(canonical(reverse(T))[0]==T,'canonical target self-converse')
        # Deterministic nontrivial vertex permutations test label invariance;
        # the key contains the complete directed matrix, never a gain kernel.
        for order in [tuple(reversed(range(10))),tuple((7*v+3)%10 for v in range(10))]:need(canonical(normalize(T,order))[0]==T,'full canonical key relabeling invariant')
        e,group,orbits=error(T);need(e==0,'every distinct boundary class has directly verified zero error')
        pattern=tuple(sorted(map(len,stable(T))));class_patterns[pattern]+=1;aut_hist[len(group)]+=1
        targets.append({'target_class_id':i,'canonical_rows':T,'stable_cells':stable(T),'error':e,'automorphism_order':len(group),'vertex_orbits':orbits,'incident_edges_from_positive_slice':key_counts[T]})
    adjacency=Counter()
    for rec in target_records:
        rec['target_class_id']=ids[rec.pop('canonical_key')];adjacency[rec['source_index'],rec['target_class_id']]+=1
    need(sum(adjacency.values())==len(target_records)==4020,'all 4020 boundary edges assigned to exact target classes')
    need(sum(key_counts.values())==4020,'boundary class incidence conservation')
    source_counts=Counter()
    for (source,target),m in adjacency.items():source_counts[source]+=m
    need(source_counts==Counter({n['source_representative_index']:n['boundary_edges'] for n in wiring['nodes']}),'per-source boundary degrees preserved')
    # This graph contains all edges touching the positive slice, not all
    # edges among its union of vertices; zero-zero edges are unenumerated.
    graph={('P',n['source_representative_index']):set() for n in wiring['nodes']}
    graph.update({('Z',i):set() for i in range(len(keys))})
    for a,b,_ in wiring['adjacency_multiplicities']:graph['P',a].add(('P',b));graph['P',b].add(('P',a))
    for a,b in adjacency:graph['P',a].add(('Z',b));graph['Z',b].add(('P',a))
    unseen=set(graph);components=[]
    while unseen:
        start=min(unseen);seen={start};queue=[start]
        while queue:
            for v in graph[queue.pop()]-seen:seen.add(v);queue.append(v)
        unseen-=seen;components.append({'positive_nodes':sorted(i for typ,i in seen if typ=='P'),'zero_error_nodes':sorted(i for typ,i in seen if typ=='Z')})
    result={'status':'PASS','scope':'complete blue edge neighborhood of the 24 positive-error nodes; no zero-zero or black edges reconstructed','base_wiring_manifest_sha256':pins['base_wiring_manifest_sha256'],'positive_nodes':24,'distinct_zero_error_boundary_nodes':len(keys),'all_neighborhood_vertices':24+len(keys),'internal_positive_edges':90,'boundary_edges':4020,'all_edges_incident_to_positive_slice':4110,'boundary_presentation_stable_patterns':[[list(k),v] for k,v in sorted(patterns.items())],'boundary_class_stable_patterns':[[list(k),v] for k,v in sorted(class_patterns.items())],'boundary_class_automorphism_orders':sorted(aut_hist.items()),'boundary_adjacency_multiplicity_histogram':sorted(Counter(adjacency.values()).items()),'boundary_target_incidence_degree_histogram':sorted(Counter(key_counts.values()).items()),'components_of_incident_edge_graph':components,'boundary_adjacency':[[a,b,m] for (a,b),m in sorted(adjacency.items())],'target_classes':targets,'boundary_edge_witnesses':target_records,'checks':dict(COUNTS)}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('boundary_adjacency','target_classes','boundary_edge_witnesses','checks','components_of_incident_edge_graph')},indent=2,sort_keys=True))
if __name__=='__main__':main()
