#!/usr/bin/env python3
"""Boundary class audit by ALL within-stable-cell permutations, not rotations.
Canonical keys retain all 100 adjacency entries, with independent labels.
"""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import permutations,product,combinations
import hashlib,json
HERE=Path(__file__).resolve().parent
N=10

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def stable(T):
    cells=[tuple(range(N))]
    while True:
        nxt=[]
        for cell in cells:
            bins={}
            for v in cell:bins.setdefault(tuple(sum(T[v][w] for w in c) for c in cells),[]).append(v)
            nxt.extend(tuple(bins[k]) for k in sorted(bins))
        if len(nxt)==len(cells):return tuple(nxt)
        cells=nxt

def normalize(T,order):return tuple(bytes(T[a][b] for b in order) for a in order)

def orders(cells):
    # Exhaustive within-cell bijections, not intrinsic cycles/rotations.
    for parts in product(*(tuple(permutations(c)) for c in cells)):
        yield tuple(v for part in parts for v in part)

def canonical(T):
    cells=stable(T)
    best=None;witness=None;count=0
    for order in orders(cells):
        candidate=normalize(T,order);count+=1
        if best is None or candidate<best:best=candidate;witness=order
    return best,witness,count

def automorphisms(T):
    cells=stable(T);source=tuple(v for c in cells for v in c);group=[]
    for target in orders(cells):
        p=[0]*N
        for a,b in zip(source,target):p[a]=b
        if all(T[i][j]==T[p[i]][p[j]] for i,j in combinations(range(N),2)):group.append(tuple(p))
    return group

def error(T,group):
    todo=set(range(N));parts=[];cell_size={v:len(c) for c in stable(T) for v in c}
    while todo:
        v=min(todo);orb={p[v] for p in group};need(orb<=todo,'exact automorphism orbits');todo-=orb;parts.append(sorted(orb))
    return sum(cell_size[o[0]]-len(o) for o in parts),parts

def decode(word):
    T=[[0]*N for _ in range(N)]
    for i in range(N-1):T[i][i+1]=1
    k=0
    for i in range(N):
        for j in range(i+2,N):T[i][j]=(word>>k)&1;T[j][i]=1-T[i][j];k+=1
    return tuple(bytes(row) for row in T)

def rows(T):return [sum(b<<j for j,b in enumerate(row)) for row in T]

def components(vertices,pairs):
    parent={v:v for v in vertices}
    def root(v):
        while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
        return v
    for a,b in pairs:parent[root(a)]=root(b)
    groups=defaultdict(list)
    for v in sorted(vertices):groups[root(v)].append(v)
    return sorted(groups.values())

def main():
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for file,digest in pins.items():need(hashlib.sha256((HERE/file).read_bytes()).hexdigest()==digest,'pinned source corruption: '+file)
    need(pins['sources/candidate-MANIFEST.json']=='141cbcb8e2696238c2db7db40e4bcd9570438d791a4c405e662427671db8f6da','exact frozen neighborhood identity')
    need(pins['sources/independent-wiring-MANIFEST.json']=='53b674adb896c5f689e38b82fcbefd3a21ea0dfe9c19e5d96dba8bb43f806c9c','exact independent wiring freeze')
    wiring=json.loads((HERE/'sources/independent-wiring-receipt.json').read_text())
    producer=json.loads((HERE/'sources/candidate-receipt.json').read_text())
    member=defaultdict(list);source=Counter();presentation_patterns=Counter();witnesses=[];search=Counter()
    for edge in wiring['boundary_edges']:
        T=decode(edge['complement_word'])
        need(all(T[i][i]==0 for i in range(N)) and all(T[i][j]+T[j][i]==1 for i,j in combinations(range(N),2)),'decoded exact tournament')
        need(all(T[i][j]!=T[N-1-i][N-1-j] for i,j in combinations(range(N),2)),'literal self-converse boundary board')
        key,order,count=canonical(T);search[count]+=1
        need(normalize(T,order)==key,'actual full-matrix canonical witness')
        member[key].append(edge['word']);source[edge['source_index']]+=1
        presentation_patterns[tuple(sorted(map(len,stable(T))))]+=1
        witnesses.append((edge,key,order))
    keys=sorted(member);ids={key:i for i,key in enumerate(keys)}
    need(len(keys)==1031,'1031 distinct boundary classes under exhaustive full-arc canonicalization')
    targets=[];aut=Counter();patterns=Counter();allgroups={}
    for key in keys:
        group=automorphisms(key);e,orbits=error(key,group)
        need(e==0,'direct actual stable-cell/aut-orbit zero error')
        pattern=tuple(sorted(map(len,stable(key))));patterns[pattern]+=1;aut[len(group)]+=1;allgroups[key]=group
        reverse=tuple(bytes(0 if i==j else 1-key[i][j] for j in range(N)) for i in range(N))
        need(canonical(reverse)[0]==key,'each target class actually self-converse')
        targets.append({'independent_target_id':ids[key],'canonical_rows':rows(key),'stable_cells':[list(c) for c in stable(key)],'automorphism_order':len(group),'vertex_orbits':orbits,'error':e,'incident_edges':len(member[key]),'source_words':sorted(member[key])})
    # Validate every producer class and coordinate witness against full matrices.
    producer_key={};pid_to_independent={}
    for rec in producer['target_classes']:
        T=tuple(bytes((r>>j)&1 for j in range(N)) for r in rec['canonical_rows'])
        key,order,_=canonical(T);pid=rec['target_class_id'];need(key in member,'producer target is an actual boundary class')
        need(key not in producer_key.values(),'producer does not split an ordinary isomorphism class')
        producer_key[pid]=key;pid_to_independent[pid]=ids[key]
        group=automorphisms(T);e,orbits=error(T,group)
        need(rec['automorphism_order']==len(group) and rec['error']==e==0,'producer exact automorphism/error data')
        need(rec['vertex_orbits']==orbits and rec['stable_cells']==[list(c) for c in stable(T)],'producer actual orbit/cell coordinates')
        need(rec['incident_edges_from_positive_slice']==len(member[key]),'producer target incidence count')
    need(len(producer_key)==len(keys),'producer class partition complete and unsplit')
    expected={edge['word']:(edge,key) for edge,key,order in witnesses}
    seen=set()
    for rec in producer['boundary_edge_witnesses']:
        word=rec['source_word'];need(word not in seen and word in expected,'producer exact boundary source word');seen.add(word)
        edge,key=expected[word];T=decode(rec['target_word']);order=rec['canonical_vertex_order'];pid=rec['target_class_id']
        need(sorted(order)==list(range(N)),'producer vertex-order permutation')
        need(rec['source_index']==edge['source_index'] and rec['target_word']==edge['complement_word'],'producer exact source/target pairing')
        need(producer_key[pid]==key,'producer full target class identification')
        raw=producer['target_classes'][pid]['canonical_rows']
        need(rows(normalize(T,order))==raw,'producer actual vertex-coordinate arc-isomorphism witness')
    need(len(seen)==4020,'all 4020 producer target witnesses verified')
    adjacency=Counter((edge['source_index'],ids[key]) for edge,key,order in witnesses)
    translated=sorted([a,pid_to_independent[b],m] for a,b,m in producer['boundary_adjacency'])
    need(translated==[[a,b,m] for (a,b),m in sorted(adjacency.items())],'entire multiplicity adjacency after label alignment')
    vertices={('P',n['source_representative_index']) for n in wiring['nodes']}|{('Z',i) for i in range(len(keys))}
    pairs=[(('P',a),('P',b)) for a,b,m in wiring['adjacency_multiplicities'] for _ in range(m)]
    pairs.extend((('P',a),('Z',b)) for (a,b),m in adjacency.items() for _ in range(m))
    comps=[{'positive_nodes':[v for t,v in c if t=='P'],'independent_zero_nodes':[v for t,v in c if t=='Z']} for c in components(vertices,pairs)]
    candidate_components=sorted((tuple(c['positive_nodes']),tuple(sorted(pid_to_independent[v] for v in c['zero_error_nodes']))) for c in producer['components_of_incident_edge_graph'])
    need(candidate_components==sorted((tuple(c['positive_nodes']),tuple(c['independent_zero_nodes'])) for c in comps),'complete incident graph components after exact label alignment')
    need(sorted((len(c['positive_nodes']),len(c['independent_zero_nodes'])) for c in comps)==[(4,238),(20,793)],'two scoped incident-edge components')
    need(len(vertices)==1055 and len(pairs)==4110,'exact incident graph vertices and edges')
    for a,b in [('boundary_class_stable_patterns',[[list(k),v] for k,v in sorted(patterns.items())]),('boundary_presentation_stable_patterns',[[list(k),v] for k,v in sorted(presentation_patterns.items())]),('boundary_class_automorphism_orders',[list(kv) for kv in sorted(aut.items())]),('boundary_adjacency_multiplicity_histogram',[list(kv) for kv in sorted(Counter(adjacency.values()).items())]),('boundary_target_incidence_degree_histogram',[list(kv) for kv in sorted(Counter(len(v) for v in member.values()).items())])]:need(producer[a]==b,'producer histogram: '+a)
    result={'status':'PASS','scope':'complete blue edges incident to the 24 positive-error order-ten self-converse nodes, with all boundary target classes identified; no zero-zero or black edges','source_pins':pins,'method':'exhaustive ALL within-actual-stable-cell permutations; complete directed matrices; exact full automorphism groups; no producer cycle/rotation restriction','positive_nodes':24,'distinct_zero_error_boundary_nodes':len(keys),'all_neighborhood_vertices':len(vertices),'all_edges_incident_to_positive_slice':len(pairs),'internal_positive_edges':90,'boundary_edges':4020,'canonical_search_sizes':sorted(search.items()),'boundary_presentation_stable_patterns':[[list(k),v] for k,v in sorted(presentation_patterns.items())],'boundary_class_stable_patterns':[[list(k),v] for k,v in sorted(patterns.items())],'boundary_class_automorphism_orders':sorted(aut.items()),'boundary_adjacency_multiplicity_histogram':sorted(Counter(adjacency.values()).items()),'boundary_target_incidence_degree_histogram':sorted(Counter(len(v) for v in member.values()).items()),'producer_to_independent_target_ids':sorted(pid_to_independent.items()),'boundary_adjacency':[[a,b,m] for (a,b),m in sorted(adjacency.items())],'components_of_incident_edge_graph':comps,'target_classes':targets,'boundary_edge_witnesses':[{'source_index':edge['source_index'],'source_word':edge['word'],'target_word':edge['complement_word'],'independent_target_id':ids[key],'independent_vertex_order':order} for edge,key,order in witnesses]}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
