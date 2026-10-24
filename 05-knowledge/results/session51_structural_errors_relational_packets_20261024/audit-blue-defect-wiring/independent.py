#!/usr/bin/env python3
"""Independent complete-HP / literal 10x10 matrix blue wiring audit.
No producer code imports, no anti-involution enumeration, no 36-bit pairing.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json
HERE=Path(__file__).resolve().parent
N=10

def need(test,reason):
    if not test: raise RuntimeError(reason)

def hp_paths(matrix):
    """Every directed HP, by ordinary adjacency DFS, without symmetry pruning."""
    successors=[tuple(j for j in range(N) if matrix[i][j]) for i in range(N)]
    def visit(path,unused):
        if not unused:
            yield path
            return
        for nxt in successors[path[-1]]:
            if nxt in unused:
                yield from visit(path+(nxt,),unused-{nxt})
    vertices=set(range(N))
    for start in range(N): yield from visit((start,),vertices-{start})

def is_blue(matrix,path):
    # Exact anti-reflection condition, including all directed pairs.
    return all(matrix[path[i]][path[j]] != matrix[path[N-1-i]][path[N-1-j]]
               for i in range(N) for j in range(i+1,N))

def normalize(matrix,path):
    return tuple(bytes(matrix[a][b] for b in path) for a in path)

def line_partner(matrix):
    """Literally reverse every non-path arc, retaining consecutive path arcs."""
    return tuple(bytes(0 if i==j else matrix[i][j] if abs(i-j)==1 else 1-matrix[i][j]
                       for j in range(N)) for i in range(N))

def cyclic(matrix):
    return sum(matrix[a][b]*matrix[b][c]*matrix[c][a]+matrix[a][c]*matrix[c][b]*matrix[b][a]
               for a,b,c in combinations(range(N),3))

def stable_cells(matrix):
    cells=[tuple(range(N))]
    while True:
        refined=[]
        for cell in cells:
            buckets={}
            for v in cell:
                sig=tuple(sum(matrix[v][w] for w in target) for target in cells)
                buckets.setdefault(sig,[]).append(v)
            refined.extend(tuple(buckets[key]) for key in sorted(buckets))
        if len(refined)==len(cells): return tuple(refined)
        cells=refined

def encoded_for_comparison(matrix):
    # This is used only to compare output to the producer's external format.
    return sum(matrix[i][j]<<k for k,(i,j) in enumerate(( (a,b) for a in range(N) for b in range(a+2,N))))

def components(vertices,edge_pairs):
    parent={v:v for v in vertices}
    def find(v):
        while parent[v]!=v:
            parent[v]=parent[parent[v]]; v=parent[v]
        return v
    for a,b in edge_pairs:
        if a in parent and b in parent: parent[find(a)]=find(b)
    groups={}
    for v in sorted(vertices):groups.setdefault(find(v),[]).append(v)
    return sorted(groups.values())

def main():
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for name,digest in pins.items():need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,'source SHA256 mismatch: '+name)
    need(pins['sources/candidate-MANIFEST.json']=='cd89cad144786c67f076c235409c90e3b4de21e0a5f81926139184068f1d676d','candidate frozen identity')
    need(pins['sources/converse-MANIFEST.json']=='30b62a8fdf93065b19f071050d9fcfe5ad08c5f79cfcfbf30afcdae53babde5a','inherited base identity')
    need(pins['sources/inherited-audit-MANIFEST.json']=='7ca02dc5c0272039b552e7cd5538bbab04923a894e04c050e08fb69de721fc7f','inherited audit identity')
    src=json.loads((HERE/'sources/converse-receipt.json').read_text())
    portrows=json.loads((HERE/'sources/blue-port-receipt.json').read_text())['representatives']
    expected={r['representative_index']:r for r in portrows}
    records={r['index']:r for r in src['representatives'] if r['self_converse']}
    need(len(records)==24,'exactly 24 inherited positive self-converse nodes')
    carrier={};counts=Counter();hps=Counter();cross={};matrices={}
    for idx,rec in sorted(records.items()):
        matrix=tuple(bytes((row>>j)&1 for j in range(N)) for row in rec['rows'])
        need(all(matrix[i][i]==0 for i in range(N)) and all(matrix[i][j]+matrix[j][i]==1 for i,j in combinations(range(N),2)),'input actual tournament')
        need(rec['error']==40 and rec['automorphism_order']==1,'inherited source rigid error 40')
        degree=sorted(map(sum,matrix))
        need(degree in [[3]*5+[6]*5,[4]*5+[5]*5],'source degrees have certified types')
        cross[idx]=degree[0]-2;matrices[idx]=matrix
        for path in hp_paths(matrix):
            hps[idx]+=1
            if not is_blue(matrix,path):continue
            board=normalize(matrix,path)
            need(all(board[i][i+1]==1 for i in range(N-1)),'normalized base path')
            need(all(board[i][j]==board[N-1-j][N-1-i] for i in range(N) for j in range(N)),'literal transpose-reflection exact blue')
            need(board not in carrier,'actual labelled blue presentation collision')
            carrier[board]=(idx,path);counts[idx]+=1
        need(hps[idx]==expected[idx]['Hamiltonian_paths'],'independent full HP count')
        need(counts[idx]==expected[idx]['blue_ports'],'independent full-HP filtered blue count')
    need(sum(hps.values())==341504 and len(carrier)==4200,'complete full HP/blue totals')
    internal=[];boundary=[];degree=Counter();outward=Counter();transport=Counter();boundary_cells=Counter()
    by_word={encoded_for_comparison(board):value for board,value in carrier.items()}
    need(len(by_word)==len(carrier),'external word encoding injective')
    for board,(idx,path) in sorted(carrier.items()):
        partner=line_partner(board)
        need(partner!=board and line_partner(partner)==board,'literal C free involution')
        need(all(partner[i][i+1]==1 for i in range(N-1)),'C keeps base path')
        need(all(partner[i][j]==partner[N-1-j][N-1-i] for i in range(N) for j in range(N)),'R and C commute on blue carrier')
        # Also check commutation on the complete normalized board, not a heuristic.
        reflect=lambda a: tuple(bytes(a[N-1-j][N-1-i] for j in range(N)) for i in range(N))
        need(reflect(line_partner(board))==line_partner(reflect(board)),'literal R C commute')
        ds=list(map(sum,board));dt=list(map(sum,partner))
        need(dt==[N-1-ds[v]+(v==0)-(v==N-1) for v in range(N)],'actual endpoint-degree identity')
        need(sum(x*x for x in dt)-sum(x*x for x in ds)==2*(ds[-1]-ds[0]+1),'squared score transport')
        td=cyclic(partner)-cyclic(board)
        need(td==ds[0]-ds[-1]-1,'literal triangle transport')
        inside=partner in carrier
        transport[(cross[idx],ds[0],ds[-1],td,inside)]+=1
        word=encoded_for_comparison(board);otherword=encoded_for_comparison(partner)
        if inside:
            oid,opath=carrier[partner];degree[idx]+=1
            need(cross[idx]==cross[oid]==2 and ds[0]==5 and ds[-1]==4,'internal necessary endpoint filter')
            need(normalize(matrices[oid],opath)==partner,'internal actual target source path')
            if word<otherword:
                internal.append({'source_index':idx,'target_index':oid,'word':word,'complement_word':otherword,'source_path':list(path),'target_path':list(opath)})
        else:
            outward[idx]+=1
            boundary.append({'source_index':idx,'word':word,'complement_word':otherword})
            boundary_cells[tuple(sorted(map(len,stable_cells(partner))))]+=1
    internal.sort(key=lambda e:e['word']);boundary.sort(key=lambda e:e['word'])
    pairs=[(e['source_index'],e['target_index']) for e in internal]
    multiplicities=Counter(tuple(sorted(pair)) for pair in pairs)
    vertices=set(records);parts=components(vertices,pairs);baseline=len(parts)
    bridges=[k for k in range(len(pairs)) if len(components(vertices,pairs[:k]+pairs[k+1:]))>baseline]
    articulations=[v for v in sorted(vertices) if len(components(vertices-{v},pairs))>baseline]
    need(len(internal)==90 and sum(degree.values())==180 and len(boundary)==4020,'edge-port-boundary totals')
    need(sum(a==b for a,b in pairs)==4,'four actual loops')
    need(sorted(len(c) for c in parts)==[1,1,1,1,20],'connected components')
    need({c[0] for c in parts if len(c)==1}=={i for i in records if cross[i]==1},'exact degree-one isolated sector')
    need(max(multiplicities.values())==3 and len(multiplicities)==70,'parallel edge multiplicities')
    need(not bridges and not articulations and len(internal)-len(vertices)+baseline==71,'induced topology')
    need(sum(v for (c,s,t,d,inside),v in transport.items() if c==2 and s==5 and t==4)==2037,'necessary endpoint condition exact insufficiency')
    for i in vertices:need(counts[i]==degree[i]+outward[i],'per-node port conservation including loop double count')
    result={'status':'PASS','scope':'only the actual induced blue multigraph on the 24 positive-error self-converse order-ten nodes; no boundary-class identification or black wiring','method':'enumerate all directed HPs by DFS; literal 10x10 normalized matrices; retain path arcs and flip all other arcs; exact carrier membership; disjoint-set topology','source_pins':pins,'positive_nodes':len(vertices),'all_Hamiltonian_paths':sum(hps.values()),'blue_ports':len(carrier),'induced_blue_edges':len(internal),'induced_blue_loop_edges':sum(a==b for a,b in pairs),'induced_nonloop_blue_edges':sum(a!=b for a,b in pairs),'internal_blue_ports':sum(degree.values()),'boundary_blue_edges':len(boundary),'adjacency_multiplicities':[[a,b,m] for (a,b),m in sorted(multiplicities.items())],'adjacency_multiplicity_histogram':sorted(Counter(multiplicities.values()).items()),'induced_cycle_rank':len(internal)-len(vertices)+baseline,'bridge_edge_indices':bridges,'articulation_node_indices':articulations,'components':parts,'induced_degree_histogram':sorted(Counter(degree[i] for i in records).items()),'endpoint_transport_histogram':[[list(k),v] for k,v in sorted(transport.items())],'boundary_stable_cell_histogram':[[list(k),v] for k,v in sorted(boundary_cells.items())],'nodes':[{'source_representative_index':i,'cross_degree':cross[i],'blue_ports':counts[i],'internal_ports':degree[i],'boundary_edges':outward[i]} for i in sorted(records)],'full_HP_counts':[[i,hps[i]] for i in sorted(records)],'internal_edges':internal,'boundary_edges':boundary}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
