#!/usr/bin/env python3
"""Independent E33 converse-repair and native-counterfeit audit.
No producer-checker imports. Explicit input paths support scratch and published layouts.
"""
import argparse,array,csv,gzip,hashlib,json,sys
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
HERE=Path(__file__).resolve().parent
CHECKS=Counter()
def need(ok,label):
    CHECKS[label]+=1
    if not ok: raise RuntimeError(label)
def tsv(path):
    with gzip.open(path,'rt') as f:
        yield from csv.DictReader(f,delimiter='\t')
def clean(c): return {k:v for k,v in c.items() if v}
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--wiring',type=Path,default=HERE.parent/'wiring');ap.add_argument('--theory',type=Path,default=HERE.parent/'theory');ap.add_argument('--write',action='store_true');args=ap.parse_args()
    wiring=args.wiring;theory=args.theory
    classes=list(tsv(wiring/'classes.tsv.gz'));need(len(classes)==6880,'class_count')
    canon=[int(r['canon']) for r in classes];conv=[int(r['converse_class']) for r in classes];node=[int(r['node_id']) for r in classes];c3=[int(r['c3']) for r in classes]
    need([int(r['class_id']) for r in classes]==list(range(6880)),'class_id_order')
    need(all(conv[conv[i]]==i for i in range(6880)),'converse_involution')
    data=[];sigs=[];pairs=list(combinations(range(8),2));alltags=set()
    for ci,code in enumerate(canon):
        A=[[0]*8 for _ in range(8)]
        for k,(u,v) in enumerate(pairs): A[u][v]=(code>>(27-k))&1;A[v][u]=1-A[u][v]
        c=[0]*8;e=[[0]*8 for _ in range(8)]
        for a,b,d in combinations(range(8),3):
            if A[a][b]==A[b][d]==A[d][a]:
                for v in (a,b,d):c[v]+=1
                for u,v in combinations((a,b,d),2):e[u][v]+=1;e[v][u]+=1
        out=[[v for v in range(8) if A[u][v]] for u in range(8)]
        need(sum(c)==3*c3[ci],'independent_triangle_counts')
        def E(p=None,q=None):
            return tuple(sorted(tuple(sorted((c[v] if p is None else c[v]%p,e[u][v] if q is None else e[u][v]%q) for v in out[u])) for u in range(8)))
        S={'E43':E(4,3),'E33':E(3,3),'E_neighbor3':E(3,None),'E_neighbor2':E(2,None),'E_edge2':E(None,2)}
        for row in S['E33']: alltags.update(row)
        data.append((A,c,e,out));sigs.append(S)
    def signature(ci,name):
        if name not in sigs[ci]:
            A,c,e,out=data[ci]
            if name=='J4':
                v4=[0]*8
                for B in combinations(range(8),4):
                    if sorted(sum(A[u][v] for v in B) for u in B)==[1,1,2,2]:
                        for u in B:v4[u]+=1
                tags=[(len(out[u]),c[u],v4[u]) for u in range(8)]
                sigs[ci][name]=tuple(sorted((tags[u],tuple(sorted(tags[v] for v in out[u]))) for u in range(8)))
            else:raise RuntimeError('unknown signature '+name)
        return sigs[ci][name]
    groups=defaultdict(list);exactgroups=defaultdict(list)
    for i in range(6880):groups[signature(i,'E33')].append(i);exactgroups[signature(i,'E_neighbor3')].append(i)
    partition=sorted(sorted(g) for g in groups.values());exactpartition=sorted(sorted(g) for g in exactgroups.values())
    need(len({signature(i,'E43') for i in range(6880)})==6880,'E43_all_ordinary_classes_separate');need(len(groups)==6876,'E33_partition_count');need(partition==exactpartition,'E33_and_exact_edge_same_partition')
    bad=[i for i in range(6880) if len(groups[signature(i,'E33')])>1]
    need(len(bad)==8,'E33_collision_class_count')
    need(all(len(groups[signature(conv[i],'E33')])==1 for i in bad),'all_ambiguous_targets_singleton_converse')
    need(len(alltags)==9,'all_nine_E33_symbols_used')
    # Independent forest certificate including parallel class edges.
    parent={};degrees=Counter()
    def find(x):
        parent.setdefault(x,x)
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    for i in range(6880):
        a=('L',signature(i,'E33'));b=('R',signature(conv[i],'E33'))
        ra,rb=find(a),find(b);need(ra!=rb,'E33_signature_multigraph_acyclic');parent[ra]=rb;degrees[a]+=1;degrees[b]+=1
    need(all(degrees[('L',signature(i,'E33'))]==1 or degrees[('R',signature(conv[i],'E33'))]==1 for i in range(6880)),'every_E33_edge_has_leaf_endpoint')
    receipt=json.loads((theory/'residual_counterfeits.json').read_text());results=receipt['results']
    modes={};needed=set();needed_merged=set();source_sets=defaultdict(set)
    for name,obj in results.items():
        for mode in ('ordinary','symmetric'):
            if mode not in obj:continue
            edits=obj[mode]['edits'];delta=Counter()
            for z in edits:
                key=tuple(z[x] for x in ('source_class','target_class','start_score','end_score','blue'))
                need(key not in delta,'witness_no_duplicate_edit');delta[key]=z['delta'];needed.add(key);needed_merged.add((node[key[0]],node[key[1]],key[4]));source_sets[name].add(key[0])
            modes[name,mode]=(edits,delta)
    native={};merged=Counter();source_margins=Counter();native_target=Counter();mass=0
    for r in tsv(wiring/'ports.tsv.gz'):
        key=tuple(int(r[x]) for x in ('source_class','target_class','start_score','end_score','blue'));i,j,a,b,blue=key;z=int(r['multiplicity']);mass+=z
        if key in needed:native[key]=z
        mk=(node[i],node[j],blue)
        if mk in needed_merged:merged[mk]+=z
        native_target[j,8-a,6-b,blue]+=z
        for name,sources in source_sets.items():
            if i in sources:source_margins[name,i,a,b,blue,signature(j,name)]+=z
    need(mass==2**20,'native_port_mass')
    target_caps={}
    for r in tsv(wiring/'endpoint_capacities.tsv.gz'):
        blue=int(r['blue'] in ('True','true','1'));key=(int(r['class_id']),int(r['start_score']),int(r['end_score']),blue)
        target_caps[key]=int(r['open'])
    need(all(native_target[k]==v for k,v in target_caps.items()),'native_target_capacities_match')
    mapping=array.array('H',gzip.decompress((wiring/'class_of_mask.u16le.gz').read_bytes()))
    if sys.byteorder!='little':mapping.byteswap()
    need(len(mapping)==2**21,'mask_map_size')
    tiles=[(u,v) for u in range(6) for v in range(u+2,8)];close=tiles.index((0,7));full=2**21-1
    refl=[tiles.index((7-v,7-u)) for u,v in tiles]
    def mirror(k):i,j,a,b,blue=k;return(conv[i],conv[j],7-b,7-a,blue)
    def reflect_mask(m):return sum(((m>>k)&1)<<refl[k] for k in range(21))
    def endpoints(m):return (7-sum((m>>k)&1 for k,(u,v) in enumerate(tiles) if v==7),sum((m>>k)&1 for k,(u,v) in enumerate(tiles) if u==0))
    outcome={}
    for (name,mode),(edits,delta) in modes.items():
        row=Counter();col=Counter();directed=Counter();undirected=Counter();degree=Counter();current=Counter()
        for z in edits:
            key=tuple(z[x] for x in ('source_class','target_class','start_score','end_score','blue'));i,j,a,b,blue=key;v=delta[key];before=native.get(key,0)
            need(before==z['before'] and before+v==z['after'],'edit_matches_native_before_after')
            need(before+v>=0,'changed_native_entries_nonnegative');need(blue==0,'blue_entries_unchanged')
            need(z['source_canon']==canon[i] and z['target_canon']==canon[j],'witness_class_canons')
            need(source_margins[name,i,a,b,blue,signature(j,name)]>0,'source_signature_row_exists')
            need(target_caps[j,8-a,6-b,blue]>0,'target_endpoint_column_exists')
            row[i,a,b,blue,signature(j,name)]+=v;col[j,8-a,6-b,blue]+=v
            directed[node[i],node[j],blue]+=v;undirected[min(node[i],node[j]),max(node[i],node[j]),blue]+=v
            degree[node[i]]+=v;degree[node[j]]+=v
            need(c3[j]-c3[i]==a-b-1,'changed_cell_triangle_law')
            current[node[i]]+=v*(c3[j]-c3[i]);current[node[j]]+=v*(c3[i]-c3[j])
            masks=z['native_masks'];need(len(masks)==len(set(masks)),'literal_masks_distinct')
            if v<0:need(len(masks)>=-v,'enough_removed_literal_masks')
            for m in masks:
                need((m>>close)&1,'literal_mask_closed')
                mm=sum(((m>>k)&1)<<refl[k] for k in range(21))
                aa=7-sum((m>>k)&1 for k,(u,w) in enumerate(tiles) if w==7);bb=sum((m>>k)&1 for k,(u,w) in enumerate(tiles) if u==0)
                need((mapping[m],mapping[m^full],aa,bb,int(mm==m))==key,'literal_mask_exact_native_cell')
                need(mapping[mm]==conv[i] and mapping[mm^full]==conv[j],'literal_mask_reflection_classes')
        need(not any(row.values()),'source_signature_margins_preserved');need(not any(col.values()),'target_endpoint_margins_preserved')
        need(not any(degree.values()),'merged_literal_degrees_preserved');need(not any(current.values()),'merged_triangle_current_preserved')
        d=clean(directed);u=clean(undirected);need(bool(d) and bool(u),'merged_projection_genuinely_changes')
        need([[*k,v] for k,v in sorted(d.items())]==results[name][mode]['directed_merged_delta'],'reported_directed_delta_exact')
        need([[*k,v] for k,v in sorted(u.items())]==results[name][mode]['undirected_merged_delta'],'reported_undirected_delta_exact')
        reflected=all(delta[mirror(k)]==v for k,v in delta.items())
        if mode=='symmetric':
            need(reflected,'counterfeit_full_native_reflection')
            need(all((merged[k]+v)%2==0 for k,v in d.items() if not k[2]),'counterfeit_directed_black_twins_even')
            need(all(k[4] or k!=mirror(k) or (native.get(k,0)+v)%2==0 for k,v in delta.items()),'counterfeit_fixed_black_cells_even')
            matched=results[name][mode]['literal_port_matching'];assignment={r['source_closed_mask']:r['counterfeit_target_open_mask'] for r in matched}
            need(len(assignment)==len(matched)==4,'literal_swap_four_distinct_sources')
            need(len(set(assignment.values()))==4,'literal_swap_distinct_open_targets')
            need(set(assignment.values())=={m^full for m in assignment},'literal_swap_preserves_open_port_inventory')
            induced=Counter()
            for r in matched:
                m=r['source_closed_mask'];old=r['native_target_open_mask'];new=r['counterfeit_target_open_mask'];a,b=endpoints(m)
                need(old==m^full,'literal_swap_native_complement')
                need((m>>close)&1 and not((new>>close)&1),'literal_swap_closed_to_open')
                need(m!=reflect_mask(m) and new!=reflect_mask(new),'literal_swap_black_to_black')
                need(endpoints(new)==(8-a,6-b),'literal_swap_endpoint_compatibility')
                need(mapping[m]==r['source_class'] and mapping[old]==r['native_target_class'] and mapping[new]==r['counterfeit_target_class'],'literal_swap_all_class_ids')
                need(signature(mapping[old],name)==signature(mapping[new],name),'literal_swap_retained_target_signature')
                need(reflect_mask(m) in assignment and assignment[reflect_mask(m)]==reflect_mask(new),'literal_swap_commutes_with_reflection')
                induced[mapping[m],mapping[old],a,b,0]-=1;induced[mapping[m],mapping[new],a,b,0]+=1
            need(clean(induced)==clean(delta),'literal_swap_induces_exact_table_counterfeit')
        outcome[name+'/'+mode]={'changed_native_cells':len(delta),'native_reflection_preserved':reflected,'undirected_merged_delta':[[*k,v] for k,v in sorted(u.items())]}
    for name in ('E_neighbor2','E_edge2','J4'):need((name,'symmetric') in modes,'required_symmetric_counterfeit_present')
    need(('E_neighbor3','ordinary') in modes,'neighbor3_ordinary_counterfeit_present')
    need(('E_neighbor3','symmetric') not in modes,'no_false_neighbor3_symmetric_counterfeit')
    # The ordinary neighbor3 witness cannot be repaired merely by adding its mirror.
    D=modes['E_neighbor3','ordinary'][1];sym=D.copy()
    for k,v in D.items():sym[mirror(k)]+=v
    rows=Counter()
    for (i,j,a,b,blue),v in sym.items():rows[i,a,b,blue,signature(j,'E_neighbor3')]+=v
    need(any(rows.values()),'naive_mirror_addition_rejected_by_signature_margins')
    # A general 2x2 marginal example: distinct joint labels do not imply uniqueness.
    X=(1,1,1,1);Y=(2,0,0,2);edges=((0,0),(0,1),(1,0),(1,1))
    for side in (0,1):need([sum(x for x,e in zip(X,edges) if e[side]==i) for i in (0,1)]==[sum(y for y,e in zip(Y,edges) if e[side]==i) for i in (0,1)],'paired_label_injectivity_not_enough')
    inputs=[wiring/n for n in ('classes.tsv.gz','ports.tsv.gz','endpoint_capacities.tsv.gz','class_of_mask.u16le.gz')]+[theory/'residual_counterfeits.json']
    out={'status':'PASS','scope':'Independent finite E33 forest/singleton certificate and literal counterfeit verification. No tournament realizability claim.',
         'E33':{'classes':6880,'distinct_signatures':len(groups),'used_symbols':len(alltags),'direct_singleton_classes':6880-len(bad),'converse_repaired_classes':len(bad),
                'collision_pairs':[[canon[i] for i in g] for g in partition if len(g)>1],
                'converse_certificates':[{'class_id':i,'canon':canon[i],'converse_class':conv[i],'converse_canon':canon[conv[i]],'converse_fiber_size':len(groups[signature(conv[i],'E33')])} for i in bad],
                'signature_multigraph_vertices':len(parent),'signature_multigraph_edges':6880,'signature_multigraph_components':len({find(v) for v in parent}),'maximum_degree':max(degrees.values())},
         'counterfeits':outcome,'conclusion':'Minimum alphabet size9 for reflection-constrained native transport in this modular family; minimum12 for ordinary class separation and unreflected transport. Both are order-eight finite statements.',
         'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'checks':dict(CHECKS)}
    dest=HERE/'reflection_audit_receipt.json'
    if args.write:dest.write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads(dest.read_text())==out,'frozen_receipt_exact_match')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
