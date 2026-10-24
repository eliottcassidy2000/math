#!/usr/bin/env python3
"""Bounded commit-pinned arithmetic review, not a native certificate replay."""
from collections import defaultdict,deque,Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse,csv,gzip,io,json,subprocess

REVISION='06306318651dfcddbdff9e9790794f11b0e39f94'
PREFIX='05-knowledge/results/balanced_native_20261023/'
CHECKS=0

def need(ok,why):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(why)


def rational_reduce(rows):
    a=[[Fraction(x) for x in row] for row in rows]
    rank=0;det=Fraction(1);sign=1
    for col in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a)) if a[i][col]),None)
        if pivot is None:continue
        if pivot!=rank:a[pivot],a[rank]=a[rank],a[pivot];sign=-sign
        p=a[rank][col];det*=p;a[rank]=[x/p for x in a[rank]]
        for i in range(rank+1,len(a)):
            q=a[i][col]
            a[i]=[x-q*y for x,y in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank,det*sign


def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--out',required=True)
    args=p.parse_args();pins={}
    def raw(path):
        data=subprocess.check_output(['git','-C',args.repo,'show',REVISION+':'+path])
        pins[path]=sha256(data).hexdigest();return data
    def read(path):return json.loads(raw(path))
    for doc in ['SIX_TRADE_LATTICE.md','NATIVE18.md','CHARGED_RETURN.md','MULTITARGET.md','six_trade.py']:
        raw(PREFIX+doc)
    inventory=read('05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json')
    opts={r['source']:r for r in inventory['records']};sources=sorted(opts)
    need(len(sources)==28,'original source universe')
    labels=[tuple(opts[s]['source_label']) for s in sources];moduli=(2,2,2,2,128)
    need(len(set(labels))==28,'labels not distinct')
    def total(indices):
        indices=list(indices)
        return tuple(sum(labels[i][j] for i in indices)%moduli[j] for j in range(5))
    moves=set();counts={}
    for k in (2,3):
        buckets=defaultdict(list)
        for subset in combinations(range(28),k):buckets[total(subset)].append(subset)
        found=set()
        for bucket in buckets.values():
            for left,right in combinations(bucket,2):
                if set(left)&set(right):continue
                vector=tuple(int(i in left)-int(i in right) for i in range(28))
                if next(x for x in vector if x)<0:vector=tuple(-x for x in vector)
                found.add(vector)
        counts[2*k]=len(found);moves|=found
    need(counts=={4:27,6:1917},'complete small-move census')
    basis=read(PREFIX+'six_trade_basis.json')
    need(basis['sources']==sources,'incoming basis source order')
    need(set(map(tuple,basis['trades']))==moves,'incoming full move set differs')
    selected=[]
    for reduced,expr in zip(basis['basis'],basis['basis_trade_combinations']):
        need(len(expr)==1 and next(iter(expr.values()))==1,'incoming basis not direct move selection')
        i=int(next(iter(expr)));v=basis['trades'][i]
        need(tuple(v) in moves and v[:27]==reduced,'basis column not actual selected move')
        need(v[27]==-sum(reduced),'lost last-coordinate reconstruction')
        selected.append(v)
    need(len(selected)==27 and len(set(map(tuple,selected)))==27,'distinct basis columns')
    matrix=[list(row) for row in zip(*[v[:27] for v in selected])]
    rank,det=rational_reduce(matrix)
    need(rank==27 and abs(det)==2048,'independent rational determinant')
    origin=labels[0]
    generators=[tuple((x-y)%m for x,y,m in zip(label,origin,moduli)) for label in labels[1:]]
    subgroup={tuple([0]*5)};todo=deque(subgroup)
    while todo:
        value=todo.popleft()
        for g in generators:
            after=tuple((a+b)%m for a,b,m in zip(value,g,moduli))
            if after not in subgroup:subgroup.add(after);todo.append(after)
    need(len(subgroup)==2048,'charge-image size based at first source')
    four=read(PREFIX+'four_move.normal.json')
    c=four['small_exact_invariant'];need(four['sources']==sources and len(c)==28,'invariant coordinates')
    four_moves=[v for v in moves if sum(x!=0 for x in v)==4]
    need(rational_reduce(four_moves)[0]==23,'four-move rank')
    for v in four_moves:need(sum(a*b for a,b in zip(c,v))==0,'small invariant not annihilator')
    native=read(PREFIX+'native18.normal.json');A=set(native['plus']);B=set(native['minus'])
    need(len(A)==len(B)==9 and not A&B,'incoming endpoint support')
    need(A=={176,1394,1457,1477,1532,1923,2087,2176,3369},'incoming positive support')
    need(B=={1425,1549,1841,1921,2185,2553,2835,2905,3000},'incoming negative support')
    f=lambda S:sum(c[sources.index(s)] for s in S)
    need((f(B),f(A))==(-3,-1),'new invariant endpoint values')
    old=lambda S:len(S&{1638,1974})
    need((old(B),old(A))==(0,0),'our earlier invariant is not blind')
    cost={s:sum(a!=b for a,b in zip(opts[s]['before']['path'],opts[s]['after']['path'])) for s in sources}
    need(sum(cost[s] for s in A|B)==45,'incoming endpoint literal clock cost')
    need(total(sources.index(s) for s in A)==total(sources.index(s) for s in B),'incoming endpoint charge')
    route=read(PREFIX+'route18.normal.json');need(route['sources']==sources,'incoming route source order')
    states=[{s for i,s in enumerate(sources) if mask>>i&1} for mask in route['route']]
    need(states[0]==B and states[-1]==A,'route endpoints')
    need(all(S-(A|B)<= {1638} for S in states),'unexpected helper source')
    for S in states:
        need(len(S)==9,'route integer margin')
        need(total(sources.index(s) for s in S)==total(sources.index(s) for s in B),'route group guard')
    support=[len(S^T) for S,T in zip(states,states[1:])]
    work=[sum(cost[s] for s in S^T) for S,T in zip(states,states[1:])]
    need(support==[6,6,6,6] and work==[14,15,15,13] and sum(work)==57,'incoming route costs')
    table=read('05-knowledge/results/session45_injective_multisource_moves_20261023/dyadic-constraints/node_labels.json')
    fiber=[r['vertices'] for r in table['fibres'] if r['label']==[1,1,0,1,18]]
    need(fiber==[[570,1965,2234]],'exact original target-charge fiber')
    domains=read('05-knowledge/results/session46_native_packet_decomposition_20261024/pinned-domains/normal.json')
    paths=[r for s in domains['domains'] for r in s['complete_pinned_paths']]
    need(len(paths)==57 and Counter(r['target'] for r in paths)=={1965:28,2234:28,1559:1},
         'inherited complete pinned target inventory')
    need(not any(r['target']==570 for r in paths),'unlicensed third target path')
    inherited='05-knowledge/results/session45_injective_multisource_moves_20261023/source-circuits/'
    menus=read(inherited+'full-search/all_menus.json')
    census=read(inherited+'receipt.normal.json')
    raw(inherited+'search_all_sources.py');raw(inherited+'check.py')
    need(len(menus)==700,'inherited complete menu count')
    need(not any(570 in r['targets'] for r in menus),'inherited grammar admits target570 menu')
    pair_sources=sorted(r['source'] for r in menus if set(r['targets'])=={1965,2234})
    need(pair_sources==sources,'inherited pair sources not same exact28')
    need(census['universe']['canonical_sources']==3528 and
         census['universe']['actual_source_paths']==1102392 and
         census['universe']['black_endpoint_arc_menus']==700,'inherited completeness scope receipt')
    original='05-knowledge/results/metagraph_order8_20261023/'
    nodes=list(csv.DictReader(io.StringIO(raw(original+'merged_nodes.tsv').decode()),delimiter='\t'))
    node=next(row for row in nodes if row['merged_id']=='570')
    need((int(node['class_a']),int(node['class_b']),int(node['fibre']),int(node['blue']))==(732,736,590,0),
         'target570 original node inventory')
    edges=list(csv.DictReader(io.StringIO(gzip.decompress(raw(original+'merged_edges.tsv.gz')).decode()),delimiter='\t'))
    incident=[row for row in edges if '570' in (row['source'],row['target'])]
    neighbors={int(row['target'] if row['source']=='570' else row['source']) for row in incident}
    need(len(incident)==len(neighbors)==213 and 570 not in neighbors,'target570 original adjacency')
    need(sum(int(row['black']) for row in incident)==590 and
         sum(int(row['blue']) for row in incident)==0,'target570 genuine black incidence')
    eligible={tuple(sorted(row['targets'])) for row in menus if set(row['targets'])<=set(fiber[0])}
    need(eligible=={(1965,2234)},'exact grammar-dependent target graph')
    # These are frozen receipt comparisons, not re-executions of their programs.
    pairs=[]
    for name in ['four_move','multi_target','native18','route18','route18_audit','six_trade']:
        need(raw(PREFIX+name+'.normal.json')==raw(PREFIX+name+'.optimized.json'),
             'incoming normal/optimized pair mismatch '+name);pairs.append(name)
    report={'status':'PASS','checks':CHECKS,'revision':REVISION,'input_sha256':pins,
            'four_and_six_move_counts':counts,'four_move_rank':23,'direct_basis_columns':27,
            'basis_determinant':int(abs(det)),'charge_image_order':len(subgroup),
            'new_invariant_values':[f(B),f(A)],'our_earlier_invariant_values':[old(B),old(A)],
            'incoming_endpoint_clock_cost':45,'incoming_route_supports':support,
            'incoming_route_clock_costs':work,'incoming_route_total_work':sum(work),
            'target_label_fiber':fiber[0],'current_pin_paths_to570':0,
            'inherited_complete_menu_scan':{'menus':700,'target570_menus':0,
                'same_pair_sources':pair_sources,'canonical_source_representatives':3528,
                'inherited_actual_paths':1102392,'full_census_rerun':False},
            'target570_original_metagraph':{'class_pair':[732,736],'fibre':590,'blue':0,
                'distinct_neighbors':213,'incident_entries':213,'black_weight':590,'loops':0},
            'pin_compatible_target_graph':{'vertices':fiber[0],'edges':[[1965,2234]],
                'edge_source_menus':28,'isolated_targets':[570]},
            'identical_receipt_pairs':pairs,
            'scope':'Independent bounded label/lattice/route arithmetic and receipt checks; no native Hamiltonian or edge-image replay'}
    Path(args.out).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','checks':CHECKS,'out':args.out}))

if __name__=='__main__':main()
