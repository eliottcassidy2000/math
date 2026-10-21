#!/usr/bin/env python3
"""Independent extension audit. Imports only the auditor's own pinned support.
No producer checker is imported or executed. Inputs remain read-only.
"""
import argparse
from collections import Counter,defaultdict
from itertools import combinations,permutations,product
from math import comb,factorial
from pathlib import Path
import hashlib,json
import core_audit_support as A

FROZEN='a6667f7a99d5fa186ba09b5b6df462c13d084599f1af2c0001aba9f9aff4c390'
BASE='a5490b17a42062db222676a44a1eb2ef1c3f0b55499b2389cc70691886f9112c'
need=A.require

def canonical(rows,weights):
    # Generic weighted tournament canonicalization, beginning with the
    # invariant outdegree rather than the producer's intrinsic-gap decoder.
    signatures=[(weights[v],rows[v].bit_count()) for v in range(len(rows))]
    unique=sorted(set(signatures));lookup={s:j for j,s in enumerate(unique)}
    colours=[lookup[s] for s in signatures]
    while True:
        cells=[[] for _ in range(max(colours)+1)]
        for v,c in enumerate(colours):cells[c].append(v)
        masks=[sum(1<<v for v in cell) for cell in cells]
        sigs=[(colours[v],tuple((rows[v]&m).bit_count() for m in masks)) for v in range(len(rows))]
        uniq=sorted(set(sigs));lookup={s:j for j,s in enumerate(uniq)}
        newer=[lookup[s] for s in sigs]
        if newer==colours:break
        colours=newer
    search=1
    for cell in cells:search*=factorial(len(cell))
    need(search<1000000,'canonical permutation universe')
    q=len(rows);width=(q+7)//8;best=None
    for pieces in product(*(permutations(cell) for cell in cells)):
        order=tuple(v for cell in pieces for v in cell)
        key=bytes([q])+bytes(weights[v] for v in order)+b''.join(sum(((rows[v]>>w)&1)<<j
                         for j,w in enumerate(order)).to_bytes(width,'little') for v in order)
        if best is None or key<best:best=key
    return best

def clones(rows):
    mask=(1<<len(rows))-1
    return [(a,b) for a,b in combinations(range(len(rows)),2)
            if not ((rows[a]^rows[b]) & (mask ^ (1<<a) ^ (1<<b)))]

def numbers(top):
    F=[1];I=[0];P=[0,0];H=[0,0]
    for k in range(1,top+1):
        F.append(F[-1]*(2*k-1));I.append(F[k]-sum(I[a]*F[k-a] for a in range(1,k)))
        if k>=2:
            conv=lambda n:sum(F[a]*F[n-a] for a in range(n+1))
            P.append((k-1)*F[k]-sum(P[j]*conv(k-j) for j in range(2,k)))
            H.append((2*k-2)*F[k-1]-sum(H[j]*conv(k-j) for j in range(2,k)))
    return F,I,P,H

def enumerate_records(k,d,strong,allowed='ABCD'):
    n=2*k;m=n-1 if strong else n+1
    for edges in A.pairings(tuple(range(n))):
        if strong and not A.cut_test(edges,n):continue
        if 'A' in allowed:
            for ds in A.words(d,m,k-1):yield 'A',edges,n,tuple(k-x for x in ds)
        if 'B' in allowed:
            for hole in range(m):
                for ds in A.words(d-k,m-1,k-1):
                    gaps=list(k-x for x in ds);gaps.insert(hole,0)
                    yield 'B',edges,n,tuple(gaps)
        if 'D' in allowed and d==2*k:
            for holes in combinations(range(m),2):
                yield 'D',edges,n,tuple(0 if j in holes else k for j in range(m))
    if 'C' in allowed and d>=k+1:
        for edges in A.stars(k):
            if strong and not A.cut_test(edges,2*k-1):continue
            for ds in A.words(d-k-1,m-1,k-1):
                yield 'C',edges,2*k-1,tuple(k-x for x in ds)

def check_layer(k,d,strong,allowed,F,I,P,H):
    fibres={};accepted=Counter();rejected=Counter();relabels=0
    for kind,edges,n,gaps in enumerate_records(k,d,strong,allowed):
        rows,weights=A.tournament(edges,n,gaps,strong)
        pair=clones(rows);actual=(0,)+gaps+(0,) if strong else gaps
        forbidden=any(b==a+1 and actual[b]==0 for a,b in edges)
        need(bool(pair)==forbidden,'complete literal omitted-gap clone criterion')
        if forbidden:
            rejected[kind]+=1;continue
        need(sum(weights)==(2*k*k+(k if strong else 3*k))-d,'literal weighted mass')
        if strong:need(A.strong_actual(rows),'literal strong connectivity')
        key=canonical(rows,weights)
        fibres[key]=fibres.get(key,'')+kind
        accepted[kind]+=1
        if accepted[kind]<=2:
            order=tuple(reversed(range(len(rows))))
            changed=tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)
            need(canonical(changed,tuple(weights[v] for v in order))==key,'canonical covariance')
            need(bool(A.homogeneous_pairs(rows))==bool(pair),'independent explicit pair predicate')
            relabels+=1
    histogram=Counter(','.join(sorted(v)) for v in fibres.values())
    owner=d==2*k
    need(set(histogram)<=({'A','B','C','D','C,C,D'} if owner else {'A','B','C'}),'exhaustive owner fibre shape')
    m=2*k-1 if strong else 2*k+1;E=I[k] if strong else F[k]
    J=(2*k-3)*I[k-1] if strong else F[k];p=P[k] if strong else (k-1)*F[k]
    h=H[k] if strong else (2*k-2)*F[k-1]
    K=(comb(2*k-4,2)*I[k-2] if k>=3 else 0) if strong else (k-1)*F[k-1]
    L=A.formula(k,d,strong,F,I,P)
    if owner:
        need(histogram['C,C,D']==h,'one exact triple fibre for every marked span-two edge')
        counts={'A':E*A.coefficient(d,m,k-1),
                'B':(m*E-J)*A.coefficient(k,m-1,k-1),
                'C':p*A.coefficient(k-1,m-1,k-1),
                'D':E*comb(m,2)-(m-1)*J+K}
        for kind in allowed:need(accepted[kind]==counts[kind],'bounded coefficient family presentation count')
        target=sum(counts[t] for t in allowed)-2*h
        need(len(fibres)==target,'owner quotient exact consequence')
        if allowed=='ABCD':need(target==L+K-2*h,'K minus two H correction')
    else:
        need(len(fibres)==sum(accepted.values())==L,'extended corridor class formula')
    return dict(k=k,delta=d,strong=strong,scope=allowed,classes=len(fibres),raw_records=sum(accepted.values()),
         accepted=dict(accepted),rejected=dict(rejected),fibre_histogram=dict(histogram),
         formal_corridor_expression=L,owner_correction=(K-2*h if owner else None),relabel_checks=relabels)

def local_rows(n,edge):
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    if edge:
        a,b=edge;rows[a]&=~(1<<b);rows[b]|=1<<a
    return tuple(rows)

def pointed(rows,port=None):
    best=None;automorphisms=0
    for mapping in permutations(range(len(rows))):
        if all(((rows[a]>>b)&1)==((rows[mapping[a]]>>mapping[b])&1)
               for a in range(len(rows)) for b in range(len(rows))) and (port is None or mapping[port]==port):
            automorphisms+=1
        if port is not None and mapping[0]!=port:continue
        value=tuple(sum(((rows[a]>>b)&1)<<j for j,b in enumerate(mapping)) for a in mapping)
        if best is None or value<best:best=value
    return best,automorphisms

def local_audit(source):
    literal=json.loads((source/'local-ports-check.json').read_text())
    full=[];mixed=[]
    for name,edge,port,degree in [('12',(0,2),4,0),('13',(0,4),2,2),('23',(2,4),0,4)]:
        rows=local_rows(5,edge);key,autos=pointed(rows,port)
        need(rows[port].bit_count()==degree,'literal full-cell distinguished degree')
        entry=next(x for x in literal['full_two_unit_cells'] if x['pair']==name)
        need(tuple(entry['rows'])==rows and entry['port']==port and entry['degree']==degree,'parse full-cell literal')
        need(tuple(entry['rooted_key'])==key,'full-cell exact pointed isomorphism key')
        full.append(dict(name=name,port=port,degree=degree,rooted_key=key,pointed_automorphisms=autos))
    need(len({r['rooted_key'] for r in full})==3,'distinct five-cell owners')
    for name,edge,port,degree in [('X13',(0,3),1,2),('X23',(1,3),0,3),('Y12',(0,2),3,0),('Y13',(0,3),2,1)]:
        rows=local_rows(4,edge);key,autos=pointed(rows,port);bare,_=pointed(rows)
        need(rows[port].bit_count()==degree,'literal mixed-cell distinguished degree')
        entry=next(x for x in literal['missing_plus_unit_cells'] if x['name']==name)
        need(tuple(entry['rows'])==rows and entry['port']==port and entry['degree']==degree,'parse mixed-cell literal')
        need(tuple(entry['rooted_key'])==key and tuple(entry['bare_key'])==bare,'mixed-cell exact isomorphism key')
        mixed.append(dict(name=name,port=port,degree=degree,rooted_key=key,bare_key=bare,pointed_automorphisms=autos))
    need(len({r['rooted_key'] for r in mixed})==4,'distinct four-cell owners')
    need(mixed[0]['bare_key']==mixed[3]['bare_key'] and mixed[0]['rooted_key']!=mixed[3]['rooted_key'],'port-erasure hostile')
    examples=json.loads((source/'extension-check.json').read_text())['three_record_collisions']
    collision_rows=[]
    for entry in examples:
        k=entry['k'];strong=entry['strong'];qs=[]
        for edge_name,gap_name,n in [('matching','matching_gaps',2*k),('p3_first','p3_gaps',2*k-1),('p3_second','p3_gaps',2*k-1)]:
            rows,weights=A.tournament(tuple(map(tuple,entry[edge_name])),n,tuple(entry[gap_name]),strong)
            supplied=entry['actual_tournaments'][len(qs)]
            need(tuple(supplied['rows'])==rows and tuple(supplied['weights'])==weights,'reconstruct literal whole collision')
            need(not clones(rows),'collision is not clone normalization')
            if strong:need(A.strong_actual(rows),'literal strong collision')
            qs.append((rows,weights))
        need(len({canonical(*q) for q in qs})==1,'literal three-record weighted isomorphism')
        start=0 if strong else 1;cell=list(range(start,start+3));source_port=start+1
        mappings=[]
        for target,(target_rows,target_weights) in enumerate(qs[1:]):
            target_port=start if target==0 else start+2;valid=[]
            for local in permutations(cell):
                mapping=list(range(len(qs[0][0])))
                for a,b in zip(cell,local):mapping[a]=b
                if mapping[source_port]!=target_port:continue
                if all(qs[0][1][v]==target_weights[mapping[v]] for v in range(len(mapping))) and all(
                    ((qs[0][0][a]>>b)&1)==((target_rows[mapping[a]]>>mapping[b])&1)
                    for a in range(len(mapping)) for b in range(len(mapping))):valid.append(mapping)
            need(len(valid)==1,'unique isomorphism fixes EVERY outside vertex and actual external port')
            mappings.append(valid[0])
        collision_rows.append(dict(k=k,strong=strong,matching_port=source_port,outside_fixed_isomorphisms=mappings))
    return dict(full_cells=full,mixed_cells=mixed,literal_collisions=collision_rows)

def main():
    p=argparse.ArgumentParser();here=Path(__file__).resolve().parent
    p.add_argument('--source',type=Path,default=here.parent/'session18-zero-heavy-extension')
    p.add_argument('--base',type=Path,default=here.parent/'session18-zero-heavy')
    p.add_argument('--repo',type=Path,default=here.parent/'math');args=p.parse_args()
    source=args.source;repo=args.repo;base=args.base
    need(A.digest(source/'manifest.json')==FROZEN,'frozen extension manifest')
    manifest=json.loads((source/'manifest.json').read_text())
    need(manifest['base_packet_manifest_sha256']==BASE==A.digest(base/'manifest.json'),'frozen base identity')
    for folder in (source,base):
        mm=json.loads((folder/'manifest.json').read_text())
        for name,sha in mm['files'].items():need(A.digest(folder/name)==sha,'frozen artifact '+name)
        for name,sha in mm['dependencies'].items():
            legacy='05-knowledge/results/charge_palettes_boundary_20261020g/zero_heavy.md'
            target=here/'frozen_base_zero_heavy.md' if folder==base and name==legacy else repo/name
            need(A.digest(target)==sha,'dependency '+name)
    need(manifest['incoming_deficit_k_commit']=='e265436678941418bc34a9eb3b697acd989524d0','incoming correction lineage')
    F,I,P,H=numbers(32);incidence=[];algebra_checks=0
    for k in range(2,8):
        all_h=all_K=strong_h=strong_K=0
        for edges in A.pairings(tuple(range(2*k))):
            h=sum(b-a==2 for a,b in edges);K=comb(sum(b-a==1 for a,b in edges),2)
            all_h+=h;all_K+=K
            if A.cut_test(edges,2*k):strong_h+=h;strong_K+=K
        need(all_h==(2*k-2)*F[k-1] and all_K==(k-1)*F[k-1],'unrestricted marked-edge incidences')
        need(strong_h==H[k] and strong_K==(comb(2*k-4,2)*I[k-2] if k>=3 else 0),'strong marked-edge incidences')
        incidence.append(dict(k=k,all_H=all_h,strong_H=strong_h,all_K=all_K,strong_K=strong_K))
    for k in range(4,33):
        for strong in (False,True):
            m=2*k-1 if strong else 2*k+1;E=I[k] if strong else F[k]
            J=(2*k-3)*I[k-1] if strong else F[k];p3=P[k] if strong else (k-1)*F[k]
            K=comb(2*k-4,2)*I[k-2] if strong else (k-1)*F[k-1]
            counts=(E*A.coefficient(2*k,m,k-1), (m*E-J)*A.coefficient(k,m-1,k-1),
                    p3*A.coefficient(k-1,m-1,k-1),E*comb(m,2)-(m-1)*J+K)
            need(sum(counts)==A.formula(k,2*k,strong,F,I,P)+K,'presentation cancellation without quotient')
            algebra_checks+=1
    local=local_audit(source)
    corridor=[]
    for k,d,strong in [(3,4,s) for s in (False,True)]+[(3,5,s) for s in (False,True)]+[(4,6,s) for s in (False,True)]+[(4,7,True)]:
        corridor.append(check_layer(k,d,strong,'ABC',F,I,P,H))
    source_corridor=json.loads((source/'extension-check.json').read_text())['layers']
    for ref in source_corridor:
        row=next(r for r in corridor if (r['k'],r['delta'],r['strong'])==(ref['k'],ref['delta'],ref['strong']))
        need(row['classes']==ref['classes'] and row['accepted']==ref['types'],'literal extended corridor receipt')
        need(sum(row['rejected'].values())==ref['clone_rejections'],'literal corridor clone count')
    owner=[]
    for k in (2,3,4):
        for strong in (False,True):owner.append(check_layer(k,2*k,strong,'ABCD' if k<=3 else 'CD',F,I,P,H))
    source_owner=json.loads((source/'first-owner-check.json').read_text())
    for ref in source_owner['layers']:
        row=next(r for r in owner if (r['k'],r['strong'])==(ref['k'],ref['strong']))
        for field in ('accepted','rejected','classes','raw_records','fibre_histogram'):need(row[field]==ref[field],'literal first-owner '+field)
        need(row['formal_corridor_expression']==ref['corridor_extrapolation'] and row['owner_correction']==ref['correction'],'literal owner arithmetic')
    for ref in json.loads((source/'owner-four-check.json').read_text())['rows']:
        row=next(r for r in owner if r['k']==4 and r['strong']==ref['strong'])
        need(row['accepted']==ref['accepted'] and row['classes']==ref['classes'],'k4 full CD universe')
        need(row['fibre_histogram']==ref['fibres'] and sum(row['rejected'].values())==ref['clone_rejections'],'k4 all owner fibres')
    for ref in source_owner['matching_incidence']:
        need(ref==next(r for r in incidence if r['k']==ref['k']),'matching incidence receipt')
    old=json.loads((repo/'05-knowledge/results/charge_palettes_boundary_20261020g/zero-heavy-mass.json').read_text())
    old2=next(r for r in old['complete_budget_compiler'] if r['k']==2)
    need(old2['all_zero_heavy_by_mass']['10']==owner[0]['classes']==114,'complete k2 compiler general')
    need(old2['strong_zero_heavy_by_mass']['6']==owner[1]['classes']==11,'complete k2 compiler strong')
    need(A.digest(source/'manifest.json')==FROZEN and A.digest(base/'manifest.json')==BASE,'frozen manifests unchanged')
    for name,sha in manifest['files'].items():need(A.digest(source/name)==sha,'extension unchanged '+name)
    print(json.dumps(dict(status='PASS independent frozen extension audit; mathematical scopes reviewed separately',
       extension_manifest_sha256=FROZEN,base_manifest_sha256=BASE,incoming_commit=manifest['incoming_deficit_k_commit'],
       source_files_authenticated=len(manifest['files']),dependency_files_authenticated=len(manifest['dependencies']),
       frozen_base_source_capture={'path':'05-knowledge/results/charge_palettes_boundary_20261020g/zero_heavy.md',
       'commit':'c4a7875652a5ea516e93a31c3c6410d2204773b2',
       'sha256':A.digest(here/'frozen_base_zero_heavy.md')},
       matching_incidence=incidence,presentation_identity_checks=algebra_checks,local_ports=local,
       extended_corridor=corridor,owner_layers=owner,
       retained_records_checked=sum(r['raw_records'] for r in corridor+owner)),indent=2,sort_keys=True))

if __name__=='__main__':main()
