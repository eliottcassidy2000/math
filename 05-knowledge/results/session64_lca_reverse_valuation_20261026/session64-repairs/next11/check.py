#!/usr/bin/env python3
"""Complete closed endpoint-four order11 finite certificate.

Uses the sealed parent Gaussian engine/old-bank words and local h-ledgers.
All224 source cells and all25 disjoint endpair subcells are checked.
"""
from pathlib import Path
from itertools import combinations
from hashlib import sha256
from collections import Counter
import importlib.util,json,sys
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'check.py'
EXPECTED_PARENT='2004c55c07a00b60eda8db0a91d5f31d6f09ee43adcc7ed7f96d2ee120011641'
if sha256(PARENT.read_bytes()).hexdigest()!=EXPECTED_PARENT:raise RuntimeError('sealed parent changed')
spec=importlib.util.spec_from_file_location('sealed_parent',PARENT);p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
gf,bank,need=p.gf,p.bank,p.need

def endpair_words(h,chord):
    n=11;a,b,x,y,c,z=0,1,2,8,9,10;D=list(range(3,8));U=list(range(4,8))
    def forward(k):
        L,R=list(range(3,k)),list(range(k,11))
        return [[a,b]+L+R+[x],[b]+L+[a]+R+[x],[b,x]+L+[a]+R]
    def backward(k):
        V=list(range(k,11));tail=list(range(3,k))
        return [[a]+V+[b,x]+tail,[x,a]+V+[b]+tail,[x,a,b]+tail+V]
    if h==2:
        if chord:return [('h2_xc',[],[[a,b,x]+D+[z,y,c],[a,x]+D+[z,y,b,c],[a,x]+D+[y,b,c,z]])]
        return [('h2_cx_3c',[(3,c,1)],[[z,a,b,x]+D+[y,c],[z,a,x]+D+[b,y,c],[a,x]+D+[b,y,c,z]]),
                ('h2_cx_c3',[(3,c,0)],[[b,x,a,3]+U+[y,c,z],[b,x,a,3,c,z]+U+[y],[a,b,x,3,c,z]+U+[y]])]
    if h==3:
        W=list(range(4,8));R=[y,c,z]
        return [('h3_forward',[(1,3,1)],forward(y)),
                ('h3_backward_4b',[(1,3,0),(1,4,0)],backward(4)),
                ('h3_backward_b4',[(1,3,0),(1,4,1)],[[a,b]+W+R+[x,3],[b]+W+[a]+R+[x,3],[b,x,3]+W+[a]+R])]
    return [('h%d_forward'%h,[(1,3,1)],forward(h)),('h%d_backward'%h,[(1,3,0)],backward(5 if h==4 else 4))]

def inspect_cell(raw,eq,words,label):
    basis=gf.eliminate(eq);need(basis is not None,('empty cell',label))
    src=[[gf.reduce_form(f,basis) for f in row] for row in raw];at=src;rawat=raw;blue=[];models=set();xor=set()
    for step,q in enumerate(words,1):
        need(sorted(q)==list(range(11)),('not HP',label,step))
        need(all(at[u][v]==1 for u,v in zip(q,q[1:])),('illegal old join',label,step))
        full=eq+gf.blue_equations(rawat,q);trial=gf.eliminate(full)
        if trial is not None:
            need(55-len(trial)==0,('non-singleton blue locus',label,step))
            free,solutions=p.enumerate_models(raw,full);need(len(solutions)==1,'singleton must mean one model')
            blue.append({'step':step,'dimension':len(free),'source_rows':list(solutions[0])});models.update(solutions)
        at=gf.k_move(at,q);rawat=gf.k_move(rawat,q);xor^=p.edge_set(q)
    need(at==gf.k_move(src,list(range(11))) and xor==p.edge_set(list(range(11))),('old target',label))
    return {'label':label,'source_dimension':55-len(basis),'old_words':words,'feasible_blue_loci':blue},models

def main():
    data=json.loads((HERE/'found.json').read_text());need(len(data)==7,'seven data sources')
    indexed={tuple(c['source_rows']):c for c in data};need(len(indexed)==7,'seven distinct sources')
    totals=Counter();failures=[];allmodels=set();partitions=[];leaves=[]
    for I in combinations(range(2,10),3):
        for opened,chord in ((True,False),(True,True),(False,False),(False,True)):
            raw,eq=p.source_equations(11,I,opened,chord);basis=gf.eliminate(eq)
            need(basis is not None,'empty complete sourcecell');totals['all_source_cells']+=1
            if opened or not {8,9}<=set(I):
                name,qs=bank.words(11,set(I),opened,chord);r,models=inspect_cell(raw,eq,qs,name)
                totals['base_cells']+=1;totals['base_failed_cells' if models else 'base_passed_cells']+=1
                if models:failures.append(dict(r,ports=list(I),middle_open=opened,x_to_c=chord))
                allmodels.update(models)
            else:
                totals['endpair_cells']+=1;h=next(iter(set(I)-{8,9}));cells=[];sizes=[];guards_all=[]
                for name,guards,qs in endpair_words(h,chord):
                    extra=[raw[u][v]^value for u,v,value in guards];g=gf.eliminate(eq+extra)
                    need(g is not None,('vacuous endpair leaf',name));sizes.append(1<<(55-len(g)));guards_all.append(extra)
                    r,models=inspect_cell(raw,eq+extra,qs,name)
                    totals['endpair_leaves']+=1;totals['endpair_failed_leaves' if models else 'endpair_passed_leaves']+=1
                    rec=dict(r,ports=list(I),middle_open=opened,x_to_c=chord,guards=guards)
                    leaves.append(rec)
                    if models:failures.append(rec)
                    allmodels.update(models)
                need(sum(sizes)==1<<(55-len(basis)),('incomplete endpair partition',I,chord))
                need(all(gf.eliminate(eq+g+h) is None for i,g in enumerate(guards_all) for h in guards_all[i+1:]),('overlapping endpair leaves',I,chord))
                partitions.append({'ports':list(I),'x_to_c':chord,'source_dimension':55-len(basis),'leaf_sizes':sizes,'disjoint_complete':True})
    need(dict(totals)=={'all_source_cells':224,'base_cells':212,'base_passed_cells':208,'endpair_cells':12,'endpair_leaves':25,'endpair_passed_leaves':23,'base_failed_cells':4,'endpair_failed_leaves':2},('counts',dict(totals)))
    need(allmodels==set(indexed),'exact failed-locus union differs from seven rows')
    need(sum(len(r['feasible_blue_loci']) for r in failures)==7,'seven singleton blue loci')
    certificates=[];allpairs=set(combinations(range(11),2))
    for A in sorted(allmodels):
        rec=indexed[A];I=rec['ports'];P=list(range(11));B=p.move(A,P);at=A;states=[p.state_record(A)];steps=[];xor=set()
        need(p.directed(A,P) and p.witness(A,P) is None and p.arc(A,10,0),'source not closed blue')
        need(A[0].bit_count()==6 and A[10].bit_count()==4,'endpoint degrees')
        need(set(p.neighbourhoods(A,0)['in'])==set(I)|{10},'incoming port set')
        need(set(p.neighbourhoods(A,10)['out'])=={0}|{10-v for v in I},'outgoing reflected port set')
        for q in rec['words']:
            need(sorted(q)==P and p.directed(at,q),'illegal repair');black=p.witness(at,q);need(black is not None,'blue repair')
            nxt=p.move(at,q);before=p.clone_data(at);after=p.clone_data(nxt)
            oldtw={tuple(e) for e in before['twin_pairs']};newtw={tuple(e) for e in after['twin_pairs']}
            steps.append({'word':q,'blackness_witness':black,'word_endpoint_ports':[p.neighbourhoods(at,q[0]),p.neighbourhoods(at,q[-1])],
                          'flipped_unordered_pairs':[list(e) for e in sorted(allpairs-p.edge_set(q))],
                          'clone_transport':{'new_twin_pairs':[list(e) for e in sorted(newtw-oldtw)],'lost_twin_pairs':[list(e) for e in sorted(oldtw-newtw)],
                                             'literal_block_intersections':[[sorted(set(u)&set(v)) for v in after['ordered_maximal_transitive_modules']] for u in before['ordered_maximal_transitive_modules']]}})
            states.append(p.state_record(nxt));xor^=p.edge_set(q);at=nxt
        need(xor==p.edge_set(P) and at==B,'literal repair target')
        certificates.append({'source_rows':list(A),'ports':I,'steps':steps,'states':states,'literal_target_rows':list(B)})
    # Hostile control: h7 two different blue stages are not one source.
    h7=[r for r in failures if r['label']=='h7_forward' and not r['x_to_c']]
    need(len(h7)==1 and len(h7[0]['feasible_blue_loci'])==2,'h7 double locus lost')
    need(h7[0]['feasible_blue_loci'][0]['source_rows']!=h7[0]['feasible_blue_loci'][1]['source_rows'],'h7 distinct sources accidentally identified')
    print(json.dumps({'status':'PASS','n':11,'scope':'complete closed marked endpoint-four replacement at order11; seven exceptional literal sources get three-black repairs; no open-path or full black-core claim',
                      'parent_check_sha256':EXPECTED_PARENT,'counts':dict(totals),'endpair_partitions':partitions,'endpair_leaves':leaves,'failed_leaves':failures,
                      'singleton_blue_loci':7,'distinct_exceptional_sources':7,'certificates':certificates},indent=2,sort_keys=True))

if __name__=='__main__':main()
