#!/usr/bin/env python3
"""Complete closed endpoint-four orders8 and9 finite certificate.

Uses the sealed parent Gaussian engine/old-bank words and local h-ledgers.
All120 source cells are checked, including entire illegal order8 leaves.
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

def endpair_words(n,h,chord):
    a,b,x,y,c,z=0,1,2,n-3,n-2,n-1;D=list(range(3,n-3));U=list(range(4,n-3))
    def forward(k):
        L,R=list(range(3,k)),list(range(k,n))
        return [[a,b]+L+R+[x],[b]+L+[a]+R+[x],[b,x]+L+[a]+R]
    def backward(k):
        V=list(range(k,n));tail=list(range(3,k))
        return [[a]+V+[b,x]+tail,[x,a]+V+[b]+tail,[x,a,b]+tail+V]
    if h==2:
        if chord:return [('h2_xc',[],[[a,b,x]+D+[z,y,c],[a,x]+D+[z,y,b,c],[a,x]+D+[y,b,c,z]])]
        return [('h2_cx_3c',[(3,c,1)],[[z,a,b,x]+D+[y,c],[z,a,x]+D+[b,y,c],[a,x]+D+[b,y,c,z]]),
                ('h2_cx_c3',[(3,c,0)],[[b,x,a,3]+U+[y,c,z],[b,x,a,3,c,z]+U+[y],[a,b,x,3,c,z]+U+[y]])]
    if h==3:
        W=list(range(4,n-3));R=[y,c,z]
        return [('h3_forward',[(1,3,1)],forward(y)),
                ('h3_backward_4b',[(1,3,0),(1,4,0)],backward(4)),
                ('h3_backward_b4',[(1,3,0),(1,4,1)],[[a,b]+W+R+[x,3],[b]+W+[a]+R+[x,3],[b,x,3]+W+[a]+R])]
    return [('h%d_forward'%h,[(1,3,1)],forward(h)),('h%d_backward'%h,[(1,3,0)],backward(5 if h==4 else 4))]


def inspect_cell(n,raw,eq,words,label):
    m=n*(n-1)//2;basis=gf.eliminate(eq);need(basis is not None,('empty cell',label))
    src=[[gf.reduce_form(f,basis) for f in row] for row in raw];at=src;rawat=raw;blue=[];models=set();xor=set();illegal=[]
    for step,q in enumerate(words,1):
        need(sorted(q)==list(range(n)),('not HP',label,step))
        illegal += [{'step':step,'join':[u,v],'forced_arc_value':at[u][v] if at[u][v]<2 else None} for u,v in zip(q,q[1:]) if at[u][v]!=1]
        full=eq+gf.blue_equations(rawat,q);trial=gf.eliminate(full)
        if trial is not None:
            need(m-len(trial)==0,('non-singleton blue locus',label,step))
            free,solutions=p.enumerate_models(raw,full);need(len(solutions)==1,'singleton must mean one model')
            blue.append({'step':step,'dimension':len(free),'source_rows':list(solutions[0])});models.update(solutions)
        at=gf.k_move(at,q);rawat=gf.k_move(rawat,q);xor^=p.edge_set(q)
    need(at==gf.k_move(src,list(range(n))) and xor==p.edge_set(list(range(n))),('old target',label))
    if illegal:
        need(n==8 and label=='h4_backward' and m-len(basis)==3,'unexpected illegal leaf')
        need(all(r['join']==[0,5] and r['forced_arc_value']==0 for r in illegal),'illegal boundary not actual false port')
        free,solutions=p.enumerate_models(raw,eq);need(len(solutions)==8,'complete illegal leaf size')
        models.update(solutions)
    return {'label':label,'source_dimension':m-len(basis),'old_words':words,'illegal_old_joins':illegal,'feasible_blue_loci':blue},models

def main():
    data=json.loads((HERE/'found.json').read_text());need(len(data)==21,'twenty-one data sources')
    indexed={tuple(c['source_rows']):c for c in data};need(len(indexed)==21,'twenty-one distinct sources')
    allmodels=set();partitions=[];failures=[];leaves=[];counts={}
    expected={8:{'all_source_cells':40,'base_cells':34,'base_passed_cells':34,'endpair_cells':6,'endpair_leaves':13,'endpair_passed_leaves':11,'endpair_failed_leaves':2},
              9:{'all_source_cells':80,'base_cells':72,'base_passed_cells':67,'base_failed_cells':5,'endpair_cells':8,'endpair_leaves':17,'endpair_passed_leaves':17}}
    for n in (8,9):
        m=n*(n-1)//2;totals=Counter()
        for I in combinations(range(2,n-1),3):
            for opened,chord in ((True,False),(True,True),(False,False),(False,True)):
                raw,eq=p.source_equations(n,I,opened,chord);basis=gf.eliminate(eq);need(basis is not None,'empty complete sourcecell');totals['all_source_cells']+=1
                if opened or not {n-3,n-2}<=set(I):
                    name,qs=bank.words(n,set(I),opened,chord);r,models=inspect_cell(n,raw,eq,qs,name)
                    totals['base_cells']+=1;totals['base_failed_cells' if models else 'base_passed_cells']+=1
                    if models:failures.append(dict(r,n=n,ports=list(I),middle_open=opened,x_to_c=chord))
                    allmodels.update(models)
                else:
                    totals['endpair_cells']+=1;h=next(iter(set(I)-{n-3,n-2}));sizes=[];guards_all=[]
                    for name,guards,qs in endpair_words(n,h,chord):
                        extra=[raw[u][v]^value for u,v,value in guards];g=gf.eliminate(eq+extra)
                        need(g is not None,('vacuous endpair leaf',name));sizes.append(1<<(m-len(g)));guards_all.append(extra)
                        r,models=inspect_cell(n,raw,eq+extra,qs,name)
                        totals['endpair_leaves']+=1;totals['endpair_failed_leaves' if models else 'endpair_passed_leaves']+=1
                        rec=dict(r,n=n,ports=list(I),middle_open=opened,x_to_c=chord,guards=guards);leaves.append(rec)
                        if models:failures.append(rec)
                        allmodels.update(models)
                    need(sum(sizes)==1<<(m-len(basis)),('incomplete endpair partition',n,I,chord))
                    need(all(gf.eliminate(eq+g+h) is None for i,g in enumerate(guards_all) for h in guards_all[i+1:]),('overlapping endpair leaves',n,I,chord))
                    partitions.append({'n':n,'ports':list(I),'x_to_c':chord,'source_dimension':m-len(basis),'leaf_sizes':sizes,'disjoint_complete':True})
        need(dict(totals)==expected[n],('counts',n,dict(totals)));counts[n]=dict(totals)
    need(allmodels==set(indexed),'exact exceptional union differs from21rows')
    need(Counter(map(len,allmodels))=={8:16,9:5},'exact exceptional order counts')
    certificates=[]
    for A in sorted(allmodels,key=lambda r:(len(r),r)):
        n=len(A);rec=indexed[A];I=rec['ports'];P=list(range(n));B=p.move(A,P);at=A;states=[p.state_record(A)];steps=[];xor=set();allpairs=set(combinations(range(n),2))
        need(p.directed(A,P) and p.witness(A,P) is None and p.arc(A,n-1,0),'source not closed blue')
        need(A[0].bit_count()==n-5 and A[n-1].bit_count()==4,'endpoint degrees')
        need(set(p.neighbourhoods(A,0)['in'])==set(I)|{n-1},'incoming port set')
        need(set(p.neighbourhoods(A,n-1)['out'])=={0}|{n-1-v for v in I},'outgoing reflected port set')
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
        certificates.append({'n':n,'source_rows':list(A),'ports':I,'steps':steps,'states':states,'literal_target_rows':list(B)})
    # Full literal order-seven obstruction; no source or HP filter.
    from itertools import permutations
    H=(6,4,8,19,103,79,15);q=tuple(range(7));T=p.move(H,q)
    need(p.directed(H,q) and p.witness(H,q) is None and p.arc(H,6,0) and H[0].bit_count()==2 and H[-1].bit_count()==4,'hostile source')
    need(T==(122,124,121,112,32,64,16),'hostile literal target')
    hp=[w for w in permutations(range(7)) if p.directed(T,w)]
    need(len(hp)==9 and all(p.witness(T,w) is None for w in hp),'target must have no black HP')
    print(json.dumps({'status':'PASS','orders':[8,9],'scope':'complete closed marked endpoint-four replacement at8 and9, together with independent higher-order packets gives sharp threshold n>=8; full black-core remains open',
                      'parent_check_sha256':EXPECTED_PARENT,'counts':counts,'endpair_partitions':partitions,'endpair_leaves':leaves,'failed_leaves':failures,
                      'distinct_exceptional_sources':21,'certificates':certificates,
                      'sharpness_hostile':{'n':7,'source_rows':H,'target_rows':T,'target_Hamiltonian_words':hp,'black_target_words':0}},indent=2,sort_keys=True))

if __name__=='__main__':main()
