#!/usr/bin/env python3
"""All labelled tournaments through order five, and structural hostile controls."""
import argparse
from collections import Counter
from itertools import combinations, permutations
import json
from pathlib import Path
import random
from engine import *


def canonical_small(a):
    return min(tuple(a[i][j] for i in p for j in p) for p in permutations(range(len(a))))


def exhaustive():
    rows=[]
    class_banks={}
    for n in range(1,6):
        count=1 << (n*(n-1)//2)
        source_aut_hist=Counter(); target_aut_hist=Counter()
        h_total=0; arc_checks=0; nonextreme=0; roots_nonunique=0; classes={}
        for mask in range(count):
            t=tournament(n,mask);r=regularize(t);z=2*n
            validate(t);validate(r)
            need(r == regularize_arcs(t),'independent construction mismatch')
            need(scores(r) == (n,)*(2*n+1),'regularity')
            need(induced(r,range(n)) == t and induced(r,range(n,2*n)) == t,'two source copies')
            signatures=[tuple(r[n+i][j] for j in range(n)) for i in range(n)]
            need(len(set(signatures)) == n,'open-in-neighborhood signatures collide')
            deg=scores(t)
            for i in range(n):
                need(common_out(r,i,n+i) == 0,'matched pair has common out-neighbor')
                need(common_out(r,z,i) == deg[i],'z U common count')
                need(common_out(r,z,n+i) == n-1-deg[i],'z V common count')
            no_extreme=all(0<d<n-1 for d in deg)
            need(universal_common_out_vertices(r) == ([z] if no_extreme else []),'canonical-root predicate')
            nonextreme+=no_extreme
            sa=brute_automorphisms(t)
            da=isomorphisms(r,r)
            root_aut=[p for p in da if p[z]==z]
            expected={diagonal_lift(p) for p in sa}
            need(set(root_aut)==expected,'root-stabilizer faithfulness')
            if no_extreme:
                need(set(da)==expected,'unmarked automorphism faithfulness on nonextreme domain')
            source_aut_hist[len(sa)]+=1;target_aut_hist[len(da)]+=1
            roots_nonunique+=len(da)>len(sa)
            # Converse-equivariance uses a fixed swap of U and V.
            swap=tuple(range(n,2*n))+tuple(range(n))+(z,)
            need(is_map(regularize(converse(t)),converse(r),swap),'converse intertwiner')
            used=set()
            for i,j in combinations(range(n),2):
                support=block(n,i,j)
                need(len(support)==4 and not (support & used),'blocks not edge-disjoint')
                used |= support
                a,b=(i,j) if t[i][j] else (j,i)
                need(directed_cycle(r,(a,b,n+a,n+b)),'source edge not a directed four-cycle')
                changed=reverse_edges(t,{frozenset((i,j))})
                need(regularize(changed)==reverse_edges(r,support),'single-edge block intertwiner')
                need(scores(reverse_edges(r,support))==scores(r),'cycle reversal changed degree')
                arc_checks+=1
            need(len(used)==4*(n*(n-1)//2),'block support count')
            paths=hamilton_paths(t)
            need(len(paths)==hamilton_path_count_dp(t),'permutation/DP path count disagreement')
            # Independent distinguished-cycle scan, without using path list.
            distinguished=[(z,)+tuple(v for i in q for v in (i,n+i))
                           for q in permutations(range(n))]
            distinguished=[c for c in distinguished if directed_cycle(r,c)]
            need(set(distinguished)=={lifted_cycle(p) for p in paths},'path-cycle bijection')
            need(len(paths)%len(sa)==0,'nonintegral ordinary path orbit mass')
            for path in paths:
                cycle=lifted_cycle(path)
                need(source_from_cycle(cycle)==path,'path recovery')
                # The lifted action on actual paths has trivial stabilizers.
                need(sum(tuple(p[i] for i in path)==path for p in sa)==1,'path action not free')
                fixed={frozenset(path[k:k+2]) for k in range(n-1)}
                support=set().union(*(block(n,*tuple(e)) for e in unordered_edges(n)-fixed)) if n>=3 else set()
                cpartner=path_complement(t,path)
                liftpartner=reverse_edges(r,support)
                need(liftpartner==regularize(cpartner),'native source C block intertwiner')
                need(directed_cycle(liftpartner,cycle),'source C failed to retain cycle')
                need(len(support)==2*(n-1)*(n-2),'lifted source C edge count')
                full=cycle_complement(r,cycle)
                need((full==liftpartner)==(n==1),'source C confused with full target noncycle C')
                need(len(changed_edges(full,liftpartner))==5*(n-1),'complement loss ledger size')
            h_total+=len(paths)
            key=canonical_small(t)
            classes.setdefault(key,t)
        representatives=list(classes.values())
        # Pairwise rooted faithfulness and target isomorphism classes are checked
        # with the generic matrix search, not by comparing source canonical keys.
        target_reps=[]; collisions=[]; nonextreme_reps=[]; merged_source_reps=[]; merged_target_reps=[]
        for i,t in enumerate(representatives):
            r=regularize(t)
            if all(0<d<n-1 for d in scores(t)):nonextreme_reps.append(t)
            if not any(isomorphic(t,v) or isomorphic(t,converse(v)) for v in merged_source_reps):
                merged_source_reps.append(t)
            if not any(isomorphic(r,v) or isomorphic(r,converse(v)) for v in merged_target_reps):
                merged_target_reps.append(r)
            matching=[k for k,v in enumerate(target_reps) if isomorphic(r,v)]
            if matching:collisions.append([i,matching[0]])
            else:target_reps.append(r)
            for j,u in enumerate(representatives):
                rooted=bool(isomorphisms(r,regularize(u),fixed={2*n:2*n},stop_after=1))
                need(rooted==(i==j),'distinct source classes have rooted-isomorphic doubles')
                rooted_dual=bool(isomorphisms(r,converse(regularize(u)),fixed={2*n:2*n},stop_after=1))
                need(rooted_dual==isomorphic(t,converse(u)),'rooted merged-class faithfulness')
        for i,t in enumerate(nonextreme_reps):
            for j,u in enumerate(nonextreme_reps):
                need(isomorphic(regularize(t),regularize(u))==(i==j),'unmarked restricted faithfulness')
        self_converse=sum(isomorphic(t,converse(t)) for t in representatives)
        nonextreme_self_converse=sum(isomorphic(t,converse(t)) for t in nonextreme_reps)
        nonextreme_merged=[]
        for t in nonextreme_reps:
            if not any(isomorphic(t,v) or isomorphic(t,converse(v)) for v in nonextreme_merged):
                nonextreme_merged.append(t)
        need(h_total==__import__('math').factorial(n)*(1 << ((n-1)*(n-2)//2)),'total labelled path incidence')
        rows.append({'n':n,'labelled_tournaments':count,'source_ordinary_classes':len(classes),
                     'source_self_converse_ordinary_classes':self_converse,
                     'source_nonextreme_ordinary_classes':len(nonextreme_reps),
                     'source_nonextreme_self_converse_classes':nonextreme_self_converse,
                     'source_nonextreme_merged_classes':len(nonextreme_merged),
                     'target_unmarked_classes':len(target_reps),'source_converse_merged_classes':len(merged_source_reps),
                     'target_unmarked_converse_merged_classes':len(merged_target_reps),
                     'labelled_nonextreme_tournaments':nonextreme,
                     'labelled_hamilton_paths':h_total,'single_arc_block_checks':arc_checks,
                     'source_automorphism_order_histogram':dict(sorted(source_aut_hist.items())),
                     'target_automorphism_order_histogram':dict(sorted(target_aut_hist.items())),
                     'targets_with_nonfixed_original_root':roots_nonunique,
                     'target_classes_lost_by_forgetting_root':len(classes)-len(target_reps)})
        class_banks[n]=representatives
    return rows,class_banks


def collision_witness(q):
    one=tournament(1,0)
    a=order_join(one,q);b=order_join(q,one)
    da=regularize(a);db=regularize(b);p=recenter_source_map(q)
    need(is_map(da,db,p),'source-to-sink recenter full-matrix identity')
    return {'source_order':len(a),'target_order':len(da),'Q':matrix_lines(q),
            'source_1_plus_Q':matrix_lines(a),'source_Q_plus_1':matrix_lines(b),
            'D_1_plus_Q':matrix_lines(da),'D_Q_plus_1':matrix_lines(db),
            'explicit_D_isomorphism_old_to_new':list(p),
            'source_scc_sizes':list(scc_sizes_in_order(a)),
            'other_source_scc_sizes':list(scc_sizes_in_order(b)),
            'source_isomorphic':isomorphic(a,b),
            'source_isomorphic_to_other_converse':isomorphic(a,converse(b)),
            'source_aut_order':len(isomorphisms(a,a)),
            'target_aut_order':len(isomorphisms(da,da)),
            'source_hamilton_paths':hamilton_path_count_dp(a)}


def recenter_controls():
    total=0
    for m in range(0,6):
        for mask in range(1 << (m*(m-1)//2)):
            q=tournament(m,mask);one=tournament(1,0)
            a=regularize(order_join(one,q));b=regularize(order_join(q,one))
            need(is_map(a,b,recenter_source_map(q)),'all-small recenter identity')
            total+=1
    c3=((0,1,0),(0,0,1),(1,0,0))
    # 0 -> 1 -> 2 -> 3 -> 0 plus 0 -> 2 and 1 -> 3.
    s4=((0,1,1,0),(0,0,1,1),(0,0,0,1),(1,0,0,0))
    need(scc_sizes_in_order(s4)==(4,),'S4 witness not strongly connected')
    ordinary=collision_witness(c3)
    merged_small=collision_witness(order_join(one,c3))
    need(not merged_small['source_isomorphic'] and not merged_small['source_isomorphic_to_other_converse'],'smaller merged witness failed')
    need(merged_small['source_scc_sizes']==[1,1,3] and merged_small['other_source_scc_sizes']==[1,3,1],'smaller SCC types')
    merged=collision_witness(order_join(c3,s4))
    need(not ordinary['source_isomorphic'] and ordinary['source_isomorphic_to_other_converse'],'ordinary witness wrong scope')
    need(not merged['source_isomorphic'] and not merged['source_isomorphic_to_other_converse'],'merged witness failed')
    need(merged['source_scc_sizes']==[1,3,4] and merged['other_source_scc_sizes']==[3,4,1],'SCC obstruction types')
    return total,ordinary,merged_small,merged


def hostile_controls():
    rejected=[]
    t=((0,1),(0,0));r=regularize(t);n=2
    # Swapping the second copy to the converse is NOT this regularization.
    bad=[list(row) for row in r]
    bad[n][n+1],bad[n+1][n]=bad[n+1][n],bad[n][n+1]
    need(scores(bad)!=(n,)*(2*n+1),'converse-copy mutant was not detected')
    rejected.append('converse-copy substitution loses regularity on transitive order two')
    bad=[list(row) for row in r]
    for i,j in combinations(range(n),2):
        for x,y in ((i,n+j),(j,n+i)):
            bad[x][y],bad[y][x]=bad[y][x],bad[x][y]
    need(scores(bad)!=(n,)*(2*n+1),'same-direction-cross mutant was not detected')
    rejected.append('same-direction cross rule loses regularity')
    bad=reverse_edges(r,{frozenset((i,n+i)) for i in range(n)})
    need(scores(bad)!=(n,)*(2*n+1),'diagonal reversal mutant was not detected')
    rejected.append('reversing only the diagonal matching loses regularity')
    need(not directed_cycle(r,(4,0,2,1,3)),'forward HP order mutant was not detected')
    rejected.append('forward rather than reversed source path is not a directed target cycle')
    need(scores(reverse_edges(r,{frozenset((0,1))}))!=scores(r),'single target edge mutant was not detected')
    rejected.append('single target arc reversal is not the degree-preserving source lift')
    need(cycle_complement(r,lifted_cycle((0,1)))!=regularize(path_complement(t,(0,1))),'full-C mutant was not detected')
    rejected.append('full target noncycle complement differs already at source order two')
    need(universal_common_out_vertices(r)==[],'unrestricted canonical root mutant was not detected')
    rejected.append('root common-neighbor criterion requires no source and no sink')
    need(len(isomorphisms(r,r))==5 and len(brute_automorphisms(t))==1,'unrestricted full-Aut mutant was not detected')
    rejected.append('unrestricted full automorphism equality fails at source order two')
    # All target Hamilton cycles versus the distinguished alternating subset.
    hc=[(4,)+p for p in permutations(range(4)) if directed_cycle(r,(4,)+p)]
    need(len(hc)>len(hamilton_paths(t)),'all-Hamilton-cycles mutant was not detected')
    rejected.append('distinguished cycle count is not all target Hamilton cycle count')
    return rejected,{'D_transitive_2_all_cycles_rooted_at_z':len(hc),'distinguished_cycles':1}


def larger_controls():
    rng=random.Random(520043)
    count=0
    for n in range(6,11):
        for _ in range(12):
            t=tournament(n,rng.getrandbits(n*(n-1)//2));r=regularize(t)
            need(scores(r)==(n,)*(2*n+1),'sample regularity')
            need(r==regularize_arcs(t),'sample independent construction')
            for i,j in combinations(range(n),2):
                partner=reverse_edges(t,{frozenset((i,j))})
                need(regularize(partner)==reverse_edges(r,block(n,i,j)),'sample block identity')
            count+=1
    return {'deterministic_samples':count,'source_orders':[6,7,8,9,10],'seed':520043,
            'claims':'sample controls only; not an exhaustive census'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output');args=parser.parse_args()
    rows,_=exhaustive()
    A={0:1,**{x['n']:x['source_ordinary_classes'] for x in rows}}
    S={0:1,**{x['n']:x['source_self_converse_ordinary_classes'] for x in rows}}
    lower_bounds=[]
    for x in rows:
        n=x['n']
        if n<2:continue
        N=A[n]-2*A[n-1]+A[n-2]
        J=S[n]-S[n-2]
        need(x['source_nonextreme_ordinary_classes']==N,'nonextreme ordinary inclusion-exclusion')
        need(x['source_nonextreme_self_converse_classes']==J,'nonextreme self-converse deletion')
        need(2*x['source_nonextreme_merged_classes']==N+J,'nonextreme converse orbit count')
        lower_bounds.append({'source_order':n,'target_regular_order':2*n+1,'ordinary_lower_bound':N,
                             'self_converse_lower_bound':J,'merged_lower_bound':(N+J)//2})
    need(all(x['source_ordinary_classes']==x['target_unmarked_classes'] for x in rows if x['n']<4),'ordinary threshold lower-order control')
    need(all(x['source_converse_merged_classes']==x['target_unmarked_converse_merged_classes'] for x in rows if x['n']<5),'merged threshold lower-order control')
    total,ordinary,merged_small,merged=recenter_controls()
    mutations,hc=hostile_controls()
    result={'status':'PASS','universe':'all labelled tournaments of source orders 1..5',
            'total_labelled_tournaments':sum(r['labelled_tournaments'] for r in rows),
            'first_ordinary_collision_source_order':next(x['n'] for x in rows if x['source_ordinary_classes']>x['target_unmarked_classes']),
            'first_merged_collision_source_order':next(x['n'] for x in rows if x['source_converse_merged_classes']>x['target_unmarked_converse_merged_classes']),
            'exhaustive_rows':rows,'no_extreme_class_count_and_regular_lower_bound_controls':lower_bounds,'all_small_recenter_checks_Q_orders_0_through_5':total,
            'ordinary_collision_order_4':ordinary,'merged_collision_order_5':merged_small,'merged_collision_order_8':merged,
            'mutation_controls_rejected':mutations,'all_vs_distinguished_cycle_control':hc,
            'larger_controls':larger_controls(),
            'scope_limits':['The Type-I regularization is known from arXiv:2605.04846v1 Lemma 4.3.',
                            'No first-novelty claim for the structural extractions.',
                            'Source order five is the first merged collision, controlled by exhaustive smaller class pairs.',
                            'No native full-target complement or port-clock equivalence.',
                            'No new general formula for A000568.']}
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(payload)
    print(payload,end='')

if __name__=='__main__':main()
