"""Independent forced-chain/path audit for the separate direct contraction."""
from collections import Counter
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path
import sys

from audit_reverse import check, partitions, raw_weight


def run(source):
    sys.path.insert(0,str(source))
    import reverse_ordinary as full_impl
    import contracted_ordinary as small_impl
    vertices=edges=forced=pathchecks=modular=0
    supportchecks=0
    samples=[]
    strict_prefix_example=None
    zero_nonforced=None
    orders=list(range(51))+[60,80,100]
    for n in orders:
        full=full_impl.build_dag(n)
        small=small_impl.build_contracted_dag(n)
        small_index={v:i for i,v in enumerate(small.vertices)}
        original_index={v:i for i,v in enumerate(full.vertices)}
        odd=full_impl.odd_factorials(n)
        # Define contraction from the actual graph, without canonical_vertex.
        # Branching availability, not merely k=0 or J=0, decides removal.
        removed={v for v,(ell,s,D) in enumerate(full.vertices)
                 if ell>0 and (s==0 or ell>s)}
        survivors=set(range(len(full.vertices)))-removed
        check({full.vertices[v] for v in survivors}==set(small.vertices),('exact survivor set',n))
        image=[]
        for v in range(len(full.vertices)):
            w=v
            while w in removed:
                row=full.edges[w]
                check(len(row)==1,('forced uniqueness',n,w))
                e=row[0]
                check(e.k==0 and e.shift==0 and full_impl.edge_factor(e,odd)==1,('unit chain edge',n,w))
                check(full.vertices[w][1:]==full.vertices[e.target][1:],('exact S/D chain preservation',n,w))
                w=e.target
            image.append(small_index[full.vertices[w]])
            check(full.potential[v]==small.potential[image[-1]],('all suffix potentials descend',n,v))
        for v in survivors:
            sv=image[v]
            expected=[(image[e.target],e.ell,e.k,e.u,e.shift) for e in full.edges[v]]
            actual=[(e.target,e.ell,e.k,e.u,e.shift) for e in small.edges[sv]]
            check(expected==actual,('each exact labelled local edge survives',n,v))
            check(full.prefix_cost[v]==small.prefix_cost[sv],('surviving prefix minimum',n,v))
        for v in removed:
            if strict_prefix_example is None and full.prefix_cost[v]>small.prefix_cost[image[v]]:
                strict_prefix_example={'n':n,'source':full.vertices[v],'canonical_target':small.vertices[image[v]],'source_P':full.prefix_cost[v],'target_P':small.prefix_cost[image[v]]}
        a,z,Y,Z,_=full_impl.exact_audit(full)
        b,w,Ys,Zs,_=full_impl.exact_audit(small)
        check((a,z)==(b,w),('exact count normalization preserved',n))
        for i,vertex in enumerate(small.vertices):
            j=original_index[vertex]
            check((Ys[i],Zs[i])==(Y[j],Z[j]),('exact survivor aggregates',n,vertex))
        for K in (1,2,3,8,17):
            residue,work,values=full_impl.normalized_modular(small,K)
            answer,_,_=full_impl.normalized_modular(full,K)
            check(residue==answer==z%(1<<K),('contracted normalized value',n,K))
            check(values==[x%(1<<K) for x in Zs],('all contracted modular coordinates',n,K))
            check(len(values)==len(small.vertices),('full contracted value-array allocation',n,K))
            modular+=1
            if n in (30,50,60,80,100) and K in (2,8):
                samples.append({'n':n,'K':K,'input_states':len(small.vertices)-1,'all_edges':sum(map(len,small.edges)),'modular_updates':work['edge_updates']})
        if n<=32:
            support={K:set() for K in (1,2,3,8)}
            full_rows=[{e.k:(j,e) for j,e in enumerate(row)} for row in full.edges]
            small_rows=[{e.k:(j,e) for j,e in enumerate(row)} for row in small.edges]
            seen_paths=set()
            for part in partitions(n):
                counts=Counter(part)
                fv=sv=0
                projected=[]
                while fv!=full.terminal:
                    ell=full.vertices[fv][0]
                    j,e=full_rows[fv][counts[ell]]
                    if fv not in removed:
                        projected.append((ell,e.k,e.u,e.shift))
                    fv=e.target
                actual=[]
                path_edges=[]
                while sv!=small.terminal:
                    ell=small.vertices[sv][0]
                    j,e=small_rows[sv][counts[ell]]
                    path_edges.append((sv,j))
                    actual.append((ell,e.k,e.u,e.shift))
                    sv=e.target
                check(actual==projected,('complete partition path contraction',n,part))
                path_signature=tuple(path_edges)
                check(path_signature not in seen_paths,('injective complete partition paths',n,part))
                seen_paths.add(path_signature)
                _,valuation=raw_weight(part)
                for K in support:
                    if valuation<small.baseline+K:
                        support[K].update(path_edges)
                pathchecks+=1
            # Count every DAG path independently, confirming surjectivity as well.
            path_count=[0]*len(small.vertices);path_count[0]=1
            for v,row in enumerate(small.edges):
                for e in row:
                    path_count[e.target]+=path_count[v]
            check(path_count[-1]==len(seen_paths),('partition-path bijection',n))
            for K,expected in support.items():
                actual={(v,j) for v,row in enumerate(small.edges) for j,e in enumerate(row)
                        if small.prefix_cost[v]+small.reduced_cost(v,e)<K}
                check(expected==actual,('exact contracted precision edge support',n,K))
                supportchecks+=1
        for v,row in enumerate(full.edges):
            if v in removed:
                continue
            for e in row:
                if zero_nonforced is None and e.k==0 and e.shift>0:
                    zero_nonforced={'n':n,'vertex':full.vertices[v],'J':e.shift,'outdegree':len(row)}
        vertices+=len(small.vertices)
        edges+=sum(map(len,small.edges))
        forced+=len(removed)
    check(strict_prefix_example is not None and zero_nonforced is not None,'contraction hostiles present')
    return {'status':'PASS','source_sha256':hashlib.sha256((source/'contracted_ordinary.py').read_bytes()).hexdigest(),
            'orders':{'consecutive':[0,50],'additional':[60,80,100]},
            'surviving_vertices_including_terminals':vertices,'surviving_edges':edges,'forced_vertices_removed':forced,
            'modular_full_coordinate_comparisons':modular,'K_values':[1,2,3,8,17],
            'partition_path_bijection_orders':[0,32],'complete_partition_paths':pathchecks,
            'exact_budget_support_equalities':supportchecks,'prefix_P_not_constant_hostile':strict_prefix_example,
            'k_zero_not_forced_hostile':zero_nonforced,'samples':samples,
            'scope':'Direct construction with only deterministic unaffordable/early-terminal T=1,J=0 skips contracted; no arbitrary k=0 deletion, defect quotient, signed excess, runtime or asymptotic claim.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parent.parent/'session64-counting')
    args=parser.parse_args()
    print(json.dumps(run(args.source),indent=2,sort_keys=True))
