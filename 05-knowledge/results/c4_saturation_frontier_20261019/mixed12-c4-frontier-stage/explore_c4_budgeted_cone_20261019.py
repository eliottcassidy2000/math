#!/usr/bin/env python3
"""Bounded c=4 exploration; numerical LP candidates require exact checking.
The canonical core census is finite; the cone search uses a fixed sample/order
budget and stops on an exact obstruction to the generator-cone screen only.
"""
from itertools import combinations, product
from pathlib import Path
from collections import Counter
from fractions import Fraction
import argparse
import hashlib
import importlib.util
import json
import random
import time


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def closure(pred):
    p=list(pred)
    for k in range(len(p)):
        for j in range(len(p)):
            if p[j]>>k&1:p[j]|=p[k]
    return tuple(p)


def common_posets(c):
    edges=list(combinations(range(c),2));seen=set()
    for mask in range(1<<len(edges)):
        p=[0]*c
        for i,(a,b) in enumerate(edges):
            if mask>>i&1:p[b]|=1<<a
        q=closure(p)
        if q not in seen:seen.add(q);yield q


def canonical_cores(c):
    full=(1<<c)-1
    for zpred in common_posets(c):
        ideals=[mask for mask in range(1<<c) if all(not(mask>>i&1) or not(zpred[i]&~mask) for i in range(c))]
        upsets=[full^mask for mask in ideals]
        for lower,D1,D2 in product(upsets,ideals,ideals):
            for upper_chain in (0,1):
                if upper_chain and D1&~D2:continue
                for ba1,ba2 in product((0,1),repeat=2):
                    if lower&D1 and not ba1:continue
                    if lower&D2 and not ba2:continue
                    if upper_chain and ba1>ba2:continue
                    p=[0,4,0]+[(zpred[i]<<3)|(4 if lower>>i&1 else 0) for i in range(c)]
                    p.extend((1|(D1<<3)|(4 if ba1 else 0),1|(D2<<3)|(4 if ba2 else 0)|((1<<(c+3)) if upper_chain else 0)))
                    yield tuple(p)


def degrees(p):
    n=len(p)
    return tuple(n-1-p[v].bit_count()-sum(bool(p[u]>>v&1) for u in range(n)) for v in range(n))


def census():
    result=[];kept=[]
    for c in range(5):
        count=survivors=0;isolated=Counter();seen=set()
        for p in canonical_cores(c):
            require(p not in seen,'Duplicate canonical core');seen.add(p)
            require(closure(p)==p,'Canonical construction not transitively closed')
            count+=1;deg=degrees(p)
            missing=sum(d>7 for d in deg);isolated[str(missing)]+=1
            if max(deg)<=7:
                survivors+=1
                if c==4:kept.append(p)
        if c<4:require(count==(7,40,357,4824)[c],'Smaller-family census mismatch')
        result.append({'common_count':c,'canonical_cores':count,'core_degree_at_most_7':survivors,'vertices_over_budget_histogram':dict(sorted(isolated.items()))})
    return result,kept


def load_verifier():
    path=Path(__file__).resolve().parents[1]/'mixed12-freeze-stage/frozen_mixed12_20261019/verify_mixed12_witness_migration_20261019.py'
    spec=importlib.util.spec_from_file_location('mixed12_verified',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def exact_cone_obstruction(rows, columns, keys, oriented, M, pred, solution):
    support=[j for j,v in enumerate(solution) if v>1e-9]
    for bound in (1000,1000000,1000000000):
        weights={j:Fraction(float(solution[j])).limit_denominator(bound) for j in support}
        weights={j:w for j,w in weights.items() if w>0}
        mass=sum(weights.values())
        weights={j:w/mass for j,w in weights.items()}
        margins=[sum(weights[j]*Fraction(int(row[j]),columns[j][0]) for j in weights) for row in rows]
        if min(margins)>0:
            probabilities=[(m+2)/3 for m in margins]
            return {'status':'EXACT_GENERATOR_CONE_OBSTRUCTION','scope':'Failure of coefficientwise weighted-comparison cone screen; not an actual exterior or conjecture counterexample.','predecessor_masks':list(pred),'core_incomparability_degrees':list(degrees(pred)),'candidate_majority_order':list(M),'core_extensions':sum(1 for _ in load_verifier().extensions(pred)),'generators':len(keys),'positive_mixture':[{'prefix':list(keys[j][0]),'suffix':list(keys[j][1]),'conditional_word_count':columns[j][0],'convex_weight':str(w),'generator_weight':str(w/columns[j][0])} for j,w in sorted(weights.items())],'forward_comparison_probabilities':[{'pair':list(pair),'probability':str(p)} for pair,p in zip(oriented,probabilities)],'minimum_probability':str(min(probabilities)),'minimum_3p_minus_2':str(min(margins))}
    return None


def search(kept, max_cores, max_orders):
    import numpy as np
    from scipy.optimize import linprog
    verifier=load_verifier();rng=random.Random(202610190734)
    # Sparse cores are the largest cones; visit a deterministic stratified sample.
    by_edges={}
    for p in kept:by_edges.setdefault(sum(x.bit_count() for x in p),[]).append(p)
    chosen=[]
    for group in by_edges.values():rng.shuffle(group)
    while len(chosen)<min(max_cores,len(kept)):
        for edge_count in sorted(by_edges):
            group=by_edges[edge_count]
            if group and len(chosen)<max_cores:chosen.append(group.pop())
    stats={'core_budget':max_cores,'order_budget':max_orders,'cores_examined':0,'orders_examined':0,'single_row_certificates':0,'two_row_certificates':0,'LP_calls':0,'nonpositive_LP_optima':0,'positive_LP_candidates':0}
    for pred in chosen:
        words=verifier.extensions(pred);n=len(pred)
        inc=[(a,b) for a,b in combinations(range(n),2) if not(pred[b]>>a&1 or pred[a]>>b&1)]
        keys,columns,rankwords=verifier.generators(pred,words,inc)
        col=np.array(columns,dtype=np.int64);E=col[:,0];forward=(3*col[:,1:]-2*E[:,None]).T;backward=(E[:,None]-3*col[:,1:]).T
        admissible=[]
        for M,ranks in rankwords:
            span=ranks[1]-ranks[0]
            if span>0 and all(abs(ranks[a]-ranks[b])<=span for a,b in inc):admissible.append((M,ranks))
        rng.shuffle(admissible);stats['cores_examined']+=1
        for M,ranks in admissible:
            if stats['orders_examined']>=max_orders:return stats,None
            stats['orders_examined']+=1
            reverse=np.array([ranks[a]>ranks[b] for a,b in inc]);rows=np.where(reverse[:,None],backward,forward)
            if np.any(rows.max(axis=1)<=0):stats['single_row_certificates']+=1;continue
            found=False
            for i in range(len(rows)):
                if np.any((rows[i]+rows[i:]).max(axis=1)<=0):found=True;break
            if found:stats['two_row_certificates']+=1;continue
            stats['LP_calls']+=1
            normalized=rows/E[None,:]
            # Maximize a common strict majority margin in a convex mixture of
            # normalized generator columns, with a free scalar margin.
            objective=np.r_[np.zeros(len(keys)),-1.]
            Aub=np.c_[-normalized,np.ones(len(rows))]
            eq=np.r_[np.ones(len(keys)),0.][None,:]
            sol=linprog(objective,A_ub=Aub,b_ub=np.zeros(len(rows)),A_eq=eq,b_eq=[1.],bounds=[(0,None)]*len(keys)+[(None,None)],method='highs')
            require(sol.success,'LP solver failure: '+str(sol.message))
            if sol.x[-1]>1e-9:
                stats['positive_LP_candidates']+=1
                oriented=[(b,a) if rev else (a,b) for (a,b),rev in zip(inc,reverse)]
                exact=exact_cone_obstruction(rows,columns,keys,oriented,M,pred,sol.x[:-1])
                if exact:return stats,exact
            else:stats['nonpositive_LP_optima']+=1
        print(json.dumps({'progress':stats}),flush=True)
    return stats,None


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--cores',type=int,default=80);parser.add_argument('--orders',type=int,default=20000);parser.add_argument('--census-only',action='store_true');args=parser.parse_args()
    counts,kept=census();root=Path(__file__).resolve().parent
    base={'scope':'Exploration only, no extension of the c<=3 theorem.','core_census':counts,'c4_degree_budget_survivors':len(kept),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (root/'c4_budgeted_core_census_20261019.json').write_text(json.dumps(base,indent=2)+'\n');print(json.dumps(base),flush=True)
    if not args.census_only:
        stats,witness=search(kept,args.cores,args.orders)
        result={**base,'bounded_search':stats,'obstruction':witness}
        (root/'c4_budgeted_cone_exploration_20261019.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
