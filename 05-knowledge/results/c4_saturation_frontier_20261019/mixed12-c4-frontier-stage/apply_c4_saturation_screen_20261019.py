#!/usr/bin/env python3
"""Bounded application of the proved zero-increment rule to c=4 cores.
Four-comparison obstructions concern only that restricted migration screen.
"""
from pathlib import Path
from itertools import combinations
import importlib.util,json,random,hashlib
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('targeted',ROOT/'explore_c4_exception_extensions_20261019.py')
T=importlib.util.module_from_spec(spec);spec.loader.exec_module(T)
F=T.F

def main():
    verifier=F.load_verifier();rng=random.Random(202610190738)
    eligible=[p for p in F.canonical_cores(4) if max(F.degrees(p))<=7 and T.patterns(p)];rng.shuffle(eligible)
    stats={'core_budget':400,'cores_examined':0,'original_generator_columns':0,'retained_generator_columns':0,'cores_with_generator_reduction':0,'orders_examined':0,'single_row_certificates':0,'two_row_certificates':0,'orders_not_screened_by_two_rows':0,'four_comparison_instances':0,'original_weight9_instances_passing':0}
    obstruction=None;first_weight9_failure=None;retention_examples=[]
    for pred in eligible[:400]:
        words=verifier.extensions(pred);inc=[(a,b) for a,b in combinations(range(9),2) if not(pred[b]>>a&1 or pred[a]>>b&1)]
        keys,columns,rankwords=verifier.generators(pred,words,inc)
        saturated={v for v,d in enumerate(F.degrees(pred)) if d==7}
        keep=[j for j,(u,v) in enumerate(keys) if not(saturated.intersection(u) or saturated.intersection(v))]
        stats['cores_examined']+=1;stats['original_generator_columns']+=len(keys);stats['retained_generator_columns']+=len(keep)
        if len(keep)<len(keys):
            stats['cores_with_generator_reduction']+=1
            if len(retention_examples)<3:retention_examples.append({'predecessor_masks':list(pred),'saturated_vertices':sorted(saturated),'original_columns':len(keys),'retained_columns':len(keep)})
        kept_keys=[keys[j] for j in keep];kept_columns=[columns[j] for j in keep]
        array=np.array(kept_columns,dtype=np.int64);E=array[:,0];forward=(3*array[:,1:]-2*E[:,None]).T;backward=(E[:,None]-3*array[:,1:]).T
        for z,a1,a2 in T.patterns(pred):
            candidate=next(((M,r) for M,r in rankwords if r[1]>r[0] and all(abs(r[a]-r[b])<=r[1]-r[0] for a,b in inc) and r[2]<r[0]<r[z]<r[1]<r[a2]<r[a1]),None)
            if candidate is None:continue
            M,r=candidate;stats['four_comparison_instances']+=1
            pairs=((0,z),(z,1),(1,a2),(a2,a1));rows4=np.array([forward[inc.index(tuple(sorted(pair)))] if pair[0]<pair[1] else backward[inc.index(tuple(sorted(pair)))] for pair in pairs])
            deficits=np.array([3,4,1,1])@rows4
            if max(deficits)<=0:stats['original_weight9_instances_passing']+=1
            elif first_weight9_failure is None:
                j=int(np.argmax(deficits));first_weight9_failure={'predecessor_masks':list(pred),'pattern':[z,a1,a2],'prefix':list(kept_keys[j][0]),'suffix':list(kept_keys[j][1]),'weighted_generator_margin':int(deficits[j]),'E':int(E[j])}
            if obstruction is None:
                norm=rows4/E[None,:];sol=linprog(np.r_[np.zeros(len(keep)),-1.],A_ub=np.c_[-norm,np.ones(4)],b_ub=np.zeros(4),A_eq=np.r_[np.ones(len(keep)),0.][None,:],b_eq=[1.],bounds=[(0,None)]*len(keep)+[(None,None)],method='highs')
                F.require(sol.success,'Four-comparison LP failure')
                if sol.x[-1]>1e-9:
                    obstruction=F.exact_cone_obstruction(rows4,kept_columns,kept_keys,pairs,M,pred,sol.x[:-1])
                    if obstruction:
                        obstruction['status']='EXACT_FOUR_COMPARISON_SCREEN_OBSTRUCTION'
                        obstruction['scope']='Only the four fixed exceptional migration comparisons fail in the saturation-pruned generator cone. This is not a full-comparison cone obstruction, an actual exterior, or a conjecture counterexample.'
                        obstruction['saturated_vertices']=sorted(saturated)
        target=[]
        for M,r in rankwords:
            span=r[1]-r[0]
            if span>0 and all(abs(r[a]-r[b])<=span for a,b in inc) and any(r[2]<r[0]<r[z]<r[1]<r[a2]<r[a1] for z,a1,a2 in T.patterns(pred)):target.append((M,r))
        rng.shuffle(target)
        for M,r in target[:300]:
            stats['orders_examined']+=1;reverse=np.array([r[a]>r[b] for a,b in inc]);rows=np.where(reverse[:,None],backward,forward)
            if np.any(rows.max(axis=1)<=0):stats['single_row_certificates']+=1;continue
            if any(np.any((rows[i]+rows[i:]).max(axis=1)<=0) for i in range(len(rows))):stats['two_row_certificates']+=1;continue
            stats['orders_not_screened_by_two_rows']+=1
    result={'scope':'Exploration only; the zero-increment lemma is proved separately. No c=4 theorem.','verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'stats':stats,'retention_examples':retention_examples,'first_original_weight9_failure':first_weight9_failure,'four_comparison_obstruction':obstruction}
    (ROOT/'c4_saturation_screen_exploration_20261019.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
