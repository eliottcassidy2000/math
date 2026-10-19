#!/usr/bin/env python3
"""Bounded targeted exploration of seven-thin c=4 extensions of the c=1 exception.
No implication for all c=4 cores. LP-derived findings are accepted only after
exact rational verification; cone obstructions are not actual posets.
"""
from itertools import combinations
from fractions import Fraction
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import random
import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frontier',ROOT/'explore_c4_budgeted_cone_20261019.py')
F=importlib.util.module_from_spec(spec);spec.loader.exec_module(F)


def patterns(p):
    if p[8]>>7&1:return []
    result=[]
    for a1,a2 in ((7,8),(8,7)):
        for z in range(3,7):
            if p[z]>>2&1 and p[a1]>>z&1 and not(p[a2]>>z&1) and p[a2]>>2&1:
                result.append((z,a1,a2))
    return result


def exact_dual(rows, candidate):
    support=[i for i,w in enumerate(candidate) if w>1e-9]
    for cap in (1000,1000000,1000000000):
        weights={i:Fraction(float(candidate[i])).limit_denominator(cap) for i in support}
        denominator=math.lcm(*(w.denominator for w in weights.values()))
        integers={i:int(w*denominator) for i,w in weights.items() if w>0}
        common=math.gcd(*integers.values())
        integers={i:w//common for i,w in integers.items()}
        values=[sum(w*int(rows[i,j]) for i,w in integers.items()) for j in range(rows.shape[1])]
        if max(values)<=0:return integers,max(values),min(values)
    return None


def main():
    verifier=F.load_verifier();rng=random.Random(202610190738)
    eligible=[p for p in F.canonical_cores(4) if max(F.degrees(p))<=7 and patterns(p)]
    rng.shuffle(eligible);selected=eligible[:400]
    stats={'eligible_exception_extension_cores':len(eligible),'core_budget':400,'order_budget':60000,'per_core_order_budget':300,'cores_examined':0,'orders_examined':0,'single_row_certificates':0,'two_row_certificates':0,'exact_LP_dual_certificates':0,'LP_calls':0,'unresolved_LP_cases':0}
    certificates=[];obstruction=None
    for pred in selected:
        words=verifier.extensions(pred);n=9
        inc=[(a,b) for a,b in combinations(range(n),2) if not(pred[b]>>a&1 or pred[a]>>b&1)]
        keys,columns,rankwords=verifier.generators(pred,words,inc)
        column=np.array(columns,dtype=np.int64);E=column[:,0];forward=(3*column[:,1:]-2*E[:,None]).T;backward=(E[:,None]-3*column[:,1:]).T
        target=[]
        for M,r in rankwords:
            span=r[1]-r[0]
            if span<=0 or any(abs(r[a]-r[b])>span for a,b in inc):continue
            if any(r[2]<r[0]<r[z]<r[1]<r[a2]<r[a1] for z,a1,a2 in patterns(pred)):target.append((M,r))
        rng.shuffle(target);stats['cores_examined']+=1
        for M,ranks in target[:300]:
            if stats['orders_examined']>=60000:break
            stats['orders_examined']+=1
            reverse=np.array([ranks[a]>ranks[b] for a,b in inc]);rows=np.where(reverse[:,None],backward,forward)
            if np.any(rows.max(axis=1)<=0):stats['single_row_certificates']+=1;continue
            found=False
            for i in range(len(rows)):
                if np.any((rows[i]+rows[i:]).max(axis=1)<=0):found=True;break
            if found:stats['two_row_certificates']+=1;continue
            stats['LP_calls']+=1
            normalized=rows/E[None,:];objective=np.r_[np.zeros(len(keys)),-1.]
            sol=linprog(objective,A_ub=np.c_[-normalized,np.ones(len(rows))],b_ub=np.zeros(len(rows)),A_eq=np.r_[np.ones(len(keys)),0.][None,:],b_eq=[1.],bounds=[(0,None)]*len(keys)+[(None,None)],method='highs')
            F.require(sol.success,'LP did not complete')
            oriented=[(b,a) if rev else (a,b) for (a,b),rev in zip(inc,reverse)]
            if sol.x[-1]>1e-9:
                obstruction=F.exact_cone_obstruction(rows,columns,keys,oriented,M,pred,sol.x[:-1])
                if obstruction:break
            else:
                exact=exact_dual(rows,-sol.ineqlin.marginals)
                if exact:
                    weights,hi,lo=exact;stats['exact_LP_dual_certificates']+=1
                    certificates.append({'predecessor_masks':list(pred),'M':list(M),'weighted_forward_pairs':[[*oriented[i],w] for i,w in weights.items()],'maximum_generator_margin':hi,'minimum_generator_margin':lo})
                    continue
            stats['unresolved_LP_cases']+=1
        if obstruction or stats['orders_examined']>=60000:break
    result={'scope':'Bounded exploration of c=4 core-degree<=7 extensions of the known exception. No c=4 theorem, and any positive cone witness is only a screen limitation.','verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'stats':stats,'exact_weighted_certificates':certificates,'obstruction':obstruction}
    (ROOT/'c4_exception_extension_exploration_20261019.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
