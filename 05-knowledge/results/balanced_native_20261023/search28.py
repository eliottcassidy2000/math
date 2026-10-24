#!/usr/bin/env python3
"""Bounded richer-supply search. A failure is NOT a support maximum."""
from itertools import combinations
from random import Random
from functools import reduce
from operator import xor
from pathlib import Path
import json, subprocess

ROOT=Path(__file__).resolve().parents[3]
PATH='05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json'
raw=subprocess.check_output(['git','show','0f9dec9c3e1d:'+PATH],cwd=ROOT)
D=json.loads(raw);options={x['source']:x for x in D['records']}
labels={s:(sum(v<<i for i,v in enumerate(x['source_label'][:4])),x['source_label'][4]) for s,x in options.items()}
ids=sorted(labels)
def total(seq):return reduce(xor,(labels[s][0] for s in seq),0),sum(labels[s][1] for s in seq)%128
def proper(A,B):
    for j in range(1,len(A)//2+1):
        bank={total(c):c for c in combinations(A,j)}
        for c in combinations(B,j):
            if total(c) in bank:return bank[total(c)],c
    return None
def main():
    rng=Random(20261023);found={};attempts={};trial_budget=20
    for k in (9,10):
        count=0;cuts={}
        for trial in range(trial_budget):
            z=ids[:];rng.shuffle(z);left=z[:14];right=z[14:];bank={}
            for a in combinations(left,k):bank.setdefault(total(a),[]).append(a)
            for b in combinations(right,k):
                for a in bank.get(total(b),[]):
                    count+=1;p=proper(a,b)
                    if p is None:
                        found[str(2*k)]={'plus':sorted(a),'minus':sorted(b),'sum':total(a),'trial':trial,
                          'owner_clock_cost':sum(options[s]['cost'] for s in a+b)}
                        break
                    cuts[len(p[0])]=cuts.get(len(p[0]),0)+1
                if str(2*k) in found:break
            if str(2*k) in found:break
        attempts[str(2*k)]={'matching_candidates':count,'proper_half_sizes':cuts,'partition_trials':trial+1}
    return dict(status='FINITE-EXACT bounded search, not a complete supply maximum',seed=20261023,
      supply=28,partition_trial_budget=trial_budget,found=found,attempts=attempts)
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
