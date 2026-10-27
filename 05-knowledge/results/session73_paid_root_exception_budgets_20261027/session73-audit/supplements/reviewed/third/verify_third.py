#!/usr/bin/env python3
import json
from pathlib import Path
from functools import lru_cache
root=Path(__file__).resolve().parent;r=json.loads((root/'receipt.json').read_text());b=json.loads((root/'input_two_step_receipt.json').read_text())['best_literal_candidate'];n=11

def check(x,msg):
 if not x:raise RuntimeError(msg)
def fas(T):
 @lru_cache(None)
 def solve(S):
  return 0 if not S else min(solve(S-{v})+sum(T[u][v] for u in S if u!=v) for v in S)
 return solve(frozenset(range(n)))
def tau(T,P):
 check(set(P)==set(range(n)) and len(P)==n,'not permutation');check(all(T[i][j] for i,j in zip(P,P[1:])),'not legal')
 pe={frozenset(p) for p in zip(P,P[1:])}
 return [[False if i==j else T[i][j] if frozenset((i,j)) in pe else not T[i][j] for j in range(n)] for i in range(n)]
def rows(T):return [sum(1<<j for j in range(n) if T[i][j]) for i in range(n)]
sq={x*x%n for x in range(1,n)};T=[[i!=j and (j-i)%n in sq for j in range(n)] for i in range(n)];values=[fas(T)]
for P in (b['P'],b['Q'],r['best_target']['R']):T=tau(T,P);values.append(fas(T))
check(values==[20,15,13,12],'FAS sequence');check(rows(T)==r['best_target']['rows'],'last rows')
check(r['rank_draws']==32 and r['path_draws']==128 and r['distinct_targets']==93 and r['final_rank_draws']==0,'domain')
check((root/'receipt.json').read_bytes()==(root/'receipt.optimized.json').read_bytes(),'optimized mismatch')
print(json.dumps({'status':'VERIFIED bounded no-witness receipt and literal best third path; no depth3 impossibility or upper-bound claim','independent_prepend_FAS':values,'additional_path':r['best_target']['R'],'whole_complement_clock':3},indent=2,sort_keys=True))
