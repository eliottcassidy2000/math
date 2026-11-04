#!/usr/bin/env python3
"""Exact prime-valuation/Farkas certificate; no solver required for replay."""
import json
S=(3,5,7,11,13)
W={2:-4,3:-9,5:-13,7:-14,11:-17,13:-18}
def factor(n):
    d={};p=2
    while p*p<=n:
        while n%p==0: d[p]=d.get(p,0)+1;n//=p
        p+=1
    if n>1:d[n]=d.get(n,0)+1
    return d

def sigma(p,a):return sum(p**j for j in range(a+1))
def column(p,a):
    d=factor(sigma(p,a));d[p]=d.get(p,0)-a
    return {q:e for q,e in d.items() if e}
def score(d):return sum(W.get(q,256)*e for q,e in d.items())
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
rows=[]
for p in S:
    for a in range(1,7):
        d=column(p,a);t=score(d)
        require(t>=0, 'dual inequality')
        rows.append({'p':p,'a':a,'sigma':sigma(p,a),'column':d,'dual_score':t})
require(score({2:1})==-4,'target strict separation')
# Positive controls: even perfect 6 and 28 have exact valuation target e_2.
for blocks in [[(2,1),(3,1)],[(2,2),(7,1)]]:
    total={}
    for p,a in blocks:
        for q,e in column(p,a).items():total[q]=total.get(q,0)+e
    total={q:e for q,e in total.items() if e}
    require(total=={2:1},'even perfect control')
    require(any(score(column(p,a))<0 for p,a in blocks),'must not exclude positive control')
# Hostile: naive integer weights almost work, but 7^3 exposes the missing debt.
naive={2:-1,3:-2,5:-3,7:-3,11:-4,13:-4}
naive_score=sum(naive.get(q,64)*e for q,e in column(7,3).items())
require(naive_score==-1,'hostile weak separator')
print(json.dumps({'status':'FINITE-EXACT certificate for 30 blocks, arbitrary nonnegative multiplicities','weights':W,'outside_weight':256,'target_score':-4,'minimum_column_score':min(r['dual_score'] for r in rows),'naive_hostile_score':naive_score,'rows':rows},indent=2,sort_keys=True))
