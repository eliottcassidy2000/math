#!/usr/bin/env python3
"""Slow exhaustive replay, approximately four minutes in the research environment.
Run alongside the main check_square_incidence.py and receipt.json.
Writes COLOR_FIRST_FAILURE_CENSUS.json; printed timings are informational only.
"""
from pathlib import Path
ns={};exec(Path('check_square_incidence.py').read_text().split('rows=[];allP=[];allX=[]')[0],ns)
import json,collections,time
rows=json.load(open('receipt.json'))['rows'];st={r['n']:r['count']>0 for r in rows};deg={r['n']:sum(max(z['edge'])==r['n'] for z in r['edge_counts']) for r in rows};out=[]
for n in range(25,37):
 t=time.time();P=rows[n-1]['paths'] if n<=32 else ns['paths'](n);st[n]=bool(P);deg[n]=len(ns['graph'](n)[n]);D=collections.defaultdict(list);H=collections.Counter()
 for p in P:
  if p[0]>p[-1]:p=p[::-1]
  word=tuple(0 if st[max(a,b)] else 1 if deg[max(a,b)]==2 else 2 for a,b in zip(p,p[1:]));D[(p[0],p[-1],word)].append(p);H[(p[0],p[-1],tuple(word.count(c) for c in range(3)))]+=1
 row=dict(n=n,path_count=len(P),word_labels=len(D),histogram_labels=len(H),max_word_fibre=max(map(len,D.values())),first_collision=next((dict(endpoints=[s,t],word=w,paths=p[:2]) for (s,t,w),p in D.items() if len(p)>1),None))
 out.append(row);print(n,len(P),len(D),round(time.time()-t,2),flush=True);json.dump(out,open('COLOR_FIRST_FAILURE_CENSUS.json','w'),indent=2)
