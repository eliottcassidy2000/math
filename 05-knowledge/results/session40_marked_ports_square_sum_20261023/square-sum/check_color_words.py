#!/usr/bin/env python3
"""Verify birth-color word finite positives and explicit order80 hostile.
Run next to receipt.json and the two witness files. No external dependencies.
"""
import json
from pathlib import Path
from math import isqrt
from collections import Counter,defaultdict
D=Path(__file__).resolve().parent

def check(p,msg):
 if not p:raise RuntimeError(msg)
def valid(p,n):return len(p)==n and set(p)==set(range(1,n+1)) and all(isqrt(a+b)**2==a+b for a,b in zip(p,p[1:]))
R=json.loads((D/'receipt.json').read_text())['rows'];W=json.loads((D/'birth_witnesses_25_80.json').read_text());H=json.loads((D/'COLOR_WORD_HOSTILE_80.json').read_text())
check(set(map(int,W))==set(range(25,81)),'birth universe')
for n,p in W.items():check(valid(p,int(n)),f'invalid birth witness {n}')
st={r['n']:r['count']>0 for r in R};st.update({int(n):True for n in W});deg={n:isqrt(2*n-1)-isqrt(n) for n in st}
def color(a,b):
 v=max(a,b)
 return 1 if st[v] else 2 if deg[v]==2 else 3
out=[];collision26=None
for r in R[24:]:
 words=Counter();hist=defaultdict(list)
 for p in r['paths']:
  if p[0]>p[-1]:p=p[::-1]
  w=tuple(color(a,b) for a,b in zip(p,p[1:]));words[(p[0],p[-1],w)]+=1;hist[(p[0],p[-1],tuple(w.count(c) for c in (1,2,3)))].append(p)
 check(max(words.values())==1,f'ordered word positive fails {r["n"]}')
 if r['n']==26:
  collision26=next(dict(endpoints=[s,t],histogram=h,paths=p) for (s,t,h),p in hist.items() if len(p)>1)
 out.append(dict(n=r['n'],paths=r['count'],endpoint_word_labels=len(words),endpoint_histogram_labels=len(hist),largest_histogram_fibre=max(map(len,hist.values()))))
check(H['n']==80 and len(H['paths'])==2,'hostile scope')
for p in H['paths']:
 check(valid(p,80),'hostile invalid');check((p[0],p[-1])==(1,2),'hostile endpoints');check(all(max(a,b)>=25 for a,b in zip(p,p[1:])),'old edge present');check(all(color(a,b)==1 for a,b in zip(p,p[1:])),'word not allfirst')
check(H['paths'][0]!=H['paths'][1],'hostiles identical')
print(json.dumps(dict(status='FINITE-EXACT positive through32; explicit counterexample at80, not claimed minimal',birth_witnesses_checked=len(W),finite_rows=out,histogram_hostile26=collision26,word_hostile80=dict(endpoints=[1,2],color_word=[1]*79,paths=H['paths'])),indent=2))
