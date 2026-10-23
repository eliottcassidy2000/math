#!/usr/bin/env python3
from pathlib import Path
from itertools import product,combinations
from math import isqrt
from collections import Counter
import json,hashlib
D=Path(__file__).resolve().parent
S=next(p for p in (D.parent/'color-holdouts',D.parent/'session42-color-holdouts') if p.exists())
checks=0
def need(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
rows={32:json.loads((S/'train32.json').read_text())}
counts=json.loads((D/'independent-counts.json').read_text())
for rec in counts:
 n=rec['N'];ps=json.loads((S/f'paths{n}.json').read_text());rows[n]=ps
 need(2*len(ps)==rec['oriented'],'independent total')
for n,ps in rows.items():
 need(len(set(map(tuple,ps)))==len(ps),'duplicates')
 for p in ps:
  need(sorted(p)==list(range(1,n+1)) and p[0]<p[-1],'vertex coverage and orientation')
  need(all(isqrt(a+b)**2==a+b for a,b in zip(p,p[1:])),'square adjacency')
# Derive Zeckendorf digits independently, retaining conventional F2=1.
f=[1,2]
while f[-1]<80:f.append(f[-1]+f[-2])
def zeck(n):
 out=[];remaining=n
 while remaining:
  k=max(i for i,v in enumerate(f) if v<=remaining)
  out.append(k+2);remaining-=f[k]
 return out
V=[]
for n in range(1,81):
 z=zeck(n);V.append((n,sum(int(c) for c in bin(n)[2:]),isqrt(n),len(z),z[-1],sum(z)))
E=sorted([(a,b) for b in range(1,81) for a in range(1,b) if isqrt(a+b)**2==a+b],key=lambda t:(t[1],t[0]));idx={e:i for i,e in enumerate(E)}
F=[]
for a,b in E:
 u,v=V[a-1],V[b-1];F.append((a,b,isqrt(a+b),u[1]+v[1],u[2]+v[2],u[3]+v[3],u[4]+v[4]))
encoded={}
for kind in ('vertex','edge'):
 encoded[kind]={}
 for n,ps in rows.items():
  encoded[kind][n]=[(bytes([p[0],p[-1]]),bytes([x-1 for x in p] if kind=='vertex' else [idx[tuple(sorted((a,b)))] for a,b in zip(p,p[1:])])) for p in ps]
ledger=json.loads((S/'ledger.json').read_text());expected={(x['kind'],x['q'],tuple(x['coefficients'])):x for x in ledger}
need(len(expected)==1549,'unique training ledger')
scored={};categories={};witnesses=0
for kind,features in [('vertex',V),('edge',F)]:
 for q in (2,3):
  unique=set();train=0;survive=0;failures=Counter()
  for coeff in product(range(q),repeat=len(features[0])):
   if not any(coeff):continue
   color=bytes(sum(c*x for c,x in zip(coeff,v))%q for v in features)
   ren={};canon=bytes(ren.setdefault(v,len(ren)) for v in color)
   if canon in unique:continue
   unique.add(canon);table=bytes(color)+bytes(256-len(color));key=(kind,q,coeff)
   labels=lambda n:[end+word.translate(table) for end,word in encoded[kind][n]]
   if len(set(labels(32)))<len(rows[32]):
    need(key not in expected,'failed training entry included');continue
   train+=1;need(key in expected,'training success omitted');item=expected[key]
   actualcounts=[];failure=0
   for n in range(32,37):
    lab=labels(n);actualcounts.append(len(set(lab)))
    if len(set(lab))<len(rows[n]) and not failure:failure=n
   need(actualcounts==item['counts32to36'],'all exact score counts')
   need(failure==item['first_failure'],'failure horizon')
   if failure:
    i,j=item['witness_indices'];lab=labels(failure)
    need(i!=j and lab[i]==lab[j] and rows[failure][i]!=rows[failure][j],'actual collision witness');witnesses+=1
   else:survive+=1
   failures[failure]+=1;scored[key]=failure
  categories[kind+str(q)]={'distinct_candidates':len(unique),'training_survivors':train,'all_holdout_survivors':survive,'first_failures':dict(sorted(failures.items()))}
need(len(scored)==1549,'all ledger scores audited')
need(sum(a['all_holdout_survivors'] for a in categories.values())==1263,'survivor total')
# Unsigned incidence map: test all labeled simple graphs through order5 over F2,F3.
graphs=0;kernelcontrols=0
for n in range(1,6):
 pairs=list(combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  edges=[p for i,p in enumerate(pairs) if mask>>i&1];adj=[[] for _ in range(n)]
  for a,b in edges:adj[a].append(b);adj[b].append(a)
  parity={};components=0;bipartite=0
  for v in range(n):
   if v in parity:continue
   components+=1;ok=True;parity[v]=0;queue=[v]
   for a in queue:
    for b in adj[a]:
     if b not in parity:parity[b]=parity[a]^1;queue.append(b)
     elif parity[b]==parity[a]:ok=False
   bipartite+=ok
  for q in (2,3):
   size=sum(all((c[a]+c[b])%q==0 for a,b in edges) for c in product(range(q),repeat=n))
   need(size==q**(components if q==2 else bipartite),'unsigned incidence kernel')
   kernelcontrols+=1
  graphs+=1
for q in (2,3):
 for n in range(1,9):
  for colors in product(range(q),repeat=n):
   word=[(a+b)%q for a,b in zip(colors,colors[1:])];recovered=[colors[0]]
   for s in word:recovered.append((s-recovered[-1])%q)
   need(tuple(recovered)==colors,'start-marked color recurrence')
source_pins={name:hashlib.sha256((S/name).read_bytes()).hexdigest() for name in ['enumerate.cpp','test.cpp','verify.py','ledger.json','carrier_hashes.json','train32.json','receipt.json']+[f'paths{n}.json' for n in range(33,37)]}
print(json.dumps({'status':'PASS','checks':checks,'independently_enumerated_carrier_counts':counts,'fixed_grammar_scores':categories,'ledger_entries':len(scored),'collision_witnesses':witnesses,'incidence_graphs':graphs,'incidence_kernel_controls':kernelcontrols,'source_pins':source_pins,'scope':'full Q33..36 census independently counted; Q32 completeness inherited; all finite scores verified, no all-size decoder claim'},indent=2))
