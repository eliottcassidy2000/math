#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,subprocess,sys,tempfile
from math import isqrt
from itertools import product
from collections import Counter
D=Path(__file__).resolve().parent

def ck(x,m):
 if not x:raise RuntimeError(m)
H=json.loads((D/'carrier_hashes.json').read_text())
if '--fresh' in sys.argv or any(not (D/f'paths{n}.json').exists() for n in range(33,37)):
 with tempfile.TemporaryDirectory() as tmp:
  exe=Path(tmp)/'enumerate';subprocess.run(['g++','-O3','-std=c++17',str(D/'enumerate.cpp'),'-o',str(exe)],check=True)
  for n in range(33,37):
   with (D/f'paths{n}.json').open('w') as out:
    run=subprocess.run([str(exe),str(n)],stdout=out,stderr=subprocess.PIPE,text=True,check=True)
   (D/f'enumerate{n}.log').write_text(run.stderr)
F=[0,1,1]
while F[-1]<10000:F.append(F[-1]+F[-2])
def zeck(n):
 out=[]
 for i in range(len(F)-1,1,-1):
  if F[i]<=n:n-=F[i];out.append(i)
 return out
E=[(a,b) for b in range(1,81) for a in range(1,b) if isqrt(a+b)**2==a+b];ei={e:i for i,e in enumerate(E)}
def vf(n):
 z=zeck(n);return (n,n.bit_count(),isqrt(n),len(z),min(z),sum(z))
def ef(e):
 a,b=e;A=vf(a);B=vf(b);return (a,b,isqrt(a+b),A[1]+B[1],A[2]+B[2],A[3]+B[3],A[4]+B[4])
rows={32:json.loads((D/'train32.json').read_text())};encoded={}
for n in range(33,37):
 raw=(D/f'paths{n}.json').read_bytes();ck(hashlib.sha256(raw).hexdigest()==H[str(n)]['sha256'],'carrier digest');rows[n]=json.loads(raw);ck(len(rows[n])==H[str(n)]['count'],'carrier count')
for n,ps in rows.items():
 ck(len(set(map(tuple,ps)))==len(ps),'duplicate')
 for p in ps:
  ck(len(p)==n and set(p)==set(range(1,n+1)) and p[0]<p[-1],'path owner/orientation')
  ck(all(isqrt(a+b)**2==a+b for a,b in zip(p,p[1:])),'nonedge')
 encoded[n]={'vertex':[(bytes((p[0],p[-1])),bytes(x-1 for x in p)) for p in ps], 'edge':[(bytes((p[0],p[-1])),bytes(ei[tuple(sorted(e))] for e in zip(p,p[1:]))) for p in ps]}
ledger=json.loads((D/'ledger.json').read_text());expected={(a['kind'],a['q'],tuple(a['coefficients'])):a for a in ledger};ck(len(expected)==1549,'ledger size')
verified={};distinct={};witnesses=[]
for kind,fs in [('vertex',[vf(n) for n in range(1,81)]),('edge',[ef(e) for e in E])]:
 for q in (2,3):
  seen=set();survivors=0
  for co in product(range(q),repeat=len(fs[0])):
   if not any(co):continue
   colors=bytes(sum(a*b for a,b in zip(co,f))%q for f in fs);ren={};canonical=bytes(ren.setdefault(v,len(ren)) for v in colors)
   if canonical in seen:continue
   seen.add(canonical);lut=colors+bytes(256-len(colors));key=(kind,q,co)
   def labels(n):return [end+w.translate(lut) for end,w in encoded[n][kind]]
   if len(set(labels(32)))!=392:
    ck(key not in expected,'false training survivor');continue
   ck(key in expected,'missing training survivor');a=expected[key];counts=[392];first=0
   for n in range(33,37):
    lab=labels(n);count=len(set(lab));counts.append(count)
    if count!=len(lab) and not first:first=n
   ck(counts==a['counts32to36'] and first==a['first_failure'],'holdout count mismatch')
   if first:
    i,j=a['witness_indices'];lab=labels(first);ck(i!=j and lab[i]==lab[j],'collision witness')
   verified[key]=first
  distinct[f'{kind}{q}']=len(seen)
summary={}
for kind in ('vertex','edge'):
 for q in (2,3):
  cohort=[a for a in ledger if a['kind']==kind and a['q']==q];freq=Counter(a['first_failure'] for a in cohort)
  summary[f'{kind}{q}']=dict(training_survivors=len(cohort),first_failure_counts={str(k):v for k,v in sorted(freq.items())},survive_all=freq[0])
  for n in (33,34,35,36):
   a=next((a for a in cohort if a['first_failure']==n),None)
   if a:
    i,j=a['witness_indices'];witnesses.append(dict(kind=kind,q=q,coefficients=a['coefficients'],first_failure=n,paths=[rows[n][i],rows[n][j]]))
ck(len(verified)==len(expected),'cohort mismatch')
(D/'witnesses.json').write_text(json.dumps(witnesses,indent=2)+'\n')
print(json.dumps(dict(status='PASS',scope='Complete heldout Q33..Q36, fixed feature grammar selected on Q32 only; no all-size claim',path_counts={str(n):len(p) for n,p in rows.items()},distinct_candidates=distinct,cohorts=summary,total_survivors=sum(v['survive_all'] for v in summary.values()),total_failures=sum(v['training_survivors']-v['survive_all'] for v in summary.values()),witness_count=len(witnesses)),indent=2))
