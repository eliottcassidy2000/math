#!/usr/bin/env python3
"""Independent standard-library exact replay. No LP or numerical dependencies."""
import gzip,json,hashlib,itertools,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def check(x,msg):
 if not x:raise RuntimeError(msg)

def catalogue(c):
 # Independently enumerate all allowed relation subsets and deduplicate closure.
 n=9+c;optional=list(itertools.combinations(range(2,n),2))
 out=set()
 for mask in range(1<<len(optional)):
  rel={(2,1)}|{(0,a) for a in range(3+c,n)}|{e for i,e in enumerate(optional) if mask>>i&1}
  while True:
   new=rel|{(a,d) for a,b in rel for t,d in rel if b==t}
   if new==rel:break
   rel=new
  out.add(tuple(sum(1<<a for a,b in rel if b==v) for v in range(n)))
 return out

def replay(p,rec):
 n=len(p);c=n-9;D=[{u for u in range(n) if p[v]>>u&1} for v in range(n)];U=[{v for v in range(n) if u in D[v]} for u in range(n)];deg=[n-1-len(D[v])-len(U[v]) for v in range(n)];check(max(deg)<=7,'degree')
 inc=[(a,b) for a,b in itertools.combinations(range(n),2) if a not in D[b] and b not in D[a]]
 seeds={0,1}|{v for v in range(n) if deg[v]==7}
 def bound(u,v):
  ma=0 if (D[u]|{u})&seeds else 7-deg[u];mb=0 if (U[v]|{v})&seeds else 7-deg[v]
  return len(D[v]-D[u])+len(U[u]-U[v])+ma+mb
 if 'analytic_cycle' in rec:
  arcs=set((a,b) for b in range(n) for a in D[b]);arcs|={(b,a) for a,b in inc for a,b in ((a,b),(b,a)) if bound(a,b)<=1};arcs|={(0,1)}|{(0,z) for z in range(3,n-6)}|{(z,1) for z in range(3,n-6)}
  check(set(map(tuple,rec['analytic_cycle']))==arcs,'cycle arcs')
  while True:
   nxt=arcs|{(a,d) for a,b in arcs for t,d in arcs if b==t}
   if nxt==arcs:break
   arcs=nxt
  check(any((v,v) in arcs for v in range(n)),'no cycle');return
 if 'analytic_balanced' in rec:
  a,b=rec['analytic_balanced'];check((min(a,b),max(a,b)) in inc,'comparable analytic witness');check(bound(a,b)<=1 and bound(b,a)<=1,'bad analytic bound');return
 # Independent recursive extension generation with mutable remaining set.
 words=[]
 def extend(w,remaining):
  if not remaining:words.append(tuple(w));return
  for v in sorted(remaining):
   if not(D[v]&remaining):extend(w+[v],remaining-{v})
 extend([],set(range(n)))
 valid=[];analytic=0
 for w in words:
  rank={v:i for i,v in enumerate(w)};s=rank[1]-rank[0]
  if s<=0 or any(abs(rank[a]-rank[b])>s for a,b in inc):continue
  oriented=[(a,b) if rank[a]<rank[b] else (b,a) for a,b in inc]
  if any(bound(a,b)<=1 for a,b in oriented):analytic+=1
  else:valid.append((w,rank))
 check(analytic==rec['analytic_directed_orders'],'analytic order count');check(len(valid)==rec['residual_orders'],'residual count')
 if not valid:check(not rec['certificates'],'extraneous certs');return
 # Count all generator conditional events independently, with only retained ones.
 sat={v for v in range(n) if deg[v]==7};gens={}
 for w in words:
  r={v:i for i,v in enumerate(w)}
  for i in range(min(r[0],r[1])+1):
   for j in range(max(r[0],r[1])+1,n+1):
    u,v=w[:i],w[j:]
    if sat.intersection(u+v):continue
    item=gens.setdefault((u,v),[0]* (1+len(rec['certificates'])))
    item[0]+=1
    for k,cert in enumerate(rec['certificates']):item[k+1]+=sum(weight*int(r[a]<r[b]) for a,b,weight in cert['terms'])
 check(len(gens)==rec['generators'],'generator count')
 for k,cert in enumerate(rec['certificates']):
  check(cert['terms'] and all(weight>0 and (min(a,b),max(a,b)) in inc for a,b,weight in cert['terms']),'invalid terms')
  total=sum(weight for a,b,weight in cert['terms']);vals=[3*data[k+1]-2*total*data[0] for data in gens.values()]
  check(max(vals)<=0 and max(vals)==cert['max_margin'],'generator inequality')
 stream=hashlib.sha256()
 for w,r in valid:
  k=next((k for k,cert in enumerate(rec['certificates']) if all(r[a]<r[b] for a,b,weight in cert['terms'])),None)
  check(k is not None,'uncovered majority order');stream.update(json.dumps([w,k]).encode())
 check(stream.hexdigest()==rec['case_digest'],'order assignment digest')

def main():
 start=time.time();raw=(ROOT/'certificates.json.gz').read_bytes();records=json.loads(gzip.decompress(raw));summary=json.loads((ROOT/'summary.json').read_text());check(summary['status']=='ALL_CERTIFIED','incomplete')
 totals=[]
 for c,expected in enumerate((96428,)):
  cat=catalogue(c);rr=[r for r in records if r['c']==c];check(len(cat)==expected==len(rr),'catalogue cardinality');check({tuple(r['pred']) for r in rr}==cat,'catalogue equality')
  for r in rr:replay(tuple(r['pred']),r)
  totals.append({'c':c,'cores':len(rr),'analytic_balanced_cores':sum('analytic_balanced' in r for r in rr),'analytic_cycle_cores':sum('analytic_cycle' in r for r in rr),'directed_screen_orders':sum(r.get('analytic_directed_orders',0) for r in rr),'exact_generator_orders':sum(r.get('residual_orders',0) for r in rr),'distinct_generator_certificates':sum(len(r.get('certificates',[])) for r in rr)})
  print(json.dumps(totals[-1]),flush=True)
 out={'status':'EXACT_REPLAY_PASSED','certificate_sha256':hashlib.sha256(raw).hexdigest(),'catalogue_and_replay':totals,'seconds':time.time()-start}
 (ROOT/'replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
