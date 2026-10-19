#!/usr/bin/env python3
"""Independent standard-library exact replay. No LP or numerical dependencies."""
import gzip,json,hashlib,itertools,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
EXPECTED_ANALYTIC_SHA256='ace74505712cc4960ac1565570e351b657ec245fcffcbb7523c825f31c9e70b4'
EXPECTED_RESIDUAL_SHA256='8a8cf27be186e9d14bd2b8f4f52faacad706208f10b6cd00abda3251fe647570'

def check(x,msg):
 if not x:raise RuntimeError(msg)

def replay(p,rec):
 n=len(p);c=n-6;D=[{u for u in range(n) if p[v]>>u&1} for v in range(n)];U=[{v for v in range(n) if u in D[v]} for u in range(n)];deg=[n-1-len(D[v])-len(U[v]) for v in range(n)];check(max(deg)<=7,'degree')
 inc=[(a,b) for a,b in itertools.combinations(range(n),2) if a not in D[b] and b not in D[a]]
 seeds={0,1}|{v for v in range(n) if deg[v]==7}
 def bound(u,v):
  ma=0 if (D[u]|{u})&seeds else 7-deg[u];mb=0 if (U[v]|{v})&seeds else 7-deg[v]
  return len(D[v]-D[u])+len(U[u]-U[v])+ma+mb
 if 'analytic_cycle' in rec:
  arcs=set((a,b) for b in range(n) for a in D[b]);arcs|={(b,a) for a,b in inc for a,b in ((a,b),(b,a)) if bound(a,b)<=1};arcs|={(0,1)}|{(0,z) for z in range(4,n-2)}|{(z,1) for z in range(4,n-2)}
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
 valid=[];analytic=0;short=0
 for w in words:
  rank={v:i for i,v in enumerate(w)};s=rank[1]-rank[0]
  if s<=0 or any(abs(rank[a]-rank[b])>s for a,b in inc):continue
  if rec.get('six_thin_theorem') and s<=3:short+=1;continue
  oriented=[(a,b) if rank[a]<rank[b] else (b,a) for a,b in inc]
  if any(bound(a,b)<=1 for a,b in oriented):analytic+=1
  else:valid.append((w,rank))
 check(short==rec.get('short_span_orders',0),'short-span count')
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


def independent_catalogue():
 # Construct successor-closed subsets from highest label downward, independent
 # of the discovery catalogue's predecessor-ideal extension direction.
 natural=0;invalid=0;out=set();up=[0]*7
 def visit(v):
  nonlocal natural,invalid
  if v<0:
   natural+=1
   q=[sum(1<<u for u in range(7) if up[u]>>w&1) for w in range(7)]
   p=[0,12]+[m<<2 for m in q]
   for a in (7,8):p[a]|=1
   n=len(p)
   check(all(not(p[v]>>u&1) or not(p[u]&~p[v]) for v in range(n) for u in range(n)),'catalogue closure')
   degrees=[n-1-p[v].bit_count()-sum(bool(m>>v&1) for m in p) for v in range(n)]
   if max(degrees)>7:invalid+=1
   else:out.add(tuple(p))
   return
  available=((1<<7)-1)^((1<<(v+1))-1);mask=available
  while True:
   if all(not(mask>>u&1) or not(up[u]&~mask) for u in range(v+1,7)):
    up[v]=mask;visit(v-1)
   if not mask:break
   mask=(mask-1)&available
 visit(6)
 check(natural==96428 and invalid==13441 and len(out)==82987,'catalogue counts')
 return out


def main():
 start=time.time()
 raw_a=(ROOT/'analytic_cover_20261019.json.gz').read_bytes()
 raw_r=(ROOT/'residual_certificates_20261019.json.gz').read_bytes()
 check(hashlib.sha256(raw_a).hexdigest()==EXPECTED_ANALYTIC_SHA256,'analytic data hash')
 check(hashlib.sha256(raw_r).hexdigest()==EXPECTED_RESIDUAL_SHA256,'residual data hash')
 cover=json.loads(gzip.decompress(raw_a));records=json.loads(gzip.decompress(raw_r))
 check(cover['status']=='ANALYTIC_COVER_COMPLETE','incomplete analytic source')
 items=cover['cores'];catalogue=independent_catalogue()
 check(len(items)==len(catalogue) and {tuple(t['pred']) for t in items}==catalogue,'independent catalogue equality')
 expected_residual={i for i,t in enumerate(items) if len(t)==1}
 check(len(records)==len(expected_residual)==31889,'residual cardinality')
 check({r['core_index'] for r in records}==expected_residual,'residual coverage')
 analytic_balanced=analytic_cycle=0
 for item in items:
  p=tuple(item['pred'])
  if 'balanced_pair' in item:
   replay(p,{'analytic_balanced':item['balanced_pair']});analytic_balanced+=1
  elif 'cycle_arcs' in item:
   replay(p,{'analytic_cycle':item['cycle_arcs']});analytic_cycle+=1
 check(analytic_balanced==47298 and analytic_cycle==3800,'analytic counts')
 certs=[];directed=orders=short=0
 for index,r in enumerate(records):
  p=tuple(r['pred']);check(p==tuple(items[r['core_index']]['pred']),'core correspondence')
  replay(p,r);directed+=r.get('analytic_directed_orders',0);orders+=r.get('residual_orders',0);short+=r.get('short_span_orders',0);certs.extend(r.get('certificates',[]))
  if (index+1)%1000==0:print(json.dumps({'replayed_residual_cores':index+1,'generator_orders':orders,'seconds':time.time()-start}),flush=True)
 check(short==0,'unexpected short-span use: c3 common neighbors force span at least four')
 summary={'status':'ALL_82987_DEGREE_VALID_CORES_INDEPENDENTLY_EXACT_REPLAYED','scope':'Seven-thin mixed(2,2), c=3 maximum-span root morphology only. No c=4 claim.','natural_nonroot_posets':96428,'degree_invalid':13441,'degree_valid_cores':82987,'analytic_balanced_cores':analytic_balanced,'analytic_cycle_cores':analytic_cycle,'residual_cores':len(records),'directed_screen_orders':directed,'generator_orders':orders,'short_span_orders':short,'distinct_certificates':len(certs),'maximum_certificate_support':max(len(c['terms']) for c in certs),'maximum_integer_weight':max(t[2] for c in certs for t in c['terms']),'analytic_sha256':hashlib.sha256(raw_a).hexdigest(),'residual_sha256':hashlib.sha256(raw_r).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.time()-start}
 (ROOT/'independent_replay_receipt_20261019.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
