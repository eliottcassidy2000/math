#!/usr/bin/env python3
"""Independent exact analytic-census and sparse-probe verification; stdlib only."""
from pathlib import Path
from itertools import combinations
import gzip,hashlib,importlib.util,json,time
ROOT=Path(__file__).resolve().parent
A_SHA='0ee6b7a5b6ae756310ec800474908e5e8e59ec0a2fb266910312dc2279cc7d52'
P_SHA='aa11834bdba584e3a39c27aec1cedc542738af4cfbf5e9bc1c45aa1c590f9aa5'
s=importlib.util.spec_from_file_location('exact_cone_replay',ROOT/'exact_cone_replay_20261019.py');R=importlib.util.module_from_spec(s);s.loader.exec_module(R)
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def decode(code):
 check(isinstance(code,int) and 0<=code<1<<28,'invalid relation code')
 q=[(code>>(v*(v-1)//2))&((1<<v)-1) for v in range(8)];p=[0,12]+[m<<2 for m in q];p[8]|=1;p[9]|=1
 check(all(not(p[v]>>u&1) or not(p[u]&~p[v]) for u in range(10) for v in range(10)),'nontransitive core')
 return tuple(p)
def code(p):return sum((p[v+2]>>2)<<(v*(v-1)//2) for v in range(8))
def structure(p):
 D=[{u for u in range(10) if p[v]>>u&1} for v in range(10)];U=[{v for v in range(10) if u in D[v]} for u in range(10)]
 degree=[9-len(D[v])-len(U[v]) for v in range(10)];check(max(degree)<=7,'degree invalid')
 seed={0,1}|{v for v in range(10) if degree[v]==7}
 alpha=[0 if (D[v]|{v})&seed else 7-degree[v] for v in range(10)];beta=[0 if (U[v]|{v})&seed else 7-degree[v] for v in range(10)]
 def bound(u,v):return len(D[v]-D[u])+len(U[u]-U[v])+alpha[u]+beta[v]
 return D,U,bound
BASE={(0,1)}|{(0,z) for z in range(4,8)}|{(z,1) for z in range(4,8)}
def verify_item(p,item):
 D,U,bound=structure(p);kind=item[1]
 if kind==0:
  a,b=item[2];check(a!=b and a not in D[b] and b not in D[a],'balanced pair comparable')
  check(bound(a,b)<=1 and bound(b,a)<=1,'bad balanced witness')
 elif kind==1:
  cycle=item[2];check(len(cycle)>=3 and cycle[0]==cycle[-1],'bad cycle form')
  for a,b in zip(cycle,cycle[1:]):
   check(a!=b and (a in D[b] or (a,b) in BASE or (a not in D[b] and b not in D[a] and bound(b,a)<=1)),'unforced cycle arc')
 else:
  check(kind==2,'unknown class');reach=list(p)
  for a,b in BASE:reach[b]|=1<<a
  for a,b in combinations(range(10),2):
   if a not in D[b] and b not in D[a]:
    if bound(a,b)<=1:reach[a]|=1<<b
    if bound(b,a)<=1:reach[b]|=1<<a
  for k in range(10):
   for v in range(10):
    if reach[v]>>k&1:reach[v]|=reach[k]
  check(not any(reach[v]>>v&1 for v in range(10)),'residual has analytic cycle')
def verify_universe(codes):
 up=[0]*8;natural=valid=0
 def visit(v):
  nonlocal natural,valid
  if v<0:
   natural+=1;comp=[up[w].bit_count()+sum(bool(m>>w&1) for m in up) for w in range(8)]
   if not all(comp[w]>=(2 if 2<=w<6 else 1) for w in range(8)):return
   valid+=1;k=0
   for a,m in enumerate(up):
    for b in range(a+1,8):
     if m>>b&1:k|=1<<(b*(b-1)//2+a)
   check(k in codes,'independent catalogue missing core');codes.remove(k);return
  available=((1<<8)-1)^((1<<(v+1))-1);mask=available
  while True:
   if all(not(mask>>u&1) or not(up[u]&~mask) for u in range(v+1,8)):
    up[v]=mask;visit(v-1)
   if not mask:break
   mask=(mask-1)&available
 visit(7)
 check(natural==2800472 and valid==1165395 and not codes,'universe cardinality/equality')
 return natural,valid

def main():
 start=time.time();a=ROOT/'analytic_classification_20261019.jsonl.gz';p=ROOT/'exact_probes_20261019.json.gz'
 check(hashlib.sha256(a.read_bytes()).hexdigest()==A_SHA and hashlib.sha256(p.read_bytes()).hexdigest()==P_SHA,'data hash')
 classes=[0,0,0];codes=set();residual=[]
 with gzip.open(a,'rt') as f:
  for i,line in enumerate(f):
   item=json.loads(line);k=item[0];check(k not in codes,'duplicate core');codes.add(k);pred=decode(k);verify_item(pred,item);classes[item[1]]+=1
   if item[1]==2:residual.append([k,i])
   if (i+1)%100000==0:print(json.dumps({'verified_analytic_records':i+1,'classes':classes,'seconds':time.time()-start}),flush=True)
 check(classes==[733492,216636,215267],'analytic class counts')
 inv=json.loads(gzip.decompress((ROOT/'residual_inventory_20261019.json.gz').read_bytes()));check(inv['classification_sha256']==A_SHA and inv['codes']==residual,'residual inventory exact equality')
 natural,valid=verify_universe(codes)
 probes=json.loads(gzip.decompress(p.read_bytes()));expected=set(range(16))|{i-1 for i in range(10000,215268,10000)}
 check(len(probes)==37 and {r['residual_index'] for r in probes}==expected,'probe schedule')
 orders=lp=0
 for r in probes:
  i=r['residual_index'];check(tuple(r['pred'])==decode(residual[i][0]) and r['valid_index']==residual[i][1],'probe core correspondence');R.replay(tuple(r['pred']),r);check(r.get('short_span_orders',0)==0,'short span');orders+=r.get('residual_orders',0);lp+=r.get('LP_calls',0)
 check(orders==366 and lp==0,'probe totals')
 out={'status':'COMPLETE_ANALYTIC_CENSUS_AND_SPARSE_PROBES_INDEPENDENTLY_REPLAYED','scope':'Complete degree/analytic catalogue and 37 exact residual probes only. The remaining residual cores are not certified by this packet.','natural_posets':natural,'degree_invalid':natural-valid,'degree_valid_cores':valid,'analytic_balanced':classes[0],'analytic_cycles':classes[1],'residual_cores':classes[2],'exact_probes':len(probes),'probe_generator_orders':orders,'probe_LP_calls':lp,'classification_sha256':A_SHA,'probe_sha256':P_SHA,'inventory_sha256':hashlib.sha256((ROOT/'residual_inventory_20261019.json.gz').read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.time()-start}
 (ROOT/'independent_census_replay_receipt_20261019.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
