"""Exact upper boundary for primitive one-slot-per-source circuits in the fixed 20-option supply."""
from itertools import combinations
from functools import reduce
from operator import xor
from pathlib import Path
import json,hashlib,sys,collections
H=Path(__file__).resolve().parent;D=json.loads((H/'pool.json').read_text());L={int(s):tuple(v) for s,v in D['labels'].items()};ids=sorted(L)
def total(S):return reduce(xor,(L[s][0] for s in S),0),sum(L[s][1] for s in S)%128
def subwitness(A,B):
 for k in range(1,len(A)//2+1):
  bag={total(C):C for C in combinations(A,k)}
  for C in combinations(B,k):
   if total(C) in bag:return [list(bag[total(C)]),list(C)]
 return None
p={s:L[s][0]|((L[s][1]%2)<<4) for s in ids};allpar=reduce(xor,p.values(),0);rows=[];tests=0;primitive=[]
for omit in combinations(ids,2):
 if p[omit[0]]^p[omit[1]]!=allpar:continue
 S=[s for s in ids if s not in omit];anchor=S[0];matches=0;witnesses=[]
 for rest in combinations(S[1:],8):
  A=(anchor,)+rest;B=tuple(s for s in S if s not in A);tests+=1
  if total(A)!=total(B):continue
  matches+=1;w=subwitness(A,B)
  if w is None:primitive.append({'plus':A,'minus':B})
  else:witnesses.append({'plus':A,'minus':B,'proper':w})
 rows.append({'omitted':omit,'balanced_relations':matches,'proper_subrelations':witnesses})
if len(ids)!=20 or allpar!=6 or tests!=72930 or sum(r['balanced_relations'] for r in rows)!=1129 or primitive:raise RuntimeError('sharp sixteen boundary consequence')
result={'twenty_source_parity':allpar,'eighteen_support_partitions_tested':tests,'eighteen_relations':sum(r['balanced_relations'] for r in rows),'primitive_eighteen':primitive,'rows':rows}
certificate=H/'max-support-certificate.json'
if '--write' in sys.argv:certificate.write_text(json.dumps(result,indent=2)+'\n')
else:
 if json.loads(certificate.read_text())!=json.loads(json.dumps(result)):raise RuntimeError('frozen maximum-support certificate')
print(json.dumps({'status':'PASS','twenty_source_parity':allpar,'partitions_tested':tests,'eighteen_relations':result['eighteen_relations'],'primitive_eighteen_count':len(primitive),'proper_subrelation_half_sizes':dict(collections.Counter(len(w['proper'][0]) for r in rows for w in r['proper_subrelations'])),'certificate_sha256':hashlib.sha256(certificate.read_bytes()).hexdigest()},indent=2))
