#!/usr/bin/env python3
"""Independent literal-projection and majority audit using original electrical labels."""
import argparse,json,hashlib,itertools
from pathlib import Path
P=argparse.ArgumentParser();P.add_argument('--repo',type=Path,required=True);P.add_argument('--packet',type=Path,required=True);P.add_argument('--out',type=Path,required=True);a=P.parse_args()
checks=0
def require(test,why):
 global checks
 checks+=1
 if not test:raise RuntimeError(why)
Epath=a.repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json';Cpath=a.packet/'candidates.json';E=json.loads(Epath.read_text())['order8'];X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4'];cases=json.loads(Cpath.read_text())['found']
expected={6:(2,2),8:(2,2),10:(8,64),12:(16,2048),14:(68,16384),16:(196,65536)}
provided={10:[0,34,249,32],12:[0,210,387,130],14:[0,646,652,644],16:[0,1806,1811,1794]}
reports=[]
for _,case in sorted(cases.items(),key=lambda z:int(z[0])):
 sources=case['plus']+case['minus'];n=len(sources);k=n//2;signs=[1]*k+[-1]*k;full=(1<<n)-1;states={4:[],128:[]}
 def balanced(mask):return (mask&((1<<k)-1)).bit_count()==(mask>>k).bit_count()
 for mask in range(1<<n):
  if not balanced(mask):continue
  selected=[j for j in range(n)if mask>>j&1]
  if any(sum(signs[j]*Y[t][sources[j]]for j in selected)%2 for t in range(4)):continue
  value=sum(signs[j]*X[sources[j]]for j in selected)
  for modulus in states:
   if value%modulus==0:states[modulus].append(mask)
 for modulus,Slist in states.items():
  S=set(Slist);require(0 in S and full in S and {full^z for z in S}==S,'constant endpoints and complement closure')
  pairprojections={}
  for i,j in itertools.combinations(range(n),2):
   projection={(z>>i&1,z>>j&1)for z in S};require(projection in ({(0,0),(1,1)},{(0,0),(0,1),(1,0),(1,1)}),'every binary projection is literal equality or unrestricted')
   pairprojections[i,j]=projection
  # Full literal pair membership test, without an equality-component shortcut.
  closure=[mask for mask in range(1<<n)if all((mask>>i&1,mask>>j&1)in projection for (i,j),projection in pairprojections.items())]
  wanted=expected[n]if modulus==4 else(2,2)
  require((len(S),len(closure))==wanted,'exact full relation and literal pairwise closure counts')
  bad=None;balanced_bad=None
  for left,right in itertools.combinations(Slist,2):
   majority=left&right
   if majority not in S:
    if bad is None:bad=[0,left,right,majority]
    if balanced(majority):balanced_bad=[0,left,right,majority];break
  require((bad is None)==(set(closure)==S),'exact Boolean closure comparison')
  if modulus==4 and n>=10:
   z,left,right,majority=provided[n]
   require(z in S and left in S and right in S,'saved hostile three assignments accepted from original electrical data')
   require(((z&left)|(left&right)|(right&z))==majority and majority not in S,'saved coordinatewise majority genuinely rejected')
  else:require(S=={0,full},'full higher bank or small H8 case is exactly the equality relation')
  if balanced_bad:
   z,left,right,majority=balanced_bad;require(all(v in S for v in (z,left,right))and balanced(majority)and majority not in S,'optional majority hostile retains target margins but fails current')
  reports.append({'slots':n,'harmonic_modulus':8 if modulus==4 else 256,'accepted_states':len(S),'pairwise_closure_states':len(closure),'equality_projection_pairs':sum(len(v)==2 for v in pairprojections.values()),'saved_majority_hostile':provided.get(n)if modulus==4 else None,'balanced_current_majority_hostile':balanced_bad})
# Direct truth-table verification of the local majority fact behind the general
# proof. Conjunction and existential projection then preserve this operation.
for signs in itertools.product((0,1),repeat=2):
 clause=[(x,y)for x,y in itertools.product((0,1),repeat=2)if (x^signs[0])or(y^signs[1])]
 for triple in itertools.product(clause,repeat=3):
  m=tuple(int(sum(t[i]for t in triple)>=2)for i in range(2));require(m in clause,'every signed binary clause is majority closed')
report={'status':'PASS','checks':checks,'reports':reports,'proof_scope':'The visible literal-bit relation cannot be an existential projection of any 2-CNF relation when a majority hostile exists. No claim excludes a nonlinear change of visible representation.','source_hashes':{'electrical':hashlib.sha256(Epath.read_bytes()).hexdigest(),'candidates':hashlib.sha256(Cpath.read_bytes()).hexdigest()}}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
