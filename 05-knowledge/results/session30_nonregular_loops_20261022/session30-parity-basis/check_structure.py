import json,numpy as np
from check import contract,basis,cores,reference,ROOT
from parity_engine import evaluate
full=json.loads((ROOT/'parity-kernels.json').read_text())
bridges=0;loops=0
for c in full:
 edges=c['edges'];V=c['V']
 for i,(a,b) in enumerate(edges):
  palette={t['indices'][i] for t in c['terms']}
  if a==b:
   loops+=1;assert palette<={0,1}
  seen={a}
  while True:
   new=seen|{v for j,(u,v) in enumerate(edges) if i!=j and u in seen}|{u for j,(u,v) in enumerate(edges) if i!=j and v in seen}
   if new==seen:break
   seen=new
  if b not in seen:
   bridges+=1;assert palette<={0,3}
   if c['weight']==1:assert palette=={0}
rng=np.random.default_rng(3030);worst=0;checks=0
for n in [5,8,12]:
 A=rng.normal(size=(n,n));M=A-A.T;M*=.35/np.linalg.norm(M,2)
 for z in [1.,.73+.17j]:
  B=basis(M,z)
  for w in [1,2]:
   refs=reference(M,w,z,True);dense=evaluate(M,w,z,True);reverse=evaluate(-M,w,z,True)
   for c,(_,ref),(_,d),(_,rev) in zip([c for c in cores if c['weight']==w],refs,dense,reverse):
    s=contract(c,B);errs=[abs(s-ref),abs(d-ref),abs(rev-ref)];worst=max(worst,*errs);checks+=1
    assert max(errs)<1e-10
out={'status':'PASS','weighted_skew_core_checks':checks,'each_compares':['sparse quotient','compressed dense','global matrix reversal'],'maximum_absolute_error':float(worst),'loop_occurrences':loops,'bridge_occurrences':bridges,'bridge_palette':'G/P only across the 31 audited tables; G only at weight one'}
(ROOT/'structure-results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
