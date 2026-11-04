#!/usr/bin/env python3
import itertools,json,pathlib,collections,hashlib
import prepare as p
P=pathlib.Path(__file__).resolve().parent
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def independent_profile(r):
 n=len(r);H=[0]*(1<<n);H[0]=1
 for s in range(1,1<<n):
  vs=[i for i in range(n) if s>>i&1]
  H[s]=sum(all(r[a]>>b&1 for a,b in zip(t,t[1:])) for t in itertools.permutations(vs))
 D=[[0]*(n+1) for _ in H];D[0][0]=1
 for s in range(1,1<<n):
  t=s
  while t:
   for k in range(1,n+1):D[s][k]+=H[t]*D[s^t][k-1]
   t=(t-1)&s
 return D[-1][1:]
def symbolic_coeff(r):
 C=[0]*1024
 for path in itertools.permutations(range(7)):
  factors=[];valid=True
  for u,v in zip(path,path[1:]):
   if u<5 and v<5:
    if not(r[u]>>v&1):valid=False;break
   elif u>=5 and v>=5:
    if (u,v)!=(5,6):valid=False;break
   elif u<5:factors.append((u+5*(v==6),1))
   else:factors.append((v+5*(u==6),0))
  if not valid:continue
  poly={0:1}
  for j,positive in factors:
   new=collections.defaultdict(int)
   for m,c in poly.items():
    if positive:new[m|(1<<j)]+=c
    else:new[m]+=c;new[m|(1<<j)]-=c
   poly=new
  for m,c in poly.items():C[m]+=c
 return C
def main():
 records=json.loads((P/'five_child_records.json').read_text());summary=json.loads((P/'summary.json').read_text())
 rs=[p.rows(c) for c in range(1024)];qs=[independent_profile(r) for r in rs];canon=[p.canonical(r) for r in rs]
 check(len(records)==12 and sorted(set(canon))==[c['code'] for c in records],'isomorphism classes')
 for c in records:
  rep=c['code'];check(c['profile']==qs[rep],'profile');check(c['labeled_count']==canon.count(rep),'multiplicity')
  check(c['modules']==p.modules(rs[rep]),'modules');check(len(c['targets'])==1<<len(c['modules']),'family completeness')
  for f,t in enumerate(c['targets']):
   fam=[m for j,m in enumerate(c['modules']) if f>>j&1];code=p.toggle(rep,fam)
   check(t==dict(family=fam,code=code,canonical_code=canon[code],profile=qs[code]),'family exact record')
   check(all(a<=b for a,b in zip(qs[rep],qs[code])),'profile monotonicity')
  for modulus in (8,16):check(c['normalized_mod'+str(modulus)]==[v*pow(qs[rep][0],-1,modulus)%modulus for v in qs[rep]],'normalization')
 controls=json.loads((P/'named_controls.json').read_text());counts=collections.Counter()
 for rec,ctl in zip(records,controls):
  for side,rows in [('source',rec['rows']),('converse',p.rows(rec['code']^1023))]:
   check(ctl[side]['rows']==rows,'control rows');co=symbolic_coeff(rows);check(co==ctl[side]['coefficients'],'complete symbolic coefficients')
   vals=co[:]
   for j in range(10):
    for m in range(1024):
     if m>>j&1:vals[m]+=vals[m^(1<<j)]
   check(vals==ctl[side]['values'],'complete basis values')
   for m,c in enumerate(co):
    if c:counts[m.bit_count()]+=1
 check(max(counts)==4,'degree four');check(len(controls)==12,'control classes')
 receipt=dict(status='PASS',method='Independent ordered subset convolution profiles and symbolic seven-vertex permutation polynomials',summary=summary,control_vectors=24,all_named_coefficients=24576,nonzero_coefficients_by_degree=dict(counts),sha256={name:hashlib.sha256((P/name).read_bytes()).hexdigest() for name in ['prepare.py','verify.py','five_child_records.json','named_controls.json','summary.json','control_summary.json','PROOF.md']})
 (P/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
