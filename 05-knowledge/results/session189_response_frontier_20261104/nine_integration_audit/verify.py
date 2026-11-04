#!/usr/bin/env python3
"""Independent frontier-nine restoration review; no producer module imports."""
import json,hashlib,itertools,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
P=ROOT.parent/'session189-frontier-nine'
if not P.is_dir():P=ROOT.parent/'frontier_nine'

def need(x,*msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def valid(Q):return all(not(r>>i&1) and all(((r>>j)&1)+((Q[j]>>i)&1)==1 for j in range(len(Q)) if i!=j) for i,r in enumerate(Q))
def modules(Q):
 n=len(Q);out=[]
 for mask in range(1,1<<n):
  if all((r&mask) in (0,mask) for i,r in enumerate(Q) if not mask>>i&1):out.append(mask)
 return out

def tree(Q,U):
 ms=modules(Q);strong=[a for a in ms if not any((a&b) and (a&~b) and (b&~a) for b in ms)]
 records=[];changed=[]
 for M in strong:
  if M.bit_count()<2:continue
  proper=[a for a in strong if a!=M and a&M==a];children=sorted(a for a in proper if not any(a!=b and a&b==a for b in proper))
  reps=[(a&-a).bit_length()-1 for a in children]
  sides=[]
  for X in (Q,U):sides.append([sum(1<<j for j,v in enumerate(reps) if X[u]>>v&1) for u in reps])
  records.append(dict(node=M,children=children,source_quotient=sides[0],target_quotient=sides[1],changed=sides[0]!=sides[1]))
  if sides[0]!=sides[1]:changed.append(M)
 highest=[a for a in changed if not any(a!=b and a&b==a for b in changed)]
 return ms,records,highest

def all_profiles(Q):
 # A separate backward endpoint/cut polynomial recurrence computes all q_h
 # for all subsets simultaneously, rather than the producer's set partitions.
 n=len(Q);need(n<=13,'profile order');width=64;bound=(1<<(n*width))-1
 incoming=[sum(1<<i for i,r in enumerate(Q) if r>>j&1) for j in range(n)]
 dp=[[0]*n for _ in range(1<<n)];total=[0]*(1<<n)
 for mask in range(1,1<<n):
  ends=mask
  while ends:
   bit=ends&-ends;ends-=bit;j=bit.bit_length()-1;prior=mask^bit
   if not prior:v=1
   else:
    v=total[prior]<<width;allowed=prior&incoming[j]
    while allowed:
     b=allowed&-allowed;allowed-=b;v+=dp[prior][b.bit_length()-1]
    v&=bound
   dp[mask][j]=v;total[mask]+=v
  need(total[mask]>>(n*width)==0,'packed carry bound')
 out=[[1]]
 for mask in range(1,1<<n):out.append([0]+[(total[mask]>>(64*h))&((1<<64)-1) for h in range(mask.bit_count())])
 return out

def extension(Q,a,c):
 n=len(Q);full=(1<<n)-1
 return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(Q)]+[(1<<(n+1))|(full^a),full^c]

def main():
 dependencies=json.loads((P/'DEPENDENCIES.json').read_text())
 for d in dependencies:need(sha(P/d['path'])==d['sha256'],'dependency pin',d['path'])
 cases=json.loads((P/'controls.json').read_text());banks=json.loads((P/'named_banks.json').read_text());profile_records=[];physical=[];metadata=[];named_count=0;actual=[]
 for case in cases:
  Q=case['source'];U=case['target'];V=case['intermediate'];need(all(valid(X) for X in (Q,U,V)),'physical tournament rows')
  ms,records,highest=tree(Q,U)
  need(ms==case['original_modules'] and records==case['strong_node_records'] and highest==case['highest_changed_nodes'],'complete strong structure')
  family=case['family'];need(all(M in ms for M in family),'original modules')
  toggled=list(Q)
  for i in range(len(Q)):
   for j in range(i+1,len(Q)):
    if sum((M>>i&1) and (M>>j&1) for M in family)%2:toggled[i]^=1<<j;toggled[j]^=1<<i
  need(toggled==U,'actual simultaneous XOR')
  need(set(highest)=={case['core'],case['other']},'two highest nodes')
  core_record=next(r for r in records if r['node']==case['core'])
  expected_core=([4,17,42,131,45,75,31,119] if case['name'].startswith('ordered') else [16,5,41,67,46,139,55,95])
  need(core_record['source_quotient']==expected_core and len(core_record['children'])==8,'specified residual prime-eight quotient')
  need(set(modules(expected_core))=={*(1<<i for i in range(8)),255},'residual primality')
  for i in range(len(Q)):
   for j in range(i+1,len(Q)):
    inside=bool((case['core']>>i&1) and (case['core']>>j&1))
    need(((V[i]>>j)&1)==(((U if inside else Q)[i]>>j)&1),'exact intermediate core-only replacement')
  a,c=case['query'];core=case['core'];other=case['other'];full=(1<<len(Q))-1
  for i in range(len(Q)):
   if not(core>>i&1):
    want=bool(Q[i]&core);need((Q[i]&core) in (0,core),'uniform root exterior')
    need(bool(a>>i&1)==bool(c>>i&1)==want,'actual exterior port bits')
  need((a&other) in (0,other) and (c&other) in (0,other),'restored node remains module')
  prof=[all_profiles(X) for X in (Q,V,U)]
  need(prof[0]==prof[1]==prof[2],'every induced subset complete profile equality')
  need(prof[0]==case['hereditary_profiles_source'] and prof[2]==case['hereditary_profiles_target'],'complete producer profile arrays')
  profile_records.append(dict(name=case['name'],source_profiles=prof[0],intermediate_profiles=prof[1],target_profiles=prof[2]))
  for role,X in [('source',Q),('intermediate',V),('target',U)]:
   bank=next(b for b in banks if b['case']==case['name'] and b['role']==role)
   need(bank['rows']==X,'named bank source rows');n=len(X);masks=[x for x in range(1<<n) if x.bit_count()<=2];need(masks==bank['physical_masks'],'complete named mask basis')
   pos={m:i for i,m in enumerate(masks)}
   subsets={m:[x for x in masks if x&m==x] for m in masks}
   for i,m in enumerate(masks):
    for j,p in enumerate(masks):
     coeff=sum((-1)**((m^x).bit_count()+(p^y).bit_count())*bank['values'][pos[x]][pos[y]] for x in subsets[m] for y in subsets[p])
     need(coeff==bank['coefficients'][i][j],'complete named Mobius coefficients')
     named_count+=1
   probes=[(a,c),(0,0),(full,full),(full,0),(0,full),(sum(1<<i for i in range(n) if i%2),sum(1<<i for i in range(n) if i%3)),(full^a,full^c),(c,a)]
   for z,(aa,cc) in enumerate(dict.fromkeys(probes)):
    value=sum(bank['coefficients'][i][j] for i,m in enumerate(masks) if m&aa==m for j,p in enumerate(masks) if p&cc==p)
    physical.append(extension(X,aa,cc));metadata.append(dict(case=case['name'],role=role,a=aa,c=cc,rows=physical[-1],polynomial_value=value))
  actual.append(dict(name=case['name'],masks=[a,c],values=case['actual_response_values']))
 with tempfile.TemporaryDirectory() as td:
  exe=str(Path(td)/'physical')
  subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'physical_backward.cpp'),'-o',exe],check=True)
  inp=str(len(physical))+'\n'+''.join(str(len(X))+'\n'+' '.join(map(str,X))+'\n' for X in physical)
  vals=list(map(int,subprocess.check_output([exe],input=inp.encode()).split()))
 need(len(vals)==len(metadata),'all independent physical values')
 for rec,value in zip(metadata,vals):need(value==rec['polynomial_value'],'literal versus complete polynomial',rec);rec['literal_value']=value
 for rec in actual:
  values=[next(r['literal_value'] for r in metadata if r['case']==rec['name'] and r['role']==role and [r['a'],r['c']]==rec['masks']) for role in ('source','intermediate','target')]
  need(values==rec['values'] and values[1]==values[2] and (values[0]-values[1])%4==2,'exact same-query restoration')
 (ROOT/'profile_recount.json').write_text(json.dumps(profile_records,sort_keys=True,separators=(',',':'))+'\n')
 (ROOT/'physical_controls.json').write_text(json.dumps(metadata,sort_keys=True,indent=2)+'\n')
 result=dict(status='PASS',scope='Independent restoration controls and mathematical review; complete named Mobius data checked, selected physical queries recounted. Prime-eight completeness is inherited from its source-pinned independent audit, not re-censused here.',cases=actual,full_named_coefficient_positions=named_count,independent_physical_queries=len(vals),maximum_physical_order=max(map(len,physical)),complete_profile_arrays=len(profile_records)*3,all_induced_profiles=sum(len(r['source_profiles'])*3 for r in profile_records),dependency_pins={d['path']:d['sha256'] for d in dependencies},producer_pins={n:sha(P/n) for n in ('PROOF.md','verify.py','controls.json','named_banks.json')},audit_data_pins={n:sha(ROOT/n) for n in ('profile_recount.json','physical_controls.json')})
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
