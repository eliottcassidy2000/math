#!/usr/bin/env python3
"""Exact scale phases, linear guarded parent banks and a timing obstruction."""
from pathlib import Path
from hashlib import sha256
from importlib.util import spec_from_file_location,module_from_spec
from collections import Counter
import json,sys
HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('supply_flow',HERE.parent/'annulus_residual_flow_20261020k/check.py')
f=module_from_spec(spec);sys.modules[spec.name]=f;spec.loader.exec_module(f)
n=f.native;o=f.old.owned
CHECKS=0
def ck(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise ValueError(msg)
def enc(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def integer(x,lo=0):
 if type(x) is not int or x<lo:raise ValueError('literal integer domain')
 return x
def primitive(u):
 integer(u,1)
 if not n.source(u):raise ValueError('primitive source')
 return u
def phase(u,a):
 primitive(u);integer(a,1)
 if a not in (1,2):raise ValueError('native first branch')
 js=[j for j in range(3) if (2**a*(3*u+1)*4**j-(2**a+3))%9==0]
 ck(len(js)==1,'unique inverse exponent modulo three')
 return js[0]
def inverse(u,a,j):
 primitive(u);integer(a,1);integer(j)
 if a not in (1,2):raise ValueError('native first branch')
 z=2**a*(3*u+1)*4**j-(2**a+3)
 if z%9:raise ValueError('inverse integrality')
 s=z//9
 if not n.source(s) or (2 if s%8==1 else 1)!=a:raise ValueError('native source branch')
 ck(n.edge(s)==(u,a+2*j),'literal primitive edge and braid clock')
 return s
def annulus_child(N,u,a):
 integer(N,1);primitive(u);integer(a,1)
 if u>=4*N:raise ValueError('unchanged old parent')
 j=phase(u,a);s=inverse(u,a,j)
 while s<4*N:s=64*s+7*(2**a+3);j+=3
 ck(4*N<=s<256*N and inverse(u,a,j)==s,'exact unique annulus child')
 return s,j
def child_residue(u,a,j,k,modulus):
 primitive(u);integer(a,1);integer(j);integer(k);integer(modulus,1)
 if a not in (1,2) or j%3!=phase(u,a):raise ValueError('authenticated native phase')
 z=(2**a*(3*u+1)*pow(4,j+3*k,9*modulus)-(2**a+3))%(9*modulus)
 ck(z%9==0,'full denominator modular division')
 return z//9
def rotate(v,k=1):
 for _ in range(k%3):v=(-v[1],v[0]-v[1])
 return v
def annulus_charge(N,s):
 integer(N,1);primitive(s)
 if not 4*N<=s<256*N:raise ValueError('new annulus source')
 b=next(i for i in range(6) if s<(8<<i)*N);r=s%3
 return ((0,0),n.U[r],tuple(-z for z in n.L[(r+2)%3]),n.U[r],n.L[r],n.U[r])[b]
def progression_count(lo,hi,residues):
 integer(lo);integer(hi)
 if hi<lo or type(residues) is not tuple or len(set(residues))!=len(residues):raise ValueError('literal progression bank')
 for r in residues:
  integer(r)
  if r>=72:raise ValueError('canonical residue')
 return sum((hi-1-r)//72-(lo-1-r)//72 for r in residues)
def linear_spec(N):
 integer(N,1);cut=(3*N+1)//2
 return (('A',cut,2*N,(9,27,63),1,1),('B',N,cut,(17,35,71),2,2))
def authenticate_pair(N,u,s,a,j,routes,family):
 integer(N,1);primitive(u);primitive(s);integer(a,1);integer(j)
 if type(family) is not str:raise ValueError('literal family label')
 ck(u<4*N<=s<256*N and s not in routes.V,'old/new domain and first-arrival guard')
 ck(inverse(u,a,j)==s,'full unchanged source formula')
 ut,uc,ud=routes.get(u);st,sc,sd=routes.get(s)
 ck((st,sc,sd)==(ut,uc+a+2*j,ud+1),'actual parent and child terminal/clock seam')
 rows=[]
 for x,t,c,d in ((u,ut,uc,ud),(s,st,sc,sd)):
  charge=n.charge(64*N,x);ck(charge==n.independent_charge(64*N,x),'independent charge')
  row=o.Row(x,x%15,x%35,t,c,d,charge);o.audit_row(N,row,routes);rows.append(row)
 value=n.beta(o.signed(rows[0]),o.signed(rows[1]))
 ck(value==(-1)**a*n.beta(rows[0].unsigned,rows[1].unsigned),'tail clock cancels only in pair sign')
 return {'family':family,'parent':u,'source':s,'branch':a,'strips':j,'terminal':ut,
         'parent_clock':uc,'child_clock':sc,'parent_depth':ud,'child_depth':sd,
         'parent_row':o.payload(rows[0]),'source_row':o.payload(rows[1]),'beta':value}
def linear_pairs(N,routes):
 out=[]
 for label,lo,hi,rs,j,weight in linear_spec(N):
  us=[u for u in range(lo,hi) if u%72 in rs]
  ck(len(us)==progression_count(lo,hi,rs),'exact full native parent progression count')
  for u in us:
   s=inverse(u,2,j);p=authenticate_pair(N,u,s,2,j,routes,label)
   ck(p['beta']==-3*weight,'uniform negative parent-child interaction')
   out.append(p)
 return tuple(out)
def slots(rows,D):
 integer(D,1);ans=[];owners=[];seen=set()
 if type(rows) is not tuple:raise ValueError('immutable physical pair bank')
 for i,p in enumerate(rows):
  for key in ('parent','source','child_depth','beta'):
   if type(p[key]) is not int:raise ValueError('literal receipt integer')
  address=(p['parent'],p['source'])
  if address in seen:raise ValueError('physical pair address reused')
  seen.add(address)
  if p['beta']>=0 or p['beta']%3:raise ValueError('negative three-unit owner')
  for t in range(-p['beta']//3):
   ans.append(f.Interval(len(ans),p['child_depth'],D));owners.append([i,t,*address])
 return tuple(ans),owners
def reject(fn):
 try:fn()
 except (ValueError,TypeError):return True
 return False
def small_controls():
 exact=0
 for N in range(1,17):
  expected={}
  for u in range(1,4*N,2):
   if n.source(u):
    for a in (1,2):
     s,j=annulus_child(N,u,a);expected[s]=(u,a+2*j)
  actual={}
  for s in range(4*N+1,256*N,2):
   if not n.source(s):continue
   ck(annulus_charge(N,s)==n.charge(64*N,s),'six-band formula versus actual B-chain')
   u,c=n.edge(s)
   if u<4*N:actual[s]=(u,c)
  ck(actual==expected,'complete one-step inverse annulus fibre')
  exact+=len(actual)
 for N in range(1,289):
  for label,lo,hi,rs,j,w in linear_spec(N):
   us=[u for u in range(lo,hi) if u%72 in rs]
   ck(len(us)==progression_count(lo,hi,rs),'floor count independently enumerated')
   if N%144==0:ck(len(us)==N//48,'exact linear count at aligned scale')
   for u in us:
    s=inverse(u,2,j)
    ck(u<4*N<=s<256*N,'growing family exact source bands')
    ck(n.beta(n.charge(N,u),n.charge(64*N,s))==-3*w,'all small family charge values')
    ck((s+35* (384 if label=='A' else 1536))%105==s%105,'literal progression joint-register period35')
 for bad in (True,1.0):
  ck(reject(lambda:annulus_child(bad,1,1)),'nonliteral scale rejected')
  ck(reject(lambda:child_residue(1,1,1,bad,105)),'nonliteral phase index rejected')
 ck(reject(lambda:annulus_child(1,7,1)),'new parent cannot masquerade as old')
 return {'complete_small_scales':16,'exact_inverse_rows':exact,'family_scales':288}
def fixed_phases(V):
 N=16384;rows=[];hist=[Counter() for _ in range(3)];nonzero=[]
 ck(max(V)<4*N,'all627 frontier parents genuinely old')
 for u in sorted(V):
  for a in (1,2):
   s,j=annulus_child(N,u,a);old=n.charge(N,u);base=n.charge(64*N,s);values=[]
   for k in range(6):
    scale=64**k;sk=scale*s+(2**a+3)*(scale-1)//9
    ck(annulus_child(scale*N,u,a)==(sk,j+3*k),'exact scale recurrence of unique inverse slot')
    vk=n.charge(64*scale*N,sk)
    ck(vk==rotate(base,7*(2**a+3)*k),'all scale phase charge law')
    ck(n.charge(scale*N,u)==old,'fixed old source charge scale invariance')
    value=(-1)**a*n.beta(old,vk);values.append(value)
    for mod in (15,35,105):ck(child_residue(u,a,j,k,mod)==sk%mod,'full denominator prime-register reader')
   ck(values[:3]==values[3:] and sum(values[:3])==0,'three-phase interaction sum is zero')
   ck(child_residue(u,a,j,6,105)==s%105,'joint-register period divides six')
   for k,v in enumerate(values[:3]):hist[k][v]+=1
   if any(values):nonzero.append([u,a,s,j,values[:3]])
   rows.append([u,a,s,j,*values[:3]])
 ck(len(rows)==1254,'fixed parent bank has exactly bounded1254 slots')
 return {'N':N,'parent_count':len(V),'slot_count':len(rows),'rows':rows,
         'phase_histograms':[sorted(h.items()) for h in hist],'first_nonzero_hostile':nonzero[0]}
def finite_completion(routes,prior):
 N=8192;D=257;linear=linear_pairs(N,routes);used={(p['parent'],p['source']) for p in linear}
 jobs=tuple(f.Interval(b['id'],b['start'],b['end']) for b in prior['removed_blocks'])
 units,unit_owners=slots(linear,D);cap=f.capacity(units,D)
 selected=f.propose(jobs,cap);lower,cuts=f.dual(jobs,cap)
 cert=f.audit_solution(jobs,cap,selected,cuts)
 ck(len(linear)==342 and len(units)==513 and len(selected)==121 and lower==38 and cuts==(1,),
    'exact linear supply and initial-time obstruction')
 supplement=[]
 for u in sorted(routes.V):
  if u>=4*N:continue
  for a in (1,2):
   s,j=annulus_child(N,u,a)
   if s in routes.V or (u,s) in used:continue
   value=(-1)**a*n.beta(n.charge(N,u),n.charge(64*N,s))
   if value!=-3:continue
   p=authenticate_pair(N,u,s,a,j,routes,'early_frontier')
   ck(p['parent_depth']==0 and p['child_depth']==1,'actual early frontier owner')
   supplement.append(p);used.add((u,s))
   if len(supplement)==lower:break
  if len(supplement)==lower:break
 ck(len(supplement)==38,'bounded early-time supplement')
 allrows=linear+tuple(supplement);allunits,allowners=slots(allrows,D);allcap=f.capacity(allunits,D)
 final=f.audit_solution(jobs,allcap,tuple(sorted(x.identity for x in jobs)),())
 segments=f.allocate(jobs,allunits,D);f.audit_allocation(jobs,allunits,D,segments)
 ledger=[]
 for d in range(D):
  C=sum(p['beta'] for p in allrows if p['child_depth']<=d)
  U=prior['residual_debit'][d];J=prior['original_bank_debit'][d]
  ck(C==-3*allcap[d] and U+C<=0 and J+C<=0,'complete source-owned scalar bank payment')
  # Coefficients of an arbitrary unchanged R: (R-C)+(J+C)=R+J.
  ck(-C+(J+C)==J and -C+(U+C)==U,'exact residual transfer, not energy gain')
  ledger.append([d,cap[d],allcap[d],U,C,J+C])
 ck(reject(lambda:slots(allrows+(allrows[0],),D)),'duplicate physical pair is rejected before slot splitting')
 return {'N':N,'linear_rows':linear,'linear_unit_owners':unit_owners,'linear_certificate':cert,
         'frontier_supplement':supplement,'final_certificate':final,'unit_owners':allowners,
         'ownership_segments':segments,'ledger_columns':['depth','linear_units','total_units','U','C','J_plus_C'],
         'ledger':ledger,'latest_linear_birth':max(p['child_depth'] for p in linear),
         'previous_completed_count':len(prior['global_certificate']['selected']),
         'remaining_count':len(jobs),'scope':'finite sibling-bank completion; full8192 R and energy not estimated'}
def run():
 controls=small_controls()
 fp=HERE.parent/'session18_native_annulus/frontier.json';fb=fp.read_bytes().replace(b'\r\n',b'\n')
 ck(sha256(fb).hexdigest()=='391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f','frozen frontier')
 front=json.loads(fb);routes=n.Routes(front['frontier'],dict(front['frontier_ranks']))
 bp=HERE.parent/'annulus_residual_flow_20261020k/receipt.json';bb=bp.read_bytes().replace(b'\r\n',b'\n')
 ck(sha256(bb).hexdigest()=='35974096cbf893fb209ce5c47d915a9837b4a048ac7579b1d0806b46e106ed4d','frozen interval bank')
 phase_receipt=fixed_phases(routes.V);finite=finite_completion(routes,json.loads(bb)['native_cases'][1])
 return {'scope':'bounded fixed-parent obstruction; guarded linear supply; source-owned finite completion; uniform153N OPEN',
         'controls':controls,'fixed_frontier_phases':phase_receipt,'case8192':finite,
         'frontier_sha256':sha256(fb).hexdigest(),'prior_sha256':sha256(bb).hexdigest(),
         'route_audit':routes.audit(),'checks':CHECKS}
if __name__=='__main__':
 result=run();data=enc(result);path=HERE/'receipt.json'
 if '--write' in sys.argv:path.write_bytes(data)
 else:ck(path.read_bytes().replace(b'\r\n',b'\n')==data,'saved exact receipt')
 print(json.dumps({'checks':CHECKS,'sha256':sha256(data).hexdigest(),
                   'linear':len(result['case8192']['linear_rows']),'supplement':len(result['case8192']['frontier_supplement'])}))
