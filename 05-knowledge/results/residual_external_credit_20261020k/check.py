#!/usr/bin/env python3
"""Targeted old/new interaction transfer. No physical energy improvement."""
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
from hashlib import sha256
import json,sys
HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('external_flow',HERE.parent/'annulus_residual_flow_20261020k/check.py')
f=module_from_spec(spec);sys.modules[spec.name]=f;spec.loader.exec_module(f)
n=f.native;o=f.old.owned
CHECKS=0
def ck(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise ValueError(msg)
def enc(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def inverse(u,a,j):
 for x,lo in ((u,1),(a,1),(j,0)):
  if type(x) is not int or x<lo:raise ValueError('literal inverse parameters')
 if a not in (1,2) or not n.source(u):raise ValueError('primitive parent and first branch')
 num=(2**a)*(3*u+1)*4**j-(2**a+3)
 if num%9:raise ValueError('native integrality')
 s=num//9
 if not n.source(s) or (2 if s%8==1 else 1)!=a:raise ValueError('native source branch')
 ck(n.edge(s)==(u,a+2*j),'native inverse formula')
 return s
def reject(fn):
 try:fn()
 except ValueError:return True
 return False
def lifted_residue(u,a,j,k,modulus):
 for x,lo in ((u,1),(a,1),(j,0),(k,0),(modulus,1)):
  if type(x) is not int or x<lo:raise ValueError('literal lifted residue inputs')
 if a not in (1,2) or not n.source(u):raise ValueError('native lifted parent')
 den=9;num=((2**a)*(3*u+1)*pow(4,j+9*k,den*modulus)-(2**a+3))%(den*modulus)
 if num%den:raise ValueError('lift integrality')
 return (num//den)%modulus
def run():
 N=4096;D=257
 bp=HERE.parent/'annulus_residual_flow_20261020k/receipt.json';bb=bp.read_bytes().replace(b'\r\n',b'\n')
 ck(sha256(bb).hexdigest()=='35974096cbf893fb209ce5c47d915a9837b4a048ac7579b1d0806b46e106ed4d','frozen flow source')
 base=json.loads(bb)['native_cases'][0];jobs=sorted(base['removed_blocks'],key=lambda b:(b['start'],b['id']))
 ck(len(jobs)==94 and min(b['start'] for b in jobs)>=1,'literal residual job universe')
 fb=(HERE.parent/'session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n',b'\n')
 ck(sha256(fb).hexdigest()=='391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f','frozen native frontier source')
 frontier=json.loads(fb)
 routes=n.Routes(frontier['frontier'],dict(frontier['frontier_ranks']))
 rows={};credits=[];examined=0
 # Only one-step native predecessors of already-grounded old frontier sources.
 # No full orbit search or broad route enumeration occurs in this discovery.
 for y in range(4*N+1,256*N,2):
  if not n.source(y):continue
  examined+=1;vy=n.charge(64*N,y)
  if y in routes.V:continue
  if vy==(0,0):continue
  u,clock=n.edge(y)
  if u not in routes.V or u>=4*N:continue
  vu=n.charge(64*N,u);sy=vy if clock%2==0 else (-vy[0],-vy[1])
  if n.beta(vu,sy)!=-3:continue
  a=2 if y%8==1 else 1;j=(clock-a)//2
  ck(inverse(u,a,j)==y,'exact symbolic inverse/source binding')
  for s in (u,y):
   t,c,d=routes.get(s);v=n.charge(64*N,s);r=o.Row(s,s%15,s%35,t,c,d,v)
   o.audit_row(N,r,routes);ck(v==n.independent_charge(64*N,s),'independent charge')
   rows[s]=o.payload(r)
  ck(tuple(rows[u][3:6])==(u,0,0) and tuple(rows[y][3:6])==(u,clock,1),'actual first-arrival endpoints')
  credits.append([u,y,u,0,clock,0,1,-3,a,j])
  if len(credits)==len(jobs):break
 ck(len(credits)==94,'bounded one-step credit compiler succeeds')
 ck(len({tuple(c[:2]) for c in credits})==94,'unique original old/new pair addresses')
 ck(all(c[0]<4*N<=c[1]<256*N for c in credits),'old/new disjoint from all-new sibling addresses')
 assignments=[[b['id'],c[0],c[1],b['start'],b['end']] for b,c in zip(jobs,credits)]
 ck(len({r[0] for r in assignments})==94,'each positive owner assigned exactly once')
 for assignment,credit in zip(assignments,credits):
  ck(assignment[1:3]==credit[:2] and assignment[3]>=max(credit[5:7]),'actual pair owner born before assigned interval')
 # Exact source ledger transfer, including all physical horizons and post-interval tail.
 pb=(HERE.parent/'annulus_block_completion_20261020j/receipt.json').read_bytes().replace(b'\r\n',b'\n')
 ck(sha256(pb).hexdigest()=='b85cfd8a95962c827d9e21d1418a8fe7f4086c6cca6203f042f547d32ad3b609','frozen physical ledger source')
 prior=json.loads(pb)
 ledger=prior['full_native_control']['full_ledger'];out=[]
 ck(len(ledger)==D,'complete physical horizon universe')
 for d,row in enumerate(ledger):
  C=sum(c[7] for c in credits if max(c[5],c[6])<=d)
  U=base['residual_debit'][d];J=base['original_bank_debit'][d]
  R=row[6];Rp=R-C;I=row[5];M=row[2]
  ck(C<=0 and U+C<=0,'all residual positive jobs paid by arrived external owners')
  ck(R+J==I==Rp+(J+C),'physical interaction ledger exactly unchanged')
  ck(Rp+(U+C)==R+U,'retaining transferred slack gives identical upper bound')
  out.append([d,C,U,R,Rp,I,M,Rp-M])
 ck(reject(lambda:inverse(87,True,1)),'bool branch rejected')
 ck(reject(lambda:inverse(87,1,1.0)),'float inverse depth rejected')
 # d=0 shows why a negative interaction cannot be credited before both sources arrive.
 ck(out[0][1]==0 and out[1][1]==-282,'no future old/new credit')
 ck(len({tuple(c[:2]) for c in credits+[credits[0]]})<95,'duplicate pair-address hostile detected')
 lifts=[]
 for k in (0,1,2,7):
  scale=64**(3*k);registers=[]
  for u,y,t,oc,c,od,nd,w,a,j in credits:
   yp=scale*y+(2**a+3)*(scale-1)//9
   ck(inverse(u,a,j+9*k)==yp,'same-parent three-lift family')
   ck(n.charge(64*scale*N,u)==n.charge(64*N,u),'old charge retained at scale ray')
   ck(n.charge(64*scale*N,yp)==n.charge(64*N,y),'new band/mod3 charge retained')
   vp=n.charge(64*scale*N,yp);vp=vp if c%2==0 else (-vp[0],-vp[1])
   ck(n.beta(n.charge(64*scale*N,u),vp)==-3,'lifted actual signed credit')
   for mod in (15,35,105):ck(lifted_residue(u,a,j,k,mod)==yp%mod,'full-denominator prime registers')
   registers.append([u,yp%15,yp%35,c+18*k])
  lifts.append({'k':k,'scale':scale,'rows':registers})
 return {'scope':'finite N4096 old/new scalar interaction transfer; no physical energy improvement; no uniform claim',
  'source_flow_sha256':sha256(bb).hexdigest(),'N':N,'examined_primitive_new_sources':examined,
  'last_scanned_source':y,'credits_columns':['old','new','terminal','oldclock','newclock','olddepth','newdepth','beta','branch','strippings'],
  'credits':credits,'source_rows':[rows[s] for s in sorted(rows)],'assignments':assignments,
  'ledger_columns':['d','C','residual_positive','R','Rprime','I','M','Rprime_minus_M'],'ledger':out,
  'maximum_repaired_remainder_excess':max(r[-1] for r in out),'scale_lift_controls':lifts,
  'route_audit':routes.audit(),'checks':CHECKS}
if __name__=='__main__':
 result=run();data=enc(result);p=HERE/'receipt.json'
 if '--write' in sys.argv:p.write_bytes(data)
 else:ck(p.read_bytes().replace(b'\r\n',b'\n')==data,'saved receipt')
 print(json.dumps({'checks':CHECKS,'sha256':sha256(data).hexdigest(),'last':result['last_scanned_source'],
 'credits':len(result['credits']),'max_remainder_excess':result['maximum_repaired_remainder_excess']}))
