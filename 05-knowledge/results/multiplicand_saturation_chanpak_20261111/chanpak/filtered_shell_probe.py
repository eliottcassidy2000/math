"""Targeted balanced-shell constructors, filtered by exact leading coefficient."""
import subprocess,json
from pathlib import Path
from head_lift_probe import ends,converse
ns={'__name__':'inherited_compiler'}
exec(subprocess.check_output(['git','show','f92e070445:05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py']).decode(),ns)
records=json.loads(subprocess.check_output(['git','show','f92e070445:05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json']))['records']
bank=[(0,),(2,0)]
for r in records:
 T=tuple(r['adjacency_outmasks']);bank.append(T)
 if r['ordinary_multiplicity']==2:bank.append(tuple(converse(T)))
eligible=[]
for C in bank:
 e=ends(C);h=[sum(x) for x in e];h[0]=1;M=(1<<len(C))-1;Q=sum(h[m]*h[M^m] for m in range(M+1));odd=Q
 while odd%2==0:odd//=2
 if 35%odd==0:eligible.append((C,Q))
print('eligible',len(eligible),'of',len(bank),flush=True)
hits=[];count=ns['count']
for index,(C,Q) in enumerate(eligible):
 sigs=[ns['signature'](C,O) for O in range(1<<len(C))]
 for s in (3,4,5):
  for r0 in range(5):
   coeff=35*2**(s+r0)
   if coeff%Q:continue
   P=coeff//Q
   for u in range(1,P+1):
    if P%u:continue
    v=P//u
    for a in range(-3,5):
     for b in range(-3,5):
      for O,sig in enumerate(sigs):
       if not all(count(sig,u*2**k-a,v*2**k-b)==1+2**s+coeff*4**k for k in (4,5,6)):continue
       start=next((j for j in range(5) if all(u*2**k-a>=0 and v*2**k-b>=0 and count(sig,u*2**k-a,v*2**k-b)==1+2**s+coeff*4**k for k in range(j,9))),None)
       hit={'core':C,'Q2':Q,'s':s,'r0':r0,'u':u,'v':v,'a':a,'b':b,'cut':O,'start_k':start,'signature':sig};hits.append(hit);print(json.dumps(hit),flush=True)
Path(__file__).with_name('filtered_shell_receipt.json').write_text(json.dumps({'eligible_cores':len(eligible),'hits':hits},indent=2)+'\n')
print('hits',len(hits),flush=True)
