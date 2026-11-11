"""Targeted two-flank search for the remaining 35-grid heads."""
import subprocess,json
from pathlib import Path
ns={'__name__':'inherited_compiler'}
exec(subprocess.check_output(['git','show','f92e070445:05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py']).decode(),ns)
C=(2,4,9,3);sigs=[ns['signature'](C,O) for O in range(16)];count=ns['count'];hits=[]
for s in (3,4,5):
 for r0 in range(0,7):
  coeff=35*2**(s+r0)
  if coeff%32:continue
  P=coeff//32
  for u in range(1,P+1):
   if P%u:continue
   v=P//u
   for a in range(-12,13):
    for b in range(-12,13):
     for O,sig in enumerate(sigs):
      ok=True
      for k in (6,7,8):
       A=u*2**k-a;B=v*2**k-b
       if count(sig,A,B)!=1+2**s+coeff*4**k:ok=False;break
      if ok:
       start=next((j for j in range(7) if all(u*2**k-a>=0 and v*2**k-b>=0 and count(sig,u*2**k-a,v*2**k-b)==1+2**s+coeff*4**k for k in range(j,10))),None)
       h={'s':s,'r0':r0,'u':u,'v':v,'a':a,'b':b,'O':O,'start_k':start};hits.append(h);print(json.dumps(h),flush=True)
Path(__file__).with_name('balanced_shell_receipt.json').write_text(json.dumps(hits,indent=2)+'\n')
print('hits',len(hits))
