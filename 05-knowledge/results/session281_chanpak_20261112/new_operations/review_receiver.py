#!/usr/bin/env python3
"""Independent reconstruction of sibling receiver seeds and carry controls."""
import json
from verify import require,direct_profile

def ordinal(parts):
 n=sum(map(len,parts));rows=[0]*n;start=0
 for part in parts:
  for i,r in enumerate(part):rows[start+i]=(r<<start)|sum(1<<j for j in range(start+len(part),n))
  start+=len(part)
 return rows

def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]

def root(D,I):
 n=len(D);r=[sum(1<<(i+1) for i in range(n) if not I>>i&1)]
 r.extend((v<<1)|((I>>i)&1) for i,v in enumerate(D));return r

def c1(T):
 n=len(T);rows=[(T[0]<<2)|2]
 rows.extend([4|sum(1<<j for j in range(3,n+2)),1|sum(1<<j for j in range(3,n+2))])
 for i in range(1,n):rows.append((T[i]&1)|((T[i]&~1)<<2))
 return rows

def g(p):
 if p==0:return 1
 v=p//2
 return g(v) if p%2 else g(v)+v+v%2

def main():
 B=[2,4,1];C=[2,4,9,3];E=[44,1,2,6,15,30]
 recipes=[([tt(3),E,tt(1)],236),([B,tt(1),C,tt(2)],408),([tt(4),E],76),([E,tt(4)],326),([tt(3),C,B],204),([tt(1),C,B,tt(2)],368)]
 checks=[]
 for i,(parts,I) in enumerate(recipes):
  raw=root(ordinal(parts),I);key,H=direct_profile(raw,0)
  require((key[0],key[3],key[5],H)==(15,1203,[15,15,15,13,15,15][i],979+2*i),'raw seed')
  new=c1(raw);k,h=direct_profile(new,0)
  require((k[0],k[3],k[5],h)==(15,4797,15,3385+2*i),'completed seed')
  checks.append({'raw_rows':raw,'H':H,'completed_H':h,'receiver':[k[0],k[3],k[5]]})
 require(g(636)==483,'g boundary');scales=[]
 for k in range(0,30,2):
  V=4792*4**k+5;A=3385+9584*(4**k-1)//3+10*k;z=(V-27)//15;P=2*z
  require((V-27)%15==0 and z%2==0 and z>0,'carry admissibility')
  t=lambda p:15*(g(p)-1)+(V-15)*p
  require(t(P+1)-t(P)==12,'twelve abutment')
  scales.append({'k':k,'N':str(A+t(P)),'receiver_defect':str(-15*2**(P.bit_length()-1))})
 require(scales[0]['N']=='3051967','first twelve bank')
 print(json.dumps({'status':'PASS','literal_raw_and_completed_seed_checks':checks,'exact_carry_controls':scales,'scope':'New seed/control computations independently reconstructed; inherited fusion law dependency, not full re-audit'},sort_keys=True,indent=2))
if __name__=='__main__':main()
