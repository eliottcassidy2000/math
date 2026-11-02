"""Finite controls for the all-order marked-boundary divisibility lemma."""
from explore import tournament,c3_blowup,H,strong,endpoint_Q_dp,extend
from phase2 import encode
from math import gcd
from functools import reduce
import json,hashlib

def check(x,msg):
 if not x:raise RuntimeError(msg)
quotients=[(3,2),(4,4)]+[(5,int(c)) for c in sorted(json.load(open('phase2_order5.json'))['class_H'],key=int)]
ledger=[];presentations=0;direct=0
for k,code in quotients:
 Q=tournament(k,code);check(strong(Q),'strong quotient');T=c3_blowup(Q);hT=H(T);rows=[]
 for r in range(0,3*k,3):
  vs=[i for i in range(3*k) if i!=r];B=tuple(sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs);n=len(B)
  check(strong(B),'strong deletion card');s,e,q=endpoint_Q_dp(B);h=sum(s)
  C=[[q[i][j]+q[j][i] for j in range(n)] for i in range(n)]
  A=[sum(C[i][j] for j in range(n) if B[i]>>j&1) for i in range(n)]
  g=reduce(gcd,A+[C[i][j] for i in range(n) for j in range(i+1,n)])
  check(g%(2*3**max(k-3,0))==0,'forced fullstar difference divisor')
  f=[h]*(1<<n)
  for S in range(1,1<<n):
   bit=S&-S;i=bit.bit_length()-1;rest=S^bit
   f[S]=f[rest]+A[i]-sum(C[i][j] for j in range(n) if rest>>j&1)
  old=sum(1<<j for j,v in enumerate(vs) if T[r]>>v&1)
  zeros=[S for S,v in enumerate(f) if v==hT]
  check(zeros==[old],'finite unique zero-response cut')
  check(all(abs(v-hT)!=2 for v in f),'finite fullstar unit isolation')
  check(all((v-hT)%(2*3**max(k-3,0))==0 for v in f),'all completion differences obey divisor')
  for S in [old,old^1,old^3]:
   check(H(extend(B,S))==f[S],'independent direct child DP control');direct+=1
  rows.append({'deleted_vertex':r,'parent_H':h,'response_gcd':g,'zero_cuts':zeros,'completion_cuts':len(f)})
  presentations+=len(f)
 row={'quotient_order':k,'quotient_code':code,'quotient_H':H(Q),'child_order':3*k,'H':hT,'forced_difference_divisor':2*3**max(k-3,0),'parent_rows':rows}
 ledger.append(row);print(json.dumps(row,sort_keys=True,separators=(',',':')),flush=True)
summary={'quotient_classes':len(ledger),'deletion_orbit_rows':sum(len(x['parent_rows']) for x in ledger),'all_cut_presentations':presentations,'direct_child_DP_controls':direct}
print(json.dumps(summary,sort_keys=True,separators=(',',':')))
print('PASS semantic_sha256='+hashlib.sha256(json.dumps([ledger,summary],sort_keys=True,separators=(',',':')).encode()).hexdigest())
