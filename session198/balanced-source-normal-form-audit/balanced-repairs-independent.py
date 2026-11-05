"""Independent one-vertex compilation and counting of newly added C1 repairs."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
root=HERE.parent
R=json.loads((root/'balanced-source-normal-form/REPAIRS.json').read_text())
prior=json.loads((root/'c1-monotone/REPAIR_WITNESSES.json').read_text())
assert all(R[k]==v for k,v in prior.items())

def birth(T,c,x):
 n=len(T);A=T.copy()+[0]
 for u in range(n):
  if x!='O'and u==c:A[n]|=1<<u
  else:A[u]|=1<<n
 return A,n if x=='J'else c

def counts(T,c):
 n=len(T);size=1<<n;P=[[0]*n for _ in range(size)];Q=[[0]*n for _ in range(size)]
 for i in range(n):P[1<<i][i]=1
 Q[1<<c][c]=1
 for S in range(1,size):
  for u in range(n):
   if not P[S][u]:continue
   rest=T[u]&~S
   while rest:
    b=rest&-rest;rest-=b;v=b.bit_length()-1
    P[S|b][v]+=P[S][u];Q[S|b][v]+=Q[S][u]
 H=[sum(p)for p in P];H[0]=1
 U=sum(H[(size-1)^S]*sum(Q[S])for S in range(1,size))
 profile=(H[-1],P[-1][c],U,H[(size-1)^(1<<c)])
 C=sum(Q[-1][u]for u in range(n)if T[u]>>c&1)
 return profile,sum(Q[-1]),C

repairs=[]
for target,entry in R.items():
 if entry['cycles']!=1:continue
 T,c=[0],0
 for typ,L in entry['record']['source_recipe']:
  word='O'if typ=='O'else'J'+(''if L==1 else''.join('Z'if b=='1'else'O'for b in bin(L-1)[:1:-1]))
  for x in word:T,c=birth(T,c,x)
 prof,s,C=counts(T,c);assert s==1 and list(prof)==entry['record']['source']
 T,c=birth(T,c,'J');prof,s,C=counts(T,c)
 assert prof[0]==int(target)and C==1
 repairs.append({'target':int(target),'order':len(T)})
assert len(repairs)==13
out={'status':'PASS','prior_C2_C3_records_unchanged':len(prior),'new_C1_repairs':repairs,'method':'compile each M_L into J plus a low-bit O/Z word, then independent subset counting'}
(HERE/'BALANCED_REPAIRS_INDEPENDENT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
