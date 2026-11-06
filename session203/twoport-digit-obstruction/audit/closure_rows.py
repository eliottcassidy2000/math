import pathlib,sys,json
from collections import Counter
p=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(p));import verify as v

def read(p,r):
 h=int(len(p)==1);s=int(len(p)==1 and p[0][0]==r);L=int(len(p)in(1,2)and any(b==r for a,b in p));return h+L-s,s
checks=0
for ci in range(6):
 for j,p in enumerate(v.P):
  variants=[]
  if ci<2:variants=[(p,ci)]
  else:
   allowed={i for i in range(2)if (ci-2)>>i&1}
   for pre in [None]+list(range(len(p))):
    if pre is not None and p[pre][1]in allowed:continue
    for suc in [None]+list(range(len(p))):
     if suc is not None and p[suc][0]not in allowed:continue
     if pre is not None and pre==suc:continue
     parts=[x for k,x in enumerate(p)if k not in (pre,suc)]+[(2 if pre is None else p[pre][0],2 if suc is None else p[suc][1])]
     variants.append((parts,2))
  out=[0,0]
  for parts,r in variants:
   q=[(a if a==r else -1,b if b==r else -1)for a,b in parts]
   if sum(a==-1 for a,b in q)>2:continue
   z=read(q,r);out=[x+y for x,y in zip(out,z)]
  assert out==[v.C[ci][k][j]for k in range(2)],(ci,j,out)
  checks+=1
print(json.dumps(dict(status='PASS',formal_closure_columns=checks,method='Independent local attachment enumeration and final endpoint readout'),indent=2))
