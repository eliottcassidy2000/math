#!/usr/bin/env python3
"""Independent adjacency/DP checks of the signed-controller model."""
import importlib.util,itertools as it,json,math,os
from pathlib import Path
dp_path=Path(os.environ.get('HALFPATH_DP_PATH',str(Path(__file__).resolve().parent.parent/'session29-halfpath'/'check.py')))
spec=importlib.util.spec_from_file_location('halfpath',dp_path)
hp=importlib.util.module_from_spec(spec);spec.loader.exec_module(hp)
def data(m,q,lab):
 es=list(it.combinations(range(m),2));ctrl={};sgn={}
 for k,(u,v) in enumerate(es):
  c=v if q>>k&1 else u;w=u if c==v else v
  ctrl[frozenset((u,v))]=c
  # lab is the required sign for traversal u then v, with +=1.
  sgn[c,w]=(1 if lab>>k&1 else -1)*(1 if c==v else -1)
 return ctrl,sgn

def factor(m,q,lab):
 ctrl,b=data(m,q,lab);F=0;U=0;H=0
 for p in it.permutations(range(m)):
  f=1;c=0
  for j in range(1,m-1):
   v,u,w=p[j],p[j-1],p[j+1]
   if ctrl[frozenset((u,v))]==v==ctrl[frozenset((v,w))]:
    c+=1;f*=1-b[v,u]*b[v,w]
  if m>1 and ctrl[frozenset(p[-2:])]==p[-1]:
   c+=1;f*=1+b[p[-1],p[-2]]
  F+=f;U+=2**c;H+=int(all(ctrl[frozenset((u,v))]==u for u,v in zip(p,p[1:])))
 return F,U,H

def tournament(m,q,lab,odd):
 ctrl,b=data(m,q,lab);n=2*m+odd;T=[0]*n
 def put(i,j):T[i]|=1<<j
 for i in range(m):
  put(2*i,2*i+1)
  if odd:put(2*i,2*m);put(2*m,2*i+1)
 for i in range(m):
  for j in range(i+1,m):
   c=ctrl[frozenset((i,j))];w=j if c==i else i
   for si in (1,-1):
    for sj in (1,-1):
     pin=si if c==i else sj
     need=-b[c,w] if c==i else b[c,w]
     vi=2*i+(si==-1);vj=2*j+(sj==-1)
     put(vi,vj) if pin==need else put(vj,vi)
 s=tuple([x^1 for x in range(2*m)]+([2*m] if odd else []))
 return T,s

def main():
 out=[];checks=0
 for m in range(1,4):
  E=m*(m-1)//2
  for q in range(1<<E):
   total=0
   for lab in range(1<<E):
    f,U,H=factor(m,q,lab);total+=f
    for odd in (0,1):
     T,s=tournament(m,q,lab,odd)
     assert hp.dp(T,s)['fixed_hp']==f
     checks+=1
    assert f>=H and f%2==1
   assert total==math.factorial(m)*(1<<E)
   assert factor(m,q,((1<<E)-1)^q)[0]==factor(m,q,0)[2]
 for m,q,lab in [(3,0,1),(4,0,9),(5,0,150),(5,2,330)]:
  f,U,H=factor(m,q,lab)
  row={'m':m,'q':q,'lab':lab,'fixed_hp':f,'relaxed':U,'controller_H':H,'tournaments':[]}
  for odd in (0,1):
   T,s=tournament(m,q,lab,odd);assert hp.dp(T,s)['fixed_hp']==f
   row['tournaments'].append({'n':len(T),'outmasks':T,'sigma':s})
   checks+=1
  out.append(row)
 result={'independent_dp_checks':checks,'witnesses':out}
 Path(__file__).with_name('independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
