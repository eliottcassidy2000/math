from module_tools import *
from itertools import combinations
import json

def multi(T,O,ports,z):
 labels=sum([[i]*(z if i in ports else 1) for i in range(len(T))],[]);out=[0]*len(labels);P=0
 for i,x in enumerate(labels):
  if O>>x&1:P|=1<<i
  for j,y in enumerate(labels):
   if (x==y and i<j) or (x!=y and T[x]>>y&1):out[i]|=1<<j
 return tuple(out),P

def reverse(T):return tuple(sum(1<<i for i in range(len(T)) if T[i]>>j&1) for j in range(len(T)))
seeds=json.loads(Path(__file__).with_name('even_births.json').read_text())['seeds'];out=[];count=0
for T,O,q in seeds:
 T=tuple(T)
 if m.data(reverse(T),O)[2]!=q:continue
 for pair in combinations([v for v in range(6) if not O>>v&1],2):
  count+=1;Z,P=multi(T,O,pair,2);h,w,u=m.data(Z,P)
  if h-w!=1:continue
  vv=m.data(reverse(Z),P)[2]
  Z3,P3=multi(T,O,pair,3);h3,w3,u3=m.data(Z3,P3);v3=m.data(reverse(Z3),P3)[2]
  z={'T':T,'O':O,'q':q,'pair':pair,'m2':(h,h-w,u,vv),'m3':(h3,h3-w3,u3,v3)};out.append(z)
  if u==vv or u-4*q in (2,4,6):print(z,flush=True)
Path(__file__).with_name('balanced_pairs.json').write_text(json.dumps({'tested':count,'unit_first_step':out},indent=2))
print('END',count,len(out))
