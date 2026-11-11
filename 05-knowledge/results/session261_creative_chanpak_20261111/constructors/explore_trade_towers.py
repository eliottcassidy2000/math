from module_tools import *
from itertools import combinations
import json
A=(2,4,9,19,39,15);B=(2,4,57,3,11,27)
def multi(T,ports,z):
 groups=[[i]* (z if i in ports else 1) for i in range(len(T))];labels=sum(groups,[]);out=[0]*len(labels)
 for i,x in enumerate(labels):
  for j,y in enumerate(labels):
   if (x==y and i<j) or (x!=y and T[x]>>y&1):out[i]|=1<<j
 return tuple(out)
D={}
for p in combinations(range(6),2):D.setdefault(covers(multi(A,p,2)),[]).append(p)
hits=[]
for q in combinations(range(6),2):
 z=covers(multi(B,q,2))
 for p in D.get(z,[]):
  x=covers(multi(A,p,3));y=covers(multi(B,q,3))
  r={'A_ports':p,'B_ports':q,'m2':z,'m3_equal':x==y,'m3A':x,'m3B':y};hits.append(r);print(r,flush=True)
Path(__file__).with_name('trade_towers.json').write_text(json.dumps(hits,indent=2))
