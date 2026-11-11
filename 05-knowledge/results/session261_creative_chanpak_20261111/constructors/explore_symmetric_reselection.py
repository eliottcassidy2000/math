from module_tools import *
from itertools import combinations
import json

def graph(T,ports,z):
 labels=sum([[i]*(z if i in ports else 1) for i in range(len(T))],[]);out=[0]*len(labels)
 for i,x in enumerate(labels):
  for j,y in enumerate(labels):
   if (x==y and i<j) or (x!=y and T[x]>>y&1):out[i]|=1<<j
 return tuple(out),labels

def pairdata(T):
 n=len(T);e=m.endpoints(T);R=tuple(sum(1<<i for i in range(n) if T[i]>>j&1) for j in range(n));s=m.endpoints(R)
 F=(1<<n)-1;K=[[0]*n for _ in range(n)];U=[0]*n;V=[0]*n
 for A in range(1<<n):
  B=F^A;ha=sum(e[A]) if A else 1;hb=sum(e[B]) if B else 1
  for j in range(n):
   U[j]+=ha*s[B][j];V[j]+=e[A][j]*hb
   if s[B][j]:
    for i in range(n):K[i][j]+=e[A][i]*s[B][j]
 h=sum(e[-1]);out=[]
 for i,j in combinations(range(n),2):
  if h-K[i][j]-K[j][i]==1:out.append({'pair':(i,j),'U':U[i]+U[j],'V':V[i]+V[j]})
 return h,out
out=[]
for ports in [(0,4),(1,3),(0,1,3,4)]:
 for z in range(1,5 if len(ports)==2 else 3):
  T,labels=graph((2,4,9,19,7),ports,z);h,rs=pairdata(T)
  d={'ports':ports,'length':z,'H':h,'balanced':[r for r in rs if r['U']==r['V']],'all_unit':rs,'labels':labels};out.append(d);print(d,flush=True)
Path(__file__).with_name('symmetric_reselection.json').write_text(json.dumps(out,indent=2))
