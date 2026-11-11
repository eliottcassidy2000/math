from module_tools import *
import json
from itertools import combinations

def substitute(S,v,M):
 old=[x for x in range(len(S)) if x!=v];n=len(old);out=[0]*(n+len(M))
 for i,x in enumerate(old):
  for j,y in enumerate(old):
   if S[x]>>y&1:out[i]|=1<<j
  if S[x]>>v&1:out[i]|=((1<<len(M))-1)<<n
 for i in range(len(M)):
  out[n+i]=M[i]<<n
  for j,x in enumerate(old):
   if S[v]>>x&1:out[n+i]|=1<<j
 return tuple(out)

def outside_data(T,k):
 n=len(T);e=m.endpoints(T);rev=tuple(sum(1<<i for i in range(n) if T[i]>>j&1) for j in range(n));s=m.endpoints(rev)
 F=(1<<n)-1;K=[[0]*k for _ in range(k)];U=[0]*k
 for A in range(1<<n):
  B=F^A
  if not B:continue
  hs=sum(e[A]) if A else 1
  for j in range(k):
   bs=s[B][j]
   if not bs:continue
   U[j]+=hs*bs
   for i in range(k):K[i][j]+=e[A][i]*bs
 h=sum(e[-1]);ret=[]
 for O in range(1,1<<k):
  z=[i for i in range(k) if O>>i&1];w=sum(K[i][j] for i in z for j in z)
  if h-w==1:ret.append((O,h,sum(U[i] for i in z)))
 return ret

A=(2,4,9,19,39,15);hosts=[];tested=0
for n in range(3,6):
 cs={canon(tour(n,b)) for b in range(1<<(n*(n-1)//2))}
 for S in sorted(cs):
  for v in range(n):
   T=substitute(S,v,A);tested+=1
   rs=outside_data(T,n-1)
   for O,h,q in rs:
    z={'quotient':S,'port':v,'rows':T,'selected_mask':O,'H':h,'U':q};hosts.append(z);print(z,flush=True)
print('END',tested,len(hosts))
Path(__file__).with_name('unit_module_hosts.json').write_text(json.dumps({'quotient_orders':[3,4,5],'tested':tested,'hosts':hosts},indent=2))
