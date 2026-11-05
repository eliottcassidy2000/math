from functools import lru_cache
from itertools import permutations
from math import prod
import json

def need(x,s):
 if not x:raise RuntimeError(s)
def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 while len(c)>1 and not c[-1]:c.pop()
 return c
def scale(a,c):return [x*c for x in a]
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def graph(m,pad=False):
 sizes=[m,3,2];groups=sum(([k]*n for k,n in enumerate(sizes)),[]);a=[0]*sum(sizes)
 for i in range(len(a)):
  for j in range(i+1,len(a)):
   f=groups[i]==groups[j] or (groups[j]-groups[i])%3==1
   a[i if f else j]|=1<<(j if f else i)
 if pad:
  n=len(a);a=[x|(1<<n)for x in a]+[0]
 return tuple(a)
def profile(a):
 n=len(a);full=(1<<n)-1;cycles=[{}for _ in a]
 for root in range(n):
  def dfs(mask,last):
   if mask.bit_count()>=3 and mask.bit_count()%2 and a[last]>>root&1:cycles[root][mask]=cycles[root].get(mask,0)+1
   q=a[last]&~mask&full&~((1<<(root+1))-1)
   while q:
    b=q&-q;q-=b;dfs(mask|b,b.bit_length()-1)
  dfs(1<<root,root)
 @lru_cache(None)
 def pack(mask):
  if not mask:return (1,)
  b=mask&-mask;r=b.bit_length()-1;out=[0]+list(pack(mask^b))
  for cyc,count in cycles[r].items():
   if cyc&mask==cyc:out=add(out,[0]+scale(pack(mask^cyc),2*count))
  return tuple(out)
 return list(pack(full))
def basis(k,parity):
 a=[0]*(parity+1)+[1]
 for j in range(1,k+1):a=mul(a,[-j*j,0,1])
 return a
def coords(p,parity):
 a=p[:];n=(len(a)-2-parity)//2;c=[0]*(n+1)
 for k in range(n,-1,-1):
  b=basis(k,parity);v=a[len(b)-1];c[k]=v
  for j,x in enumerate(b):a[j]-=v*x
 need(not any(a),'basis remainder');return c
def det(M):
 out=0
 for p in permutations(range(len(M))):
  sign=(-1)**sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))
  out+=sign*prod(M[i][p[i]]for i in range(len(p)))
 return out
P={m:profile(graph(m))for m in range(7)}
need(graph(1)==(14,60,56,48,33,1),'root rows')
need(graph(1,True)==(78,124,120,112,97,65,0),'padded rows')
for m in (0,1,3,5):need(profile(graph(m,True))==[0]+P[m],'literal padding')
columns=[P[1],add(P[3],scale(P[1],-1)),add(add(P[5],scale(P[3],-5)),scale(P[1],4))]
r=[coords(p,1)for p in columns] # parity=1 denotes x^2 start
q=[coords([0]+p,0)for p in columns]
M=[[r[j][i]for j in range(3)]for i in range(3)]
PM=[[q[j][i]for j in range(3)]for i in range(3)]
need(M==[[23,486,2880],[17,424,6220],[1,49,1571]],'physical minor')
need(PM==[[23,486,2880],[91,2182,27760],[26,865,20359]],'padded minor')
need(det(M)==-468310 and det(PM)==-16859160,'two negative determinants')
need(all(x>=0 for c in r+q for x in c),'positive columns')
print(json.dumps({'status':'PASS','physical_profiles':P,'unpadded_columns':r,'padded_columns':q,'M':M,'PM':PM,'det':det(M),'padded_det':det(PM),'physical_profile_count':11},sort_keys=True,indent=2))
