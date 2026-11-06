"""Exact unordered path-forest seed vectors for two retained active ports."""
from collections import Counter
from functools import lru_cache
BULK=-1

def path_tables(a,ports):
 n=len(a);N=1<<n;lab={v:i for i,v in enumerate(ports)};types=[Counter()for _ in range(N)];starts=[]
 for first in range(n):
  dp=[[0]*n for _ in range(N)];dp[1<<first][first]=1
  for mask in range(1,N):
   for last in range(n):
    z=dp[mask][last]
    if not z:continue
    types[mask][(lab.get(first,BULK),lab.get(last,BULK))]+=z
    choices=a[last]&~mask
    while choices:
     bit=choices&-choices;choices-=bit;j=bit.bit_length()-1;dp[mask|bit][j]+=z
  starts.append(sum(dp[N-1]))
 return types,starts

def forest(a,ports=(0,1),p=2):
 a=tuple(a);ports=tuple(ports);types,starts=path_tables(a,ports);maxparts=len(ports)+p
 @lru_cache(None)
 def F(mask):
  if not mask:return {():1}
  least=mask&-mask;sub=mask;out=Counter()
  while sub:
   if sub&least:
    for old,count in F(mask^sub).items():
     if len(old)>=maxparts:continue
     used=sum(x==BULK for x,y in old)
     for part,mul in types[sub].items():
      if used+int(part[0]==BULK)>p:continue
      out[tuple(sorted(old+(part,)))]+=count*mul
   sub=(sub-1)&mask
  return dict(out)
 return F((1<<len(a))-1),starts

def read(v,root=0):
 h=s=L=0
 for pattern,count in v.items():
  if len(pattern)==1:
   h+=count
   if pattern[0][0]==root:s+=count
  if len(pattern)in(1,2)and any(e==root for _,e in pattern):L+=count
 return h,s,L,h+L-s

if __name__=='__main__':
 import sys,json,time
 sys.path.insert(0,'session203-positive-profile-defects');from counting import labelled,nine
 t=time.time();checks=0
 for n in range(2,6):
  for a in labelled(n):
   v,starts=forest(a);h,s,e,L,R,d,C,b,T=nine(a)
   if read(v)!=(h,s,L,h+L-s)or starts[0]!=s:raise RuntimeError((a,read(v),nine(a)))
   checks+=1
 print(json.dumps(dict(status='PASS',labeled_graphs=checks,elapsed=time.time()-t)))
