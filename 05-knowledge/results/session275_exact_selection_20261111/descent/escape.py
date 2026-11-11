"""A literal reusable seam escapes the normalized-21 gadget content."""
import json
from pathlib import Path

T=(2,508,433,165,9,17,189,49,249)

def rooted(rows,a):
 n=len(rows); full=(1<<n)-1
 ends=[[0]*n for _ in range(full+1)];start=[[0]*n for _ in range(full+1)];h=[0]*(full+1);h[0]=1
 for v in range(n):ends[1<<v][v]=1
 start[1<<a][a]=1
 for s in range(1,full+1):
  for v in range(n):
   if not(s>>v&1):continue
   q=s^(1<<v)
   if q:
    ends[s][v]=sum(ends[q][u] for u in range(n) if q>>u&1 and rows[u]>>v&1)
    start[s][v]=sum(start[q][u] for u in range(n) if q>>u&1 and rows[u]>>v&1)
  h[s]=sum(ends[s])
 U=sum(sum(start[s])*h[full^s] for s in range(full+1))
 return (h[full],ends[full][a],U,h[full^(1<<a)])

def fan(rows,a,k):
 n=len(rows);m=((1<<k)-1)<<n
 return tuple(r if i==a else r|m for i,r in enumerate(rows))+tuple((1<<a)|sum(1<<(n+j) for j in range(i+1,k)) for i in range(k))

def require(t,msg):
 if not t:raise RuntimeError(msg)

menus=[]
for a in range(9):
 D,E,U,d=rooted(T,a);N=D+U-E
 menus.append({'root':a,'packet':[D,E,U,d],'singleton_output':N,'mod29':N%29})
 require(rooted(fan(T,a,1),a)[0]==N,'singleton law')
require(rooted(T,7)==(609,108,1066,241),'source packet')
controls=[]
for k in range(1,6):
 p=rooted(fan(T,7,k),7); expected=(825*2**k-83,241,825*2**k+241,241)
 require(p==expected,'full fan packet')
 controls.append({'k':k,'packet':p})
require(pow(2,28,29)==1 and all(pow(2,d,29)!=1 for d in (1,2,4,7,14)),'order')
require((2**28-1)%29==0 and (2**28-1)%(29**2)!=0,'simple lift')
require([k for k in range(1,29) if (825*2**k-83)%29==0]==[26],'bad class')
# Finite arithmetic controls of the all-level unique-lift proof.
classes=[];k=26;period=28
for r in range(1,5):
 mod=29**r
 require((825*pow(2,k,mod)-83)%mod==0,'root')
 classes.append({'r':r,'residue':k,'period':period})
 lifts=[k+j*period for j in range(29) if (825*pow(2,k+j*period,mod*29)-83)%(mod*29)==0]
 require(len(lifts)==1,'unique lift')
 k=lifts[0];period*=29
print(json.dumps({'root_menu':menus,'fan_controls':controls,'bad_exponent_classes':classes},indent=2))
