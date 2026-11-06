from itertools import permutations,product
from collections import defaultdict
from functools import lru_cache
import json
Z=(0,0,0,0);ONE={Z:1}
X=[{tuple(int(i==j) for i in range(4)):1}for j in range(4)]
def add(*args):
 d=defaultdict(int)
 for c,p in args:
  for m,v in p.items():d[m]+=c*v
 return {m:c for m,c in d.items()if c}
def mul(p,q):
 d=defaultdict(int)
 for m,c in p.items():
  for n,e in q.items():d[tuple(a+b for a,b in zip(m,n))]+=c*e
 return {m:c for m,c in d.items()if c}
def tensor(core,outs,incs,kind='h'):
 n=len(core);total=defaultdict(int)
 for order in permutations(range(n)):
  if kind=='s' and order[0]!=0:continue
  if kind=='e' and order[-1]!=0:continue
  p=ONE
  for gap in range(n+1):
   left=order[gap-1]if gap else None;right=order[gap]if gap<n else None
   f={Z:1}if left is None or right is None or core[left]>>right&1 else{}
   if not((kind=='s'and gap==0)or(kind=='e'and gap==n)):
    for a,b in product(range(2),repeat=2):
     if(left is None or outs[a]>>left&1)and(right is None or incs[b]>>right&1):f=add((1,f),(1,X[2*a+b]))
   p=mul(p,f)
   if not p:break
  for m,c in p.items():total[m]+=c
 return dict(total)
def profiles():
 a=[102,0,10,83,135,30,182,47];v=5;keep=[i for i in range(8)if i!=v]
 c=[sum(1<<j for j,y in enumerate(keep)if a[x]>>y&1)for x in keep]
 out=sum(1<<i for i,x in enumerate(keep)if a[x]>>v&1);inc=sum(1<<i for i,x in enumerate(keep)if a[v]>>x&1)
 outs=[out&~8,out|8];incs=[inc|8,inc&~8]
 h=tensor(c,outs,incs);s=tensor(c,outs,incs,'s');e=tensor(c,outs,incs,'e')
 cr=c[:]
 for i in range(1,len(c)):cr[i]|=1
 R=tensor(cr,outs,[x|1 for x in incs]);cl=c[:];cl[0]=(1<<len(c))-2
 L=tensor(cl,[x|1 for x in outs],incs)
 dc=[x>>1 for x in c[1:]];do=[x>>1 for x in outs];di=[x>>1 for x in incs];d=tensor(dc,do,di)
 n=len(dc);bc=[x|(1<<n)for x in dc]+[(1<<n)-1];beta=tensor(bc,[x|(1<<n)for x in do],[x|(1<<n)for x in di])
 return dict(h=h,s=s,e=e,L=L,R=R,d=d,beta=beta,B=add((1,h),(1,L),(-1,s)),V=add((1,R),(1,beta)))
# Coordinates Q,A,B,Y: x=A-1-Y, y=Y, z=Q-A-B+1+Y, w=B-1-Y.
subs=[add((1,X[1]),(-1,ONE),(-1,X[3])),X[3],add((1,X[0]),(-1,X[1]),(-1,X[2]),(1,ONE),(1,X[3])),add((1,X[2]),(-1,ONE),(-1,X[3]))]
def transform(p):
 d={}
 for m,c in p.items():
  q=ONE
  for i,power in enumerate(m):
   for j in range(power):q=mul(q,subs[i])
  d=add((1,d),(c,q))
 return d
@lru_cache(None)
def evaluate_monomial(m,mask,word):
 if not any(m):return int(mask==0)
 i=next(i for i,e in enumerate(m)if e);mm=list(m);mm[i]-=1;mm=tuple(mm);ans=0;sub=mask
 while True:
  if not sub:ok=(i!=3)
  else:
   first=(sub&-sub).bit_length()-1;last=sub.bit_length()-1;a=word>>first&1;b=word>>last&1
   ok=(i==0)or(i==1 and a==0)or(i==2 and b==1)or(i==3 and a==0 and b==1)
  if ok:ans+=evaluate_monomial(mm,mask^sub,word)
  if not sub:break
  sub=(sub-1)&mask
 return ans
ps=profiles();report={}
for name,p in ps.items():
 q=transform(p);values={}
 for word in [10,14,17,21]:
  vals=defaultdict(int)
  for m,c in q.items():vals[(m[0]+m[2])*(m[0]+m[1])]+=c*evaluate_monomial(m,31,word)
  values[word]={k:v for k,v in sorted(vals.items())if v}
 report[name]={'degree':max(map(sum,p))if p else 0,'eigenvalues':values,'tensor':[[list(m),c]for m,c in sorted(p.items())],'diagonal_tensor':[[list(m),c]for m,c in sorted(q.items())]}
for name in ['B','e','V','d']:
 vals=report[name]['eigenvalues'];base=vals[10]
 if any(vals[w]!=base for w in [14,17,21]):raise RuntimeError(('row failure',name,vals))
for name in ['h','s','L','R','beta']:
 vals=report[name]['eigenvalues'];base=vals[10]
 for w in [14,17,21]:
  diff={k:vals[w].get(k,0)-base.get(k,0)for k in set(vals[w])|set(base)}
  report[name].setdefault('differences',{})[w]={k:v for k,v in sorted(diff.items())if v}
print(json.dumps(report,indent=2))
