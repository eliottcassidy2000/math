import itertools,json
from pathlib import Path
A=[102,0,10,83,135,30,182,47]; B=[142,0,34,22,71,155,47,94]
def strip(a,v):
 keep=[i for i in range(len(a)) if i!=v]
 core=[sum(1<<j for j,y in enumerate(keep) if a[x]>>y&1) for x in keep]
 out=sum(1<<i for i,x in enumerate(keep) if a[x]>>v&1)
 inc=sum(1<<i for i,x in enumerate(keep) if a[v]>>x&1)
 return core,out,inc

def tensor(a,out,inc,kind='h'):
 n=len(a); total=[0]*(n+2)
 for p in itertools.permutations(range(n)):
  if kind=='s' and p[0]!=0:continue
  if kind=='e' and p[-1]!=0:continue
  q=[1]
  for j in range(n+1):
   l=p[j-1] if j else -1;r=p[j] if j<n else -1
   empty=int(l<0 or r<0 or a[l]>>r&1)
   filled=int((l<0 or out>>l&1) and (r<0 or inc>>r&1))
   if (kind=='s' and j==0) or(kind=='e' and j==n):filled=0
   z=[0]*(len(q)+1)
   for k,c in enumerate(q):z[k]+=c*empty;z[k+1]+=c*filled
   q=z
  for k,c in enumerate(q):total[k]+=c
 return total

def add(*args):
 z=[0]*12
 for factor,p in args:
  for k,c in enumerate(p):z[k]+=factor*c
 return z

def profiles(a,v):
 c,o,i=strip(a,v); h=tensor(c,o,i);s=tensor(c,o,i,'s');e=tensor(c,o,i,'e')
 cr=c.copy()
 for j in range(1,len(c)):cr[j]|=1
 R=tensor(cr,o,i|1)
 cl=c.copy();cl[0]=(1<<len(c))-2
 L=tensor(cl,o|1,i)
 dcore=[x>>1 for x in c[1:]];d=tensor(dcore,o>>1,i>>1)
 k=len(dcore);bc=[x|(1<<k) for x in dcore]+[(1<<k)-1]
 beta=tensor(bc,(o>>1)|(1<<k),(i>>1)|(1<<k))
 return {'h':h,'s':s,'e':e,'L':L,'R':R,'d':d,'beta':beta,'B':add((1,h),(1,L),(-1,s)),'V':add((1,R),(1,beta))}
p=profiles(A,5);q=profiles(B,3);diff={key:add((1,q[key]),(-1,p[key])) for key in p}
for key in ['B','e','V','d']:
 if any(diff[key][1:]):raise RuntimeError((key,diff[key]))
print(json.dumps({'A':p,'B':q,'delta':diff,'scope':'Exact gap tensors; coefficients multiply ordered nonempty path-cover counts of arbitrary module X.'},indent=2))
