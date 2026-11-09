#!/usr/bin/env python3
"""Independent labelled exhaustive audit; no packet code or data dependencies."""
import itertools,json
from collections import Counter

def require(condition, message):
 if not condition:
  raise RuntimeError(message)

def adj(n,code):
 a=[0]*n
 for b,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if code>>b&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a

def paths(a):
 n=len(a); f=[[0]*n for _ in range(1<<n)]
 for v in range(n):f[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   t=a[v]&~s
   while t:
    b=t&-t;t-=b;f[s|b][b.bit_length()-1]+=f[s][v]
 return f

def cut(a):
 n=len(a);return min(sum((a[v]&s).bit_count() for v in range(n) if not s>>v&1) for s in range(1,(1<<n)-1))

def mod(a):
 n=len(a)
 return any(all((a[v]&s) in (0,s) for v in range(n) if not s>>v&1) for s in range(1,(1<<n)-1) if s.bit_count()>1)

def attach(a,x):
 n=len(a); b=a.copy()+[0]
 for v in range(n):
  if x>>v&1:b[v]|=1<<n
  else:b[n]|=1<<v
 return b

def extension_form(a,f):
 n=len(a);full=(1<<n)-1
 ra=[sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n)];g=paths(ra)
 w=[[0]*n for _ in range(n)]
 for s in range(1,full):
  for i in range(n):
   if s>>i&1:
    for j in range(n):
     if not s>>j&1:w[i][j]+=f[s][i]*g[full^s][j]
 l=[f[full][i]-g[full][i]+sum(w[i]) for i in range(n)]
 return l,w

def run():
 spectra={0:{1},1:{1}};res={};checks=0
 for n in range(2,7):
  spec=set();cats=Counter();floors={};neutral=0
  for c in range(1<<(n*(n-1)//2)):
   a=adj(n,c);f=paths(a);h=sum(f[-1]);spec.add(h);la=cut(a)
   delta=min(min(v.bit_count(),n-1-v.bit_count()) for v in a)
   for k in range(1,n):
    if n<=4*k+1:require((la>=k)==(delta>=k), f"Degree threshold failed: n={n}, code={c}, k={k}")
    if la>=k:floors[k]=min(floors.get(k,h),h)
   fc=h not in spectra[n-1];ext=any(sum(f[((1<<n)-1)^(1<<v)]) not in spectra[n-2] for v in range(n))
   cat='FC' if fc else 'EXT' if ext else 'SUM' if la==0 else 'SUB' if mod(a) else 'REST';cats[cat]+=1
   require(cat!="REST" or la<=2, f"Grammar failed: n={n}, code={c}")
   if n<=5:
    l,w=extension_form(a,f)
    for x in range(1<<n):
     actual=sum(paths(attach(a,x))[-1]); pred=h+sum(l[i] for i in range(n) if x>>i&1)-sum(w[i][j]+w[j][i] for i in range(n) for j in range(i+1,n) if x>>i&1 and x>>j&1)
     require(actual==pred, f"Attachment failed: n={n}, code={c}, x={x}");checks+=1
  spectra[n]=spec;res[n]={'labelled_count':1<<(n*(n-1)//2),'spectrum':sorted(spec),'categories':dict(cats),'floors':floors}
 # Sharp strong examples for k=2,3,4; two regular blocks, one reverse arc.
 sharp=[]
 for k in range(2,5):
  q=2*k+1;n=2*q;a=[0]*n
  for i in range(n):
   for j in range(n):
    if i//q==j//q:
     if 1<=(j-i)%q<=k:a[i]|=1<<j
    elif i<q:a[i]|=1<<j
  a[0]^=1<<q;a[q]^=1
  delta=min(min(v.bit_count(),n-1-v.bit_count()) for v in a);la=cut(a)
  require(delta==k and la==1, f"Sharp witness failed: k={k}")
  sharp.append({'k':k,'n':n,'min_semidegree':delta,'lambda':la,'adjacency_outmasks':a})
 print(json.dumps({'exhaustive':res,'attachment_checks':checks,'sharp_examples':sharp},indent=2))
if __name__=='__main__':run()
