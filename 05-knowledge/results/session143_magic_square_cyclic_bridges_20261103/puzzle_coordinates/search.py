from itertools import product
from math import lcm
from time import monotonic
pts=[(1,a,b) for a in range(3) for b in range(3)]+[(0,1,a) for a in range(3)]+[(0,0,1)]
N=len(pts); ident=tuple(range(N))
def dot(p,q):return sum(x*y for x,y in zip(p,q))%3
def line(a,b):
 x,y=pts[a],pts[b]
 n=((x[1]*y[2]-x[2]*y[1])%3,(x[2]*y[0]-x[0]*y[2])%3,(x[0]*y[1]-x[1]*y[0])%3)
 return tuple(i for i,p in enumerate(pts) if dot(n,p)==0)
def move(a,b):
 if a==b:return ident
 c,d=[x for x in line(a,b) if x not in [a,b]]
 p=list(ident);p[a],p[b]=b,a;p[c],p[d]=d,c;return tuple(p)
moves={(a,b):move(a,b) for a in range(N) for b in range(N)}
def comp(p,q):return tuple(p[x] for x in q)
def power(p,k):
 v=ident
 for _ in range(k):v=comp(p,v)
 return v
def cycles(p):
 seen=set();out=[]
 for i in range(N):
  if i not in seen:
   c=[];j=i
   while j not in seen:c.append(j);seen.add(j);j=p[j]
   out.append(tuple(c))
 return tuple(out)
def order(p):return lcm(*(len(c) for c in cycles(p)))
def permword(path):
 p=ident
 for a,b in zip(path,path[1:]):p=comp(moves[a,b],p)
 return p
g=None
for k in range(1,5):
 for inner in product(range(1,N),repeat=k):
  path=(0,)+inner+(0,)
  if any(a==b for a,b in zip(path,path[1:])):continue
  p=permword(path)
  if order(p)==5:g=p;gp=path;break
 if g:break
print('gpath',gp,'gcycles',cycles(g),'g',g,flush=True)
g2=power(g,2)
# Equivalent h*g=g^2*h avoids constructing the inverse.
start=monotonic();tested=0;h=None
for k in range(1,6):
 for inner in product(range(1,N),repeat=k):
  path=(0,)+inner+(0,)
  if any(a==b for a,b in zip(path,path[1:])):continue
  p=permword(path);tested+=1
  if comp(p,g)==comp(g2,p) and order(p)==4:h=p;hp=path;break
 if h:break
 print('conjugator exhausted innerlength',k,'tested',tested,'seconds',monotonic()-start,flush=True)
print('hpath',hp if h else None,'hcycles',cycles(h) if h else None,'tested',tested,'seconds',monotonic()-start,flush=True)
import json
open(str(__import__('pathlib').Path(__file__).with_name('search_result.json')),'w').write(json.dumps({'points':pts,'g_word':gp,'g_permutation':g,'g_cycles':cycles(g),'h_word':hp if h else None,'h_permutation':h,'h_cycles':cycles(h) if h else None,'conjugator_tested':tested,'search_bound_inner_vertices':5},indent=2)+'\n')
