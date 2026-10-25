import itertools,json,math
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from edge_actions import pairs,relabel,require

@lru_cache(None)
def info(n,x):
 ps=pairs(n);total=0;even=True;aut=0
 for p in itertools.permutations(range(n)):
  eps=(-1)**sum(p[i]>p[j] for k,(i,j) in enumerate(ps) if x>>k&1)
  total+=eps
  if relabel(x,n,p)==x:
   aut+=1
   if eps<0:even=False
 return even,F(total,math.factorial(n)),aut

def comps(n,x):
 ps=pairs(n);adj=[set() for _ in range(n)]
 for k,(i,j) in enumerate(ps):
  if x>>k&1:adj[i].add(j);adj[j].add(i)
 out=[];todo=set(range(n))
 while todo:
  start=min(todo);seen={start};stack=[start]
  while stack:
   for v in adj[stack.pop()]-seen:seen.add(v);stack.append(v)
  todo-=seen;vertices=sorted(seen);c=0
  for k,(i,j) in enumerate(pairs(len(vertices))):
   if vertices[j] in adj[vertices[i]]:c|=1<<k
  out.append((len(vertices),c))
 return out

def mul(a,b,N):return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N+1)]
def euler(exponents,N,sign=1):
 # Each connected type contributes1/(1-sign*z^n).
 a=[1]+[0]*N
 for n,k in enumerate(exponents):
  if n==0:continue
  for _ in range(k):
   factor=[0]*(N+1)
   for j in range(N//n+1):factor[n*j]=sign**j
   a=mul(a,factor,N)
 return a

def main():
 N=6;connected={'all':[0]*(N+1),'even_edges':[0]*(N+1),'odd_edges':[0]*(N+1),'good':[0]*(N+1),'zero_even':[0]*(N+1)}
 counts={'t':[1],'c':[1],'r':[1]};graphs=[];records=[]
 for n in range(1,N+1):
  seen=set();t=c=r=0;gcount=0
  for x in range(1<<(n*(n-1)//2)):
   if x in seen:continue
   orbit={relabel(x,n,p) for p in itertools.permutations(range(n))};seen.update(orbit);gcount+=1
   even,alpha,aut=info(n,x);cs=comps(n,x)
   require(even==all(info(k,y)[0] for k,y in cs),'componentwise signed automorphism')
   prod=F(1)
   for k,y in cs:prod*=info(k,y)[1]
   require(alpha==prod,'alpha multiplicativity')
   if not even:continue
   t+=1;par=x.bit_count()%2
   if par==0:c+=1
   if alpha:
    require(par==0,'nonzero alpha forces total even edges');r+=1
    require(all(y.bit_count()%2==0 for k,y in cs),'nonzero components evenedge')
   if len(cs)==1:
    connected['all'][n]+=1;connected['odd_edges' if par else 'even_edges'][n]+=1
    if alpha:connected['good'][n]+=1
    elif par==0:connected['zero_even'][n]+=1
   records.append({'n':n,'mask':x,'alpha':str(alpha),'edges':x.bit_count(),'components':cs,'aut':aut})
  graphs.append(gcount);counts['t'].append(t);counts['c'].append(c);counts['r'].append(r)
 T=euler(connected['all'],N);A=euler(connected['even_edges'],N);R=euler(connected['good'],N);Z=euler(connected['zero_even'],N)
 Oplus=euler(connected['odd_edges'],N);Ominus=euler(connected['odd_edges'],N,-1)
 Oeven=[(a+b)//2 for a,b in zip(Oplus,Ominus)];D=mul(A,Ominus,N);C=[(t+d)//2 for t,d in zip(T,D)]
 require(T==counts['t'] and C==counts['c'] and R==counts['r'],'Euler coefficient direct counts')
 require(C==mul(mul(R,Z,N),Oeven,N),'rank kernel product factorization')
 K=[c-r for c,r in zip(C,R)];Kfactor=mul(R,[a-(i==0) for i,a in enumerate(mul(Z,Oeven,N))],N)
 require(K==Kfactor,'kernel series')
 require((T,C,R)==([1,1,1,2,4,12,56],[1,1,1,2,3,10,34],[1,1,1,2,3,10,32]),'known exact finite controls')
 require(connected['good']==[0,1,0,1,1,7,21],'good connected counts')
 require(K[6]>=2*C[0],'two-mode n6 lower bound boundary')
 for x in [255,2015]:
  require(info(6,x)[0] and info(6,x)[1]==0 and x.bit_count()%2==0 and len(comps(6,x))==1,'connected zero basis generator')
 require((255).bit_count()!=(2015).bit_count(),'distinct component types')
 # Eight-vertex two-claw hostile: no need enumerate8!; component theorem gives alpha0.
 claw=7 # on4vertices:01,02,03
 require(info(4,claw)[0] and info(4,claw)[1]==0,'oddedge Royle-even claw')
 out={'max_n':N,'graph_type_counts_n1_to_n6':graphs,'connected':connected,'T':T,'D':D,'C':C,'R':R,'K':K,'Z':Z,'O_even':Oeven,'records':records}
 Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
if __name__=='__main__':main()
