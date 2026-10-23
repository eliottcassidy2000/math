#!/usr/bin/env python3
"""Exact finite controls, standard library only. No claims from bare assertion."""
import itertools as it, math, json, random
checks=0
stats={}
def need(ok,why):
 global checks
 checks+=1
 if not ok: raise RuntimeError(why)
def shift(B,m,n,a,b):
 return tuple(B[((i-a)%m)*n+(j-b)%n] for i in range(m) for j in range(n))
def hnf(H,m,n):
 h=next(a for a in range(1,m+1) if (a%m,0) in H)
 v=next(b for b in range(1,n+1) if any((a,b%n) in H for a in range(m)))
 s=next(a for a in range(h) if (a,v%n) in H)
 return h,s,v
def member(a,b,h,s,v):return b%v==0 and (v*a-s*b)%(h*v)==0
rng=random.Random(42001)
blocks=0
for m,n in [(1,1),(1,2),(2,2),(2,3),(3,2),(3,3),(3,4),(4,3),(4,4)]:
 masks=range(1<<(m*n)) if m*n<=12 else sorted(set([0,(1<<16)-1]+[rng.randrange(1<<16) for _ in range(1024)]))
 for mask in masks:
  B=tuple(mask>>k&1 for k in range(m*n));images={(a,b):shift(B,m,n,a,b) for a in range(m) for b in range(n)}
  H={p for p,C in images.items() if C==B};h,s,v=hnf(H,m,n)
  need(h>0 and v>0 and 0<=s<h and m%h==0 and n%v==0 and (s*n//v)%h==0,'period conditions')
  need(len(H)==m*n//(h*v),'HNF count')
  need(sum(b==0 for a,b in H)==m//h,'horizontal kernel')
  need(sum(a==0 for a,b in H)==n*math.gcd(h,s)//(h*v),'vertical kernel')
  for (a0,b0),target in images.items():
   actual={p for p,C in images.items() if C==target}
   predicted={(a,b) for a in range(m) for b in range(n) if member(a-a0,b-b0,h,s,v)}
   need(actual==predicted,'offset relation')
  blocks+=1
stats['binary_blocks']=blocks
# Actual regular nine-cell source; translation x+1 respects the hidden modules.
def arc(x,y):
 a,u=x%3,x//3;b,v=y%3,y//3
 return x!=y and ((b-a)%3==1 if a!=b else (v-u)%3==1)
need(all(sum(arc(x,y) for y in range(9))==4 for x in range(9)),'R9 regular')
need(all(arc(x,y)==arc((x+1)%9,(y+1)%9) for x in range(9) for y in range(9)),'C9 actual automorphism')
B=tuple(int(x%3==y%3) for x in range(9) for y in range(9))
H={(a,b) for a in range(9) for b in range(9) if shift(B,9,9,a,b)==B}
need(hnf(H,9,9)==(3,1,1) and len(H)==27,'R9 cyclic kernel')
# Enumerate all automorphisms with outer rotation and three independent inner rotations.
autos=[]
for a in range(3):
 for u in it.product(range(3),repeat=3):
  p=tuple((x%3+a)%3+3*((x//3+u[x%3])%3) for x in range(9))
  need(all(arc(x,y)==arc(p[x],p[y]) for x in range(9) for y in range(9)),'R9 full chart family')
  autos.append(p)
pairs=sum(all(B[x*9+y]==B[p[x]*9+q[y]] for x in range(9) for y in range(9)) for p in autos for q in autos)
need(pairs==2187,'full relation count')
# Noncyclic C3^2 chart enumeration: graph of digit addition is not affine cyclic.
f=lambda x:(x%3+1)%3+3*(x//3)
R={(f(x),x) for x in range(9)}
o=next(iter(R));difference={((a-o[0])%9,(b-o[1])%9) for a,b in R}
closed=all(((a+c)%9,(b+d)%9) in difference for a,b in difference for c,d in difference)
need(not closed,'noncyclic hostile must fail coset closure')
need((f(0),f(1),f(2))==(1,2,0),'noncyclic explicit hostile')
stats['R9_full_chart_pairs']=pairs
# Random actual blocks translated under shared cyclic charts. Half include independently
# generated targets, so incompatibility is real and not filtered out.
networks=0;inconsistent=0
for k in range(500):
 ms=rng.choice([(3,3,3),(3,4,3),(3,6,3),(4,4,4)])
 edges=[(0,1),(1,2),(0,2)];reference=[rng.randrange(m) for m in ms]
 relations=[];equations=[]
 for i,j in edges:
  m,n=ms[i],ms[j];block=tuple(rng.randrange(2) for _ in range(m*n))
  # Include many nontrivial kernels rather than nearly always rigid random blocks.
  if k%3==0:block=tuple(int(x%math.gcd(m,n)==y%math.gcd(m,n)) for x in range(m) for y in range(n))
  a,b=(reference[i],reference[j]) if k%2==0 else (rng.randrange(m),rng.randrange(n))
  target=shift(block,m,n,a,b)
  H={(x,y) for x in range(m) for y in range(n) if shift(block,m,n,x,y)==block}
  relations.append((i,j,{(x,y) for x in range(m) for y in range(n) if shift(block,m,n,x,y)==target}))
  equations.append((i,j,a,b,*hnf(H,m,n)))
 direct={x for x in it.product(*(range(m) for m in ms)) if all((x[i],x[j]) in R for i,j,R in relations)}
 lattice={x for x in it.product(*(range(m) for m in ms)) if all(member(x[i]-a,x[j]-b,h,s,v) for i,j,a,b,h,s,v in equations)}
 need(direct==lattice,'shared network congruence')
 if direct:
  x0=next(iter(direct))
  kernel={x for x in it.product(*(range(m) for m in ms)) if all(member(x[i],x[j],h,s,v) for i,j,a,b,h,s,v in equations)}
  need(direct=={tuple((x0[i]+x[i])%ms[i] for i in range(3)) for x in kernel},'global coset')
 else:inconsistent+=1
 networks+=1
stats.update(networks=networks,inconsistent_networks=inconsistent)
# Composite C9 network: comparison count and joint class quotient via explicit block states.
for edges in [[(0,1),(1,2)],[(0,1),(1,2),(0,2)],[(0,1),(1,2),(2,3),(3,0)]]:
 v=1+max(max(e) for e in edges)
 kernel=sum(all((x[i]-x[j])%3==0 for i,j in edges) for x in it.product(range(9),repeat=v))
 need(kernel==3**(v+1),'composite kernel size')
 image={tuple((x[i]-x[j])%3 for i,j in edges) for x in it.product(range(3),repeat=v)}
 need(3**len(edges)//len(image)==3**(len(edges)-v+1),'composite alignment quotient')
# Visible-kernel base lower bounds and subgroup-intersection hostile.
identity=tuple(range(9))
full_kernel=[p for p in autos if all(p[x]%3==x%3 for x in range(9))]
cyclic_kernel=[tuple((x+a)%9 for x in range(9)) for a in (0,3,6)]
def base_size(group):
 for k in range(10):
  for marks in it.combinations(range(9),k):
   if all(p==identity or any(p[x]!=x for x in marks) for p in group):return k
need(len(full_kernel)==27 and base_size(full_kernel)==3,'full visible kernel base')
need(len(cyclic_kernel)==3 and base_size(cyclic_kernel)==1,'cyclic visible kernel base')
K1={(a,0) for a in range(3)};K2={(0,b) for b in range(3)}
need(len(K1)==len(K2)==3 and K1&K2=={(0,0)},'kernel sizes lose intersection')
stats['visible_kernel_base_sizes']={'full_R9':3,'cyclic_R9':1}
stats['checks']=checks
print(json.dumps({'status':'PASS','scope':'supplied cyclic/full-group chart carriers, not unrestricted native wiring','tests':stats},indent=2))
