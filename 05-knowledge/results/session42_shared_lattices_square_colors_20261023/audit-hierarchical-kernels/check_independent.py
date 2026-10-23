"""Independent lexicographic tournament, generators, kernels and owner bases."""
from itertools import product,combinations
import json,math
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def arc(u,v):
 for a,b in zip(u,v):
  if a!=b:return (b-a)%3==1
 return False
def generator(words,prefix):
 d=len(prefix)
 return {w:w[:d]+((w[d]+1)%3,)+w[d+1:] if w[:d]==prefix else w for w in words}
reports=[]
for k in range(1,5):
 words=list(product(range(3),repeat=k));q=len(words)
 need(all(sum(arc(u,v) for v in words)==(q-1)//2 for u in words),'internal regularity')
 gens=[]
 for d in range(k):
  for prefix in product(range(3),repeat=d):
   g=generator(words,prefix);gens.append((d,g))
   need(all(arc(u,v)==arc(g[u],g[v]) for u in words for v in words),'generator automorphism')
 need(len(gens)==(q-1)//2,'tree generator nodes')
 for l in range(k+1):
  kernel_gens=[g for d,g in gens if d>=l]
  need(len(kernel_gens)==(3**k-3**l)//2,'kernel internal-node exponent')
  for g in kernel_gens:
   need(all((u[:l]==v[:l])==(g[u][:l]==v[:l]) for u in words for v in words),'one-sided actual block preservation')
  if l<k:
   marks={w for w in words if w[-1]==0}
   need(len(marks)==3**(k-1),'one mark per bottom triple')
   need(all(any(g[w]!=w for w in marks) for g in kernel_gens),'each nontrivial node generator meets marks')
  reports.append({'k':k,'l':l,'vertices':q,'declared_kernel_order':str(3**((3**k-3**l)//2)),'kernel_base':0 if l==k else 3**(k-1),'informative':l>0})
# Complete W2 family: outer cyclic rotation and three independent child rotations.
words=list(product(range(3),repeat=2));autos=[]
for a in range(3):
 for bs in product(range(3),repeat=3):
  autos.append({(u,v):((u+a)%3,(v+bs[u])%3) for u,v in words})
need(len({tuple(g[w] for w in words) for g in autos})==81,'W2 size')
for l,wanted in [(0,81),(1,27),(2,1)]:
 kernel=[g for g in autos if all(u[:l]==g[u][:l] for u in words)]
 need(len(kernel)==wanted,'exact W2 kernel')
 found=None
 for b in range(4):
  for marks in combinations(words,b):
   survivors=[g for g in kernel if all(g[w]==w for w in marks)]
   if len(survivors)==1:found=b;break
  if found is not None:break
 need(found==(0 if l==2 else 3),'exact W2 base')
# Complete outer C3 network, three named cells with equality/complement crossblocks.
for k in (2,3):
 words=list(product(range(3),repeat=k));q=len(words)
 for l in range(1,k):
  vertices=list(product(range(3),words))
  def global_arc(x,y):
   i,u=x;j,v=y
   if i==j:return arc(u,v)
   equal=u[:l]==v[:l]
   return equal if (j-i)%3==1 else not equal
  need(all(sum(global_arc(x,y) for y in vertices)==(3*q-1)//2 for x in vertices),'global regularity')
  for i in range(3):
   for d in range(l,k):
    for prefix in product(range(3),repeat=d):
     g=generator(words,prefix)
     def image(x):return (i,g[x[1]]) if x[0]==i else x
     need(all(global_arc(x,y)==global_arc(image(x),image(y)) for x in vertices for y in vertices),'independent actual network kernel')
print(json.dumps({'status':'PASS','checks':checks,'scope':'declared tree-rotation subgroups; full unrestricted automorphism group not enumerated','kernel_controls':reports,'exact_W2_automorphisms_in_declared_family':81},indent=2))
