"""Bounded independent algebra and frozen-state arithmetic; no native transitions."""
import json, hashlib, itertools, pathlib, sys
root=pathlib.Path(sys.argv[1]); count=0

def need(x, label):
 global count
 count+=1
 if not x: raise ValueError(label)
def q(v): return v[0]*v[0]-v[0]*v[1]+v[1]*v[1]
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def d(v): return q(add(v[:2],v[2:]))-q(v[:2])
def k(v,w): return d(add(v,w))-d(v)-d(w)
def beta(v,w): return q(add(v,w))-q(v)-q(w)
def proj(v,p,r): return tuple(2*p*p*v[i]+(p*p+r*r)*v[i+2] for i in (0,1)),tuple(2*p*p*v[i]+(p*p-r*r)*v[i+2] for i in (0,1))
# A homogeneous quadratic identity is determined by basis and pair-sum values.
basis=[tuple(int(i==j) for i in range(4)) for j in range(4)]
for p,r in [(1,1),(1,2),(2,1),(3,7),(11,2)]:
 for v in basis+[add(x,y) for x,y in itertools.combinations(basis,2)]:
  g,f=proj(v,p,r); need(q(g)-q(f)==(2*p*r)**2*d(v),'quadratic coefficients')
 for v,w in itertools.product(basis,repeat=2):
  g,f=proj(v,p,r); gg,ff=proj(w,p,r)
  need(beta(g,gg)-beta(f,ff)==(2*p*r)**2*k(v,w),'bilinear coefficients')
# Independent vertex construction: support normals perpendicular to each edge,
# with both infinitesimal tangential perturbations. No angular sorting.
def maxnorm(vs):
 if not any(any(v) for v in vs): return 0
 best=0
 for a,b in set(vs):
  if a==b==0: continue
  for side in [-1,1]:
   x=y=0
   for c,e in vs:
    z=-b*c+a*e
    if z==0: z=side*(a*c+b*e)
    s=1 if z>=0 else -1
    x+=s*c; y+=s*e
   best=max(best,q((x,y)))
 return best
rows=json.loads((root/'05-knowledge/results/future_tail_blocks_20261021b/receipt.json').read_text())['cases']
out=[]
for row in rows:
 ports=[tuple(r[1:]) for r in row['arrived_ports']]+[(0,0,0,0)]
 vs=[tuple(r[1:5]) for r in row['pending']]
 need(all(v[:2]==(0,0) for v in vs),'zero aggregate old channel')
 base=sum(d(v) for v in ports)
 linear=sum(max(abs(k(c,v)) for c in ports) for v in vs)
 angular=maxnorm([v[2:] for v in vs]); bound=base+linear+angular
 expected={7:1153,4096:551844,8192:1118881,32768:4142214}[row['N']]
 need(bound==expected,'independent linked lambda1 bound')
 u=base+sum(max(0,d(v)+max(abs(k(c,v)) for c in ports)) for v in vs)+sum(abs(k(v,w)) for v,w in itertools.combinations(vs,2))
 need(u==row['envelope'][-1],'independent absolute envelope')
 out.append({'N':row['N'],'D':row['D'],'base':base,'linear':linear,'angular':angular,'lambda1_bound':bound,'U':u,'earlier_cut_predicate_cited':row['early_cuts_pass']})
print(json.dumps({'status':'PASS; independent algebra and frozen-state arithmetic only','checks':count,'cases':out},sort_keys=True,indent=2))
