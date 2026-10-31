"""Independent forward fiber closure, linear operators, modular no-go patterns."""
from pathlib import Path
from collections import deque
from itertools import product
import json,sys,hashlib
SOURCE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(SOURCE));import digital as d
raw=json.loads((SOURCE/'actual_seed_roots.json').read_text());roots=[(r['state'][0],r['state'][1],r['state'][2],r['state'][0]-r['state'][3]) for r in raw]
def step(s,b):
 h,F,q,g=s
 if b==0:return h,F+q+g,2*q,0
 if b==1:return h,F,2*q+h,h
 j=int(b[1]);return 3*h,3*F+(9-3*j)*q+(3-j)*g+(2*h if j in (1,2) else 0),12*q+3*j*h,j*h

def p(h,q):
 if q==0 or q%2 and q<h:return 0
 if q%2:return p(h,(q-h)//2)
 return p(h,q//2)+q//2+h*((q//2)%2)
fibers=0;boundaries=[]
expected=json.loads((SOURCE/'fiber_boundaries.json').read_text())
for h in sorted({s[0] for s in roots}):
 sr=[s for s in roots if s[0]==h];L=max(2*s[2]+h for s in sr);limit=2*L+h
 initial={step(s,b) for s in sr for b in (0,1)};seen=set(initial);queue=deque(initial);values={}
 while queue:
  s=queue.popleft();values.setdefault(s[2],set()).add(s[1])
  for b in (0,1):
   z=step(s,b)
   if z[2]<=limit and z not in seen:seen.add(z);queue.append(z)
 for q in range(limit+1):
  if values.get(q,set())!=d.fiber(h,q,roots):raise RuntimeError(('forward fiber',h,q))
  if p(h,q)!=d.phi(h,q):raise RuntimeError('phi')
  fibers+=1
 boundary=[]
 for q in range(L+1,limit+1):
  parent=(q-h)//2 if q%2 else q//2
  if parent<=L and values.get(q):boundary.append({'q':q,'offsets':sorted(x-p(h,q) for x in values[q])})
 item=next(r for r in expected if r['h']==h)
 if boundary!=item['boundary_nodes'] or L!=item['maximum_source_q']:raise RuntimeError('frozen boundary')
 global_offsets=sorted({v[1]-p(h,v[2]) for v in initial})
 signatures={tuple(x['offsets']) for x in boundary}
 def consecutive(xs):
  best=run=0;last=None
  for x in xs:
   run=run+1 if last is not None and x==last+2 else 1;best=max(best,run);last=x
  return best
 if item['global_offsets']!=global_offsets or item['eventual_boundary_nodes']!=len(boundary) or item['eventual_distinct_fibers']!=len(signatures) or item['maximum_eventual_fiber_size']!=max(map(len,signatures)) or item['maximum_eventual_consecutive_odd_width']!=max(map(consecutive,signatures)):raise RuntimeError('fiber summary')
 # Each listed cylinder's descendants introduce no new offset beyond its ancestor.
 for node in boundary:
  q=node['q'];K=node['offsets']
  for bits in product((0,1),repeat=5):
   t=q
   for b in bits:t=2*t+h*b
   if sorted(F-p(h,t) for F in d.fiber(h,t,roots))!=K:raise RuntimeError('frozen cylinder')
 boundaries.append(dict(h=h,nodes=len(boundary),maximum_fiber=max(len(v['offsets']) for v in boundary)))
linear=0
for j in range(4):
 for i in range(4):
  s=tuple(int(k==i) for k in range(4));h,F,q,g=s
  shadow=step(step(s,j//2),j%2);actual=step(s,'C'+str(j));delta=tuple(a-3*b for a,b in zip(actual,shadow))
  want=((0,0,0,0),(0,3*q-g+2*h,0,-2*h),(0,-3*q+g-4*h,0,2*h),(0,0,0,0))[j]
  if delta!=want:raise RuntimeError('operator defect')
  linear+=1
# Matrix mod3 proofs on all basis vectors: constants kill; C-S-C kills output;
# three consecutive nonconstant C kill whole vector. Zero h,q,g remains zero under S.
modular=0
for i in range(4):
 e=tuple(int(k==i) for k in range(4))
 for j in (0,3):
  if any(v%3 for v in step(e,'C'+str(j))):raise RuntimeError('constant triangle')
  modular+=1
 for a,b in product((1,2),repeat=2):
  for bit in (0,1):
   z=step(step(step(e,'C'+str(a)),bit),'C'+str(b))
   if any(v%3 for v in z):raise RuntimeError('separated triangles')
   modular+=1
 for word in product((1,2),repeat=3):
  z=e
  for j in word:z=step(z,'C'+str(j))
  if any(v%3 for v in z):raise RuntimeError('three triangles')
  modular+=1
 for b in (0,1):
  z=step((0,1,0,0),b)
  if z!=(0,1,0,0):raise RuntimeError('zero envelope')
print(json.dumps(dict(status='PASS',forward_fiber_nodes=fibers,boundary_families=boundaries,linear_basis_identities=linear,mod3_basis_products=modular,source_sha256=hashlib.sha256((SOURCE/'digital.py').read_bytes()).hexdigest()),sort_keys=True))
