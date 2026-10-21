import itertools,json
from pathlib import Path

def req(x,msg):
 if not x: raise RuntimeError(msg)

def rotate(t): return t[1:]+t[:1]
def canon(t): return min(t,rotate(t),rotate(rotate(t)))
rows=[]
for q in range(1,10):
 for f in range(q%2,q+1,2):
  invol=list(range(q))
  for j in range(f,q,2): invol[j],invol[j+1]=j+1,j
  reps={canon(t) for t in itertools.product(range(q),repeat=3) if len(set(t))>1}
  fixed=sum(canon(tuple(invol[x] for x in reversed(t)))==t for t in reps)
  req(len(reps)==(q**3-q)//3,'q recursion failed')
  req(fixed==f*(q-1),'f recursion failed')
  rows.append({'q':q,'f':f,'q_next':len(reps),'f_next':fixed})
# Exhaust all cyclic-invariant F3-subspaces in dimensions 3 and 5,
# then test the necessary marker-coverage containment for all reflections.
def add(x,y): return tuple((a+b)%3 for a,b in zip(x,y))
def scale(a,x): return tuple(a*b%3 for b in x)
def extend(K,v): return frozenset(add(k,scale(a,v)) for k in K for a in range(3))
linear=[]
for m in (3,5):
 zero=(0,)*m; V=list(itertools.product(range(3),repeat=m)); init=frozenset([zero]); spaces={init}; todo=[init]
 while todo:
  K=todo.pop()
  for v in V:
   if v in K: continue
   L=K; w=v
   for j in range(m): L=extend(L,w); w=rotate(w)
   if L not in spaces: spaces.add(L); todo.append(L)
 witnesses=[]
 for K in spaces:
  if len(K)==3**m: continue
  for shift in range(m):
   # Reflection j maps to shift-j; its fixed locus always contains one basis vector.
   image={tuple((x[j]+x[(shift-j)%m])%3 for j in range(m)) for x in V}
   req(not image<=K,'Proper transitive cyclic invariant kernel covered every marker')
   witnesses.append({'kernel_size':len(K),'reflection':shift,'outside_image':list(min(image-K))})
 linear.append({'m':m,'invariant_subspaces':len(spaces),'proper_kernel_reflection_tests':len(witnesses),'witnesses':witnesses})
# Removing transitivity genuinely breaks the kernel implication: reflection swaps
# coordinates 0,1 and fixes 2; its image is a proper nontransitive-invariant plane.
boundary=frozenset((a,a,b) for a in range(3) for b in range(3))
image=frozenset(((x[0]+x[1])%3,(x[0]+x[1])%3,2*x[2]%3) for x in itertools.product(range(3),repeat=3))
req(image==boundary and len(boundary)==9,'Transitivity hostile failed')
Path(__file__).with_name('hostile_receipt.json').write_text(json.dumps({'status':'PASS','palette_recursions':rows,'linear_marker_controls':linear,'without_transitivity_hostile':{'dimension':3,'proper_kernel_size':9,'image_contained':True}},indent=2)+'\n')
print('PASS',len(rows),'palette recursions;',[(r['m'],r['invariant_subspaces']) for r in linear])
