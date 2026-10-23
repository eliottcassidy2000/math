import itertools as I,math,json,hashlib,pathlib

def need(b,m):
 if not b:raise RuntimeError(m)
here=pathlib.Path(__file__).resolve().parent
source=next((p for p in (here.parent/'session43-native-ports',here.parent/'native-ports') if (p/'MANIFEST.json').is_file()),None)
need(source is not None,'cannot locate native-ports sibling packet')
manifest=json.loads((source/'MANIFEST.json').read_text())
for row in manifest['files']:
 b=(source/row['path']).read_bytes();need(len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],'manifest '+row['path'])
# Independent adjacency-set enumeration, no imports from the subject checker.
censuses=[]
for n in [3,4]:
 edges=list(I.combinations(range(n),2));perms=list(I.permutations(range(n)));classes={}
 for bits in I.product([0,1],repeat=len(edges)):
  arcs=frozenset((b,a) if bit else (a,b) for (a,b),bit in zip(edges,bits))
  key=min(tuple(int((p[a],p[b]) not in arcs) for a,b in edges) for p in perms)
  if key in classes:continue
  paths=[p for p in perms if all((a,b) in arcs for a,b in zip(p,p[1:]))]
  aut=[p for p in perms if frozenset((p[a],p[b]) for a,b in arcs)==arcs]
  orbit={tuple(int((p[a],p[b]) not in arcs) for a,b in edges if b>a+1) for p in paths}
  need(len(paths)==len(aut)*len(orbit),'free observer fibre')
  classes[key]=(len(paths),len(aut),len(orbit))
 need(sum(v[2] for v in classes.values())==2**((n-1)*(n-2)//2),'observer fibre sum')
 censuses.append({'n':n,'classes':len(classes),'fibre_triples':sorted(classes.values())})
# Independent fixed-point and two-involution budgets through n7.
budgets=[]
for n in range(3,8):
 E=[(a,b) for a,b in I.combinations(range(n),2) if b>a+1];m=len(E);mask=(1<<m)-1
 perm=[E.index((n-1-b,n-1-a)) for a,b in E]
 fixed=0;blue=0;black=0;orbits=set()
 for x in range(1<<m):
  r=sum(((x>>k)&1)<<perm[k] for k in range(m));need(r!=x^mask,'FR has no fixed points')
  fixed+=x==r
  if x<(x^mask):
   if r==x:blue+=1
   else:black+=1
  orbits.add(min(x,r,x^mask,r^mask))
 expected=2**((m+(n-1)//2)//2)
 need(fixed==expected and blue==fixed//2 and black==((1<<m)-fixed)//2 and len(orbits)==((1<<m)+fixed)//4,'budgets')
 budgets.append({'n':n,'fixed_R':fixed,'blue':blue,'black':black,'FR_orbits':len(orbits)})
# Actual six-vertex tournament, independent arcs as set.
arcs={(i*3+u,i*3+(u+1)%3) for i in range(2) for u in range(3)}
for u,v in I.product(range(3),repeat=2):arcs.add((u,v+3) if (v-u)%3 else (v+3,u))
paths=[p for p in I.permutations(range(6)) if all(e in arcs for e in zip(p,p[1:]))]
valid=[]
for a,b in I.product(range(3),repeat=2):
 p=tuple((a+j)%3 for j in range(3))+tuple(3+(b+j)%3 for j in range(3))
 if all(e in arcs for e in zip(p,p[1:])):valid.append((a,b))
need(len(valid)==6 and set((b-a)%3 for a,b in valid)=={0,1},'six path relation')
need([sum(a==v for a,b in arcs) for v in range(6)]==[3]*3+[2]*3,'almost regular')
# All 216 palette triples with three actual nine-vertex spanning templates.
palettes=[frozenset(a) for r in [1,2] for a in I.combinations(range(3),r)]
templates=[tuple((i,u) for u in range(3) for i in range(3)),tuple((i,u) for i in range(3) for u in range(3)),((0,0),(1,0),(2,0),(0,1),(2,1),(1,1),(0,2),(1,2),(2,2))]
trials=0;valid_total=0;rejected_branches=0
for Ds in I.product(palettes,repeat=3):
 blocks=dict(zip([(0,1),(0,2),(1,2)],Ds))
 def arc(i,u,j,v):
  if i==j:return (v-u)%3==1
  if i<j:return (v-u)%3 in blocks[i,j]
  return (u-v)%3 not in blocks[j,i]
 for temp in templates:
  actual=[]
  for aa in I.product(range(3),repeat=3):
   path=[(i,(u+aa[i])%3) for i,u in temp]
   if all(arc(i,u,j,v) for (i,u),(j,v) in zip(path,path[1:])):actual.append(aa)
  # Independently integrate a_0=0: enumerate two potentials, not edge-branch tuples.
  zero_root=[]
  for b,c in I.product(range(3),repeat=2):
   aa=(0,b,c)
   if all(arc(i,(u+aa[i])%3,j,(v+aa[j])%3) for (i,u),(j,v) in zip(temp,temp[1:])):zero_root.append(aa)
  need(len(actual)==3*len(zero_root),'connected kernel lift')
  # Record every edge difference; all actual assignments must have zero signed cycle sums.
  signatures={tuple((aa[j]-aa[i])%3 for (i,u),(j,v) in zip(temp,temp[1:]) if i!=j) for aa in actual}
  need(len(actual)==3*len(signatures),'unique difference branch with three lifts')
  occurrences=[(i,u,j,v) for (i,u),(j,v) in zip(temp,temp[1:]) if i!=j]
  internal=all(arc(i,u,j,v) for (i,u),(j,v) in zip(temp,temp[1:]) if i==j)
  possible_signatures={tuple((aa[j]-aa[i])%3 for i,u,j,v in occurrences) for aa in I.product(range(3),repeat=3)}
  choices=[[d for d in range(3) if arc(i,u,j,(v+d)%3)] for i,u,j,v in occurrences]
  feasible=0
  for sig in I.product(*choices):
   if sig in possible_signatures:feasible+=1
   else:rejected_branches+=1
  need(len(actual)==(3*feasible if internal else 0),'all branch choices versus shared phase feasibility')
  trials+=1;valid_total+=len(actual)
# General p^c, including composite odd modulus: fixed difference branches, forest plus cycle.
forest=[]
for p in [3,5,9]:
 for edges,c in [([(0,1,1)],3), ([(0,1,1),(1,2,2),(2,0,-3)],2), ([(0,1,1),(1,2,2),(2,3,1)],1)]:
  good=[aa for aa in I.product(range(p),repeat=4) if all((aa[j]-aa[i]-d)%p==0 for i,j,d in edges)]
  need(len(good)==p**c,'general component lifts')
  # One actual pinned owner in first component removes exactly one factor p.
  need(sum(aa[0]==0 for aa in good)==p**(c-1),'pinned owner lift')
  forest.append({'p':p,'components':c,'lifts':len(good)})
print(json.dumps({'status':'PASS','manifest_hashes_verified':len(manifest['files']),'censuses':censuses,'involution_budgets':budgets,'six_vertex_total_H':len(paths),'six_vertex_valid_charts':valid,'nine_vertex_palette_template_cases':trials,'nine_vertex_accepted_charts':valid_total,'independently_rejected_cycle_branches':rejected_branches,'composite_and_partial_forest_controls':forest},indent=2,sort_keys=True))
