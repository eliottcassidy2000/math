from itertools import product
import json

def req(c,m):
    if not c:raise RuntimeError(m)
# Independent exact hyperplane check.
vectors=list(product(range(3),repeat=4)); normals=[]
for v in vectors:
 if not any(v):continue
 if next(x for x in v if x)==1:normals.append(v)
req(len(normals)==40,'normalization')
for a in normals:
 K=[v for v in vectors if sum(x*y for x,y in zip(a,v))%3==0]
 req(len(K)==27 and any(sum(x!=0 for x in v)>=3 for v in K),'hyperplane')
# All subgroups of the standard Sylow3 subgroup C3 wr C3 on9 points.
perms=[]
for a in product(range(3),repeat=3):
 for t in range(3):perms.append(tuple(3*((i+t)%3)+(j+a[i])%3 for i in range(3) for j in range(3)))
lookup={g:i for i,g in enumerate(perms)};one=lookup[tuple(range(9))]
mul=[[lookup[tuple(h[g[x]] for x in range(9))] for h in perms] for g in perms]
def close(gens):
 S={one};todo=[one]
 for x in todo:
  for y in gens:
   z=mul[x][y]
   if z not in S:S.add(z);todo.append(z)
 return frozenset(S)
subgroups={frozenset([one])};todo=list(subgroups)
for H in todo:
 for g in range(81):
  if g in H:continue
  K=close(list(H)+[g])
  if K not in subgroups:subgroups.add(K);todo.append(K)
z=lookup[tuple(3*i+(j+1)%3 for i in range(3) for j in range(3))]
transitive=cyclic=elementary=0
for H in subgroups:
 if {perms[g][0] for g in H}!=set(range(9)):continue
 transitive+=1
 has9=False
 for g in H:
  p=perms[g];x=p[0];l=1
  while x!=0:x=p[x];l+=1
  if l==9:has9=True
 if has9:cyclic+=1;continue
 req(z in H,'diagonal missing')
 ok=False
 for g in H:
  if perms[g][0]//3==0:continue
  K=close([g,z])
  if len(K)==9 and {perms[h][0] for h in K}==set(range(9)) and mul[g][z]==mul[z][g] and mul[mul[g][g]][g]==one:ok=True
 req(ok,'regular elementary group missing');elementary+=1
print(json.dumps(dict(hyperplanes=40,subgroups=len(subgroups),transitive_subgroups=transitive,containing_regular_C9=cyclic,remaining_regular_C3xC3=elementary,status='PASS'),sort_keys=True,indent=2))
