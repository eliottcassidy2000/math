#!/usr/bin/env python3
"""Exact W(D4) fixed-set Euler classifier; geometric and retention sidecars required."""
from itertools import permutations,product
from collections import deque
import json,hashlib
GATES=0
def check(ok,label):
 global GATES
 GATES+=1
 if not ok:raise RuntimeError(label)

# Exact even signed permutations of four coordinates, acting on eight signed letters.
letters=(1,2,3,4,-1,-2,-3,-4)
idx={x:i for i,x in enumerate(letters)}
elements=[]
for p in permutations((1,2,3,4)):
 for signs in product((-1,1),repeat=4):
  if signs[0]*signs[1]*signs[2]*signs[3]!=1:continue
  im=tuple(p[i]*signs[i] for i in range(4));images=im+tuple(-x for x in im)
  elements.append(tuple(idx[x] for x in images))
check(len(elements)==192 and len(set(elements))==192,'complete even signed permutation universe')
elements.sort();lookup={g:i for i,g in enumerate(elements)}
def mulp(a,b):return tuple(a[b[i]] for i in range(8))
mult=[[lookup[mulp(a,b)] for b in elements] for a in elements]
one=lookup[tuple(range(8))];neg=lookup[tuple(idx[-x] for x in letters)]
def reflection(i,j,negative=False):
 ims=list(letters)
 for s in (-1,1):
  ims[idx[s*i]]=s*j*(-1 if negative else 1)
  ims[idx[s*j]]=s*i*(-1 if negative else 1)
 return lookup[tuple(idx[x] for x in ims)]
a=reflection(1,2);b=reflection(2,3);c=reflection(3,4);d=reflection(3,4,True)
inverses=[next(j for j in range(192) if mult[g][j]==one) for g in range(192)]
def closure(gens):
 out={one};todo=[one]
 while todo:
  x=todo.pop()
  for h in gens:
   y=mult[x][h]
   if y not in out:out.add(y);todo.append(y)
 return frozenset(out)
check(closure((a,b,c,d))==frozenset(range(192)),'four declared reflections generate the entire group')
check(neg!=one and mult[neg][neg]==one,'nontrivial central involution')
for g in range(192):check(mult[neg][g]==mult[g][neg],'central minus identity')
for leaf in [a,c,d]:
 check(mult[leaf][leaf]==one and mult[b][b]==one,'simple reflections involutive')
 check(mult[mult[leaf][b]][leaf]==mult[mult[b][leaf]][b],'D4 braid relation')
for i,j in [(a,c),(a,d),(c,d)]:check(mult[i][j]==mult[j][i],'D4 commuting leaves')
base=closure((neg,b));seen={base:(neg,b)};queue=deque([base])
while queue:
 h=queue.popleft();gens=seen[h]
 # One representative per right coset suffices: <H,g h>=<H,g>.
 covered=set(h)
 for g in range(192):
  if g in covered:continue
  covered.update(mult[g][x] for x in h)
  k=closure(gens+(g,))
  if k not in seen:seen[k]=gens+(g,);queue.append(k)
check(len(seen)==26,'complete upward closure has twenty-six subgroups')
for hh,gg in seen.items():
 check(closure(gg)==hh and base<=hh,'explicit subgroup and retained base')
 for g in range(192):check(closure(gg+(g,)) in seen,'every one-element overgroup is retained')
e=mult[mult[b][c]][inverses[b]]
node0=mult[mult[inverses[e]][a]][e]
pairs=[(b,c),(b,d),(a,e),(node0,b),(a,d),(c,d)]
rows=[]
for h in seen:
 cosets=[];assignment={}
 for g in range(192):
  if g in assignment:continue
  coset=frozenset(mult[g][x] for x in h);check(len(coset)==len(h),'full left coset');j=len(cosets);cosets.append(coset)
  for x in coset:assignment[x]=j
 reps=[min(k) for k in cosets];deg=len(cosets)
 acts={q:[assignment[mult[q][x]] for x in reps] for q in {a,b,c,d,e,node0}}
 fixed=lambda qs:sum(all(acts[q][i]==i for q in qs) for i in range(deg))
 check(deg*len(h)==192 and len(assignment)==192,'all literal cosets exhausted')
 check(all(assignment[mult[neg][x]]==assignment[x] for x in range(192)),'central involution acts trivially on this coset set')
 aa=fixed((b,));check(aa>0,'reflection fixed-point supplier retained');nn=[fixed(p) for p in pairs];chi=3*deg-8*aa+sum(nn)
 rows.append({'orderH':len(h),'degree':deg,'a':aa,'pair_counts':nn,'chi':chi,'gens':seen[h]})
passrows=[r for r in rows if r['chi']==1 and r['a']>0]
check(len(passrows)==7,'all seven labelled survivor subgroups')
check(sorted(r['degree'] for r in passrows)==[1,4,4,4,4,4,4],'only degree one or four')
for rr in passrows:
 if rr['degree']==4:
  check(rr['a']==2 and rr['pair_counts'][:3]==[1,1,1],'complete degree-four cusp counts')
  check(sorted(rr['pair_counts'][3:])==[0,0,2],'complete degree-four node counts')
# Independent direct natural eight-letter control, with free central involution.
fixed8=lambda qs:sum(all(elements[q][i]==i for q in qs) for i in range(8))
counts8=[fixed8(pair) for pair in pairs];a8=fixed8((b,));chi8=24-8*a8+sum(counts8)
check(a8==4 and counts8==[2,2,2,0,0,4] and chi8==2,'natural eight-letter parity hostile')
check(all(elements[neg][i]!=i for i in range(8)),'hostile central action is free')
profiles=sorted({(r['degree'],r['a'],tuple(r['pair_counts']),r['chi']) for r in rows})
print('Literal W(D4) full-fixed-set Euler classifier PASS; geometric consumer conditional')
print('Universe:192 even signed permutations;26 overgroups of central minusI and one reflection')
print('Euler-one survivors: trivial degree1 and six labelled degree4 stabilizers')
print('Hostile: natural transitive degree8 has Euler2; central-free actions have even Euler')
print('Always-active gates:',GATES)
print('Semantic SHA256:',hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest())
