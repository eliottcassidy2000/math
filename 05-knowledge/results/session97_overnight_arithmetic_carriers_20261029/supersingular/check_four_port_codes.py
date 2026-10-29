#!/usr/bin/env python3
"""A four-port hostile among eight cyclic Lagrangian-derived length32 codes."""
import contextlib,io,runpy
from pathlib import Path
from itertools import combinations
from collections import Counter
with contextlib.redirect_stdout(io.StringIO()):
 d=runpy.run_path(str(Path(__file__).with_name('check_cyclic_frobenius.py')))
CHECKS=0
def check(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
rank=d['rank'];basis_for=d['code_basis']
gens=[g for n,g,w in d['weight_rows'] if n==31]
bases={g:basis_for(g,31) for g in gens}
def image(v,a):return sum(((v>>i)&1)<<((a*i)%31) for i in range(31))
orbits=[];unseen=set(gens)
while unseen:
 g=min(unseen);orbit=set()
 for a in range(1,31):
  transformed=[image(v,a) for v in bases[g]]
  targets=[h for h in gens if rank(transformed+bases[h])==15]
  check(len(targets)==1,'unique multiplier image')
  orbit.add(targets[0])
 check(orbit<=unseen,'disjoint group orbits')
 unseen-=orbit;orbits.append(orbit)
check(sorted(map(len,orbits))==[2,6],'normalizer orbits2+6')
rows=[];profiles={};full_rows=set()
for g in gens:
 words=[0]
 for v in bases[g]:words += [w^v for w in words]
 extended=words+[w^((1<<32)-1) for w in words]
 octads=[w for w in extended if w.bit_count()==8]
 check(len(octads)==620,'620 octads')
 counts={t:Counter() for t in range(1,5)}
 for w in octads:
  support=[i for i in range(32) if w>>i&1]
  for t in range(1,5):counts[t].update(combinations(support,t))
 for t,value,number in [(1,155,32),(2,35,496),(3,7,4960)]:
  check(len(counts[t])==number and set(counts[t].values())=={value},'all named low-port counts agree')
 full_rows.add(frozenset(counts[4].items()))
 hist=Counter(counts[4].values());hist[0]=35960-len(counts[4])
 profile=tuple(sorted(hist.items()));profiles[g]=profile
 expected0={0:11160,1:12400,2:7440,3:3720,4:1240}
 expected1={0:0,1:34720,7:1240}
 check(dict(profile) in [expected0,expected1],'four-port profile')
 hole=next((s for s in combinations(range(32),4) if s not in counts[4]),None)
 rows.append((g,profile,hole))
for o in orbits:check(len({profiles[g] for g in o})==1,'profile constant within orbit')
check(len(set(profiles.values()))==2,'two arbitrary-coordinate isomorphism classes')
check(len(full_rows)==8,'four-port table faithfully separates all eight marked codes')
print('CHECKS',CHECKS)
print('Multiplier orbits:',sorted(sorted(o) for o in orbits))
print('Every named point/pair/triple lies in155/35/7 octads, respectively.')
for row in rows:print(row)
print('Exactly two S32-equivalence classes: multipliers prove within-orbit equivalence; four-port profiles separate orbits.')
print('Full weighted four-port tables are distinct for all eight marked codes; three-port tables coincide.')
