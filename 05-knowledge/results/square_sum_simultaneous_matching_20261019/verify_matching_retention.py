#!/usr/bin/env python3
"""Check table cut support and finite algebra controls; stdlib, no assertions."""
import hashlib,json,argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=Path(__file__).resolve().parent.parent);args=p.parse_args()
base=args.results/'seeded_eventual_cycle_theorem_20261019';checks=0;hashes={}
def check(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def read(path):
 raw=path.read_bytes();hashes[str(path.relative_to(args.results))]=hashlib.sha256(raw).hexdigest();return json.loads(raw)
expected={(f,p,d)for f in ('B','Aplus','Aminus')for p in(0,1)for d in range(24,73)}
splits={('B',0,48):0,('B',1,48):0,('Aminus',0,48):0,('Aminus',1,48):0,('Aplus',0,50):-2,('Aplus',1,46):-2};seen=set()
for path in sorted((base/'certified-pair-transitions').glob('*.json')):
 t=read(path);state=tuple(t[x]for x in('family','parity','residue'));check(state in expected and state not in seen,'state');seen.add(state)
 cut=[] if state not in splits else [splits[state]];check(t['split_prefix_blocks']==cut,'allowed prefix cut support')
 check(t['certificate_type']==('whole-block' if not cut else 'single-prefix-split'),'certificate type')
 check(len(t['words'])==2,'two outputs')
 for j,word in enumerate(t['words']):
  blocks={};heads=set();single=set()
  for atom in word:
   kind=atom[0];check(kind in('single','block','head','tail'),'typed atom')
   if kind=='single':check(len(atom)==2 and atom[1] not in single,'singleton');single.add(atom[1]);continue
   c=atom[1];check(-24<=c<=24,'offset')
   if kind=='head':check(len(atom)==2 and c in cut and c not in heads,'head');heads.add(c);continue
   check(len(atom)==4 and atom[3] in(-1,1),'orientation schema')
   check(c not in blocks and kind==('tail' if c in cut else 'block'),'exact block or tail');blocks[c]=kind
   check(atom[2]==int(49+(-1)**(state[1]+1)*c<=state[2]+j),'selected input')
  check(set(blocks)==set(range(-24,25)) and heads==set(cut) and single==set(range(1,25)),'partition')
check(seen==expected,'complete pair table')
final=read(args.results/'seeded_final_split_table_20261019/final_split_glue_certificates_20261019.json');seen=set()
for t in final:
 state=tuple(t[x]for x in('parity','residue','order0','order1'));check(state not in seen,'final duplicate');seen.add(state)
 check(t['cutparity']==1,'2 odd, 7 even');segments=set();singles=set()
 for atom in t['word']:
  check(atom[0] in('single','segment'),'final typed atom')
  if atom[0]=='single':check(len(atom)==2 and atom[1] not in singles,'final single');singles.add(atom[1]);continue
  check(len(atom)==5,'segment schema');_,c,i,s,o=atom
  check(-24<=c<=24 and s in(0,1) and o in(-1,1) and (c,s) not in segments,'segment');segments.add((c,s))
  check(i==int(49+(-1)**(state[0]+1)*c<=state[1]),'final selected input')
 check(segments=={(c,s)for c in range(-24,25)for s in(0,1)} and singles==set(range(1,25)),'final partition')
check(seen=={(p,d,a,b)for p in(0,1)for d in range(24,73)for a in(-1,1)for b in(-1,1)},'all final states')
def L(k):
 q=49**k;h=(49**(k-1)-1)//2
 return {a for a in range((q+1)//2,(3*q-1)//2+1)if abs(a-q)>h}
check(L(1)==set(range(25,74))-{49},'first level')
for k in(1,2):
 for sign in(-1,1):check({49*a+sign*c for a in L(k)for c in range(-24,25)}==L(k+1),'interval recursion both parity signs')
for k in range(1,10):
 q=49**k;check((107*q-1)//2==49*((107*(q//49)-1)//2)+24,'threshold recurrence')
 check(49*((49**(k-1)-1)//2)+24==(q-1)//2,'gap recurrence')
 # Both signs and both traversal orientations preserve each undirected lifted edge.
 for sigma in(-1,1):
  for c in range(-24,25):
   a=(q+1)//2;b=4*q-a;edge=(49*a+sigma*c,49*b-sigma*c)
   check(sum(edge)==4*49*q and set(edge)==set(edge[::-1]),'sum and reversal invariance')
print(json.dumps({'status':'passed','checks':checks,'pair_states':294,'final_states':392,'split_states':[{ 'family':f,'parity':p,'residue':d,'offset':c}for(f,p,d),c in sorted(splits.items())],'scope':'Cut-support and edge-retention audit; Hamiltonian validity remains supplied by dependency universal checkers. Segment semantics are the original-position prefix and 2–7 cuts defined in the theorem.','dependency_sha256':hashes},indent=2,sort_keys=True))
