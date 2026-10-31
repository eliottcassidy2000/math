#!/usr/bin/env python3
"""Universal parity check with a separate source colour for each offset pair."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'seeded_eventual_cycle_theorem_20261019'/'certified-pair-transitions'
def need(x,m):
 if not x:raise RuntimeError(m)

def profile(word,m):
 at=0;heads={};tails={};singletons={};cut_positions={}
 for atom in word:
  kind,c=atom[:2]
  if kind=='single':
   need(c not in singletons,'duplicate singleton');singletons[c]=at%2
   if c in (2,7):cut_positions[c]=at
   at+=1
  elif kind=='head':
   need(c not in heads,'duplicate head');heads[c]=at%2;at+=1
  else:
   _,c,i,o=atom;L=m+i
   need(c not in tails,'duplicate tail')
   shift=(at-(kind=='tail') if o==1 else at+L-1)%2
   tails[c]=shift
   if kind=='block':
    need(c not in heads,'duplicate whole head');heads[c]=shift
   at+=L-(kind=='tail')
 return heads,tails,singletons,cut_positions,at

states=offset_checks=0
for family in ('Aplus','Aminus'):
 for parity in (0,1):
  m=42+parity
  for d in range(24,73):
   record=json.loads((OLD/f'{family}_p{parity}_d{d}.json').read_text())
   a,b=[profile(w,m) for w in record['words']]
   need(a[:3]==b[:3],'per-offset head/tail/singleton colour equality')
   need(set(a[0])==set(range(-24,25))==set(a[1]),'all offsets individually retained')
   need(a[4]==49*m+d and b[4]==49*m+d+1,'typed output lengths')
   need(set(a[2])==set(range(1,25)),'singleton labels')
   for p in (a,b):need(abs(p[3][2]-p[3][7])==1 and p[3][2]%2==0,'fixed final cut')
   # Positions are affine sums of atom lengths; no source permutation enters.
   states+=1;offset_checks+=49
print(json.dumps({'status':'passed','states':states,'separate_offset_profiles':offset_checks,
                  'scope':'All complete A interface states; parity depends only on input-index parity, and named cut positions only on atom lengths, not source code.'},indent=2,sort_keys=True))
