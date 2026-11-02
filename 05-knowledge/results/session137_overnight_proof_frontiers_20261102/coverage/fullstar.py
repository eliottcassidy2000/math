"""Full one-vertex-replacement unit atlas via exact deletion-card incidence."""
from phase2 import canon,encode,check,intervals
from explore import tournament,H,strong
from itertools import combinations,permutations
from collections import defaultdict
from math import gcd
from functools import reduce
import json,hashlib

def induced(a,vs):
 return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)

def analyse(n):
 atlas=json.load(open('phase2_order%d.json'%n));hs={int(c):h for c,h in atlas['class_H'].items()};decks={};parents=defaultdict(set)
 for c in hs:
  a=tournament(n,c);deck=[]
  for r in range(n):
   b=induced(a,[i for i in range(n) if i!=r]);card=canon(n-1,encode(b));deck.append(card);parents[card].add(c)
  decks[c]=deck
 p={c:c for c in hs}
 def root(c):
  while p[c]!=c:p[c]=p[p[c]];c=p[c]
  return c
 forest=[]
 for card,cs in sorted(parents.items()):
  for c,d in combinations(sorted(cs),2):
   if abs(hs[c]-hs[d])<=2 and root(c)!=root(d):
    p[root(c)]=root(d);forest.append((c,d,card))
 comps=defaultdict(list)
 for c in hs:comps[root(c)].append(c)
 comp=sorted([{'classes':sorted(cs),'image':sorted({hs[c] for c in cs})} for cs in comps.values()],key=lambda d:(d['image'],d['classes']))
 report={'order':n,'fullstar_components':[{'classes':len(c['classes']),'image':c['image']} for c in comp],'global_intervals':intervals(set(hs.values()))}
 print(json.dumps(report,sort_keys=True,separators=(',',':')))
 data={'report':report,'components':comp,'decks':decks,'class_H':hs,'forest':forest}
 with open('fullstar_order%d.json'%n,'w') as f:json.dump(data,f,sort_keys=True,separators=(',',':'));f.write('\n')
 return data

ledger=[analyse(n) for n in [6,7]]
# A tiny consequence certificate, independent of all graph connectivity searches.
a=(6,112,74,35,13,21,57);check(H(a)==147 and strong(a),'H147 hostile');code=encode(a);n=7;pairs=list(combinations(range(n),2));cards=ledger[1]['decks'][code]
matching145=[c for c,h in ledger[1]['class_H'].items() if h==145]
othercards=sorted({card for c in matching145 for card in ledger[1]['decks'][c]})
check(not set(cards)&set(othercards),'H147/H145 deletion decks disjoint')
replacements=[]
for r in range(n):
 es=[k for k,ij in enumerate(pairs) if r in ij]
 for mask in range(1<<len(es)):
  d=code^sum(1<<k for i,k in enumerate(es) if mask>>i&1);b=tournament(n,d)
  if strong(b):replacements.append(H(b))
check(145 not in replacements and 149 not in replacements,'all 448 literal vertex replacements miss adjacent levels')
check(set(h for h in replacements if abs(h-147)<=2)=={147},'fullstar unit isolation')
three=[]
for size in [1,2,3]:
 for ks in combinations(range(len(pairs)),size):
  d=code^sum(1<<k for k in ks);b=tournament(n,d)
  if H(b)==145 and strong(b):
   check(size==3,'no one/two arbitrary-arc repair');three.append({'flipped_pairs':[pairs[k] for k in ks],'code':d,'masks':b})
check(len(three)==4,'all minimum three-arc repairs')
card_rows=[]
for r in range(n):
 base=induced(a,[i for i in range(n) if i!=r]);hbase=H(base)
 # Literal completion spectra, not the primary response-carrier formula.
 image=[]
 for S in range(64):
  child=tuple(base[i]|((1<<6) if not(S>>i&1) else 0) for i in range(6))+(S,)
  image.append(H(child))
 d=reduce(gcd,[h-hbase for h in image]);check(d==2 and strong(base),'all deletion cards strong and lattice-saturated')
 card_rows.append({'deleted_vertex':r,'H':hbase,'ear_gcd':d})
aut=sum(all(bool(a[i]>>j&1)==bool(a[p[i]]>>p[j]&1) for i in range(n) for j in range(i+1,n)) for p in permutations(range(n)))
check(aut==1,'rigid hostile')
cert={'automorphisms':aut,'deletion_card_rows':card_rows,'H147_code':code,'H147_masks':a,'H147_card_codes':cards,'H145_class_codes':matching145,'H145_card_union':othercards,'literal star presentations':448,'strong replacement H values':sorted(set(replacements)),'minimum Hamming repairs':three}
print(json.dumps({'deletion_deck_hostile':cert},sort_keys=True,separators=(',',':')))
with open('deck_barrier.json','w') as f:json.dump(cert,f,sort_keys=True,separators=(',',':'));f.write('\n')
ledger.append(cert)
print('PASS semantic_sha256='+hashlib.sha256(json.dumps(ledger,sort_keys=True,separators=(',',':')).encode()).hexdigest())
