#!/usr/bin/env python3
"""Finite set-system lemma checks only, not a poset census."""
from itertools import combinations_with_replacement
import json
counts={}
for universe,rows in [(2,2),(3,2),(3,3),(4,4)]:
 sets=[s for s in range(1<<universe) if s.bit_count()>=2];checked=0
 for family in combinations_with_replacement(sets,rows):
  if not any((a&~b).bit_count()<=1 and (b&~a).bit_count()<=1 for i,a in enumerate(family) for b in family[i+1:]):raise RuntimeError('Near-pair lemma failed')
  checked+=1
 counts[f'{universe} elements, {rows} rows']=checked
# Five two-element subsets on universes of size at most five cannot be disjoint.
for universe in range(2,6):
 sets=[s for s in range(1<<universe) if s.bit_count()==2];checked=0
 for family in combinations_with_replacement(sets,5):
  if not any(a&b for i,a in enumerate(family) for b in family[i+1:]):raise RuntimeError('Disjointness bound failed')
  checked+=1
 counts[f'{universe} elements, five 2-sets']=checked
print(json.dumps({'status':'EXACT_SET_SYSTEM_CONTROLS_PASSED','scope':'Only the elementary finite set-system lemmas used in the symbolic height-two seven-thin theorem. No poset census.','families_checked':counts},indent=2))
