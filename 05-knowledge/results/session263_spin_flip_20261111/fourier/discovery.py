#!/usr/bin/env python3
"""Optional bounded replay of the exploratory eight-filter discovery."""
from itertools import combinations
from random import Random
from filters import split,residues,F8
rng=Random(263092)
for proposal in range(1000):
 bits=rng.getrandbits(28);a=[0]*8
 for k,(i,j) in enumerate(combinations(range(8),2)):
  if bits>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 p=split(tuple(a))
 if len(set(residues(p,8)))==1:
  assert proposal==795 and tuple(a)==F8
  print({'proposal':proposal,'rows':a,'K':p});break
else:raise AssertionError('Expected positive constructor missing')
