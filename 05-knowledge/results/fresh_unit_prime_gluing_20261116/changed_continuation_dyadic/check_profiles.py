from collections import defaultdict
from pathlib import Path
f=Path(__file__).with_name('profiles.txt')
profiles={}
odd=defaultdict(set)
for l in f.read_text().splitlines():
 q,g,e,o,b=map(int,l.split())
 profiles[q,g,e,o]=b
 u=e+o
 u//=u&-u
 for i in range(9):
  if b>>i&1: odd[9*q+2*i,g].add(u)
G=(1,3,5,9,15,45)
assert len(profiles)==1175
assert max(b.bit_count() for b in profiles.values())==4
assert [(key,m) for key,b in profiles.items() for m in (7,56,448) if b&m==m]==[((163,15,1107,1122),448)]
assert profiles[79,5,545,550]==28
assert profiles[163,15,1107,1122]==450
def intersect(hs,g):return set.intersection(*(odd[h,g] for h in hs))
assert not any(intersect((1471,1481,1483),g) for g in G)
assert {g:intersect((719,727),g) for g in G if intersect((719,727),g)}=={3:{1137},9:{845,1051}}
assert not any(intersect((713,715,719,721,725,727),g) for g in (3,9))
print('Terminal profile windows and signed-zero screening checked.')
