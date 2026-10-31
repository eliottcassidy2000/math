from transfer import *
from itertools import product
import json
checks=0;seen=set()
for n in range(1,9):
 for w in product(*(range(1<<min(i,2)) for i in range(n))):
  h=count(w,2);ind=dp(w,2)
  if h!=ind:raise ValueError((w,h,ind))
  if h==47:raise ValueError('47')
  seen.add(h);checks+=1
witnesses=0;got=set()
from pathlib import Path
for line in Path(__file__).with_name('width3_witnesses.tsv').read_text().splitlines():
 h,*w=map(int,line.split())
 if h!=dp(w,3) or h!=count(w,3):raise ValueError((h,w))
 got.add(h);witnesses+=1
if got!={x for x in range(1,1654,2) if x not in (7,21)}:raise ValueError('coverage')
print(json.dumps({'width2_independent_checks':checks,'47_absent_through_n':8,'width3_independent_witnesses':witnesses,'width3_allowed_prefix':1653},sort_keys=True))
