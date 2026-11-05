"""Independent graph-DP check of the finite part of the 67 obstruction."""
from apex_transfer import *
from census import dpstats
from itertools import product
import json
raw,t,o=build_cycle_dfa();p,t,o=minimize_dfa(t,o)
minima={};small_values=set();checked=0;inventory={}
for m in range(1,8):
 inv=Counter()
 for w in product(range(4),repeat=m):
  state=p[0]
  for a in w:state=t[state][a]
  c=o[state]
  if c not in (1,2):continue
  h,gc=dpstats(adjacency(w))
  assert gc==c and counts(w)==(h,c)
  checked+=1;inv[c]+=1
  key=(m,c)
  if key not in minima or h<minima[key][0]:minima[key]=(h,w)
  if m<=6:small_values.add(h)
 inventory[m]=dict(inv)
assert minima[(7,1)][0]==75 and minima[(7,2)][0]==115
assert 67 not in small_values
assert counts((1,2,3,0,1))==dpstats(adjacency((1,2,3,0,1)))==(67,4)
print('independent low-cycle graph checks',checked)
print('minima',minima)
print('67 absent at base<=6; base7 minima 75 and 115; deletion lemma excludes every larger base')
json.dump({'checks':checked,'inventory':inventory,'minima':{str(k):v for k,v in minima.items()},'67_absent_through_base_6':True},open(__file__.replace('.py','.json'),'w'),indent=2)
