#!/usr/bin/env python3
"""Independent large-exponent and finite critical-cycle controls."""
import json
import word_projected_search as independent
fs=[(2,2),(3,20),(5,20)]
b=4*3**20*5**20
# Use the supplied known-prime factorization; never trial-divide this large b.
independent.factor=lambda candidate: fs if candidate==b else (_ for _ in ()).throw(ValueError('Unexpected denominator'))
N=10
pairs=[(0,k) for k in range(N)]+[(k,k-1) for k in range(1,N)]
r=independent.search(b,pairs)
if r['cycles'] or r['unsigned_pair_count']!=6859 or r['primitive_triangle_candidates']!=16:
 raise ValueError('Large-exponent mismatch')
# This independently reconstructed orbit search does not use the author's
# real-escape guard, and therefore checks that its new guard did not lose cycles.
critical=[]
for p,c0,m in ((3,2,2),(5,1,3)):
 levels=[]
 for k in range(1,5):
  modulus=p**k;count=0
  for c in range(c0,modulus,p):
   image=[(x*x+c)%modulus for x in range(modulus)]
   seen=set();cycles=[]
   for root in range(modulus):
    if root in seen:continue
    path=[];index={};x=root
    while x not in seen and x not in index:
     index[x]=len(path);path.append(x);x=image[x]
    if x in index:cycles.append(path[index[x]:])
    seen.update(path)
   if len(cycles)!=1 or len(cycles[0])!=m:
    raise ValueError(('critical-cycle mismatch',p,k,c,cycles))
   count+=1
  levels.append(dict(depth=k,modulus=modulus,parameter_residue_classes=count))
 critical.append(dict(p=p,c_mod_p=c0,only_exact_period=m,levels=levels))
print(json.dumps(dict(status='PASS',large_exponent_independent_search_without_escape_guard=r,critical_cycle_complete_finite_graph_controls=critical),indent=2))
