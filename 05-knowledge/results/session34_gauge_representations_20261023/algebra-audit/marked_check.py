#!/usr/bin/env python3
"""Direct labelled Hamiltonian path states, without normal-form action assumptions."""
import itertools as it,collections,json,math
from pathlib import Path
out=[]
for n in range(1,7):
 E=list(it.combinations(range(n),2));ix={e:k for k,e in enumerate(E)};N=len(E);counts=collections.Counter();total=0
 for path in it.permutations(range(n)):
  forced={tuple(sorted((a,b))):(a<b) for a,b in zip(path,path[1:])}
  free=[e for e in E if e not in forced]
  # Unique g such that g(reverse(path))=path.
  g=[0]*n
  for a,b in zip(reversed(path),path):g[a]=b
  for bits in range(1<<len(free)):
   arcs={e:v for e,v in forced.items()};arcs.update({e:bool(bits>>k&1) for k,e in enumerate(free)})
   fixed=True
   for (a,b),forward in arcs.items():
    ga,gb=g[a],g[b];image_forward=(not forward) if ga<gb else forward
    if arcs[tuple(sorted((ga,gb)))]!=image_forward:fixed=False;break
   if fixed:counts[tuple(g)]+=1
   total+=1
 m=(n-1)*(n-2)//2;q=(m+(n-1)//2)//2;z=2**(n//2)*math.factorial(n//2)
 assert total==math.factorial(n)*2**m
 assert set(counts.values())=={z*2**q}
 assert len(counts)*z==math.factorial(n)
 out.append(dict(n=n,marked_states=total,supported_permutations=len(counts),twisted_character_per_supported_permutation=z*2**q,chord_orbits=q))
Path(__file__).with_name('marked_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS direct marked states n<=6')
