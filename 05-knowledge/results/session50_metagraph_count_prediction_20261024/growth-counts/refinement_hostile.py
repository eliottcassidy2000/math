#!/usr/bin/env python3
"""Independent all-7! verification of the minimal finite refinement gap."""
from itertools import permutations
import json
T=(74,84,97,22,37,11,56)
n=len(T)
def require(ok,message):
    if not ok: raise RuntimeError(message)
require(all(T[v] & (1<<v)==0 for v in range(n)), 'no loops')
require(all(((T[v]>>w)&1)+((T[w]>>v)&1)==1 for v in range(n) for w in range(v)), 'tournament')
require([r.bit_count() for r in T]==[3]*7, 'regular hence one stable color cell')
group=[]
for g in permutations(range(n)):
    if all(((T[v]>>w)&1)==((T[g[v]]>>g[w])&1) for v in range(n) for w in range(n)):
        group.append(g)
orbits=sorted({tuple(sorted({g[v] for g in group})) for v in range(n)})
require(len(group)==3, 'automorphism order three')
require(sorted(map(len,orbits))==[1,3,3], 'three vertex orbits of sizes one three three')
require(all(g[6]==6 for g in group), 'marked root six fixed')
require(7*len(orbits)-sum(map(len,orbits))==14, 'refinement trace excess fourteen')
print(json.dumps({'rows':T,'tested_permutations':5040,'automorphisms':group,
                  'vertex_orbits':orbits,'stable_refinement_cell_size':7,
                  'root':6,'root_orbit_size':1,'trace_excess':14},sort_keys=True))
inflations=[]
for k in range(1,7):
    m=n*k
    rows=tuple(sum(1<<(w*k+j) for w in range(n) for j in range(k)
                   if ((v==w and i<j) or (v!=w and (T[v]>>w)&1)))
               for v in range(n) for i in range(k))
    require([row.bit_count() for row in rows]==[3*k+k-1-i for v in range(n) for i in range(k)],
            'inflation degrees')
    for v in range(n):
        for i in range(k):
            for j in range(k):
                count=sum((rows[v*k+i]>>(w*k+j))&1 for w in range(n))
                require(count==3+(i<j), 'stable clone-level neighbor counts')
    # Consecutive clones are exactly homogeneous pairs. Their graph has
    # canonical chain components of size k, verifying the inherited sidecar.
    observed=[]
    for a in range(m):
        for b in range(a+1,m):
            homogeneous=all(((rows[a]>>c)&1)==((rows[b]>>c)&1) for c in range(m) if c not in (a,b))
            if homogeneous: observed.append((a,b))
    expected=[(v*k+i,v*k+i+1) for v in range(n) for i in range(k-1)]
    require(observed==expected,'canonical inflated clone chains')
    inflations.append({'multiplier':k,'order':m,'stable_cells':k,'cell_size':7,
                       'fixed_root_clones':k,'homogeneous_pairs':len(observed)})
print(json.dumps({'inflation_controls':inflations},sort_keys=True))
print('PASS independent regular-seven hostile and inherited inflation controls')
