#!/usr/bin/env python3
"""Independent permutation audit of the seven-vertex refinement hostile."""
from itertools import permutations
import json
from pathlib import Path

T=(74,84,97,22,37,11,56)
N=len(T)

def check(ok,message):
    if not ok:
        raise RuntimeError(message)

def ref_cells(t):
    c=tuple(sum((row>>v)&1 for v in range(len(t))) for row in t)
    while True:
        names=sorted(set(c))
        sig=tuple((c[v],tuple(sum((t[v]>>w)&1 for w in range(len(t)) if c[w]==a)
                                for a in names)) for v in range(len(t)))
        index={s:i for i,s in enumerate(sorted(set(sig)))}
        d=tuple(index[s] for s in sig)
        if len(set(d))==len(set(c)):
            return sorted(sorted(v for v,x in enumerate(d) if x==a) for a in sorted(set(d)))
        c=d

for v in range(N):
    check(((T[v]>>v)&1)==0,'loop-free')
    for w in range(v):
        check(((T[v]>>w)&1)+((T[w]>>v)&1)==1,'exactly one arc')
aut=[]
for p in permutations(range(N)):
    if all(((T[v]>>w)&1)==((T[p[v]]>>p[w])&1) for v in range(N) for w in range(N)):
        aut.append(p)
orbits=sorted({tuple(sorted({p[v] for p in aut})) for v in range(N)})
check(len(aut)==3,'full permutation group order three')
check(sorted(map(len,orbits))==[1,3,3],'orbits are one plus three plus three')
check(ref_cells(T)==[list(range(N))],'regular tournament refinement remains one cell')
check({p[6] for p in aut}=={6},'reported marked root is fixed')
# Discrete partitions are equivariant as unlabelled partitions, but individual
# singleton cells are not invariant. Hence the theorem needs orbit-coarsening.
cycle=(2,4,1)
rot=(1,2,0)
check(all(((cycle[v]>>w)&1)==((cycle[rot[v]]>>rot[w])&1)
          for v in range(3) for w in range(3)),'C3 rotation is automorphism')
check(rot[0]!=0,'C3 singleton is not setwise invariant')
inflation_controls=[]
for k in range(1,7):
    vertices=[(v,i) for v in range(N) for i in range(k)]
    rows=tuple(sum(1<<b for b,(w,j) in enumerate(vertices)
                   if (v==w and i<j) or (v!=w and ((T[v]>>w)&1))) for v,i in vertices)
    cells=ref_cells(rows)
    expected=sorted([v*k+i for v in range(N)] for i in range(k))
    check(cells==expected,'inflated stable cells exactly clone positions')
    for i in range(k):
        orbit={p[6]*k+i for p in aut}
        check(orbit=={6*k+i},'fixed-base clones remain fixed under lifted base group')
    # The inherited maximal-transitive-module theorem supplies group completeness;
    # directly verify its distinguishing homogeneous-pair components here.
    observed=[]
    for a in range(len(rows)):
        for b in range(a):
            if all(((rows[a]>>c)&1)==((rows[b]>>c)&1)
                   for c in range(len(rows)) if c not in (a,b)):
                observed.append((b,a))
    check(set(observed)=={(v*k+i,v*k+i+1) for v in range(N) for i in range(k-1)},
          'homogeneous-pair edges recover exactly the clone chains')
    inflation_controls.append({'k':k,'order':N*k,'stable_cells':len(cells),'cell_size':7,
                               'fixed_base_clones':k,'trace_excess_from_three_base_orbits':14*k})
result={'status':'PASS','inflation_controls':inflation_controls,'rows':T,'automorphisms':aut,'orbits':orbits,
        'stable_refinement_cells':ref_cells(T),'marked_root':6,
        'root_orbit_size':1,'root_refinement_size':7,
        'all_rooted_orbit_states_trace':sum(map(len,orbits)),
        'all_rooted_orbit_states_refinement_bound':N*len(orbits),
        'trace_overestimate':N*len(orbits)-sum(map(len,orbits)),
        'partition_wording_control':'The discrete partition is equivariant, not cellwise Aut-invariant.'}
Path(__file__).with_name('refinement-witness.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
